import cupy as cp
import time
import math
import numpy as np
import pandas as pd
from csv import reader
from scipy import special
from scipy.optimize import minimize
from scipy.optimize import Bounds

import matplotlib.pyplot as plt
from scipy.stats import chi2
import gc
import inspect

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import minimize
from scipy.special import betaln
import time

import tools

# Kernel Functions (in C++ Codes)

# matrix multiplication 
dot_product_bb = cp.ReductionKernel(
    'T x, T y',        
    'T z',         
    'x * y',         
    'a + b',          
    'z = a',              
    '0',                   
    'dot_product_bb'          
)

# Kernel function for sigmoid function
sigmoid2_bb = cp.ElementwiseKernel(
    'float32 x',
    'float32 mu',
    '''
    mu = 1 / (1 + exp(- x ));
    ''',
    'sigmoid2_bb'
)

# kernel function for the difference of two digamma function 
# Sicne the digamma function can be represented as a infinite sum
# the difference can be represented as the finite summation for each gene i and each cell j
e_step1_bb = cp.ElementwiseKernel(
    'float32 e_mu, float64 e_phi, uint16 e_y, uint16 e_n',
    'float32 v2',
    '''
    double temp1_1 = 0.0;
    double temp1_2 = 0.0;

    double x_small_1 = e_mu / e_phi;
    double x_small_2 = (1.0 - e_mu) / e_phi;

    for (int i = 0; i < e_y; i++) {
        temp1_1 += 1.0 / (i + x_small_1);
    }
    for (int i = 0; i < (e_n - e_y); i++) {
        temp1_2 += 1.0 / (i + x_small_2);
    }
    v2 = (float)((temp1_1 - temp1_2) * e_mu * (1.0 - e_mu) / e_phi);
    ''',
    'e_step1_bb'
)


# kernel function for the difference of two Trigamma function 
# Sicne the trigamma function can be represented as a infinite sum
# the difference can be represented as the finite summation for each gene i and each cell j
e_step2_bb = cp.ElementwiseKernel(
    'float32 e_mu, float64 e_phi, uint16 e_y, uint16 e_n',
    'float32 v2',
    '''
    double temp1_1 = 0.0;
    double temp1_2 = 0.0;
    double temp2_1 = 0.0;
    double temp2_2 = 0.0;

    double x_small_1 = e_mu / e_phi;
    double x_small_2 = (1.0 - e_mu) / e_phi;

    for (int i = 0; i < e_y; i++) {
        temp1_1 += 1.0 / (i + x_small_1);
    }
    for (int i = 0; i < (e_n - e_y); i++) {
        temp1_2 += 1.0 / (i + x_small_2);
    }
    double temp1 = temp1_1 - temp1_2;

    for (int i = 0; i < e_y; i++) {
        temp2_1 -= 1.0 / ((i + x_small_1) * (i + x_small_1));
    }
    for (int i = 0; i < (e_n - e_y); i++) {
        temp2_2 -= 1.0 / ((i + x_small_2) * (i + x_small_2));
    }
    double temp2 = temp2_1 + temp2_2;
    v2 = (temp2 * e_mu * e_mu * (1.0 - e_mu) * (1.0 - e_mu) / (e_phi * e_phi) +
          temp1 * (1.0 - 2.0 * e_mu) * e_mu * (1.0 - e_mu) / e_phi);
    v2 = (float)(v2);
    ''',
    'e_step2_bb'
)

# Kernel function for the difference of two beta function
beta2_kernel_bb = cp.ElementwiseKernel(
    'float32 mu, uint16 y, float32 phi, uint16 n',
    'float32 bt2',
    '''
    float a1 = mu / phi;
    float a2 = (1 - mu) / phi;
    if (n < 0.1){
      bt2 = 0;
    }
    else if ( y == n){
      bt2 = lgamma(a1+y) - lgamma(a1+a2+n) - lgamma(a1) + lgamma(a1+a2);
    }
    else{
      bt2 = lgamma(a1+y) + lgamma(a2+n-y) - lgamma(a1+a2+n) - lgamma(a1) - lgamma(a2) + lgamma(a1+a2);
    }
    ''',
    'beta2_kernel_bb'
)

# kernel function of Q-function, the difference of two beta functions
m_step1_bb = cp.ElementwiseKernel(
    'float32 linpred, float32 rand, float32 prec, float32 zmat, float32 phi, uint16 y, uint16 n',
    'float32 beta2',
    '''
    float mu = 1 / (1 + exp(-(linpred + rand + zmat / prec)));
    mu = (mu > 1.0f - 1e-5f) ? (1.0f - 1e-5f) : ((mu < 1e-5f) ? 1e-5f : mu);

    float alpha_1 = mu / phi;
    float alpha_2 = (1 - mu) / phi;

    if (n < 0.1) {
        beta2 = 0;
    } else if (y == n) {
        beta2 = lgamma(alpha_1 + y) - lgamma(alpha_1 + alpha_2 + n) - lgamma(alpha_1) + lgamma(alpha_1 + alpha_2);
    } else {
        beta2 = lgamma(alpha_1 + y) + lgamma(alpha_2 + n - y) - lgamma(alpha_1 + alpha_2 + n) - lgamma(alpha_1) - lgamma(alpha_2) + lgamma(alpha_1 + alpha_2);
    }
    ''',
    'm_step1_bb'
)

# kernel function for the derivative of \beta_0 and \beta_1 
m_step2_bb = cp.ElementwiseKernel(
    'float32 linpred, float32 rand, float32 prec, float32 zmat, float32 phi, uint16 y, uint16 n',
    'float32 temp_b',
    '''
    // Calculate mu using the sigmoid function
    float mu = 1 / (1 + exp(-(linpred + rand + zmat / prec)));
    mu = (mu > 1.0f - 1e-5f) ? (1.0f - 1e-5f) : ((mu < 1e-5f) ? 1e-5f : mu); // Clipping to avoid extreme values

    // Calculate alpha_1 and alpha_2 using the input phi
    float alpha_1 = mu / phi;
    float alpha_2 = (1 - mu) / phi;

    // Compute temp_b
    if (n < 0.1) {
        temp_b = 0;
    } else {
        float temp_1 = 0.0;
        for (int i = 0; i < y; i++) {
            temp_1 += 1 / (i + alpha_1);
        }
        float temp_2 = 0.0;
        for (int i = 0; i < n - y; i++) {
            temp_2 += 1 / (i + alpha_2);
        }
        temp_b = (temp_1 - temp_2) * mu * (1 - mu);
    }
    ''',
    'm_step2_bb'
)

# Kernel function for the derivative of phi
m_step3_bb = cp.ElementwiseKernel(
    'float32 linpred, float32 rand, float32 prec, float32 zmat, float32 phi, uint16 y, uint16 n',
    'float32 temp_phi',
    '''
    // Calculate mu using the sigmoid function
    float mu = 1 / (1 + exp(-(linpred + rand + zmat / prec)));
    mu = (mu > 1.0f - 1e-5f) ? (1.0f - 1e-5f) : ((mu < 1e-5f) ? 1e-5f : mu); // Clipping to avoid extreme values

    // Calculate alpha_1 and alpha_2 using the input phi
    float alpha_1 = mu / phi;
    float alpha_2 = (1 - mu) / phi;

    // Compute temp_phi
    if (n < 0.1) {
        temp_phi = 0;
    } else {
        float temp_1 = 0.0;
        for (int i = 0; i < y; i++) {
            temp_1 += 1 / (i + alpha_1);
        }
        float temp_2 = 0.0;
        for (int i = 0; i < n - y; i++) {
            temp_2 += 1 / (i + alpha_2);
        }
        float temp_3 = 0.0;
        float phi_inv = 1 / phi;
        for (int i = 0; i < n; i++) {
            temp_3 += 1 / (i + phi_inv);
        }
        temp_phi = -temp_1 * mu - temp_2 * (1 - mu) + temp_3;
    }
    ''',
    'm_step3_bb'
)



def lbfgs_mul_bb(fun, x0, maxk=20, step_factor=0.55, sigma=0.4, *args, **kwargs):
  '''
  # Tengfei's Optimizer
  # This optimizer is actually the BFGS algorithm. I rewrite the algorithm in in Cupy and extend the algothrm to 
  # matrix operation so that the optimizer could optmize multiple functions at the same time 
  # This optimizer might not perfect and very robust compared with the optimizer in Scipy, but its performance on Cuomo data and 
  # OneK1K data are pretty good, and it could optimize multiple genes at the same time which is an unique feature
  #Paramters: 
  # fun: the function to be optimized
  # x0: initial values for the parameters of the fun
  # maxk: the maximum number of FBGS iterations
  # step_factor: a paramter of line search, it controls the speed to decrease step length to find best step length
  # sigma: a paramter of line search, it controls the Armijo rule 
  '''

  # prepare some initial data
  num_q, num_param = x0.shape
  H0 = cp.eye(num_param).reshape(1, num_param, num_param)  # use an identity matrix as the initial matrix for Hessian matrix
  H0 = cp.tile(H0, (num_q, 1, 1))  # increase the Hessian matrix to three dimentions, since we need to optimze multiple functions at the same time
  H = H0.copy()
  x = 0
  val_new = 0

  k = 1  # the first direction
  func, gk = fun(x0, *args, **kwargs)
  gk = gk.reshape(gk.shape[0], 1, gk.shape[1])  # inital gradient
  dk = -cp.sum(H0 * gk, axis=2)  # initial direction

  # Start BFGS iteration
  while k <= maxk:

    # We first need to determine the step length
    step_find = cp.zeros((num_q))
    step_length = cp.ones((num_q, 1))
    val, gk = fun(x0, *args, **kwargs)
    line_iter = 0
    oldf, old_grad = val, gk

    # Line Search 
    while line_iter < 15: # we use 15 iterations to determine the step length
      newf, new_grad = fun((x0 + step_length * dk).astype(cp.float32), *args, **kwargs)
      addtion_diff = sigma * step_length * cp.sum(gk * dk, axis=1, keepdims=True)
      addtion_diff = addtion_diff.reshape(-1)  
      step_find = step_find.reshape(-1)
      
      # Armijo rule 
      step_find = cp.where(newf < oldf + addtion_diff, 2, step_find)
      step_find = cp.where((newf < (oldf + addtion_diff) * 1.2) & (newf > oldf + addtion_diff), 1, step_find)
      step_find = step_find.reshape(-1, 1)
      
      # for genes statisfying Armijo rule, we do not change the step length 
      # for genes not statisfying Armijo rule, we continue to descrease the step length
      step_length = cp.where(step_find == 1, step_length * step_factor, step_length)
      step_length = cp.where(step_find == 0, step_length * (step_factor * step_factor * step_factor), step_length)
      line_iter += 1
      
      num_converge = cp.sum(step_find == 2).item()
      if num_converge > (num_q * 0.99): # If we have found best step length for 99% genes, we stop our work
        break
    
    # if we fail to find appropraite step length for some genes, we need to decrease the step length to small enough
    step_length = cp.where(step_find != 2, step_length * 0.1, step_length) 

    # usde the step length to update the parameters    
    x = (x0 + step_length * dk).astype(cp.float32)
    x[:,:-1] = cp.clip(x[:,:-1], -100, 100)
    x[:,-1] = cp.clip(x[:,-1], 0.0001, 100)

    val_new, gk_new = fun(x, *args, **kwargs)
    sk = (x - x0).reshape(num_q, num_param, 1)
    yk = (gk_new - gk).reshape(num_q, num_param, 1)
    factor = cp.sum(sk * yk, axis=1, keepdims=True) + 1e-8

    # Update the Hessian Matrix
    matrix1 = cp.eye(num_param).reshape(1, num_param, num_param) - sk * yk.transpose(0, 2, 1) / factor
    matrix2 = cp.eye(num_param).reshape(1, num_param, num_param) - yk * sk.transpose(0, 2, 1) / factor
    matrix3 = sk * sk.transpose(0, 2, 1) / factor
    H = matrix1 @ H @ matrix2 + matrix3
    gk_new = gk_new.reshape(gk_new.shape[0], gk_new.shape[1], 1) 
    new_dk = -cp.matmul(H, gk_new).squeeze(-1)

    # Updarte the direction for parameters
    dominator = cp.sum(sk * yk, axis=1)
    dk = cp.where(dominator > 1e-10, new_dk, dk)  
    num_bad = cp.sum(cp.where(dominator <= 0, 1, 0))
    k += 1
    x0 = x
  return val_new, x

# Q-function for BB model
def cu_q_bb(cu_param, cu_x_q, cu_y_q, cu_n_q, cu_zmat_q, cu_ghq_weights_q, cu_randint_q, cu_randint_prec_q, cu_num_genes_q):
    '''
    # Q function for BB model: num_gene is the number of genes, num_cells is the number of cells
    # cu_param: this is a matrix (num_gene x 3), represents the parameters for num_gene genes; the first two columns are b0 and b1, the third column is phi
    # cu_x_q: 2 x num_cells, the covariates with first columns as all 1. If it is the null hypothesis, X only has one column with all 1
    # cu_y_q: num_gene x num_cells, the alternative count data
    # cu_n_q: num_gene x num_cells, the total read count data
    # cu_zmat_q: num_gene x num_cells x 3, it is used to approxiamate the integral in the Q-function (the z_m term in the formula of Q function)
    # cu_ghq_weights_q: 1 x 3, it is used to approxiamate the integral in the Q-function
    # cu_randint_q: num_cells x num_individuals, it includes all randint effects for all individuals and all genes
    # cu_randint_prec_q: num_cells x num_individuals, it is used Gauss-Hermite quadrature in the formula  (the \sigma_i term in the formula of Q function)
    # cu_num_genes_q: the number of genes
    '''
    # First, we need to change dimensions of some parameters
    cu_b_q = cu_param[:,:-1]
    cu_b_q = cu_b_q.reshape(cu_b_q.shape[0],1,cu_b_q.shape[1]) # cu_b shape[num_gene,1,2]
    cu_phi_q = cu_param[:,-1]
    cu_phi_q = cp.abs(cu_phi_q.reshape(cu_phi_q.shape[0],1,1)) # cu_phi_shape: [num_gene,1,1], we take the absolute values  
    cu_y_q = cu_y_q.reshape(cu_y_q.shape[0],cu_y_q.shape[1],1)
    cu_n_q = cu_n_q.reshape(cu_n_q.shape[0],cu_n_q.shape[1],1)

    cu_linpred_q = dot_product_bb(cu_x_q, cu_b_q, axis = 2)
    cu_linpred_q = (cu_linpred_q).reshape(cu_num_genes_q,cu_n_q.shape[1],1)

    # Compute the Q function results
    cu_beta2_chunk = cp.zeros((cu_num_genes_q, 3))
    for i in range(3):
      cu_beta2 = m_step1_bb(cu_linpred_q[:,:,0], cu_randint_q[:,:,0], cu_randint_prec_q[:,:,0], cu_zmat_q[:,:,i], cu_phi_q[:,:,0], cu_y_q[:,:,0], cu_n_q[:,:,0])
      cu_beta2 = cp.sum(cu_beta2,axis=1)
      cu_beta2_chunk[:,i] = cu_beta2 * cu_ghq_weights_q[:,i]
      del cu_beta2
    cu_result_new = - cp.sum(cu_beta2_chunk, axis=1) # Q function results with shape [num_gene, 1]

    #######################################################################################
    # Compute the derivatives of \beta0 and \beta1
    cu_b_chunk = cp.zeros((cu_num_genes_q, cu_x_q.shape[2], 3))
    for i in range(0,cu_x_q.shape[2]):
      for j in range(0,3):
        cu_b_new = m_step2_bb(cu_linpred_q[:,:,0], cu_randint_q[:,:,0], cu_randint_prec_q[:,:,0], cu_zmat_q[:,:,j], cu_phi_q[:,:,0], cu_y_q[:,:,0], cu_n_q[:,:,0])
        cu_b_new *= cu_x_q[:,:,i]
        cu_b_new = cp.sum(cu_b_new, axis=1) / cu_phi_q.reshape(-1)
        cu_b_chunk[:,i,j] = cu_b_new * cu_ghq_weights_q[:,j]
        del cu_b_new
    cu_b_temp = cp.sum(cu_b_chunk, axis = 2)
    
    # compute the derivatives of \phi
    cu_phi_chunk = cp.zeros((cu_num_genes_q, 3))
    for i in range(3):
      cu_phi_new = m_step3_bb(cu_linpred_q[:,:,0], cu_randint_q[:,:,0], cu_randint_prec_q[:,:,0], cu_zmat_q[:,:,i], cu_phi_q[:,:,0], cu_y_q[:,:,0], cu_n_q[:,:,0])
      cu_phi_chunk[:,i] = cp.sum(cu_phi_new, axis=1) / np.square(cu_phi_q.reshape(-1))
      del cu_phi_new
    cu_phi_temp = cp.dot(cu_phi_chunk, cu_ghq_weights_q.reshape(-1,1))

    # combine derivatives and return results
    cu_deri = - cp.concatenate((cu_b_temp, cu_phi_temp),axis=1)
    return cu_result_new, cu_deri


def cu_bb_lkl_agq_fixvar(cu_b_phi, cu_sigm2, cu_x, cu_randint, cu_id_data, cu_y, cu_n, cu_num_genes):
  # This function is used to compute the log-likelihood
  cu_b = cu_b_phi[:, :-1]
  cu_b = cu_b.reshape(cu_b.shape[0], 1, cu_b.shape[1]) # shape[num_genes,1,2]
  cu_phi = cu_b_phi[:,-1]
  cu_phi = cu_phi.reshape(-1,1)
  cu_randint = cp.zeros_like(cu_randint)
  cu_linpred = dot_product_bb(cu_x, cu_b, axis = 2)

  for i in range(3):
    cu_df = cu_compute_randint_deriv_bb(cu_linpred, cu_randint, cu_id_data, cu_y, cu_n, cu_phi)
    cu_randint = cu_randint - 0.9 * (cu_df[:cu_num_genes,:] - cu_randint / cu_sigm2) / (cu_df[cu_num_genes:,:] - 1 / cu_sigm2)
  cp.cuda.Stream.null.synchronize()
  end_llkl_1 = time.time()

  cu_df = cu_compute_randint_deriv_bb(cu_linpred, cu_randint, cu_id_data, cu_y, cu_n, cu_phi, cu_inlcude_driv0=True)
  cu_randint_prec = cp.abs(cu_df[cu_num_genes:cu_num_genes*2,:] - 1 / cu_sigm2)
  cu_fval = cu_df[cu_num_genes*2:,:] - cu_randint * cu_randint / cu_sigm2 / 2

  llkl = cp.sum(cu_fval - 0.5 * cp.log(cu_randint_prec),axis=1) - (int(cu_randint.shape[1]) * cp.log(cu_sigm2)/2).reshape(-1)
  return llkl

def cu_compute_randint_deriv_bb(cu_linpred, cu_e_randint, cu_id_data, cu_e_y, cu_e_n, cu_e_phi, cu_inlcude_driv0 = False):
  # Compute first and second derivatives of
  # sum(lbeta(mu/phi+y,(1-mu)/phi+n-y)-lbeta(mu/phi,(1-mu)/phi)) wrt random intercept
  lin = (cp.matmul(cu_e_randint, cu_id_data)).astype(cp.float32)
  lin += cu_linpred
  cu_e_mu = sigmoid2_bb(lin)
  del lin # clear the used memeory by lin

  cu_e_mu = cp.clip(cu_e_mu,0.0000001, 0.9999999) # bound the range of cu_e_mu

  # temporary transform the data type to cp.float64
  cu_e_phi = cu_e_phi.astype(cp.float64)
  cu_v1_new = e_step1_bb(cu_e_mu, cu_e_phi, cu_e_y, cu_e_n)
  cu_v1 = cp.matmul(cu_v1_new, cu_id_data.T)
  del cu_v1_new # clear the used memeory

  cu_v2_new = e_step2_bb(cu_e_mu, cu_e_phi, cu_e_y, cu_e_n)
  cu_v2 = cp.matmul(cu_v2_new, cu_id_data.T)
  del cu_v2_new # clear the used memeory

  if cu_inlcude_driv0:
    cu_e_phi = cu_e_phi.astype(cp.float32) 
    cu_v3 = beta2_kernel_bb(cu_e_mu, cu_e_y, cu_e_phi, cu_e_n)
    cu_v3 = cp.matmul(cu_v3, cu_id_data.T)
    cu_v =  cp.concatenate((cu_v1,cu_v2,cu_v3),axis=0)
  else:
    cu_v =  cp.concatenate((cu_v1,cu_v2),axis=0)
  return cu_v.astype(cp.float32) 

def cu_vem_bb(cu_b_phi, cu_sigm2, cu_x, cu_randint, cu_randint_prec, cu_id_data, cu_y, cu_n, cu_ghq_weights, cu_ghq_nodes, cu_n_lap, cu_num_genes, cu_iteration, null):
  # Get currenr parameters b0, b1, phi)
  cu_b = cu_b_phi[:, :-1]
  cu_b = cu_b.reshape(cu_b.shape[0], 1, cu_b.shape[1]) # shape[num_genes,1,2]
  cu_phi = cu_b_phi[:,-1]
  cu_phi = cu_phi.reshape(-1,1)
  cu_linpred = dot_product_bb(cu_x, cu_b, axis = 2)

  ############
  # E - step #
  ############
  for i in range(2):
    cu_df = cu_compute_randint_deriv_bb(cu_linpred, cu_randint, cu_id_data, cu_y, cu_n, cu_phi)
    cu_randint = cu_randint - 0.9 * (cu_df[:cu_num_genes,:] - cu_randint / cu_sigm2) / (cu_df[cu_num_genes:,:] - 1 / cu_sigm2)
  cu_df = cu_compute_randint_deriv_bb(cu_linpred, cu_randint, cu_id_data, cu_y, cu_n, cu_phi)
  cu_randint_prec = cp.abs(cu_df[cu_num_genes:,:] - 1 / cu_sigm2)
  del cu_linpred # clear the used memeory by cu_linpred

  ############
  # M - step #
  ############
  cu_zmat = (cu_ghq_nodes.reshape(1,1,-1)).astype(cp.float32)
  cu_sigm2_test = cu_randint * cu_randint + cp.reciprocal(cu_randint_prec)
  rand_indicator = (cu_randint != 0.0).astype(cp.float32)
  cu_sigm2 = cp.sum(cu_sigm2_test * rand_indicator, axis = 1) / cp.sum(rand_indicator, axis=1)
  cu_sigm2 = cu_sigm2.reshape(-1,1)

  cu_randint_vector = (cp.matmul(cu_randint, cu_id_data)).reshape(cu_num_genes, cu_n.shape[1], 1) 
  cu_precision_vector = (cp.matmul(cp.sqrt(cu_randint_prec), cu_id_data)).reshape(cu_num_genes, cu_n.shape[1], 1) 
  cu_randint_vector = cu_randint_vector.astype(cp.float32)
  cu_precision_vector = cu_precision_vector.astype(cp.float32)


  q_func, bphi = lbfgs_mul_bb(cu_q_bb,cu_b_phi, maxk=10, step_factor=0.55, sigma=0.4,
                               cu_x_q=cu_x, cu_y_q=cu_y, cu_n_q=cu_n, cu_zmat_q=cu_zmat, cu_ghq_weights_q=cu_ghq_weights,
                               cu_randint_q=cu_randint_vector, cu_randint_prec_q=cu_precision_vector, cu_num_genes_q=cu_num_genes)
  cu_b_phi = cp.asarray((bphi))
  cu_b_phi = cu_b_phi.astype(cp.float32)

  # After one EM iteration, compute the updated log-likelihood
  llkl = cu_bb_lkl_agq_fixvar(cu_b_phi, cu_sigm2, cu_x, cu_randint, cu_id_data, cu_y, cu_n, cu_num_genes)
  llkl = cp.asnumpy(llkl)
  # print(f"llkl: {llkl}")
  return cu_b_phi, cu_sigm2, llkl, cu_randint, cu_randint_prec   # return all results

def VEM_bb(gene_index, num_iteration, min_iter, ynxid_data, cu_id_data, cu_n_lap, initial_param, null, cum_time):
  '''
  # gene_index: the index of gene you want to analyze by DAESC-GPU
  # num_iteration: the maximum number of EM iteration
  # min_iter: the minmum number of EM iteration for all genes
  # ynxid_data: the count sequence, it is a matrix with genes as rows and cells as columns. For example,
  #.  if the ynxid_data is the a matrix of (2a+2) x b, the first a rows is the alternative count matrix Y, the (a+1)th to (2a)th rows
  #.  is the count matrix of the total read count. The (2a+1)th row is the covariates X, the last row (2a+2)th row is the indicator of 
  #.  individuals
  # cu_id_data: one-hot matrix for dornors, row is the cells, columns is the individuals
  # cu_n_lap: 
  # initial_param: initial estimates for parameters
  # null: (True or False) whether to use null model
  # cum_time: the cummulative time of all EM iterations
  '''

  # prepare some global data
  total_genes = int((ynxid_data.shape[0] - 1) / 2)
  cu_ghq_weights = cp.array([[0.1666667, 0.6666667, 0.1666667]]).astype(cp.float32).reshape(1,3)
  cu_ghq_nodes = cp.array([[-1.732051, -2.045201e-16, 1.732051]]).astype(cp.float32).reshape(1,3)

  # prepare the initial parameters and count data, convert them into GPU
  initial_param = initial_param[gene_index,:]
  if null == False: # under H1
    cu_param_result = cp.asarray(initial_param[:,[0,1,3]]).astype(cp.float32)
    cu_sigm2_result = cp.asarray(initial_param[:,[2]]).astype(cp.float32)
  else: # under H0
    cu_param_result = cp.asarray(initial_param[:,[4,7]]).astype(cp.float32)
    cu_sigm2_result = cp.asarray(initial_param[:,[6]]).astype(cp.float32)

  cu_randint_result = cp.zeros((len(gene_index),cu_id_data.shape[0]), dtype=cp.float32)
  cu_randint_prec_result = cp.ones_like(cu_randint_result,dtype = cp.float32) * 0

  cu_y_initial = (cp.asarray(ynxid_data[gene_index,:])).astype(cp.uint16)
  cu_n_initial = (cp.asarray(ynxid_data[ np.array(gene_index) + total_genes,:])).astype(cp.uint16)
  cu_id_data = (cu_id_data).astype(cp.uint8)

  if null == False:
    cu_x_initial = cp.zeros((1,cu_y_initial.shape[1],2))
    cu_x_initial[:,:,0] = cp.zeros((1,cu_y_initial.shape[1])) + 1
    cu_x_initial[:,:,1] = cp.asarray(ynxid_data[-1,:])
  else:
    cu_x_initial = cp.zeros((1,cu_y_initial.shape[1],1))
    cu_x_initial[:,:,0] = cp.zeros((1,cu_y_initial.shape[1])) + 1
  cu_x_initial = cu_x_initial.astype(cp.float32)

  gene_index = [i for i in range(len(gene_index))]
  
  # prepare a dict to record the log-likelihood and time
  llkl_record_dict = {}
  for ele in gene_index:
    llkl_record_dict[ele] = [0]
  used_time = []

  # Start EM algorithm
  for i in range(num_iteration):
    if len(gene_index) > 0:
      print()
      print("Iteration " + str(i) + " start: " + str(len(gene_index)) + " genes")
      start_em = time.time()

      # extract the data (Y, N, X, and parameters) for unfitted genes
      cu_y = cu_y_initial[gene_index,:]
      cu_n = cu_n_initial[gene_index,:]
      cu_x = cu_x_initial[:,:]

      cu_param = cu_param_result[gene_index,:] # b0, b1, phi
      cu_sigm2 = cu_sigm2_result[gene_index,:]
      cu_randint = cu_randint_result[gene_index,:]
      cu_randint_prec = cu_randint_result[gene_index,:]

      # Start one EM iteration for unfitted data
      cu_param, cu_sigm2, llkl, cu_randint, cu_randint_prec  = cu_vem_bb(cu_param, cu_sigm2,
                                        cu_x, cu_randint, cu_randint_prec, cu_id_data, cu_y, cu_n,
                                        cu_ghq_weights, cu_ghq_nodes, cu_n_lap, len(gene_index),i, null)

      # Record results (updated parameters and llkl)
      # select which genes have already fitted very well
      new_gene_index = []
      for gene in gene_index:
        (llkl_record_dict[gene])[0] += 1
        if llkl[gene_index.index(gene)] > (llkl_record_dict[gene])[-1] - 0.001 or i < 1:
          cu_param_result[gene_index,:] = cu_param
          cu_sigm2_result[gene_index,:] = cu_sigm2
          cu_randint_result[gene_index,:] = cu_randint
          cu_randint_prec_result[gene_index,:] = cu_randint_prec
          (llkl_record_dict[gene]).append(llkl[gene_index.index(gene)])
        if i < min_iter or len(llkl_record_dict[gene]) < 4 :
          new_gene_index.append(gene)
        elif abs(((llkl_record_dict[gene])[-1] - (llkl_record_dict[gene])[-2]) / (llkl_record_dict[gene])[-1]) > 1e-7: 
          # select which genes need more EM iterations
          # another rule: (llkl_record_dict[gene])[-1] - (llkl_record_dict[gene])[-2] > - 0.01)
          new_gene_index.append(gene)
      gene_index = new_gene_index
      end_em = time.time()
      cum_time += end_em - start_em
      print(f"Iteration {i}, cummulative time = {cum_time / 60} min")

  return cu_param_result, cu_sigm2_result, llkl_record_dict, cum_time

# Functions for initial estimates


def daesc_bb_gpu(ynidx_data,
                 part =1, num_iteration=50, min_iter=20):  
  ################################################################
  # make one-hot matrix
  labels = ynx_data_cuomo[-2, :].astype(int)
  num_classes = np.max(labels) + 1
  id_data = np.eye(num_classes)[labels]
  id_data = id_data.T
  cu_id_data = cp.asarray(id_data)

  # make initial values
  ynidx_data = ynidx_data.T
  num_gene = (ynidx_data.shape[1] - 2) // 2
  results = []
  start = time.time()
  print("Start to compute the initial parameters: ")
  for i in range(num_gene):
    if i % 500 == 0:
      print(f"Finish estimating {i} genes for initial estimates")
    total_count = ynidx_data[:, i + num_gene]
    index = np.nonzero(total_count)[0]
    test_y = (ynidx_data[index, i]).astype(np.float32)
    test_n = (ynidx_data[index, i + num_gene]).astype(np.float32)
    test_subj = (ynidx_data[index, 2 * num_gene]).astype(np.float32)
    test_x = ynidx_data[index, 2 * num_gene + 1]
    ones_column = np.ones((test_x.shape[0], 1))
    test_x = np.column_stack((ones_column, test_x)).astype(np.float32)

    # get initial estimates for phi
    phi_hat = optimize_multiple_starts(test_y, test_n, num_starts=1)

    # get initial estimtae for \beta_0 and \beta_1
    y_ratio = test_y / test_n
    X = sm.add_constant(test_x)
    glm_model = sm.GLM(y_ratio, X, family=sm.families.Binomial())
    glm_result = glm_model.fit()
    glm_coefficients = glm_result.params

    # fit initial parameters for null hypothsis
    X_intercept_only = np.ones((test_x.shape[0], 1))
    glm_model_intercept = sm.GLM(y_ratio, X_intercept_only, family=sm.families.Binomial())
    glm_result_intercept = glm_model_intercept.fit()
    new_intercept = glm_result_intercept.params[0]

    results.append({
        "Intercept": glm_coefficients[0], "Coef_x1": glm_coefficients[1], "Fixed_value": 0.05, "Estimated_Phi": phi_hat,          
        "New_Intercept": new_intercept, "Zero_Value": 0,"Repeated_Fixed_value": 0.05, "Repeated_Estimated_Phi": phi_hat})
  # record all results
  results_df = pd.DataFrame(results)
  initial_bb_param = results_df.to_numpy()
  end = time.time()
  print(f"Total time for initial estimates: {(end - start) / 60} minutes")

  #########################################################################
  # Begin to train the model
  gene_index = list(range(0, num_gene))
  parts = split_array(gene_index, 1)
  all_results = []
  cum_time = end - start
  ynidx_data = ynidx_data.T

  for i, part in enumerate(parts):
    print(f"{i}-th part begins, in total {len(part)} genes")
    print("First, under H1 model: ")
    my_param_result, my_sigm2_result, my_llkl_record_dict, cum_time = VEM_bb(
        part, num_iteration=num_iteration,min_iter=min_iter, ynxid_data=ynidx_data, cu_id_data=cu_id_data,
        cu_n_lap=2, initial_param=initial_bb_param, null=False, cum_time=cum_time)

    # record results
    my_param_result = my_param_result.get()
    my_sigm2_result = my_sigm2_result.get()
    llkl_np = np.array([[my_llkl_record_dict[i][-1]] for i in range(len(part))])

    clear_memory() # clear memory

    print("Second, under H0 model: ")
    my_param_null, my_sigm2_null, my_llkl_null_dict, cum_time = VEM_bb(
        part, num_iteration=num_iteration,min_iter=min_iter, ynxid_data=ynidx_data, cu_id_data=cu_id_data,
        cu_n_lap=2, initial_param=initial_bb_param, null=True, cum_time=cum_time)

    # record results under H0
    my_param_null = my_param_null.get()
    my_sigm2_null = my_sigm2_null.get()
    llkl_null_np = np.array([[my_llkl_null_dict[i][-1]] for i in range(len(part))])

    clear_memory()

    # combine all results and return results
    for idx in range(len(part)):
      combined_result = np.hstack((
            my_param_result[idx,0],my_param_result[idx,1],my_param_result[idx,2],my_sigm2_result[idx,0],llkl_np[idx],  # results under H1
            my_param_null[idx,0],my_param_null[idx,1],my_sigm2_null[idx,0],llkl_null_np[idx])) # results under H0
      all_results.append(combined_result)

  final_results = np.array(all_results)
  final_results = np.array(final_results)
  column_names = ["b0", "b1", "phi", "sigma2","llkl","b0_null","phi_null","sigma_null","llkl_null"]
  df = pd.DataFrame(final_results, columns=column_names)
  return df





