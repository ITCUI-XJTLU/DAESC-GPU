"""
daesc_mix_scipy.py
------------------
Simplified DAESC-Mix: replaces the custom batched CuPy BFGS optimizer
with scipy L-BFGS-B.  Processes one gene at a time (no GPU batching).

Key CPU <-> GPU transfer rules (the silent error source):
  • Every call to scipy.minimize receives a 1-D float64 *numpy* array.
  • Inside objective(), immediately move to GPU with
      cp.asarray(x_np).reshape(1, num_params).astype(cp.float32)
  • cu_q_mix runs entirely on GPU and returns CuPy arrays.
  • Before returning to scipy:  .get().flatten().astype(np.float64)
  • After scipy.minimize(), result.x is numpy → back to GPU with
      cp.asarray(result.x).reshape(1, num_params).astype(cp.float32)

Parameter layout of initial_mix_param (columns):
    0  : b0        (H1 intercept)
    1  : b1        (H1 covariate coef)
    2  : sigma2    (H1 initial)
    3  : phi       (H1 initial)
    4  : z0 = 0.9  (H1 mixture weight, phase-1)
    5  : z1 = 0.1  (H1 mixture weight, phase-2)
    6  : b0_null   (H0 intercept)
    7  : 0         (placeholder, not used)
    8  : sigma2_null
    9  : phi_null
    10 : z0_null = 0.9
    11 : z1_null = 0.1

Function names carry the _scipy suffix to avoid collisions with originals.
"""

import cupy as cp
import numpy as np
import pandas as pd
import statsmodels.api as sm
import time
from scipy.optimize import minimize

# ── Re-use all GPU kernels and helper utilities from the original module ───
from daesc_mix_gpu import (
    dot_product_mix,        # batched dot-product kernel
    cu_q_mix,               # Q-function + its gradient (Mix model)
    cu_logBB_bysubj,        # log-BB likelihood aggregated by individual (used in E-step)
    cu_compute_randint_deriv_mix,  # 1st & 2nd derivatives of rand-effect log-posterior
    cu_mix_lkl_agq,         # full marginal log-likelihood (Laplace approx)
    optimize_multiple_starts,
    clear_memory,
    split_array,
)


# ═══════════════════════════════════════════════════════════════════════════
#  Single-gene EM step
# ═══════════════════════════════════════════════════════════════════════════

def cu_vem_mix_scipy(
    cu_b_phi,           # [1, num_params]          float32 GPU  (3 params H1, 2 params H0)
    cu_sigm2,           # [1, 1]                   float32 GPU
    cu_p,               # [1, 2]                   float32 GPU  (mixture weights)
    cu_x,               # [1, num_cells, num_cov]  float32 GPU
    cu_randint,         # [1, num_individuals]      float32 GPU
    cu_randint_prec,    # [1, num_individuals]      float32 GPU
    cu_id_data,         # [num_individuals, num_cells]  uint8 GPU
    cu_y,               # [1, num_cells]            uint16 GPU
    cu_n,               # [1, num_cells]            uint16 GPU
    cu_ghq_weights,     # [1, 3]                    float32 GPU
    cu_ghq_nodes,       # [1, 3]                    float32 GPU
    cu_n_lap,           # kept for API parity (not used directly here)
    null,               # bool: True → H0 model
    max_optim = 10,     # max L-BFGS-B iterations per M-step call
):
    """
    One full EM iteration (E-step + M-step) for a *single* gene.
    The M-step optimises the Q-function using scipy L-BFGS-B.

    All heavy computation stays on GPU; only the parameter vector
    travels to/from CPU to satisfy the scipy.minimize API.

    Returns
    -------
    cu_b_phi        [1, num_params]  float32 GPU  (updated)
    cu_sigm2        [1, 1]           float32 GPU  (updated)
    cu_p            [1, 2]           float32 GPU  (updated mixture weights)
    llkl_np         numpy scalar array, shape [1]
    cu_randint      [1, num_individuals]  float32 GPU  (updated)
    cu_randint_prec [1, num_individuals]  float32 GPU  (updated)
    """
    cu_num_genes    = 1
    num_individuals = cu_id_data.shape[0]
    num_cells       = cu_n.shape[1]

    # ── Shared intermediate objects (stay on GPU throughout) ─────────────────

    # 3-point Gauss-Hermite quadrature abscissas in cell-broadcast form: [1, 1, 3]
    cu_zmat = cu_ghq_nodes.reshape(1, 1, -1).astype(cp.float32)

    # float32 cast of the one-hot donor matrix for matmul
    cu_id_f32 = cu_id_data.astype(cp.float32)

    # ── Linear predictor:  X @ b  →  [1, num_cells] ─────────────────────────
    cu_b       = cu_b_phi[:, :-1].reshape(1, 1, -1)       # [1, 1, num_cov]
    cu_phi     = cu_b_phi[:, -1].reshape(-1, 1)            # [1, 1]
    cu_linpred = dot_product_mix(cu_x, cu_b, axis=2)       # [1, num_cells]

    # ── Expand random effects to cell level ───────────────────────────────────
    # cu_randint        : [1, I]  →  [1, C, 1]  (MAP estimate per cell)
    # cu_randint_prec   : [1, I]  →  [1, C, 1]  (posterior precision, used
    #                                             as GHQ σ in the quadrature)
    def _expand(rand, prec):
        rv = (cp.matmul(rand,              cu_id_f32)
              .reshape(1, num_cells, 1)
              .astype(cp.float32))
        pv = (cp.matmul(cp.sqrt(prec),     cu_id_f32)
              .reshape(1, num_cells, 1)
              .astype(cp.float32))
        return rv, pv

    cu_rv, cu_pv = _expand(cu_randint, cu_randint_prec)

    # ── E-step: 2 Newton updates for the random-intercept MAP estimates ───────
    # Mirrors the original cu_vem_mix loop exactly, but for a single gene.
    for _ in range(2):
        # ① Mixture weights  cu_wt : [1, I, 2]
        #    log p(phase=0 | data, params) ∝ log z_0 + log-BB(phase-0 linpred)
        #    log p(phase=1 | data, params) ∝ log z_1 + log-BB(phase-1 linpred)
        #    cu_logBB_bysubj aggregates per-cell log-BB values to per-individual.
        #
        # NOTE: cu_logBB_bysubj has a legacy `cu_sigm2` argument that is not
        # used inside the kernel; we pass cu_sigm2 for API consistency.
        log_w0 = cp.log(cu_p[:, [0]]) + cu_logBB_bysubj(
            cu_linpred,  cu_sigm2, cu_phi,
            cu_y, cu_n, cu_id_data,
            cu_rv, cu_pv, cu_zmat, cu_ghq_weights,
        )  # [1, I]
        log_w1 = cp.log(cu_p[:, [1]]) + cu_logBB_bysubj(
            -cu_linpred, cu_sigm2, cu_phi,
            cu_y, cu_n, cu_id_data,
            cu_rv, cu_pv, cu_zmat, cu_ghq_weights,
        )  # [1, I]

        cu_wt = cp.stack((log_w0, log_w1), axis=-1)                 # [1, I, 2]
        cu_wt = cp.exp(cu_wt - cp.max(cu_wt, axis=2, keepdims=True))
        cu_wt = cu_wt / cp.sum(cu_wt, axis=2, keepdims=True)        # normalise

        # ② 1st and 2nd derivatives of the log-posterior wrt random intercept u
        #    cu_dfX has shape [2, I]:
        #      rows [0:1] → ∂ log p / ∂ u   (gradient)
        #      rows [1:2] → ∂² log p / ∂ u²  (negative Hessian)
        cu_df0 = cu_compute_randint_deriv_mix(
            cu_linpred,  cu_randint, cu_id_data, cu_y, cu_n, cu_phi)
        cu_df1 = cu_compute_randint_deriv_mix(
            -cu_linpred, cu_randint, cu_id_data, cu_y, cu_n, cu_phi)

        # Mixture-weighted gradient and Hessian of the full conditional
        wt0, wt1 = cu_wt[:, :, 0], cu_wt[:, :, 1]         # each [1, I]
        mixed_grad = (wt0 * cu_df0[:cu_num_genes, :]
                      + wt1 * cu_df1[:cu_num_genes, :]
                      - cu_randint / cu_sigm2)              # [1, I]
        mixed_hess = (wt0 * cu_df0[cu_num_genes:, :]
                      + wt1 * cu_df1[cu_num_genes:, :]
                      - 1.0 / cu_sigm2)                    # [1, I]  (negative)

        # Newton step (step size 0.9 for stability, same as original)
        cu_randint = cu_randint - 0.9 * mixed_grad / mixed_hess

        # ③ Recompute derivatives at updated u to get the posterior precision
        cu_df0 = cu_compute_randint_deriv_mix(
            cu_linpred,  cu_randint, cu_id_data, cu_y, cu_n, cu_phi)
        cu_df1 = cu_compute_randint_deriv_mix(
            -cu_linpred, cu_randint, cu_id_data, cu_y, cu_n, cu_phi)
        cu_randint_prec = cp.abs(
            wt0 * cu_df0[cu_num_genes:, :]
            + wt1 * cu_df1[cu_num_genes:, :]
            - 1.0 / cu_sigm2
        )                                                   # [1, I]  (positive)

        # ④ Refresh cell-level expansions with updated u and prec
        cu_rv, cu_pv = _expand(cu_randint, cu_randint_prec)

    # cu_wt, cu_rv, cu_pv now hold the final E-step values
    del cu_linpred   # free GPU memory; will be recomputed inside cu_q_mix

    # ── M-step 1: update σ² (Laplace approximation) ──────────────────────────
    # E[u²] ≈ û² + 1/precision
    cu_sig_test    = cu_randint * cu_randint + cp.reciprocal(cu_randint_prec)
    rand_indicator = (cu_randint != 0.0).astype(cp.float32)        # [1, I]
    denom_sig      = cp.sum(rand_indicator, axis=1)
    denom_sig      = cp.where(denom_sig == 0, 1.0, denom_sig)      # guard /0
    cu_sigm2       = (cp.sum(cu_sig_test * rand_indicator, axis=1)
                      / denom_sig).reshape(-1, 1)                   # [1, 1]

    # ── M-step 2: update mixture proportions p ───────────────────────────────
    # Weighted average of cu_wt over individuals (ignore zero-count donors)
    rand_ind_3d = rand_indicator.reshape(1, num_individuals, 1)    # [1, I, 1]
    cu_p = (cp.sum(cu_wt * rand_ind_3d, axis=1)
            / cp.sum(rand_ind_3d, axis=1))                          # [1, 2]
    cu_p = cp.maximum(cu_p, 0.001)
    cu_p = cu_p / cp.sum(cu_p, axis=1, keepdims=True)

    # ── M-step 3: optimise Q(b, φ) via scipy L-BFGS-B ───────────────────────
    num_params = cu_b_phi.shape[1]

    # Parameter bounds: β are unconstrained; φ must be positive
    bounds = [(-100.0, 100.0)] * (num_params - 1) + [(1e-4, 100.0)]

    # Pull current parameters from GPU → CPU as 1-D float64 (scipy requirement)
    #
    # CRITICAL CPU ← GPU transfer #1
    x0_np = cp.asnumpy(cu_b_phi).flatten().astype(np.float64)

    # Cast data arrays to required dtypes (do this once, outside objective)
    cu_y_u16  = cu_y.astype(cp.uint16)
    cu_n_u16  = cu_n.astype(cp.uint16)
    cu_wt_f32 = cu_wt.astype(cp.float32)
    cu_id_u8  = cu_id_data.astype(cp.uint8)

    # All GPU tensors are captured as closure constants.
    # They must NOT be reassigned inside objective() — scipy calls objective()
    # many times with different x_np but the same GPU data.
    _cu_wt  = cu_wt_f32          # [1, I, 2]  mixture weights (E-step output)
    _cu_x   = cu_x               # [1, C, cov]
    _cu_y   = cu_y_u16           # [1, C]
    _cu_n   = cu_n_u16           # [1, C]
    _cu_z   = cu_zmat            # [1, 1, 3]  GHQ abscissas
    _cu_ghw = cu_ghq_weights     # [1, 3]     GHQ weights
    _cu_rv  = cu_rv              # [1, C, 1]  cell-level random effects
    _cu_pv  = cu_pv              # [1, C, 1]  cell-level GHQ σ
    _cu_id  = cu_id_u8           # [I, C]

    def objective(x_np):
        """
        Thin wrapper: CPU → GPU → compute Q-function → CPU.

        Parameters
        ----------
        x_np  : 1-D float64 numpy array  (parameter vector from scipy)

        Returns
        -------
        (val_cpu, grad_cpu) : (float64 scalar, 1-D float64 numpy array)
        """
        # ① CPU → GPU
        #    Reshape to [1, num_params] and cast to float32 (kernel precision)
        x_cp = cp.asarray(x_np, dtype=cp.float32).reshape(1, num_params)

        # ② GPU computation: Q-function value + gradient
        #    cu_q_mix returns (val [1], grad [1, num_params]), both CuPy
        val_cp, grad_cp = cu_q_mix(
            x_cp,
            cu_wt_q           = _cu_wt,
            cu_x_q            = _cu_x,
            cu_y_q            = _cu_y,
            cu_n_q            = _cu_n,
            cu_zmat_q         = _cu_z,
            cu_ghq_weights_q  = _cu_ghw,
            cu_randint_q      = _cu_rv,
            cu_randint_prec_q = _cu_pv,
            cu_id_data_q      = _cu_id,
            cu_num_genes_q    = 1,
            cu_iteration      = 0,      # iteration arg is unused in the kernel
        )

        # Ensure all GPU ops are finished before pulling data back
        cp.cuda.Stream.null.synchronize()

        # ③ GPU → CPU
        #    CRITICAL: .get() copies device memory to host; astype ensures float64
        val_cpu  = float(val_cp.get()[0])
        grad_cpu = grad_cp.get().flatten().astype(np.float64)

        return val_cpu, grad_cpu

    result = minimize(
        objective,
        x0_np,
        method  = 'L-BFGS-B',
        jac     = True,       # objective returns (f, g) together
        bounds  = bounds,
        options = {
            'maxiter' : max_optim,
            'ftol'    : 1e-12,
            'gtol'    : 1e-7,
            'maxfun'  : max_optim * 20,
        },
    )

    # ④ CPU → GPU: push optimised parameter vector back to device
    #    CRITICAL: must be float32 to match kernel precision
    cu_b_phi = (cp.asarray(result.x, dtype=cp.float32)
                .reshape(1, num_params))

    # ── Compute marginal log-likelihood at updated parameters ─────────────────
    # cu_mix_lkl_agq internally uses a Laplace approximation
    # (3 Newton steps from u=0) — same approach as the original.
    llkl    = cu_mix_lkl_agq(
        cu_b_phi, cu_sigm2, cu_p, cu_x,
        cu_randint, cu_id_u8, cu_y_u16, cu_n_u16, cu_num_genes,
    )
    llkl_np = cp.asnumpy(llkl)    # shape [1], numpy float32

    return cu_b_phi, cu_sigm2, cu_p, llkl_np, cu_randint, cu_randint_prec


# ═══════════════════════════════════════════════════════════════════════════
#  EM outer loop (iterates over genes one-by-one)
# ═══════════════════════════════════════════════════════════════════════════

def VEM_mix_scipy(
    gene_index,       # list of gene indices into ynxid_data rows
    num_iteration,    # max EM iterations
    min_iter,         # min EM iterations before convergence check
    ynxid_data,       # (2G+2) × C  numpy array
    cu_id_data,       # [num_individuals, num_cells]  CuPy
    cu_n_lap,         # kept for API parity
    initial_param,    # G × 12  numpy  (initial parameter table)
    null,             # bool: True → H0 model
    cum_time,         # cumulative elapsed seconds (pass-through)
    max_optim,        # max L-BFGS-B inner iterations per M-step
):
    """
    EM loop for DAESC-Mix (scipy version).
    Each gene is processed independently; no GPU batching.

    Returns
    -------
    cu_param_result   : [G, num_params]  CuPy
    cu_sigm2_result   : [G, 1]           CuPy
    cu_p_result       : [G, 2]           CuPy
    llkl_record_dict  : dict[gene_local_idx → [iter_count, llkl_1, llkl_2, ...]]
    cum_time          : updated float (seconds)
    """
    total_genes = int((ynxid_data.shape[0] - 1) / 2)

    # 3-point Gauss-Hermite quadrature constants (shared across all genes)
    cu_ghq_weights = cp.array([[0.1666667, 0.6666667, 0.1666667]],
                               dtype=cp.float32)    # [1, 3]
    cu_ghq_nodes   = cp.array([[-1.732051, -2.045201e-16, 1.732051]],
                               dtype=cp.float32)    # [1, 3]

    # ── Select initial parameter columns ──────────────────────────────────────
    # See module docstring for full column layout.
    subset = initial_param[gene_index, :]
    if not null:
        # H1: (b0, b1, phi), sigma2, (z0, z1)
        cu_param_result = cp.asarray(subset[:, [0, 1, 3]], dtype=cp.float32)  # [G, 3]
        cu_sigm2_result = cp.asarray(subset[:, [2]],       dtype=cp.float32)  # [G, 1]
        cu_p_result     = cp.asarray(subset[:, [4, 5]],    dtype=cp.float32)  # [G, 2]
    else:
        # H0: (b0_null, phi_null), sigma2_null, (z0_null, z1_null)
        cu_param_result = cp.asarray(subset[:, [6, 9]],    dtype=cp.float32)  # [G, 2]
        cu_sigm2_result = cp.asarray(subset[:, [8]],       dtype=cp.float32)  # [G, 1]
        cu_p_result     = cp.asarray(subset[:, [10, 11]],  dtype=cp.float32)  # [G, 2]

    num_genes_total = len(gene_index)
    num_individuals = cu_id_data.shape[0]

    # Random intercepts start at 0; precisions initialised to 1/σ² (same as original)
    cu_randint_result      = cp.zeros(
        (num_genes_total, num_individuals), dtype=cp.float32)
    cu_randint_prec_result = (
        cp.ones((num_genes_total, num_individuals), dtype=cp.float32)
        / cu_sigm2_result                          # broadcast [G,1] → [G, I]
    )

    # Count arrays: GPU, uint16
    cu_y_initial = cp.asarray(
        ynxid_data[gene_index, :], dtype=cp.uint16)                       # [G, C]
    cu_n_initial = cp.asarray(
        ynxid_data[np.array(gene_index) + total_genes, :], dtype=cp.uint16)  # [G, C]
    cu_id_data   = cu_id_data.astype(cp.uint8)                            # [I, C]

    # Covariate matrix: shape [1, num_cells, num_cov] (shared across genes)
    num_cells = cu_y_initial.shape[1]
    if not null:
        cu_x_initial          = cp.zeros((1, num_cells, 2), dtype=cp.float32)
        cu_x_initial[:, :, 0] = 1.0                                   # intercept
        cu_x_initial[:, :, 1] = cp.asarray(ynxid_data[-1, :])         # covariate X
    else:
        cu_x_initial          = cp.zeros((1, num_cells, 1), dtype=cp.float32)
        cu_x_initial[:, :, 0] = 1.0                                   # intercept only

    # Local (0-based) gene indices for this call
    local_index = list(range(num_genes_total))

    # llkl_record_dict[g] = [iteration_counter, llkl_at_iter_1, llkl_at_iter_2, ...]
    llkl_record_dict = {g: [0] for g in local_index}

    # ── EM iterations ─────────────────────────────────────────────────────────
    for i in range(num_iteration):
        if not local_index:
            break
        print(f"\nIteration {i}  |  active genes: {len(local_index)}")
        t0 = time.time()
        still_active = []

        for g in local_index:
            # Slice single-gene tensors (first dim = 1)
            cu_y_g            = cu_y_initial[g:g+1, :]            # [1, C]
            cu_n_g            = cu_n_initial[g:g+1, :]            # [1, C]
            cu_b_phi_g        = cu_param_result[g:g+1, :]         # [1, P]
            cu_sigm2_g        = cu_sigm2_result[g:g+1, :]         # [1, 1]
            cu_p_g            = cu_p_result[g:g+1, :]             # [1, 2]
            cu_randint_g      = cu_randint_result[g:g+1, :]       # [1, I]
            cu_randint_prec_g = cu_randint_prec_result[g:g+1, :]  # [1, I]

            # One full EM step (E-step + M-step) for this gene
            (cu_b_phi_g,
             cu_sigm2_g,
             cu_p_g,
             llkl_g,
             cu_randint_g,
             cu_randint_prec_g) = cu_vem_mix_scipy(
                cu_b_phi_g, cu_sigm2_g, cu_p_g,
                cu_x_initial,                       # shared covariate matrix
                cu_randint_g, cu_randint_prec_g,
                cu_id_data,
                cu_y_g, cu_n_g,
                cu_ghq_weights, cu_ghq_nodes,
                cu_n_lap, null, max_optim,
            )

            # Write updated values back into the result arrays
            cu_param_result[g:g+1, :]       = cu_b_phi_g
            cu_sigm2_result[g:g+1, :]       = cu_sigm2_g
            cu_p_result[g:g+1, :]           = cu_p_g
            cu_randint_result[g:g+1, :]     = cu_randint_g
            cu_randint_prec_result[g:g+1, :] = cu_randint_prec_g

            # Log-likelihood bookkeeping
            llkl_val = float(llkl_g[0])
            llkl_record_dict[g][0] += 1
            llkl_record_dict[g].append(llkl_val)

            # Convergence check (same criterion as original VEM_mix)
            history = llkl_record_dict[g]
            if i < min_iter or len(history) < 4:
                still_active.append(g)
            elif (abs((history[-1] - history[-2]) / (history[-1] + 1e-12)) > 1e-7
                  and history[-1] - history[-2] > -0.001):
                still_active.append(g)

        local_index = still_active
        elapsed = time.time() - t0
        cum_time += elapsed
        print(f"  cumulative time = {cum_time / 60:.2f} min")

    return cu_param_result, cu_sigm2_result, cu_p_result, llkl_record_dict, cum_time


# ═══════════════════════════════════════════════════════════════════════════
#  Top-level entry point  (mirrors daesc_mix_gpu.daesc_mix_gpu API)
# ═══════════════════════════════════════════════════════════════════════════

def daesc_mix_scipy(
    ynidx_data,           # (2G+2) × C  numpy array  (same format as daesc_mix_gpu)
    num_iteration = 50,
    min_iter      = 20,
    max_optim     = 10,   # max L-BFGS-B inner iterations per M-step
):
    """
    DAESC-Mix with scipy L-BFGS-B optimiser (single-gene processing).

    Parameters
    ----------
    ynidx_data    : numpy array, shape (2G + 2) × C
                    Rows 0..G-1       : alternative allele counts  Y
                    Rows G..2G-1      : total counts               N
                    Row  2G           : donor/individual IDs        (int)
                    Row  2G+1         : covariate                   X
    num_iteration : max EM iterations
    min_iter      : min EM iterations before early stopping
    max_optim     : max L-BFGS-B iterations per M-step

    Returns
    -------
    pandas.DataFrame  with columns:
        b0, b1, phi, sigma2, z0, z1, llkl,
        b0_null, phi_null, sigma_null, z0_null, z1_null, llkl_null
    """
    # ── One-hot donor matrix ──────────────────────────────────────────────────
    labels      = ynidx_data[-2, :].astype(int)
    num_classes = int(np.max(labels)) + 1
    id_data     = np.eye(num_classes)[labels].T            # [I, C]
    cu_id_data  = cp.asarray(id_data)

    # ── Initial parameter estimates ───────────────────────────────────────────
    ynidx_T  = ynidx_data.T                                # [C, 2G+2]
    num_gene = (ynidx_T.shape[1] - 2) // 2
    results  = []
    t_start  = time.time()
    print("Computing initial parameters …")

    for i in range(num_gene):
        if i % 500 == 0:
            print(f"  {i}/{num_gene} genes")

        total_count = ynidx_T[:, i + num_gene]
        idx         = np.nonzero(total_count)[0]

        test_y = ynidx_T[idx, i].astype(np.float32)
        test_n = ynidx_T[idx, i + num_gene].astype(np.float32)
        test_x = ynidx_T[idx, 2 * num_gene + 1]
        ones   = np.ones((len(idx), 1))
        test_x = np.column_stack([ones, test_x]).astype(np.float32)

        # Phi estimate from beta-binomial marginal
        phi_hat = optimize_multiple_starts(test_y, test_n, num_starts=1)

        y_ratio = test_y / np.maximum(test_n, 1e-9)

        # H1 GLM (intercept + covariate)
        X_h1    = sm.add_constant(test_x)
        glm_h1  = sm.GLM(y_ratio, X_h1, family=sm.families.Binomial()).fit()
        coef_h1 = glm_h1.params

        # H0 GLM (intercept only)
        X_h0    = np.ones((len(idx), 1))
        glm_h0  = sm.GLM(y_ratio, X_h0, family=sm.families.Binomial()).fit()
        coef_h0 = glm_h0.params[0]

        results.append({
            "Intercept"           : coef_h1[0],   # col 0: b0 H1
            "Coef_x1"             : coef_h1[1],   # col 1: b1 H1
            "Fixed_value"         : 0.05,          # col 2: sigma2 H1
            "Estimated_Phi"       : phi_hat,       # col 3: phi H1
            # cols 4, 5 inserted below (z0=0.9, z1=0.1 for H1)
            "New_Intercept"       : coef_h0,       # col 6: b0 H0
            "Zero_Value"          : 0.0,           # col 7: placeholder
            "Repeated_Fixed_value": 0.05,          # col 8: sigma2 H0
            "Repeated_Phi"        : phi_hat,       # col 9: phi H0
            # cols 10, 11 inserted below (z0=0.9, z1=0.1 for H0)
        })

    results_df = pd.DataFrame(results)
    results_df.insert(4,  "H1_9", 0.9)     # col 4:  z0 H1
    results_df.insert(5,  "H1_1", 0.1)     # col 5:  z1 H1
    results_df.insert(10, "H0_9", 0.9)     # col 10: z0 H0
    results_df.insert(11, "H_1",  0.1)     # col 11: z1 H0
    initial_mix_param = results_df.to_numpy()   # shape [G, 12]

    t_end    = time.time()
    cum_time = t_end - t_start
    print(f"Initial estimates done in {cum_time / 60:.2f} min")

    gene_index = list(range(num_gene))

    # ── H1 model ──────────────────────────────────────────────────────────────
    print("\n── H1 model ──")
    param_h1, sigm2_h1, p_h1, llkl_h1_dict, cum_time = VEM_mix_scipy(
        gene_index,
        num_iteration = num_iteration,
        min_iter      = min_iter,
        ynxid_data    = ynidx_data,
        cu_id_data    = cu_id_data,
        cu_n_lap      = 2,
        initial_param = initial_mix_param,
        null          = False,
        cum_time      = cum_time,
        max_optim     = max_optim,
    )
    param_h1_np = param_h1.get()                             # [G, 3]
    sigm2_h1_np = sigm2_h1.get()                             # [G, 1]
    p_h1_np     = p_h1.get()                                 # [G, 2]
    llkl_h1_np  = np.array([llkl_h1_dict[i][-1] for i in range(num_gene)])

    clear_memory()

    # ── H0 model ──────────────────────────────────────────────────────────────
    print("\n── H0 model ──")
    param_h0, sigm2_h0, p_h0, llkl_h0_dict, cum_time = VEM_mix_scipy(
        gene_index,
        num_iteration = num_iteration,
        min_iter      = min_iter,
        ynxid_data    = ynidx_data,
        cu_id_data    = cu_id_data,
        cu_n_lap      = 2,
        initial_param = initial_mix_param,
        null          = True,
        cum_time      = cum_time,
        max_optim     = max_optim,
    )
    param_h0_np = param_h0.get()                             # [G, 2]
    sigm2_h0_np = sigm2_h0.get()                             # [G, 1]
    p_h0_np     = p_h0.get()                                 # [G, 2]
    llkl_h0_np  = np.array([llkl_h0_dict[i][-1] for i in range(num_gene)])

    clear_memory()

    # ── Assemble output DataFrame ──────────────────────────────────────────────
    rows = []
    for idx in range(num_gene):
        rows.append([
            # H1 results
            param_h1_np[idx, 0],    # b0
            param_h1_np[idx, 1],    # b1
            param_h1_np[idx, 2],    # phi
            sigm2_h1_np[idx, 0],    # sigma2
            p_h1_np[idx, 0],        # z0
            p_h1_np[idx, 1],        # z1
            llkl_h1_np[idx],        # llkl
            # H0 results
            param_h0_np[idx, 0],    # b0_null
            param_h0_np[idx, 1],    # phi_null
            sigm2_h0_np[idx, 0],    # sigma_null
            p_h0_np[idx, 0],        # z0_null
            p_h0_np[idx, 1],        # z1_null
            llkl_h0_np[idx],        # llkl_null
        ])

    column_names = [
        "b0", "b1", "phi", "sigma2", "z0", "z1", "llkl",
        "b0_null", "phi_null", "sigma_null", "z0_null", "z1_null", "llkl_null",
    ]
    return pd.DataFrame(rows, columns=column_names)









