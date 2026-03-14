"""
daesc_bb_scipy.py
-----------------
Simplified DAESC-BB: replaces the custom batched CuPy BFGS optimizer
with scipy L-BFGS-B.  Processes one gene at a time (no batching).

Key CPU <-> GPU transfer rules (the silent error source):
  • Every call to scipy.minimize receives a 1-D float64 *numpy* array.
  • Inside objective(), we immediately move that array to GPU with
    cp.asarray(x_np).reshape(1, num_params).astype(cp.float32).
  • The Q function (cu_q_bb) runs entirely on GPU and returns CuPy arrays.
  • Before returning to scipy we call .get().flatten().astype(np.float64)
    to pull the scalar loss and the gradient back to CPU.
  • After scipy.minimize(), result.x is numpy → convert back to GPU with
    cp.asarray(result.x).reshape(1, num_params).astype(cp.float32).

Function names use the suffix  _scipy  to avoid collisions with originals.
"""

import cupy as cp
import numpy as np
import pandas as pd
import statsmodels.api as sm
import time
from scipy.optimize import minimize

# ── re-use *all* GPU kernels and helper utilities from the original module ──
from daesc_bb_gpu import (
    dot_product_bb,
    cu_q_bb,
    cu_bb_lkl_agq_fixvar,
    cu_compute_randint_deriv_bb,
    optimize_multiple_starts,
    clear_memory,
    split_array,
)


# ═══════════════════════════════════════════════════════════════════════════
#  Core single-gene EM step
# ═══════════════════════════════════════════════════════════════════════════

def cu_vem_bb_scipy(
    cu_b_phi,          # [1, num_params]   float32 GPU  (num_params=3 H1, 2 H0)
    cu_sigm2,          # [1, 1]            float32 GPU
    cu_x,              # [1, num_cells, num_cov]  float32 GPU  (shared across genes)
    cu_randint,        # [1, num_individuals]     float32 GPU
    cu_randint_prec,   # [1, num_individuals]     float32 GPU
    cu_id_data,        # [num_individuals, num_cells]  uint8 GPU
    cu_y,              # [1, num_cells]    uint16 GPU
    cu_n,              # [1, num_cells]    uint16 GPU
    cu_ghq_weights,    # [1, 3]            float32 GPU
    cu_ghq_nodes,      # [1, 3]            float32 GPU
    cu_n_lap,          # not used here (kept for API parity)
    null,              # bool (kept for API parity)
    max_optim=10,      # max L-BFGS-B iterations
):
    """
    One full EM iteration (E-step + M-step) for a *single* gene using
    scipy L-BFGS-B for the M-step Q-function optimisation.

    Returns
    -------
    cu_b_phi          updated [1, num_params]
    cu_sigm2          updated [1, 1]
    llkl_np           numpy scalar array shape [1]  (log-likelihood)
    cu_randint        updated [1, num_individuals]
    cu_randint_prec   updated [1, num_individuals]
    """
    cu_num_genes = 1  # this function is intentionally single-gene

    # ── E-step ──────────────────────────────────────────────────────────────
    # Linear predictor:  X @ b  →  [1, num_cells]
    cu_b      = cu_b_phi[:, :-1].reshape(1, 1, -1)       # [1, 1, num_cov]
    cu_phi    = cu_b_phi[:, -1].reshape(-1, 1)            # [1, 1]
    cu_linpred = dot_product_bb(cu_x, cu_b, axis=2)       # [1, num_cells]

    # 2 Newton updates for random intercepts
    # cu_compute_randint_deriv_bb returns [2, num_individuals]:
    #   rows [0:1]  →  first  derivative  (∂ log-lik / ∂ u)
    #   rows [1:2]  →  second derivative  (∂² log-lik / ∂ u²)
    for _ in range(2):
        cu_df = cu_compute_randint_deriv_bb(
            cu_linpred, cu_randint, cu_id_data, cu_y, cu_n, cu_phi
        )
        numerator   = cu_df[:cu_num_genes, :] - cu_randint / cu_sigm2   # [1, num_indiv]
        denominator = cu_df[cu_num_genes:,  :] - 1.0 / cu_sigm2         # [1, num_indiv]
        cu_randint  = cu_randint - 0.9 * numerator / denominator

    # Final precision = |second deriv of full conditional|
    cu_df = cu_compute_randint_deriv_bb(
        cu_linpred, cu_randint, cu_id_data, cu_y, cu_n, cu_phi
    )
    cu_randint_prec = cp.abs(cu_df[cu_num_genes:, :] - 1.0 / cu_sigm2)  # [1, num_indiv]
    del cu_linpred

    # ── M-step: update σ² ───────────────────────────────────────────────────
    # Laplace approximation:  E[u²] ≈ û² + 1/precision
    cu_zmat        = cu_ghq_nodes.reshape(1, 1, -1).astype(cp.float32)   # [1, 1, 3]
    cu_sigm2_test  = cu_randint * cu_randint + cp.reciprocal(cu_randint_prec)
    rand_indicator = (cu_randint != 0.0).astype(cp.float32)              # [1, num_indiv]
    denom_sig      = cp.sum(rand_indicator, axis=1)
    denom_sig      = cp.where(denom_sig == 0, 1.0, denom_sig)            # guard /0
    cu_sigm2       = (
        cp.sum(cu_sigm2_test * rand_indicator, axis=1) / denom_sig
    ).reshape(-1, 1)                                                      # [1, 1]

    # Expand random effects from individual-level → cell-level
    # [1, num_individuals] @ [num_individuals, num_cells]  →  [1, num_cells]
    # then reshape to [1, num_cells, 1] for use inside cu_q_bb
    cu_id_f32           = cu_id_data.astype(cp.float32)
    cu_randint_vector   = (
        cp.matmul(cu_randint, cu_id_f32)
        .reshape(1, cu_n.shape[1], 1)
        .astype(cp.float32)
    )                                                                     # [1, num_cells, 1]
    cu_precision_vector = (
        cp.matmul(cp.sqrt(cu_randint_prec), cu_id_f32)
        .reshape(1, cu_n.shape[1], 1)
        .astype(cp.float32)
    )                                                                     # [1, num_cells, 1]

    # ── M-step: optimise Q(b, φ) with scipy L-BFGS-B ────────────────────────
    num_params = cu_b_phi.shape[1]

    # Bounds: β parameters are unconstrained; φ must be positive
    bounds = [(-100.0, 100.0)] * (num_params - 1) + [(1e-4, 100.0)]

    # Initial point: pull current params from GPU → CPU (1-D float64)
    x0_np = cp.asnumpy(cu_b_phi).flatten().astype(np.float64)

    # Freeze all GPU tensors as closure variables so the objective is a
    # pure function of the 1-D CPU parameter vector only
    _cu_x   = cu_x                  # [1, num_cells, num_cov]
    _cu_y   = cu_y.astype(cp.uint16)
    _cu_n   = cu_n.astype(cp.uint16)
    _cu_z   = cu_zmat               # [1, 1, 3]
    _cu_ghw = cu_ghq_weights        # [1, 3]
    _cu_rv  = cu_randint_vector     # [1, num_cells, 1]
    _cu_pv  = cu_precision_vector   # [1, num_cells, 1]

    def objective(x_np):
        """
        CPU → GPU → CPU wrapper around cu_q_bb.

        x_np   : 1-D float64 numpy  (from scipy)
        returns: (scalar float64, 1-D float64 numpy)  (value, gradient)
        """
        # ① CPU → GPU: reshape to [1, num_params] as cu_q_bb expects
        x_cp = cp.asarray(x_np, dtype=cp.float32).reshape(1, num_params)

        # ② GPU computation (all intermediate tensors stay on GPU)
        val, grad = cu_q_bb(
            x_cp,
            cu_x_q            = _cu_x,
            cu_y_q            = _cu_y,
            cu_n_q            = _cu_n,
            cu_zmat_q         = _cu_z,
            cu_ghq_weights_q  = _cu_ghw,
            cu_randint_q      = _cu_rv,
            cu_randint_prec_q = _cu_pv,
            cu_num_genes_q    = 1,
        )
        # val  : [1]            CuPy float32
        # grad : [1, num_params] CuPy float32

        cp.cuda.Stream.null.synchronize()   # ensure GPU work is done

        # ③ GPU → CPU
        val_cpu  = float(val.get()[0])                          # Python float
        grad_cpu = grad.get().flatten().astype(np.float64)      # 1-D numpy

        return val_cpu, grad_cpu

    result = minimize(
        objective,
        x0_np,
        method  = 'L-BFGS-B',
        jac     = True,          # objective returns (f, g) jointly
        bounds  = bounds,
        options = {
            'maxiter' : max_optim,
            'ftol'    : 1e-12,
            'gtol'    : 1e-7,
            'maxfun'  : max_optim * 20,
        },
    )

    # ④ CPU → GPU: move optimised parameters back
    cu_b_phi = cp.asarray(result.x, dtype=cp.float32).reshape(1, num_params)

    # ── Compute log-likelihood at updated parameters ─────────────────────────
    llkl    = cu_bb_lkl_agq_fixvar(
        cu_b_phi, cu_sigm2, cu_x, cu_randint, cu_id_data, cu_y, cu_n, cu_num_genes
    )
    llkl_np = cp.asnumpy(llkl)   # [1] numpy

    return cu_b_phi, cu_sigm2, llkl_np, cu_randint, cu_randint_prec


# ═══════════════════════════════════════════════════════════════════════════
#  EM outer loop (iterates over genes one-by-one)
# ═══════════════════════════════════════════════════════════════════════════

def VEM_bb_scipy(
    gene_index,       # list of integer indices (into ynxid_data rows)
    num_iteration,    # max EM iterations
    min_iter,         # min EM iterations before convergence check
    ynxid_data,       # (2*G+2) × C  numpy  (genes × cells layout)
    cu_id_data,       # [num_individuals, num_cells]  CuPy
    cu_n_lap,         # kept for API parity
    initial_param,    # G × 8  numpy  (initial parameter table)
    null,             # bool: True → H0 model, False → H1 model
    cum_time,         # cumulative elapsed seconds (pass-through)
    max_optim,        # max L-BFGS-B inner iterations per M-step
):
    """
    EM loop for DAESC-BB (scipy version).
    Processes each gene independently; no GPU batching.
    """
    total_genes = int((ynxid_data.shape[0] - 1) / 2)

    # Gauss-Hermite quadrature nodes & weights (3-point)
    cu_ghq_weights = cp.array([[0.1666667, 0.6666667, 0.1666667]],
                               dtype=cp.float32)                  # [1, 3]
    cu_ghq_nodes   = cp.array([[-1.732051, -2.045201e-16, 1.732051]],
                               dtype=cp.float32)                  # [1, 3]

    # ── Select columns of initial_param ─────────────────────────────────────
    # initial_param columns (see daesc_bb_gpu for full description):
    #   0: b0,  1: b1,  2: sigma2,  3: phi  (H1)
    #   4: b0_null,  6: sigma2_null,  7: phi_null  (H0)
    subset = initial_param[gene_index, :]
    if not null:
        cu_param_result = cp.asarray(subset[:, [0, 1, 3]], dtype=cp.float32)  # [G, 3]
        cu_sigm2_result = cp.asarray(subset[:, [2]],       dtype=cp.float32)  # [G, 1]
    else:
        cu_param_result = cp.asarray(subset[:, [4, 7]], dtype=cp.float32)     # [G, 2]
        cu_sigm2_result = cp.asarray(subset[:, [6]],    dtype=cp.float32)     # [G, 1]

    num_genes_total = len(gene_index)
    num_individuals = cu_id_data.shape[0]

    # Random intercepts and precisions: start at zero
    cu_randint_result      = cp.zeros((num_genes_total, num_individuals),
                                       dtype=cp.float32)
    cu_randint_prec_result = cp.zeros_like(cu_randint_result)

    # Count data → GPU
    cu_y_initial = cp.asarray(
        ynxid_data[gene_index, :], dtype=cp.uint16
    )                                                                   # [G, C]
    cu_n_initial = cp.asarray(
        ynxid_data[np.array(gene_index) + total_genes, :], dtype=cp.uint16
    )                                                                   # [G, C]
    cu_id_data   = cu_id_data.astype(cp.uint8)

    # Covariate matrix X: shape [1, num_cells, num_cov]
    # (same X for every gene; broadcasting handles multi-gene in cu_q_bb)
    num_cells = cu_y_initial.shape[1]
    if not null:
        cu_x_initial          = cp.zeros((1, num_cells, 2), dtype=cp.float32)
        cu_x_initial[:, :, 0] = 1.0                                    # intercept
        cu_x_initial[:, :, 1] = cp.asarray(ynxid_data[-1, :])          # covariate
    else:
        cu_x_initial          = cp.zeros((1, num_cells, 1), dtype=cp.float32)
        cu_x_initial[:, :, 0] = 1.0

    # Local (0-based) gene indices used inside this function
    local_index = list(range(num_genes_total))

    # Log-likelihood history: dict[local_idx] = [iteration_counter, llkl_1, llkl_2, ...]
    llkl_record_dict = {g: [0] for g in local_index}

    # ── EM iterations ────────────────────────────────────────────────────────
    for i in range(num_iteration):
        if not local_index:
            break
        print(f"\nIteration {i}  |  active genes: {len(local_index)}")
        t0 = time.time()

        still_active = []

        for g in local_index:
            # Extract single-gene slices  →  first dim is always 1
            cu_y_g             = cu_y_initial[g:g+1, :]             # [1, C]
            cu_n_g             = cu_n_initial[g:g+1, :]             # [1, C]
            cu_b_phi_g         = cu_param_result[g:g+1, :]          # [1, P]
            cu_sigm2_g         = cu_sigm2_result[g:g+1, :]          # [1, 1]
            cu_randint_g       = cu_randint_result[g:g+1, :]        # [1, I]
            cu_randint_prec_g  = cu_randint_prec_result[g:g+1, :]   # [1, I]

            # One EM step (single gene)
            (cu_b_phi_g,
             cu_sigm2_g,
             llkl_g,
             cu_randint_g,
             cu_randint_prec_g) = cu_vem_bb_scipy(
                cu_b_phi_g, cu_sigm2_g,
                cu_x_initial,           # [1, C, num_cov] – shared
                cu_randint_g, cu_randint_prec_g,
                cu_id_data,
                cu_y_g, cu_n_g,
                cu_ghq_weights, cu_ghq_nodes,
                cu_n_lap, null, max_optim,
            )

            # Write results back into the result arrays
            cu_param_result[g:g+1, :]       = cu_b_phi_g
            cu_sigm2_result[g:g+1, :]       = cu_sigm2_g
            cu_randint_result[g:g+1, :]     = cu_randint_g
            cu_randint_prec_result[g:g+1, :] = cu_randint_prec_g

            # Log-likelihood bookkeeping
            llkl_val = float(llkl_g[0])
            llkl_record_dict[g][0] += 1          # iteration counter
            llkl_record_dict[g].append(llkl_val)

            # Convergence check (same criterion as original)
            history = llkl_record_dict[g]
            if i < min_iter or len(history) < 4:
                still_active.append(g)
            elif abs((history[-1] - history[-2]) / (history[-1] + 1e-12)) > 1e-7:
                still_active.append(g)

        local_index = still_active
        elapsed = time.time() - t0
        cum_time += elapsed
        print(f"  cumulative time = {cum_time / 60:.2f} min")

    return cu_param_result, cu_sigm2_result, llkl_record_dict, cum_time


# ═══════════════════════════════════════════════════════════════════════════
#  Top-level entry point  (mirrors daesc_bb_gpu API)
# ═══════════════════════════════════════════════════════════════════════════

def daesc_bb_scipy(
    ynidx_data,           # (2G+2) × C  numpy  (same format as daesc_bb_gpu)
    num_iteration = 50,
    min_iter      = 20,
    max_optim     = 10,   # max L-BFGS-B inner iterations per M-step
):
    """
    DAESC-BB with scipy L-BFGS-B optimiser (single-gene processing).

    Parameters
    ----------
    ynidx_data    : numpy array, shape (2G + 2) × C
                    Rows 0..G-1       : alternative allele counts  Y
                    Rows G..2G-1      : total counts               N
                    Row  2G           : donor/individual IDs        (int)
                    Row  2G+1         : covariate                   X
    num_iteration : max EM iterations
    min_iter      : min EM iterations before early stopping
    max_optim     : max L-BFGS-B iterations in M-step

    Returns
    -------
    pandas.DataFrame  with columns:
        b0, b1, phi, sigma2, llkl,
        b0_null, phi_null, sigma_null, llkl_null
    """
    # ── One-hot donor matrix ─────────────────────────────────────────────────
    labels     = ynidx_data[-2, :].astype(int)
    num_classes = int(np.max(labels)) + 1
    id_data     = np.eye(num_classes)[labels].T                  # [I, C]
    cu_id_data  = cp.asarray(id_data)

    # ── Initial parameter estimates ──────────────────────────────────────────
    ynidx_T  = ynidx_data.T                                      # [C, 2G+2]
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
        test_x = np.column_stack([np.ones(len(idx)), test_x]).astype(np.float32)

        phi_hat = optimize_multiple_starts(test_y, test_n, num_starts=1)

        y_ratio = test_y / np.maximum(test_n, 1e-9)

        # H1 GLM
        X_h1     = sm.add_constant(test_x)
        glm_h1   = sm.GLM(y_ratio, X_h1, family=sm.families.Binomial()).fit()
        coef_h1  = glm_h1.params

        # H0 GLM (intercept only)
        X_h0     = np.ones((len(idx), 1))
        glm_h0   = sm.GLM(y_ratio, X_h0, family=sm.families.Binomial()).fit()
        coef_h0  = glm_h0.params[0]

        results.append({
            "Intercept"           : coef_h1[0],
            "Coef_x1"             : coef_h1[1],
            "Fixed_value"         : 0.05,        # sigma2 initial (H1)
            "Estimated_Phi"       : phi_hat,
            "New_Intercept"       : coef_h0,
            "Zero_Value"          : 0.0,
            "Repeated_Fixed_value": 0.05,         # sigma2 initial (H0)
            "Repeated_Phi"        : phi_hat,
        })

    results_df       = pd.DataFrame(results)
    initial_bb_param = results_df.to_numpy()     # G × 8
    t_end            = time.time()
    print(f"Initial estimates done in {(t_end - t_start) / 60:.2f} min")

    gene_index = list(range(num_gene))
    cum_time   = t_end - t_start

    # ── H1 model ─────────────────────────────────────────────────────────────
    print("\n── H1 model ──")
    param_h1, sigm2_h1, llkl_h1_dict, cum_time = VEM_bb_scipy(
        gene_index,
        num_iteration = num_iteration,
        min_iter      = min_iter,
        ynxid_data    = ynidx_data,
        cu_id_data    = cu_id_data,
        cu_n_lap      = 2,
        initial_param = initial_bb_param,
        null          = False,
        cum_time      = cum_time,
        max_optim     = max_optim,
    )
    param_h1_np  = param_h1.get()                                # [G, 3]
    sigm2_h1_np  = sigm2_h1.get()                                # [G, 1]
    llkl_h1_np   = np.array([llkl_h1_dict[i][-1] for i in range(num_gene)])

    clear_memory()

    # ── H0 model ─────────────────────────────────────────────────────────────
    print("\n── H0 model ──")
    param_h0, sigm2_h0, llkl_h0_dict, cum_time = VEM_bb_scipy(
        gene_index,
        num_iteration = num_iteration,
        min_iter      = min_iter,
        ynxid_data    = ynidx_data,
        cu_id_data    = cu_id_data,
        cu_n_lap      = 2,
        initial_param = initial_bb_param,
        null          = True,
        cum_time      = cum_time,
        max_optim     = max_optim,
    )
    param_h0_np  = param_h0.get()                                # [G, 2]
    sigm2_h0_np  = sigm2_h0.get()                                # [G, 1]
    llkl_h0_np   = np.array([llkl_h0_dict[i][-1] for i in range(num_gene)])

    clear_memory()

    # ── Assemble output DataFrame ─────────────────────────────────────────────
    rows = []
    for idx in range(num_gene):
        rows.append([
            param_h1_np[idx, 0],   # b0
            param_h1_np[idx, 1],   # b1
            param_h1_np[idx, 2],   # phi
            sigm2_h1_np[idx, 0],   # sigma2
            llkl_h1_np[idx],       # llkl
            param_h0_np[idx, 0],   # b0_null
            param_h0_np[idx, 1],   # phi_null
            sigm2_h0_np[idx, 0],   # sigma_null
            llkl_h0_np[idx],       # llkl_null
        ])

    column_names = ["b0", "b1", "phi", "sigma2", "llkl",
                    "b0_null", "phi_null", "sigma_null", "llkl_null"]
    return pd.DataFrame(rows, columns=column_names)
































