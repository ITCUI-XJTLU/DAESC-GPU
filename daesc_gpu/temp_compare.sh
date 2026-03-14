"""
compare_models.py
─────────────────
Compare daesc_bb_gpu  vs  daesc_bb_scipy   (custom CuPy BFGS  vs  scipy L-BFGS-B)
    AND
    daesc_mix_gpu vs  daesc_mix_scipy  (same optimiser swap for the Mix model)

Usage
-----
    python compare_models.py

The script prints per-column statistics for both model families and saves
a combined scatter-plot PDF:  daesc_comparison.pdf
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats as stats

# ── model modules ─────────────────────────────────────────────────────────
import daesc_bb_gpu
import daesc_bb_scipy
import daesc_mix_gpu
import daesc_mix_scipy


# ═══════════════════════════════════════════════════════════════════════════
# 1.  Load data  &  subset to first 10 genes for quick validation
# ═══════════════════════════════════════════════════════════════════════════
print("Loading data …")
ynx_data_cuomo = pd.read_parquet(
    "/content/drive/My Drive/UW/Stat_Gene/Data/"
    "Paper_output_data/Cuomo_data/Data/ynxid_cuomo.parquet"
)
ynx_data_cuomo = ynx_data_cuomo.to_numpy().astype(np.float32)
print(f"Full data shape: {ynx_data_cuomo.shape}")

# ── Subset to first K genes ───────────────────────────────────────────────
# Layout of ynx_data_cuomo  (rows × cells):
#   rows   0 .. G-1       →  Y  (alternative counts),  G = 4102 genes
#   rows   G .. 2G-1      →  N  (total counts)
#   row    2G             →  donor / individual IDs
#   row    2G+1           →  covariate X
K       = 10
G_full  = (ynx_data_cuomo.shape[0] - 2) // 2   # = 4102

Y_sub   = ynx_data_cuomo[0        :K,            :]
N_sub   = ynx_data_cuomo[G_full   :G_full + K,   :]
id_row  = ynx_data_cuomo[2*G_full :2*G_full + 1, :]
x_row   = ynx_data_cuomo[2*G_full + 1:,          :]

ynx_data_cuomo = np.vstack([Y_sub, N_sub, id_row, x_row])   # [22, C]
print(f"Subset shape (first {K} genes): {ynx_data_cuomo.shape}")


# ═══════════════════════════════════════════════════════════════════════════
# 2.  Shared hyper-parameters
# ═══════════════════════════════════════════════════════════════════════════
ITER     = 50
MIN_ITER = 20
MAX_OPT  = 10


# ═══════════════════════════════════════════════════════════════════════════
# 3.  Run BB models
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "═"*65)
print("BB model  ──  original (custom CuPy BFGS)")
print("═"*65)
res_bb_orig = daesc_bb_gpu.daesc_bb_gpu(
    ynx_data_cuomo,
    num_iteration = ITER,
    min_iter      = MIN_ITER,
    max_optim     = MAX_OPT,
)

print("\n" + "═"*65)
print("BB model  ──  scipy L-BFGS-B")
print("═"*65)
res_bb_scipy = daesc_bb_scipy.daesc_bb_scipy(
    ynx_data_cuomo,
    num_iteration = ITER,
    min_iter      = MIN_ITER,
    max_optim     = MAX_OPT,
)


# ═══════════════════════════════════════════════════════════════════════════
# 4.  Run Mix models
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "═"*65)
print("Mix model  ──  original (custom CuPy BFGS)")
print("═"*65)
res_mix_orig = daesc_mix_gpu.daesc_mix_gpu(
    ynx_data_cuomo,
    num_iteration = ITER,
    min_iter      = MIN_ITER,
    max_optim     = MAX_OPT,
)

print("\n" + "═"*65)
print("Mix model  ──  scipy L-BFGS-B")
print("═"*65)
res_mix_scipy = daesc_mix_scipy.daesc_mix_scipy(
    ynx_data_cuomo,
    num_iteration = ITER,
    min_iter      = MIN_ITER,
    max_optim     = MAX_OPT,
)


# ═══════════════════════════════════════════════════════════════════════════
# 5.  Numerical comparison helper
# ═══════════════════════════════════════════════════════════════════════════

def print_comparison_table(res_orig, res_scipy, model_name):
    """Print Pearson r / MAD / RMSE for every shared column."""
    cols = res_orig.columns.tolist()
    print(f"\n{'═'*65}")
    print(f"  {model_name}  —  original vs scipy")
    print(f"{'─'*65}")
    print(f"{'Column':<16} {'Pearson r':>10} {'MAD':>12} {'RMSE':>12}")
    print(f"{'─'*65}")
    for col in cols:
        ov = res_orig[col].values
        sv = res_scipy[col].values
        ok = np.isfinite(ov) & np.isfinite(sv)
        if ok.sum() < 2:
            print(f"{col:<16}  (insufficient finite values)")
            continue
        r, _  = stats.pearsonr(ov[ok], sv[ok])
        mad   = np.mean(np.abs(ov[ok] - sv[ok]))
        rmse  = np.sqrt(np.mean((ov[ok] - sv[ok]) ** 2))
        print(f"{col:<16} {r:>10.6f} {mad:>12.6f} {rmse:>12.6f}")
    print(f"{'═'*65}")


def print_llkl_diff(res_orig, res_scipy, model_name, suffixes=("llkl", "llkl_null")):
    """Compare log-likelihoods and report how often scipy is strictly better."""
    for col, label in zip(suffixes, ("H1", "H0")):
        ov = res_orig[col].values
        sv = res_scipy[col].values
        ok = np.isfinite(ov) & np.isfinite(sv)
        d  = sv[ok] - ov[ok]
        n_better = int((d > 0).sum())
        print(f"\n  {model_name}  Log-lik diff ({label})  scipy − original:")
        print(f"    mean={d.mean():.6f}  std={d.std():.6f}"
              f"  min={d.min():.6f}  max={d.max():.6f}")
        print(f"    scipy strictly better: {n_better}/{ok.sum()} "
              f"({100*n_better/ok.sum():.1f}%)")


# ── Print tables ──────────────────────────────────────────────────────────
print_comparison_table(res_bb_orig,  res_bb_scipy,  "DAESC-BB")
print_comparison_table(res_mix_orig, res_mix_scipy, "DAESC-Mix")

print("\nFirst 5 rows  ──  BB original:");  print(res_bb_orig.head())
print("\nFirst 5 rows  ──  BB scipy:");     print(res_bb_scipy.head())
print("\nFirst 5 rows  ──  Mix original:"); print(res_mix_orig.head())
print("\nFirst 5 rows  ──  Mix scipy:");    print(res_mix_scipy.head())

# ── Log-likelihood comparison ─────────────────────────────────────────────
print_llkl_diff(res_bb_orig,  res_bb_scipy,  "DAESC-BB")
print_llkl_diff(res_mix_orig, res_mix_scipy, "DAESC-Mix")


# ═══════════════════════════════════════════════════════════════════════════
# 6.  Scatter plots  (BB + Mix side by side in one PDF)
# ═══════════════════════════════════════════════════════════════════════════

def add_scatter_block(fig, gs, row_offset, res_orig, res_scipy, label_orig, label_scipy, block_title):
    """
    Add one block of scatter-plot panels (one per column) starting at
    row `row_offset` within the GridSpec `gs`.
    Returns the number of rows consumed.
    """
    cols  = res_orig.columns.tolist()
    ncols = 3
    nrows = int(np.ceil(len(cols) / ncols))

    for k, col in enumerate(cols):
        r = row_offset + k // ncols
        c = k % ncols
        ax = fig.add_subplot(gs[r, c])

        ov = res_orig[col].values
        sv = res_scipy[col].values
        ok = np.isfinite(ov) & np.isfinite(sv)

        ax.scatter(ov[ok], sv[ok], s=6, alpha=0.5,
                   color="steelblue", rasterized=True)

        lo = min(ov[ok].min(), sv[ok].min())
        hi = max(ov[ok].max(), sv[ok].max())
        ax.plot([lo, hi], [lo, hi], "r--", lw=0.8)

        if ok.sum() >= 2:
            r_val, _ = stats.pearsonr(ov[ok], sv[ok])
            ax.set_title(f"{block_title} · {col}  (r={r_val:.4f})", fontsize=8)
        else:
            ax.set_title(f"{block_title} · {col}  (n/a)", fontsize=8)

        ax.set_xlabel(label_orig,  fontsize=7)
        ax.set_ylabel(label_scipy, fontsize=7)
        ax.tick_params(labelsize=6)

    return nrows


# Determine grid dimensions
bb_cols  = res_bb_orig.columns.tolist()
mix_cols = res_mix_orig.columns.tolist()
NCOLS    = 3
bb_rows  = int(np.ceil(len(bb_cols)  / NCOLS))
mix_rows = int(np.ceil(len(mix_cols) / NCOLS))
total_rows = bb_rows + mix_rows + 2   # +2 for section title rows

fig = plt.figure(figsize=(5 * NCOLS, 4 * total_rows))
gs  = gridspec.GridSpec(total_rows, NCOLS, figure=fig,
                        hspace=0.55, wspace=0.38)

# ── BB block ──────────────────────────────────────────────────────────────
# Section title (spans all columns)
ax_title_bb = fig.add_subplot(gs[0, :])
ax_title_bb.axis("off")
ax_title_bb.text(0.5, 0.5, "DAESC-BB: Original (CuPy BFGS) vs Scipy L-BFGS-B",
                 ha="center", va="center", fontsize=13, fontweight="bold")

add_scatter_block(fig, gs,
                  row_offset=1,
                  res_orig=res_bb_orig,
                  res_scipy=res_bb_scipy,
                  label_orig="Original (CuPy BFGS)",
                  label_scipy="Scipy L-BFGS-B",
                  block_title="BB")

# ── Mix block ─────────────────────────────────────────────────────────────
mix_start = 1 + bb_rows
ax_title_mix = fig.add_subplot(gs[mix_start, :])
ax_title_mix.axis("off")
ax_title_mix.text(0.5, 0.5, "DAESC-Mix: Original (CuPy BFGS) vs Scipy L-BFGS-B",
                  ha="center", va="center", fontsize=13, fontweight="bold")

add_scatter_block(fig, gs,
                  row_offset=mix_start + 1,
                  res_orig=res_mix_orig,
                  res_scipy=res_mix_scipy,
                  label_orig="Original (CuPy BFGS)",
                  label_scipy="Scipy L-BFGS-B",
                  block_title="Mix")

plt.savefig("daesc_comparison.pdf", bbox_inches="tight", dpi=150)
print("\nScatter plots saved to  daesc_comparison.pdf")