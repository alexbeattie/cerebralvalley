"""Matplotlib figures for the ISM map and the AlphaGenome cross-check.

Headless by design (Agg backend): the demo writes PNGs to `artifacts/` on a server with no
display. Kept deliberately plain -- these are diagnostic plots for a judge to eyeball, not a
polished figure panel.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # no display needed; must be set before pyplot import

import matplotlib.pyplot as plt
import numpy as np

from .encoding import BASES


def plot_ism_heatmap(
    matrix: np.ndarray,
    seq: str,
    out_path: str,
    *,
    title: str = "ISM mutation-effect map",
    highlight: tuple[int, int] | None = None,
) -> str:
    """Heatmap of the (L, 4) ISM delta matrix (bases on y, position on x). Saves a PNG.

    A diverging colormap centered at 0 shows activity-increasing vs -decreasing mutations;
    `highlight` (start, end) draws a box around the known/interesting span.
    """
    length = matrix.shape[0]
    vmax = float(np.abs(matrix).max()) or 1.0
    fig, ax = plt.subplots(figsize=(min(0.06 * length + 2, 20), 2.6))
    im = ax.imshow(
        matrix.T,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_yticks(range(4))
    ax.set_yticklabels(list(BASES))
    ax.set_xlabel("position (bp)")
    ax.set_ylabel("alt base")
    ax.set_title(title)
    if highlight is not None:
        start, end = highlight
        ax.add_patch(
            plt.Rectangle((start - 0.5, -0.5), end - start, 4, fill=False, edgecolor="black", lw=1.5)
        )
    fig.colorbar(im, ax=ax, label="Δ predicted activity")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return out_path


def plot_importance_track(
    track: np.ndarray,
    out_path: str,
    *,
    title: str = "Per-base ISM importance",
    highlight: tuple[int, int] | None = None,
) -> str:
    """Line plot of the per-position importance track. Saves a PNG."""
    fig, ax = plt.subplots(figsize=(min(0.06 * track.shape[0] + 2, 20), 2.6))
    ax.plot(track, color="#1f77b4", lw=1.0)
    ax.set_xlabel("position (bp)")
    ax.set_ylabel("mean |Δ activity|")
    ax.set_title(title)
    if highlight is not None:
        start, end = highlight
        ax.axvspan(start - 0.5, end - 0.5, color="orange", alpha=0.3, label="planted motif")
        ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return out_path


def plot_comparison_scatter(
    our_deltas: np.ndarray,
    ag_scores: np.ndarray,
    out_path: str,
    *,
    labels: list[str] | None = None,
    correlation: float | None = None,
) -> str:
    """Scatter of our-model ISM delta vs AlphaGenome magnitude for the cross-check. Saves PNG."""
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.scatter(our_deltas, ag_scores, color="#2ca02c", s=40, zorder=3)
    if labels is not None:
        for x, y, lab in zip(our_deltas, ag_scores, labels):
            ax.annotate(lab, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axhline(0, color="gray", lw=0.6)
    ax.axvline(0, color="gray", lw=0.6)
    ax.set_xlabel("regmodel |ISM Δ| (ours)")
    ax.set_ylabel("AlphaGenome magnitude")
    sub = f"  (Pearson r={correlation:.2f})" if correlation is not None and np.isfinite(correlation) else ""
    ax.set_title(f"regmodel vs AlphaGenome{sub}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return out_path
