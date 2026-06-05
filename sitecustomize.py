"""Repository-wide runtime customizations for the CHO metabolomics pipeline.

Python imports ``sitecustomize`` automatically when it is present on ``sys.path``.
The pipeline launches each step with the repository root as the working directory,
so this file is loaded before the plotting scripts import matplotlib.

Purpose for v1.0
----------------
The calculation scripts were validated during development while still carrying
legacy figure labels such as Figure 10, 15, and 21. For the v1.0 release we keep
those validated computational scripts stable and standardize only the user-facing
figure labels at render time.
"""

FIGURE_TITLE_REPLACEMENTS = {
    "Figure 5. Exchange Rate Comparison": "Figure 4. Exchange Rate Comparison",
    "Figure 10. Clone/group central-metabolism and mAb-pathway flux heatmap": "Figure 5. Clone/group central-metabolism and mAb-pathway flux heatmap",
    "Figure 10B. Row-wise relative flux pattern across clones/groups": "Figure 6. Row-wise relative flux pattern across clones/groups",
    "Figure 11. High-vs-Low group flux differences across central metabolism and mAb pathways": "Figure 7. High-vs-Low group flux differences across central metabolism and mAb pathways",
    "Figure 21. Genome-scale/internal FVA High-vs-Low separation": "Figure 9. Full/internal/all FVA High-vs-Low separation",
    "Figure 21B. Genome-scale/internal FVA flexibility delta": "Figure 10. Full/internal/all FVA flexibility delta",
    "CHO Fed-Batch Clone Flux-Comparison Summary": "Figure 11. CHO Fed-Batch Clone Flux-Comparison Summary",
    "Figure 13.": "Supplementary Figure 1.",
}


def _standardize_title(text):
    if not isinstance(text, str):
        return text
    out = text
    for old, new in FIGURE_TITLE_REPLACEMENTS.items():
        out = out.replace(old, new)
    return out


try:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    _orig_axes_set_title = Axes.set_title
    _orig_figure_suptitle = Figure.suptitle

    def _patched_axes_set_title(self, label, *args, **kwargs):
        return _orig_axes_set_title(self, _standardize_title(label), *args, **kwargs)

    def _patched_figure_suptitle(self, t, *args, **kwargs):
        return _orig_figure_suptitle(self, _standardize_title(t), *args, **kwargs)

    Axes.set_title = _patched_axes_set_title
    Figure.suptitle = _patched_figure_suptitle
except Exception:
    # Never let title normalization interfere with numerical pipeline execution.
    pass
