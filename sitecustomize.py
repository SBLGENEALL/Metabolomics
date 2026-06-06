"""Repository-wide runtime customizations for the CHO metabolomics pipeline.

Python imports ``sitecustomize`` automatically when it is present on ``sys.path``.
The pipeline launches each step with the repository root as the working directory,
so this file is loaded before the plotting scripts import matplotlib.

Purpose
-------
Reserve a repository-wide title normalization hook. Current step scripts use the
final release titles directly, so no runtime replacements are required.
"""

FIGURE_TITLE_REPLACEMENTS = {}


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
