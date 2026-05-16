#!/usr/bin/env python
"""Build a standalone Escher HTML overlay for iCHO3K exchange fluxes.

This script is optional. It requires COBRApy and Escher. If they are not
installed, it exits cleanly and leaves the CSV/JSON overlay files available.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = ROOT / "models" / "iCHO3K" / "Model" / "iCHO3K_cho_prod_generic_unblocked.json"
DEFAULT_REACTION_DATA = ROOT / "results" / "interactive_run" / "icho3k_inputs" / "escher_reaction_data_mean_flux.json"
DEFAULT_OUTPUT = ROOT / "results" / "interactive_run" / "icho3k_inputs" / "escher_flux_overlay.html"


def patch_cobra_cache() -> None:
    """Keep COBRApy cache inside the project when default user cache is blocked."""
    try:
        import appdirs
    except ImportError:
        return
    cache_dir = ROOT / ".cache" / "cobrapy"
    cache_dir.mkdir(parents=True, exist_ok=True)
    appdirs.user_cache_dir = lambda appname=None, appauthor=None, **kwargs: str(cache_dir)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--map-json", type=Path, help="Optional Escher map JSON compatible with iCHO3K reaction IDs")
    parser.add_argument("--reaction-data", type=Path, default=DEFAULT_REACTION_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    patch_cobra_cache()
    try:
        import cobra
        import escher
    except ImportError:
        print("COBRApy/Escher not installed; skipped Escher HTML generation.")
        return

    if not args.model.exists():
        raise SystemExit(f"Model not found: {args.model}")
    if not args.reaction_data.exists():
        raise SystemExit(f"Reaction data not found: {args.reaction_data}")

    model = cobra.io.load_json_model(str(args.model))
    reaction_data = json.loads(args.reaction_data.read_text(encoding="utf-8"))
    builder_kwargs = {"model": model, "reaction_data": reaction_data}
    if args.map_json:
        if not args.map_json.exists():
            raise SystemExit(f"Escher map JSON not found: {args.map_json}")
        builder_kwargs["map_json"] = args.map_json.read_text(encoding="utf-8")
    builder = escher.Builder(**builder_kwargs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    builder.save_html(str(args.output))
    print(f"Wrote Escher overlay HTML to {args.output}")


if __name__ == "__main__":
    main()
