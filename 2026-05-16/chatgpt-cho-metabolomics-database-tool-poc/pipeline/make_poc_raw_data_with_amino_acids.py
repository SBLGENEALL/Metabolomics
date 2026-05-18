#!/usr/bin/env python
"""Create a plausible amino-acid augmented CHO raw-data workbook for POC runs."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd


AA_BASELINE_MM = {
    "Ala": 1.1,
    "Arg": 1.0,
    "Asn": 2.2,
    "Asp": 0.55,
    "Cys": 0.18,
    "Gly": 1.0,
    "His": 0.45,
    "Ile": 0.85,
    "Leu": 1.25,
    "Lys": 1.05,
    "Met": 0.35,
    "Phe": 0.45,
    "Pro": 0.8,
    "Ser": 1.7,
    "Thr": 0.95,
    "Trp": 0.12,
    "Tyr": 0.28,
    "Val": 1.05,
}

AA_DAILY_CONSUMPTION_MM = {
    "Ala": -0.02,
    "Arg": 0.09,
    "Asn": 0.18,
    "Asp": 0.04,
    "Cys": 0.018,
    "Gly": 0.05,
    "His": 0.035,
    "Ile": 0.075,
    "Leu": 0.11,
    "Lys": 0.10,
    "Met": 0.032,
    "Phe": 0.04,
    "Pro": 0.045,
    "Ser": 0.13,
    "Thr": 0.08,
    "Trp": 0.012,
    "Tyr": 0.022,
    "Val": 0.09,
}

FEED4_AA_MM = {
    "Ala": 1.0,
    "Arg": 8.0,
    "Asn": 12.0,
    "Asp": 2.0,
    "Cys": 1.0,
    "Gly": 4.0,
    "His": 3.0,
    "Ile": 7.0,
    "Leu": 9.0,
    "Lys": 8.0,
    "Met": 3.0,
    "Phe": 4.0,
    "Pro": 5.0,
    "Ser": 10.0,
    "Thr": 7.0,
    "Trp": 1.0,
    "Tyr": 2.0,
    "Val": 8.0,
}

CELLBOOST_AA_MM = {
    "Ala": 0.8,
    "Arg": 6.0,
    "Asn": 9.0,
    "Asp": 1.5,
    "Cys": 0.7,
    "Gly": 3.0,
    "His": 2.0,
    "Ile": 5.5,
    "Leu": 7.0,
    "Lys": 6.0,
    "Met": 2.0,
    "Phe": 3.0,
    "Pro": 4.0,
    "Ser": 7.5,
    "Thr": 5.5,
    "Trp": 0.7,
    "Tyr": 1.5,
    "Val": 6.0,
}


def parse_day(value) -> float:
    text = str(value)
    digits = "".join(ch for ch in text if ch.isdigit() or ch == ".")
    return float(digits) if digits else float("nan")


def clone_factor(sample_id: str) -> float:
    digits = "".join(ch for ch in str(sample_id) if ch.isdigit())
    idx = int(digits) if digits else 1
    return 0.88 + 0.035 * idx


def deterministic_noise(sample_id: str, day: float, aa: str) -> float:
    seed = sum(ord(ch) for ch in f"{sample_id}-{day}-{aa}")
    return ((seed % 17) - 8) / 100.0


def add_amino_acid_scenario(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    first_day = out.groupby("Sample ID")["DAY"].transform(lambda s: min(parse_day(v) for v in s))
    for aa, baseline in AA_BASELINE_MM.items():
        values = []
        for _, row in out.iterrows():
            day = parse_day(row["DAY"])
            elapsed = max(0.0, day - parse_day(first_day.loc[row.name]))
            factor = clone_factor(row["Sample ID"])
            feed4_umol = float(row.get("Feed4 mL", 0) or 0) * FEED4_AA_MM[aa]
            cellboost_umol = float(row.get("CellBoost mL", 0) or 0) * CELLBOOST_AA_MM[aa]
            volume = float(row.get("Culture Volume mL", 30) or 30)
            feed_contribution = (feed4_umol + cellboost_umol) / max(volume, 1.0)
            consumption = AA_DAILY_CONSUMPTION_MM[aa] * elapsed * factor
            value = baseline + feed_contribution - consumption + deterministic_noise(row["Sample ID"], day, aa)
            values.append(max(0.01, round(value, 4)))
        out[aa] = values
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sheet", default=None)
    args = parser.parse_args()

    xl = pd.ExcelFile(args.input)
    sheet = args.sheet or xl.sheet_names[0]
    raw = pd.read_excel(args.input, sheet_name=sheet)
    augmented = add_amino_acid_scenario(raw)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        augmented.to_excel(writer, sheet_name=sheet, index=False)
        guide = pd.DataFrame(
            {
                "note": [
                    "POC-only simulated amino acid values.",
                    "Column names follow standard 3-letter amino acid abbreviations.",
                    "Only amino-acid measurement columns are added; feed amino-acid composition columns are not added.",
                ]
            }
        )
        guide.to_excel(writer, sheet_name="AA_POC_Guide", index=False)
    print(f"Wrote amino-acid augmented POC workbook to {args.output}")


if __name__ == "__main__":
    main()
