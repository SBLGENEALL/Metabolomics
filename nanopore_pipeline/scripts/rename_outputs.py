#!/usr/bin/env python3
"""
rename_outputs.py

Bulk-rename/organize already-existing analysis outputs (from this pipeline,
EPI2ME, or any other tool) so their filenames match the sample names in
samplesheet.csv instead of raw barcode IDs.

For each row "barcode_id,reference_fasta,sample_name" in the samplesheet,
every file or folder under --results-dir whose name contains the
barcode_id is renamed by replacing that barcode_id with sample_name
(both in the path and, optionally, inside the file's content for
text-based formats such as .vcf, .fasta and .txt headers).

Usage:
    python scripts/rename_outputs.py \
        --samplesheet samplesheet.csv \
        --results-dir results \
        [--rename-content] [--dry-run]
"""
import argparse
import csv
import os
import shutil


def load_samplesheet(path):
    mapping = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            barcode = row["barcode_id"].strip()
            sample = row["sample_name"].strip()
            if barcode and sample:
                mapping.append((barcode, sample))
    # Rename longer barcode IDs first to avoid partial-match collisions
    mapping.sort(key=lambda x: len(x[0]), reverse=True)
    return mapping


def rename_paths(results_dir, mapping, dry_run):
    # Walk bottom-up so files are renamed before their parent directories
    for root, dirs, files in os.walk(results_dir, topdown=False):
        for name in files + dirs:
            old_path = os.path.join(root, name)
            new_name = name
            for barcode, sample in mapping:
                if barcode in new_name:
                    new_name = new_name.replace(barcode, sample)
            if new_name != name:
                new_path = os.path.join(root, new_name)
                print(f"{old_path}  ->  {new_path}")
                if not dry_run:
                    if os.path.exists(new_path):
                        raise FileExistsError(f"Target already exists: {new_path}")
                    shutil.move(old_path, new_path)


def rename_content(results_dir, mapping, dry_run, extensions):
    for root, _dirs, files in os.walk(results_dir):
        for name in files:
            if not any(name.endswith(ext) for ext in extensions):
                continue
            path = os.path.join(root, name)
            try:
                with open(path, "r", errors="ignore") as f:
                    text = f.read()
            except (UnicodeDecodeError, OSError):
                continue
            new_text = text
            for barcode, sample in mapping:
                new_text = new_text.replace(barcode, sample)
            if new_text != text:
                print(f"updating contents: {path}")
                if not dry_run:
                    with open(path, "w") as f:
                        f.write(new_text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", default="samplesheet.csv")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument(
        "--rename-content",
        action="store_true",
        help="Also replace barcode IDs with sample names inside text files "
             "(.fasta, .vcf, .txt, .md headers)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the changes that would be made without touching anything",
    )
    args = parser.parse_args()

    mapping = load_samplesheet(args.samplesheet)
    if not mapping:
        print("No barcode/sample_name pairs found in samplesheet.")
        return

    print("Renaming using:")
    for barcode, sample in mapping:
        print(f"  {barcode} -> {sample}")
    print()

    rename_paths(args.results_dir, mapping, args.dry_run)

    if args.rename_content:
        rename_content(
            args.results_dir, mapping, args.dry_run,
            extensions=(".fasta", ".fa", ".vcf", ".txt", ".md"),
        )

    if args.dry_run:
        print("\n(dry run - no files were changed)")


if __name__ == "__main__":
    main()
