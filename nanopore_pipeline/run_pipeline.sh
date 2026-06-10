#!/usr/bin/env bash
#
# Oxford Nanopore (MinION) reference-mapping pipeline
#
# Folder-name based matching: for every experiment folder that exists
# under BOTH references_dir and data_dir (same folder name, e.g.
# "20260610_pUC19_test"), the pipeline maps that experiment's reads
# against the matching reference(s). Output files are named after the
# reference FASTA (the real vector/sample name), not the barcode ID,
# so no manual renaming is needed afterwards.
#
# Layout:
#   references/<date>_<experiment_name>/<reference_name>.fasta
#   data/raw/<date>_<experiment_name>/...        (fastq.gz, optionally
#                                                   inside barcodeXX/ subfolders)
#   results/<date>_<experiment_name>/<reference_name>/...
#
# Usage:
#   ./run_pipeline.sh [config.yaml]
#
# Requires: minimap2, samtools, bcftools (or medaka if variant_caller=medaka)
#           NanoFilt (optional, only used if min_read_length/quality > 0)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-$SCRIPT_DIR/config.yaml}"

# ---- minimal flat YAML reader (key: value, no nesting) -------------------
read_cfg() {
    local key="$1" default="$2"
    local val
    val=$(grep -E "^${key}:" "$CONFIG" | head -n1 | sed -E "s/^${key}:[[:space:]]*//; s/[[:space:]]*(#.*)?$//")
    if [[ -z "$val" ]]; then echo "$default"; else echo "$val"; fi
}

REF_ROOT="$SCRIPT_DIR/$(read_cfg references_dir references)"
DATA_ROOT="$SCRIPT_DIR/$(read_cfg data_dir data/raw)"
RESULTS_ROOT="$SCRIPT_DIR/$(read_cfg results_dir results)"
PRESET="$(read_cfg minimap2_preset map-ont)"
THREADS="$(read_cfg threads 4)"
MIN_LEN="$(read_cfg min_read_length 0)"
MIN_QUAL="$(read_cfg min_read_quality 0)"
VARIANT_CALLER="$(read_cfg variant_caller bcftools)"

echo "==================================================================="
echo " Nanopore reference-mapping pipeline"
echo "   references  : $REF_ROOT"
echo "   data        : $DATA_ROOT"
echo "   results     : $RESULTS_ROOT"
echo "   preset      : $PRESET   threads: $THREADS"
echo "==================================================================="

mkdir -p "$RESULTS_ROOT"

shopt -s nullglob
for EXP_DIR in "$REF_ROOT"/*/; do
    EXP_NAME="$(basename "$EXP_DIR")"
    DATA_EXP_DIR="$DATA_ROOT/$EXP_NAME"

    echo ""
    echo "==================================================================="
    echo "Experiment: $EXP_NAME"
    echo "==================================================================="

    if [[ ! -d "$DATA_EXP_DIR" ]]; then
        echo "  [SKIP] no matching data folder: $DATA_EXP_DIR" >&2
        continue
    fi

    # Collect reference fasta files for this experiment
    REF_FILES=("$EXP_DIR"*.fasta "$EXP_DIR"*.fa "$EXP_DIR"*.fna)
    if [[ ${#REF_FILES[@]} -eq 0 ]]; then
        echo "  [SKIP] no reference fasta found in $EXP_DIR" >&2
        continue
    fi

    # Collect all fastq(.gz) files for this experiment, recursively
    # (handles both flat layout and barcodeXX/ subfolders)
    FASTQ_FILES=()
    while IFS= read -r -d '' f; do FASTQ_FILES+=("$f"); done \
        < <(find "$DATA_EXP_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fastq" \) -print0)

    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "  [SKIP] no fastq files found under $DATA_EXP_DIR" >&2
        continue
    fi

    EXP_RESULTS_DIR="$RESULTS_ROOT/$EXP_NAME"
    mkdir -p "$EXP_RESULTS_DIR"

    # 1. Merge all reads for this experiment once (shared across references)
    MERGED_FASTQ="$EXP_RESULTS_DIR/${EXP_NAME}.merged.fastq.gz"
    if [[ ! -s "$MERGED_FASTQ" ]]; then
        echo "  [1/5] Merging $(printf '%d' ${#FASTQ_FILES[@]}) read file(s) -> $(basename "$MERGED_FASTQ")"
        : > "${MERGED_FASTQ%.gz}.tmp"
        for f in "${FASTQ_FILES[@]}"; do
            if [[ "$f" == *.gz ]]; then
                zcat "$f" >> "${MERGED_FASTQ%.gz}.tmp"
            else
                cat "$f" >> "${MERGED_FASTQ%.gz}.tmp"
            fi
        done
        gzip -c "${MERGED_FASTQ%.gz}.tmp" > "$MERGED_FASTQ"
        rm -f "${MERGED_FASTQ%.gz}.tmp"
    else
        echo "  [1/5] Merged reads already exist, skipping"
    fi

    # 2. Optional QC filtering with NanoFilt
    READS_FOR_MAPPING="$MERGED_FASTQ"
    if [[ "$MIN_LEN" -gt 0 || "$MIN_QUAL" -gt 0 ]] && command -v NanoFilt >/dev/null 2>&1; then
        FILTERED_FASTQ="$EXP_RESULTS_DIR/${EXP_NAME}.filtered.fastq.gz"
        echo "  [2/5] QC filtering (length>=$MIN_LEN, quality>=$MIN_QUAL) -> $(basename "$FILTERED_FASTQ")"
        zcat "$MERGED_FASTQ" | NanoFilt -l "$MIN_LEN" -q "$MIN_QUAL" | gzip > "$FILTERED_FASTQ"
        READS_FOR_MAPPING="$FILTERED_FASTQ"
    else
        echo "  [2/5] QC filtering skipped"
    fi

    # 3-5. Map against every reference fasta found in this experiment folder
    for REF_PATH in "${REF_FILES[@]}"; do
        REF_FILE="$(basename "$REF_PATH")"
        REF_NAME="${REF_FILE%.*}"
        SAMPLE_DIR="$EXP_RESULTS_DIR/$REF_NAME"
        mkdir -p "$SAMPLE_DIR"

        echo ""
        echo "  --- Reference: $REF_FILE -> results/$EXP_NAME/$REF_NAME/ ---"

        # 3. Map to reference with minimap2, sort with samtools
        SORTED_BAM="$SAMPLE_DIR/${REF_NAME}.sorted.bam"
        echo "  [3/5] Mapping -> $(basename "$SORTED_BAM")"
        minimap2 -ax "$PRESET" -t "$THREADS" "$REF_PATH" "$READS_FOR_MAPPING" \
            | samtools sort -@ "$THREADS" -o "$SORTED_BAM" -
        samtools index "$SORTED_BAM"
        samtools flagstat "$SORTED_BAM" > "$SAMPLE_DIR/${REF_NAME}.flagstat.txt"
        samtools depth -a "$SORTED_BAM" > "$SAMPLE_DIR/${REF_NAME}.depth.txt"

        # 4. Consensus sequence
        CONSENSUS_FASTA="$SAMPLE_DIR/${REF_NAME}.consensus.fasta"
        echo "  [4/5] Building consensus -> $(basename "$CONSENSUS_FASTA")"
        samtools consensus -a -f fasta "$SORTED_BAM" > "$CONSENSUS_FASTA"
        sed -i "1s/.*/>${REF_NAME}/" "$CONSENSUS_FASTA"

        # 5. Variant calling
        VCF="$SAMPLE_DIR/${REF_NAME}.vcf.gz"
        echo "  [5/5] Calling variants ($VARIANT_CALLER) -> $(basename "$VCF")"
        if [[ "$VARIANT_CALLER" == "medaka" ]] && command -v medaka_haploid_variant >/dev/null 2>&1; then
            medaka_haploid_variant -i "$READS_FOR_MAPPING" -r "$REF_PATH" -o "$SAMPLE_DIR/medaka"
            cp "$SAMPLE_DIR/medaka/medaka.annotated.vcf" "$SAMPLE_DIR/${REF_NAME}.vcf" 2>/dev/null || true
            bgzip -f "$SAMPLE_DIR/${REF_NAME}.vcf"
        else
            bcftools mpileup -f "$REF_PATH" "$SORTED_BAM" 2>/dev/null \
                | bcftools call -mv -Oz -o "$VCF"
            bcftools index -f "$VCF"
        fi

        # Per-sample summary report
        REPORT="$SAMPLE_DIR/${REF_NAME}_report.md"
        MAPPED=$(grep "mapped (" "$SAMPLE_DIR/${REF_NAME}.flagstat.txt" | head -1)
        MEAN_DEPTH=$(awk '{sum+=$3; n++} END {if (n>0) printf "%.2f", sum/n; else print "0"}' "$SAMPLE_DIR/${REF_NAME}.depth.txt")
        {
            echo "# Report: $REF_NAME"
            echo ""
            echo "- Experiment: $EXP_NAME"
            echo "- Reference: $REF_FILE"
            echo "- Mapping: $MAPPED"
            echo "- Mean depth: ${MEAN_DEPTH}x"
            echo "- Consensus: ${REF_NAME}.consensus.fasta"
            echo "- Variants: ${REF_NAME}.vcf.gz"
        } > "$REPORT"

        echo "  Done: $SAMPLE_DIR"
    done
done

echo ""
echo "Pipeline finished. Results in: $RESULTS_ROOT"
