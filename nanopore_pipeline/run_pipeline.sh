#!/usr/bin/env bash
#
# Oxford Nanopore (MinION) reference-mapping pipeline
#
# For every row in samplesheet.csv:
#   barcode_id, reference_fasta, sample_name
#
# the pipeline maps that barcode's reads against the matching reference
# and writes ALL outputs under results/<sample_name>/, already named
# after <sample_name> instead of the raw barcode ID — so no manual
# renaming step is needed afterwards.
#
# Usage:
#   ./run_pipeline.sh [samplesheet.csv] [config.yaml]
#
# Requires: minimap2, samtools, bcftools (or medaka if variant_caller=medaka)
#           NanoFilt (optional, only used if min_read_length/quality > 0)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLESHEET="${1:-$SCRIPT_DIR/samplesheet.csv}"
CONFIG="${2:-$SCRIPT_DIR/config.yaml}"

# ---- minimal flat YAML reader (key: value, no nesting) -------------------
read_cfg() {
    local key="$1" default="$2"
    local val
    val=$(grep -E "^${key}:" "$CONFIG" | head -n1 | sed -E "s/^${key}:[[:space:]]*//; s/[[:space:]]*(#.*)?$//")
    if [[ -z "$val" ]]; then echo "$default"; else echo "$val"; fi
}

FASTQ_DIR="$SCRIPT_DIR/$(read_cfg fastq_dir data/raw/fastq_pass)"
REF_DIR="$SCRIPT_DIR/$(read_cfg reference_dir references)"
RESULTS_DIR="$SCRIPT_DIR/$(read_cfg results_dir results)"
PRESET="$(read_cfg minimap2_preset map-ont)"
THREADS="$(read_cfg threads 4)"
MIN_LEN="$(read_cfg min_read_length 0)"
MIN_QUAL="$(read_cfg min_read_quality 0)"
VARIANT_CALLER="$(read_cfg variant_caller bcftools)"

echo "==================================================================="
echo " Nanopore reference-mapping pipeline"
echo "   samplesheet : $SAMPLESHEET"
echo "   fastq_dir   : $FASTQ_DIR"
echo "   reference   : $REF_DIR"
echo "   results_dir : $RESULTS_DIR"
echo "   preset      : $PRESET   threads: $THREADS"
echo "==================================================================="

mkdir -p "$RESULTS_DIR"

# Skip header line, ignore blank lines / lines starting with '#'
tail -n +2 "$SAMPLESHEET" | grep -v '^[[:space:]]*$' | grep -v '^[[:space:]]*#' | \
while IFS=',' read -r BARCODE_ID REF_FASTA SAMPLE_NAME; do
    BARCODE_ID="$(echo "$BARCODE_ID" | xargs)"
    REF_FASTA="$(echo "$REF_FASTA" | xargs)"
    SAMPLE_NAME="$(echo "$SAMPLE_NAME" | xargs)"

    SAMPLE_DIR="$RESULTS_DIR/$SAMPLE_NAME"
    REF_PATH="$REF_DIR/$REF_FASTA"
    BARCODE_DIR="$FASTQ_DIR/$BARCODE_ID"

    echo ""
    echo "------------------------------------------------------------"
    echo "Sample: $SAMPLE_NAME  (barcode: $BARCODE_ID, ref: $REF_FASTA)"
    echo "------------------------------------------------------------"

    if [[ ! -d "$BARCODE_DIR" ]]; then
        echo "  [SKIP] fastq folder not found: $BARCODE_DIR" >&2
        continue
    fi
    if [[ ! -f "$REF_PATH" ]]; then
        echo "  [SKIP] reference not found: $REF_PATH" >&2
        continue
    fi

    mkdir -p "$SAMPLE_DIR"

    # 1. Concatenate all fastq(.gz) files for this barcode
    MERGED_FASTQ="$SAMPLE_DIR/${SAMPLE_NAME}.fastq.gz"
    if [[ ! -s "$MERGED_FASTQ" ]]; then
        echo "  [1/5] Merging reads -> $(basename "$MERGED_FASTQ")"
        shopt -s nullglob
        FILES=("$BARCODE_DIR"/*.fastq.gz "$BARCODE_DIR"/*.fastq)
        shopt -u nullglob
        if [[ ${#FILES[@]} -eq 0 ]]; then
            echo "  [SKIP] no fastq files in $BARCODE_DIR" >&2
            continue
        fi
        : > "${MERGED_FASTQ%.gz}.tmp"
        for f in "${FILES[@]}"; do
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
        FILTERED_FASTQ="$SAMPLE_DIR/${SAMPLE_NAME}.filtered.fastq.gz"
        echo "  [2/5] QC filtering (length>=$MIN_LEN, quality>=$MIN_QUAL) -> $(basename "$FILTERED_FASTQ")"
        zcat "$MERGED_FASTQ" | NanoFilt -l "$MIN_LEN" -q "$MIN_QUAL" | gzip > "$FILTERED_FASTQ"
        READS_FOR_MAPPING="$FILTERED_FASTQ"
    else
        echo "  [2/5] QC filtering skipped"
    fi

    # 3. Map to reference with minimap2, sort with samtools
    SORTED_BAM="$SAMPLE_DIR/${SAMPLE_NAME}.sorted.bam"
    echo "  [3/5] Mapping to $REF_FASTA -> $(basename "$SORTED_BAM")"
    minimap2 -ax "$PRESET" -t "$THREADS" "$REF_PATH" "$READS_FOR_MAPPING" \
        | samtools sort -@ "$THREADS" -o "$SORTED_BAM" -
    samtools index "$SORTED_BAM"
    samtools flagstat "$SORTED_BAM" > "$SAMPLE_DIR/${SAMPLE_NAME}.flagstat.txt"
    samtools depth -a "$SORTED_BAM" > "$SAMPLE_DIR/${SAMPLE_NAME}.depth.txt"

    # 4. Consensus sequence
    CONSENSUS_FASTA="$SAMPLE_DIR/${SAMPLE_NAME}.consensus.fasta"
    echo "  [4/5] Building consensus -> $(basename "$CONSENSUS_FASTA")"
    samtools consensus -a -f fasta "$SORTED_BAM" > "$CONSENSUS_FASTA"
    sed -i "1s/.*/>${SAMPLE_NAME}/" "$CONSENSUS_FASTA"

    # 5. Variant calling
    VCF="$SAMPLE_DIR/${SAMPLE_NAME}.vcf.gz"
    echo "  [5/5] Calling variants ($VARIANT_CALLER) -> $(basename "$VCF")"
    if [[ "$VARIANT_CALLER" == "medaka" ]] && command -v medaka_haploid_variant >/dev/null 2>&1; then
        medaka_haploid_variant -i "$READS_FOR_MAPPING" -r "$REF_PATH" -o "$SAMPLE_DIR/medaka"
        cp "$SAMPLE_DIR/medaka/medaka.annotated.vcf" "$SAMPLE_DIR/${SAMPLE_NAME}.vcf" 2>/dev/null || true
        bgzip -f "$SAMPLE_DIR/${SAMPLE_NAME}.vcf"
    else
        bcftools mpileup -f "$REF_PATH" "$SORTED_BAM" 2>/dev/null \
            | bcftools call -mv -Oz -o "$VCF"
        bcftools index -f "$VCF"
    fi

    # Per-sample summary report
    REPORT="$SAMPLE_DIR/${SAMPLE_NAME}_report.md"
    MAPPED=$(grep "mapped (" "$SAMPLE_DIR/${SAMPLE_NAME}.flagstat.txt" | head -1)
    MEAN_DEPTH=$(awk '{sum+=$3; n++} END {if (n>0) printf "%.2f", sum/n; else print "0"}' "$SAMPLE_DIR/${SAMPLE_NAME}.depth.txt")
    {
        echo "# Report: $SAMPLE_NAME"
        echo ""
        echo "- Barcode: $BARCODE_ID"
        echo "- Reference: $REF_FASTA"
        echo "- Mapping: $MAPPED"
        echo "- Mean depth: ${MEAN_DEPTH}x"
        echo "- Consensus: ${SAMPLE_NAME}.consensus.fasta"
        echo "- Variants: ${SAMPLE_NAME}.vcf.gz"
    } > "$REPORT"

    echo "  Done: $SAMPLE_DIR"
done

echo ""
echo "Pipeline finished. Results in: $RESULTS_DIR"
