#!/usr/bin/env bash
# Usage: bash parse_scores_li.sh [output_file] [score_dir] [n_jobs]
# Requires: gawk, h5dump, GNU parallel (or xargs)
# Resumable: writes to plain CSV throughout; gzips only at the end.
# On resume, validates and trims any corrupted trailing lines before continuing.

set -euo pipefail

OUTPUT="${1:-results_li.csv.gz}"
SCORE_DIR="${2:-.}"
JOBS="${3:-4}"

# Plain CSV used during processing — kept until gzip succeeds
CSV_WORK="${OUTPUT%.gz}"   # e.g. results_li.csv

HEADER='"dataset","ref","preprocessing_mixRNA","feature_selection_mixRNA","preprocessing_RNA","feature_selection_RNA","preprocessing_scRNA","feature_selection_scRNA","deconvolution_rna","preprocessing_mixMET","feature_selection_mixMET","preprocessing_MET","feature_selection_MET","deconvolution_met","late_integration","aid","aid_norm","aitchison","aitchison_norm","jsd","jsd_norm","mae","mae_norm","pearson_col","pearson_col_norm","pearson_row","pearson_row_norm","pearson_tot","pearson_tot_norm","rmse","rmse_norm","score_aggreg","sdid","sdid_norm","spearman_col","spearman_col_norm","spearman_row","spearman_row_norm","spearman_tot","spearman_tot_norm"'

# Expected number of comma-separated columns in a valid data row
EXPECTED_COLS=40

# ── process_file ──────────────────────────────────────────────────────────────
process_file() {
    local file="$1"
    local base
    base=$(basename "$file")

    local name="${base#score-li-}"
    name="${name%.h5}"

    local pattern='^([^_]+)_([^_]+)_mixRNA_([^_]+)_([^_]+)_RNA_([^_]+)_([^_]+)_scRNA_([^_]+)_([^_]+)_([^_]+)_mixMET_([^_]+)_([^_]+)_MET_([^_]+)_([^_]+)_([^_]+)_([^_]+)$'

    if [[ ! "$name" =~ $pattern ]]; then
        echo "SKIP (no match): $base" >&2
        return 0
    fi

    local dataset="${BASH_REMATCH[1]}"    ref="${BASH_REMATCH[2]}"
    local pp_mixRNA="${BASH_REMATCH[3]}"  fs_mixRNA="${BASH_REMATCH[4]}"
    local pp_RNA="${BASH_REMATCH[5]}"     fs_RNA="${BASH_REMATCH[6]}"
    local pp_scRNA="${BASH_REMATCH[7]}"   fs_scRNA="${BASH_REMATCH[8]}"
    local deconv_rna="${BASH_REMATCH[9]}"
    local pp_mixMET="${BASH_REMATCH[10]}" fs_mixMET="${BASH_REMATCH[11]}"
    local pp_MET="${BASH_REMATCH[12]}"    fs_MET="${BASH_REMATCH[13]}"
    local deconv_met="${BASH_REMATCH[14]}" late_int="${BASH_REMATCH[15]}"

    local METRIC_COLS="aid aid_norm aitchison aitchison_norm jsd jsd_norm mae mae_norm pearson_col pearson_col_norm pearson_row pearson_row_norm pearson_tot pearson_tot_norm rmse rmse_norm score_aggreg sdid sdid_norm spearman_col spearman_col_norm spearman_row spearman_row_norm spearman_tot spearman_tot_norm"

    local metrics
    metrics=$(h5dump -m "%.15g" "$file" 2>/dev/null | gawk -v cols="$METRIC_COLS" '
    BEGIN {
        split(cols, col_order, " ")
        current_group = ""
        in_data = 0
        reading = 0
    }
    /GROUP / {
        match($0, /GROUP "([a-z_]+)"/, arr)
        if (arr[1] != "") {
            current_group = arr[1]
            in_data = 0
            reading = 0
        }
    }
    /DATASET "data"/ { in_data = 1 }
    in_data && /DATA \{/ { reading = 1 }
    reading && /\(0\):/ {
        match($0, /\(0\): *(.+)$/, arr)
        val = arr[1]
        gsub(/[[:space:]]+$/, "", val)
        if (val == "-nan" || val == "nan" || val == "NaN" || val == "-NaN" || val == "inf" || val == "-inf") val = "NA"
        metrics[current_group] = val
        reading = 0
        in_data = 0
    }
    END {
        n = length(col_order)
        for (i = 1; i <= n; i++) {
            col = col_order[i]
            printf "%s", (col in metrics ? metrics[col] : "NA")
            if (i < n) printf ","
        }
        printf "\n"
    }')

    if [[ -z "$metrics" ]]; then
        echo "SKIP (h5dump failed or empty): $base" >&2
        return 0
    fi

    printf '"%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s",%s\n' \
        "$dataset" "$ref" \
        "$pp_mixRNA" "$fs_mixRNA" \
        "$pp_RNA" "$fs_RNA" \
        "$pp_scRNA" "$fs_scRNA" \
        "$deconv_rna" \
        "$pp_mixMET" "$fs_mixMET" \
        "$pp_MET" "$fs_MET" \
        "$deconv_met" "$late_int" \
        "$metrics"
}

export -f process_file

# ── find all files ────────────────────────────────────────────────────────────
echo "Scanning $SCORE_DIR for score-li-*.h5 ..." >&2
mapfile -t ALL_FILES < <(find "$SCORE_DIR" -maxdepth 1 -name 'score-li-*.h5' | sort)
TOTAL=${#ALL_FILES[@]}
echo "Found $TOTAL files total." >&2

if [[ $TOTAL -eq 0 ]]; then
    echo "No files found. Exiting." >&2
    exit 1
fi

# ── resume: validate and load existing CSV ────────────────────────────────────
declare -A DONE

if [[ -f "$CSV_WORK" ]]; then
    echo "Existing CSV found: $CSV_WORK — validating trailing lines..." >&2

    # Count total lines before any trimming
    TOTAL_LINES=$(wc -l < "$CSV_WORK")

    # Scan from the bottom up and remove any incomplete/corrupt lines.
    # A valid data line must have exactly EXPECTED_COLS comma-separated fields.
    # The header line (line 1) is always kept.
    LINES_REMOVED=0
    while [[ "$(wc -l < "$CSV_WORK")" -gt 1 ]]; do
        last_line=$(tail -1 "$CSV_WORK")
        col_count=$(echo "$last_line" | gawk -F',' '{print NF}')

        if [[ "$col_count" -eq "$EXPECTED_COLS" ]]; then
            break   # Last line looks valid — stop trimming
        else
            echo "  Removing corrupt/incomplete line ($col_count cols): ${last_line:0:80}..." >&2
            # Remove the last line in-place
            sed -i '$ d' "$CSV_WORK"
            (( LINES_REMOVED++ ))
        fi
    done

    if [[ $LINES_REMOVED -gt 0 ]]; then
        echo "  Trimmed $LINES_REMOVED incomplete line(s) from the end of the CSV." >&2
    else
        echo "  Trailing lines look valid — no trimming needed." >&2
    fi

    # Build DONE set from the now-validated CSV (skip header on line 1)
    while IFS= read -r stem; do
        DONE["$stem"]=1
    done < <(tail -n +2 "$CSV_WORK" | gawk -F',' '
    {
        for (i=1; i<=NF; i++) gsub(/"/, "", $i)
        printf "score-li-%s_%s_mixRNA_%s_%s_RNA_%s_%s_scRNA_%s_%s_%s_mixMET_%s_%s_MET_%s_%s_%s_%s\n",
            $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15
    }')

    ALREADY=${#DONE[@]}
    echo "Already processed: $ALREADY valid rows — resuming." >&2

elif [[ -f "$OUTPUT" ]]; then
    # No plain CSV but a .gz exists — decompress it to resume from
    echo "No plain CSV found but $OUTPUT exists — decompressing to resume..." >&2
    if ! gzip -t "$OUTPUT" 2>/dev/null; then
        echo "ERROR: $OUTPUT is corrupted and cannot be resumed. Remove it and restart." >&2
        exit 1
    fi
    gzip -dk "$OUTPUT"   # produces $CSV_WORK, keeps $OUTPUT intact
    # Recurse: re-run the validation block above by restarting the script
    exec "$0" "$@"
fi

# ── filter to unprocessed files ───────────────────────────────────────────────
FILES=()
for f in "${ALL_FILES[@]}"; do
    stem=$(basename "$f" .h5)
    if [[ -z "${DONE[$stem]+x}" ]]; then
        FILES+=("$f")
    fi
done

TODO=${#FILES[@]}
echo "To process: $TODO files — using $JOBS parallel jobs..." >&2

if [[ $TODO -eq 0 ]]; then
    echo "Nothing to do. Output is already complete." >&2
    # If only CSV exists (no gz yet), still gzip it
    if [[ -f "$CSV_WORK" && ! -f "$OUTPUT" ]]; then
        echo "Gzipping $CSV_WORK → $OUTPUT ..." >&2
        gzip -k "$CSV_WORK"
        rm -f "$CSV_WORK"
        echo "Done → $OUTPUT" >&2
    fi
    exit 0
fi

# ── write header if starting fresh ───────────────────────────────────────────
if [[ ! -f "$CSV_WORK" ]]; then
    echo "$HEADER" > "$CSV_WORK"
fi

# ── process new files, append directly to CSV ────────────────────────────────
if command -v parallel &>/dev/null; then
    printf '%s\n' "${FILES[@]}" \
        | parallel -j "$JOBS" --env PATH --bar --line-buffer process_file {} \
        >> "$CSV_WORK"
else
    echo "GNU parallel not found, falling back to xargs -P $JOBS" >&2
    printf '%s\0' "${FILES[@]}" \
        | xargs -0 -P "$JOBS" -I{} bash -c 'process_file "$@"' _ {} \
        >> "$CSV_WORK"
fi

# ── gzip at the end, keep CSV until gz is confirmed valid ────────────────────
echo "Processing complete — gzipping $CSV_WORK → $OUTPUT ..." >&2
gzip -k "$CSV_WORK"           # -k keeps the CSV alongside the .gz

# Verify the gz is valid before removing the CSV
if gzip -t "$OUTPUT" 2>/dev/null; then
    echo "Gzip verified OK — removing plain CSV." >&2
    rm -f "$CSV_WORK"
    echo "Done → $OUTPUT" >&2
else
    echo "ERROR: gzip verification failed — keeping $CSV_WORK as backup." >&2
    rm -f "$OUTPUT"
    exit 1
fi