#!/usr/bin/env bash
# Usage: bash parse_scores_li.sh [output_file] [score_dir] [n_jobs]
# Requires: gawk, h5dump, GNU parallel (or xargs)
# Resumable: if output file exists, already-processed files are skipped.

set -euo pipefail

OUTPUT="${1:-results_li.csv.gz}"
SCORE_DIR="${2:-.}"
JOBS="${3:-4}"

HEADER='"dataset","ref","preprocessing_mixRNA","feature_selection_mixRNA","preprocessing_RNA","feature_selection_RNA","preprocessing_scRNA","feature_selection_scRNA","deconvolution_rna","preprocessing_mixMET","feature_selection_mixMET","preprocessing_MET","feature_selection_MET","deconvolution_met","late_integration","aid","aid_norm","aitchison","aitchison_norm","jsd","jsd_norm","mae","mae_norm","pearson_col","pearson_col_norm","pearson_row","pearson_row_norm","pearson_tot","pearson_tot_norm","rmse","rmse_norm","score_aggreg","sdid","sdid_norm","spearman_col","spearman_col_norm","spearman_row","spearman_row_norm","spearman_tot","spearman_tot_norm"'

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

# ── resume: find already-done files ──────────────────────────────────────────
TEMP_NEW=$(mktemp /tmp/results_li_new.XXXXXX.csv.gz)
trap 'rm -f "$TEMP_NEW"' EXIT

if [[ -f "$OUTPUT" ]]; then
    echo "Existing output found: $OUTPUT — checking already-processed files..." >&2

    # Reconstruct expected basenames from the CSV columns and match back to filenames.
    # Each row encodes the full pipeline, so we rebuild the stem and check against ALL_FILES.
    # Simpler and more robust: store processed stems in a hash set via awk.
    declare -A DONE
    while IFS= read -r stem; do
        DONE["$stem"]=1
    done < <(zcat "$OUTPUT" | tail -n +2 | gawk -F',' '
    {
        # Strip quotes from all string fields
        for (i=1; i<=NF; i++) gsub(/"/, "", $i)
        # Reconstruct the filename stem from the 15 metadata columns
        # score-li-{1}_{2}_mixRNA_{3}_{4}_RNA_{5}_{6}_scRNA_{7}_{8}_{9}_mixMET_{10}_{11}_MET_{12}_{13}_{14}_{15}
        printf "score-li-%s_%s_mixRNA_%s_%s_RNA_%s_%s_scRNA_%s_%s_%s_mixMET_%s_%s_MET_%s_%s_%s_%s\n",
            $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15
    }')

    ALREADY=${#DONE[@]}
    echo "Already processed: $ALREADY files." >&2

    # Filter to only unprocessed files
    FILES=()
    for f in "${ALL_FILES[@]}"; do
        stem=$(basename "$f" .h5)
        if [[ -z "${DONE[$stem]+x}" ]]; then
            FILES+=("$f")
        fi
    done
else
    FILES=("${ALL_FILES[@]}")
fi

TODO=${#FILES[@]}
echo "To process: $TODO files — using $JOBS parallel jobs..." >&2

if [[ $TODO -eq 0 ]]; then
    echo "Nothing to do. Output is already complete." >&2
    exit 0
fi

# ── process new files ─────────────────────────────────────────────────────────
{
    # Write header only if starting fresh
    if [[ ! -f "$OUTPUT" ]]; then
        echo "$HEADER"
    fi

    if command -v parallel &>/dev/null; then
        printf '%s\n' "${FILES[@]}" \
            | parallel -j "$JOBS" --env PATH --bar --line-buffer process_file {}
    else
        echo "GNU parallel not found, falling back to xargs -P $JOBS" >&2
        printf '%s\0' "${FILES[@]}" \
            | xargs -0 -P "$JOBS" -I{} bash -c 'process_file "$@"' _ {}
    fi
} | gzip -c > "$TEMP_NEW"

# ── merge old + new ───────────────────────────────────────────────────────────
if [[ -f "$OUTPUT" ]]; then
    echo "Merging with existing output..." >&2
    TEMP_MERGE=$(mktemp /tmp/results_li_merge.XXXXXX.csv.gz)
    trap 'rm -f "$TEMP_NEW" "$TEMP_MERGE"' EXIT
    # Concatenate: old (with header) + new (without header, skip first line)
    { zcat "$OUTPUT"; zcat "$TEMP_NEW" | tail -n +2; } | gzip -c > "$TEMP_MERGE"
    mv "$TEMP_MERGE" "$OUTPUT"
else
    mv "$TEMP_NEW" "$OUTPUT"
fi

echo "Done → $OUTPUT" >&2