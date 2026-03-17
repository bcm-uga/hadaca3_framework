#!/usr/bin/env bash
# Usage: bash parse_scores_li.sh [output_file] [score_dir] [n_jobs]
# Requires: gawk, h5dump, GNU parallel (or xargs)
#
# Parses files matching:
#   score-li-{dataset}_{ref}_mixRNA_{pp_mRNA}_{fs_mRNA}_RNA_{pp_RNA}_{fs_RNA}
#             _scRNA_{pp_scRNA}_{fs_scRNA}_{deconv_rna}
#             _mixMET_{pp_mixMET}_{fs_mixMET}_MET_{pp_MET}_{fs_MET}_{deconv_met}
#             _{late_integration}.h5

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

    # 15 components:
    # {dataset}_{ref}_mixRNA_{pp_mRNA}_{fs_mRNA}_RNA_{pp_RNA}_{fs_RNA}
    # _scRNA_{pp_scRNA}_{fs_scRNA}_{deconv_rna}
    # _mixMET_{pp_mixMET}_{fs_mixMET}_MET_{pp_MET}_{fs_MET}_{deconv_met}_{late_int}
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
    metrics=$(h5dump "$file" 2>/dev/null | gawk -v cols="$METRIC_COLS" '
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

# ── find files ────────────────────────────────────────────────────────────────
echo "Scanning $SCORE_DIR for score-li-*.h5 ..." >&2
mapfile -t FILES < <(find "$SCORE_DIR" -maxdepth 1 -name 'score-li-*.h5'  | sort)
TOTAL=${#FILES[@]}
echo "Found $TOTAL files — processing with $JOBS parallel jobs..." >&2

if [[ $TOTAL -eq 0 ]]; then
    echo "No files found. Exiting." >&2
    exit 1
fi

{
    echo "$HEADER"
    if command -v parallel &>/dev/null; then
        printf '%s\n' "${FILES[@]}" \
            | parallel -j "$JOBS" --env PATH --bar --line-buffer process_file {}
    else
        echo "GNU parallel not found, falling back to xargs -P $JOBS" >&2
        printf '%s\0' "${FILES[@]}" \
            | xargs -0 -P "$JOBS" -I{} bash -c 'process_file "$@"' _ {}
    fi
} | gzip -c > "$OUTPUT"

echo "Done → $OUTPUT" >&2