#!/usr/bin/env bash
# Usage: bash parse_scores_ei.sh [output_file] [score_dir] [n_jobs]
# Default output: results_ei.csv.gz
# Default input dir: current directory
#
# Parses files matching:
#   score-ei-{dataset}_{ref}_mixRNA_{pp_mRNA}_{fs_mRNA}_RNA_{pp_RNA}_{fs_RNA}
#             _scRNA_{pp_scRNA}_{fs_scRNA}_mixMET_{pp_mixMET}_{fs_mixMET}
#             _MET_{pp_MET}_{fs_MET}_{early_int}_{deconv}.h5

set -euo pipefail

OUTPUT="${1:-results_ei.csv.gz}"
SCORE_DIR="${2:-.}"
JOBS="${3:-4}"



HEADER='"dataset","ref","preprocessing_mixRNA","feature_selection_mixRNA","preprocessing_RNA","feature_selection_RNA","preprocessing_scRNA","feature_selection_scRNA","preprocessing_mixMET","feature_selection_mixMET","preprocessing_MET","feature_selection_MET","early_integration","deconvolution","aid","aid_norm","aitchison","aitchison_norm","jsd","jsd_norm","mae","mae_norm","pearson_col","pearson_col_norm","pearson_row","pearson_row_norm","pearson_tot","pearson_tot_norm","rmse","rmse_norm","score_aggreg","sdid","sdid_norm","spearman_col","spearman_col_norm","spearman_row","spearman_row_norm","spearman_tot","spearman_tot_norm"'
 

# ── process one file ──────────────────────────────────────────────────────────
process_file() {
    local file="$1"
    local METRIC_COLS="aid aid_norm aitchison aitchison_norm jsd jsd_norm mae mae_norm pearson_col pearson_col_norm pearson_row pearson_row_norm pearson_tot pearson_tot_norm rmse rmse_norm score_aggreg sdid sdid_norm spearman_col spearman_col_norm spearman_row spearman_row_norm spearman_tot spearman_tot_norm"

    # echo $PATH >&2 
    # echo $(which h5dump) >&2
    # echo PWD =  $(pwd) >&2
    # echo $(ls -rthl $1 ) >&2
    # echo $(h5dump $file ) >&2


    # echo COLS =  $METRIC_COLS >&2
    # echo $(which gawk)

    local base
    base=$(basename "$file")

    local name="${base#score-ei-}"
    name="${name%.h5}"

    # 14 components:
    # {dataset}_{ref}_mixRNA_{pp_mRNA}_{fs_mRNA}_RNA_{pp_RNA}_{fs_RNA}
    # _scRNA_{pp_scRNA}_{fs_scRNA}_mixMET_{pp_mixMET}_{fs_mixMET}
    # _MET_{pp_MET}_{fs_MET}_{early_int}_{deconv}
    local pattern='^([^_]+)_([^_]+)_mixRNA_([^_]+)_([^_]+)_RNA_([^_]+)_([^_]+)_scRNA_([^_]+)_([^_]+)_mixMET_([^_]+)_([^_]+)_MET_([^_]+)_([^_]+)_([^_]+)_([^_]+)$'

    if [[ ! "$name" =~ $pattern ]]; then
        echo "SKIP (no match): $base" >&2
        return 0
    fi

    local dataset="${BASH_REMATCH[1]}"  ref="${BASH_REMATCH[2]}"
    local pp_mixRNA="${BASH_REMATCH[3]}"  fs_mixRNA="${BASH_REMATCH[4]}"
    local pp_RNA="${BASH_REMATCH[5]}"     fs_RNA="${BASH_REMATCH[6]}"
    local pp_scRNA="${BASH_REMATCH[7]}"   fs_scRNA="${BASH_REMATCH[8]}"
    local pp_mixMET="${BASH_REMATCH[9]}"  fs_mixMET="${BASH_REMATCH[10]}"
    local pp_MET="${BASH_REMATCH[11]}"    fs_MET="${BASH_REMATCH[12]}"
    local early_int="${BASH_REMATCH[13]}" deconv="${BASH_REMATCH[14]}"

    # Dump HDF5 once, extract all scalar metrics with awk
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

    printf '"%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s",%s\n' \
        "$dataset" "$ref" \
        "$pp_mixRNA" "$fs_mixRNA" \
        "$pp_RNA" "$fs_RNA" \
        "$pp_scRNA" "$fs_scRNA" \
        "$pp_mixMET" "$fs_mixMET" \
        "$pp_MET" "$fs_MET" \
        "$early_int" "$deconv" \
        "$metrics"
}

export -f process_file



# ── find files ────────────────────────────────────────────────────────────────
echo "Scanning $SCORE_DIR for score-ei-*.h5 ..." >&2
# mapfile -t FILES < <(find "$SCORE_DIR" -maxdepth 1 -name 'score-ei-*.h5' -type f | sort)
mapfile -t FILES < <(find "$SCORE_DIR" -maxdepth 1 -name 'score-ei-*.h5'  | sort)
TOTAL=${#FILES[@]}
echo "Found $TOTAL files — processing with $JOBS parallel jobs..." >&2

if [[ $TOTAL -eq 0 ]]; then
    echo "No files found. Exiting." >&2
    exit 1
fi

# ── parallel processing → gzip stream ────────────────────────────────────────
{
    echo "$HEADER"
    if command -v parallel &>/dev/null; then
        # GNU parallel: best for 1M+ files, progress bar included
        printf '%s\n' "${FILES[@]}" \
            | parallel -j "$JOBS" --env PATH --bar --line-buffer  process_file {}
            
    else
        # Fallback: xargs -P (no progress bar)
        echo "GNU parallel not found, falling back to xargs -P $JOBS" >&2
        printf '%s\0' "${FILES[@]}" \
            | xargs -0 -P "$JOBS" -I{} bash -c 'process_file "$@"' _ {}
    fi
} | gzip -c > "$OUTPUT"

# process_file ${FILES[1]}

# echo "len=${#FILES[@]}"
# echo "FILES[0]='${FILES[0]}'"
# echo "FILES[1]='${FILES[1]}'"

# # Et ce que printf envoie réellement
# printf '%s\n' "${FILES[1]}" | cat -A

echo "Done → $OUTPUT" >&2