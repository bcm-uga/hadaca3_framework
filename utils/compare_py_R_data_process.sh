#!/bin/bash

set -e

echo "=== Running compatibility check between data_processing.R and data_processing.py ==="

# dir_path="~/projects/hadaca3_framework/"
dir_path="/home/nicolashomberg/projects/hadaca3_framework/"


R_SCRIPT=$dir_path"utils/data_processing.R"
PY_SCRIPT=$dir_path"utils/data_processing.py"

# INPUT1=$dir_path"data/mixes1_insilicodirichletCopule_pdac.h5"
# INPUT2=$dir_path"data/ref.h5"

INPUT1=$dir_path"data/mixes1_insilicodirichletCopule_pdac.h5"
INPUT2=$dir_path"data/ref.h5"

TMP_DIR="tmp"
rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"



# Output files
R_OUT1="$TMP_DIR/r_output_mixes1.h5"
R_OUT2="$TMP_DIR/r_output_ref.h5"
PY_OUT1="$TMP_DIR/py_output_mixes1.h5"
PY_OUT2="$TMP_DIR/py_output_ref.h5"


echo "--- Running Python script on $INPUT1 and $INPUT2"
python3 -c "from data_processing import read_hdf5, write_global_hdf5; \
            d1 = read_hdf5('$INPUT1'); write_global_hdf5('$PY_OUT1', d1); \
            d2 = read_hdf5('$INPUT2'); write_global_hdf5('$PY_OUT2', d2);"

echo "--- Running R script on $INPUT1 and $INPUT2"
Rscript -e "source('$R_SCRIPT'); \
            data1 <- read_hdf5('$INPUT1'); write_global_hdf5('$R_OUT1', data1); \
            data2 <- read_hdf5('$INPUT2'); write_global_hdf5('$R_OUT2', data2);"



echo "--- Comparing outputs using h5diff"

diffs=0
for file in ref mixes1; do
  R_FILE="$TMP_DIR/r_output_${file}.h5"
  PY_FILE="$TMP_DIR/py_output_${file}.h5"

  echo "\n\nComparing $R_FILE vs $PY_FILE"
  h5diff -v "$R_FILE" "$PY_FILE"
  # if h5diff "$R_FILE" "$PY_FILE"; then
  #   echo "✅ $file: No differences"
  # else
  #   echo "❌ $file: Differences found"
  #   diffs=$((diffs + 1))
  # fi
done

# if [ "$diffs" -eq 0 ]; then
#   echo "✅ All HDF5 files are compatible!"
# else
#   echo "❌ $diffs mismatches found in HDF5 compatibility check."
#   exit 1
# fi
