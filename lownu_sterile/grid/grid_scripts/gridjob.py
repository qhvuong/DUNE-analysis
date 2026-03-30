import os

# Define variants and their destination directories
variants = {
    "default":        {"extra_flag": "", "suffix": "freeUt42"},
    "fixedUt42":      {"extra_flag": "--fixedUt42", "suffix": "fixedUt42"},
    "originalTgt":    {"extra_flag": "--originalTgt", "suffix": "originalTgt"},
    "FCsample_stat":  {"extra_flag": "--FCsample \"stat\"", "suffix": "FCstat"},
    "FCsample_flux":  {"extra_flag": "--FCsample \"flux\"", "suffix": "FCflux"},
    "FCsample_det":   {"extra_flag": "--FCsample \"det\"", "suffix": "FCdet"},
    "FCsample_sig":   {"extra_flag": "--FCsample \"sig\"", "suffix": "FCsig"},
}

template = """#!/usr/bin/env bash
export IFDH_CP_UNLINK_ON_ERROR=1
export IFDH_CP_MAXRETRIES=1
export IFDH_DEBUG=0

ifdh_mkdir_p() {{
  local dir=$1
  local force=$2
  if [ `ifdh ls $dir 0 $force | wc -l` -gt 0 ]
  then
      :
  else
      ifdh_mkdir_p `dirname $dir` $force
      ifdh mkdir $dir $force
  fi
}}

FITSPERJOB=$1

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup ifdhc
setup root v6_22_08d -q e20:p392:prof
setup jobsub_client
setup cigetcert

ifdh cp /pnfs/dune/persistent/users/qvuong/FCgrid/FCgrid.tar.gz FCgrid.tar.gz
tar -xzf FCgrid.tar.gz
mv FCgrid/* .

RESULTS_FILE="job_output_${{PROCESS}}.txt"
LOGFILE="job_log_${{PROCESS}}.txt"

if [ ! -f "$RESULTS_FILE" ]; then
  echo -e "universe\\tbf0\\tbf1\\tbf2\\tbf3\\tbfchi2\\ttgtchi2\\tratio" > "$RESULTS_FILE"
fi

for RUN in $(seq 0 $(($FITSPERJOB - 1)) ); do
  POINT=$(($PROCESS * $FITSPERJOB + $RUN))

  {{
    echo "========== POINT $POINT =========="
    echo "Start time: $(date)"
    echo "Command: ./DoTemplateFit --u ${{POINT}} {extra}"
    START=$(date +%s)

    RESULT=$(./DoTemplateFit --u ${{POINT}} {extra} 2>&1)
    STATUS=$?
    END=$(date +%s)

    if [ $STATUS -eq 0 ]; then
      echo "$RESULT" >> "$RESULTS_FILE"
    else
      echo "❌ Fit failed for POINT $POINT with exit code $STATUS" >> "$RESULTS_FILE"
      echo -e "${{POINT}}\\tNaN\\tNaN\\tNaN\\tNaN\\tNaN\\tNaN\\tNaN" >> "$RESULTS_FILE"
    fi

    echo "Exit code: $STATUS"
    echo "End time: $(date)"
    echo "Elapsed time: $((END - START)) seconds"
    echo
  }} >> "$LOGFILE" 2>&1
done

ifdh mkdir -p /pnfs/dune/scratch/users/qvuong/output/FCgrid/{suffix}/outFiles/
ifdh mkdir -p /pnfs/dune/scratch/users/qvuong/output/FCgrid/{suffix}/logs/

ifdh cp -D "$RESULTS_FILE" /pnfs/dune/scratch/users/qvuong/output/FCgrid/{suffix}/outFiles/
ifdh cp -D "$LOGFILE"      /pnfs/dune/scratch/users/qvuong/output/FCgrid/{suffix}/logs/
"""

# Generate one script per variant
for variant_name, variant in variants.items():
    extra_flag = variant["extra_flag"]
    suffix = variant["suffix"]

    script_content = template.format(
        extra=extra_flag,
        suffix=suffix
    )

    filename = f"/pnfs/dune/persistent/users/qvuong/FCgrid/submit_grid_{suffix}.sh"
    with open(filename, "w") as f:
        f.write(script_content)
    os.chmod(filename, 0o755)
    print(f"✅ Created {filename}")
