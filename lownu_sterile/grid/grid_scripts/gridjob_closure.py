import os

# Define your two sets of tgtPar configs
configs = {
    "physicalBoundary": [
        {"name": "smallDm2", "tgtPar": [0.35, 0.45, 0.2, 0.5]},
        {"name": "medDm2",   "tgtPar": [0.35, 0.45, 0.2, 5.]},
        {"name": "largeDm2", "tgtPar": [0.35, 0.45, 0.2, 90.]},
    ],
    "Asimov": [
        {"name": "smallDm2", "tgtPar": [0.05, 0.05, 0., 0.5]},
        {"name": "medDm2",   "tgtPar": [0.03, 0.01, 0., 5.]},
        {"name": "largeDm2", "tgtPar": [0.04, 0.01, 0., 90.]},
    ],
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
    echo "Command: ./DoTemplateFit --u ${{POINT}} --tgtPar {tgtPar}"
    START=$(date +%s)

    RESULT=$(./DoTemplateFit --u ${{POINT}} --tgtPar {tgtPar} 2>&1)
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

ifdh cp -D "$RESULTS_FILE" /pnfs/dune/scratch/users/qvuong/output/FCgrid/{group}/{name}/outFiles/
ifdh cp -D "$LOGFILE"    /pnfs/dune/scratch/users/qvuong/output/FCgrid/{group}/{name}/logs/
"""

# Generate one script per config per group
for group, setups in configs.items():
    for config in setups:
        tgtPar_str = " ".join(str(x) for x in config["tgtPar"])
        script_content = template.format(
            tgtPar=tgtPar_str,
            group=group,
            name=config["name"]
        )
        filename = f"/pnfs/dune/persistent/users/qvuong/FCgrid/submit_grid_{group}_{config['name']}.sh"
        with open(filename, "w") as f:
            f.write(script_content)
        os.chmod(filename, 0o755)
        print(f"✅ Created {filename}")