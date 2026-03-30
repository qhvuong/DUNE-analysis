#!/bin/bash
set -euo pipefail

# if [ "$#" -ne 5 ]; then
#   echo "Usage: $0 ue42 um42 ut42 dm2 nPoints" >&2
#   exit 1
# fi

# Expect 5 args
if [ "$#" -ne 5 ]; then
  echo "Usage: $0 ue42 um42 ut42 dm2 nPoints" >&2
  exit 1
fi

UE42=$1; UM42=$2; UT42=$3; DM2=$4; NPOINTS=$5

# Ensure binary exists
if ! command -v ./DoTemplateFit >/dev/null 2>&1; then
  echo "Error: ./DoTemplateFit not found or not executable." >&2
  exit 2
fi

# Avoid oversubscription if DoTemplateFit uses OpenMP/threads
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}

# Kill children on Ctrl-C
trap 'echo "Stopping..."; pkill -P $$ || true; wait || true; exit 130' INT TERM

# Timestamped logs (optional)
ts=$(date +%Y%m%d_%H%M%S)

# Launch in background
./DoTemplateFit --scanMode ue42_um42 --tgtPar "$UE42" "$UM42" "$UT42" "$DM2" \
  --nPoints "$NPOINTS" --job 0 --njobs 1 > "logs/log_ue_um_${ts}.txt" 2>&1 &

./DoTemplateFit --scanMode ue42_dm2  --tgtPar "$UE42" "$UM42" "$UT42" "$DM2" \
  --nPoints "$NPOINTS" --job 0 --njobs 1 > "logs/log_ue_dm_${ts}.txt" 2>&1 &

./DoTemplateFit --scanMode um42_dm2  --tgtPar "$UE42" "$UM42" "$UT42" "$DM2" \
  --nPoints "$NPOINTS" --job 0 --njobs 1 > "logs/log_um_dm_${ts}.txt" 2>&1 &

./DoTemplateFit --scanMode ut42_dm2  --tgtPar "$UE42" "$UM42" "$UT42" "$DM2" \
  --nPoints "$NPOINTS" --job 0 --njobs 1 > "logs/log_ut_dm_${ts}.txt" 2>&1 &

wait
echo "All scans done: tgtPar=($UE42 $UM42 $UT42 $DM2), nPoints=$NPOINTS"
