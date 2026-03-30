#!/bin/bash
set -euo pipefail

usage() {
  echo "Usage: $0 ue42 um42 ut42 dm2 nPoints [noNuE] [fixedL]"
  echo "  noNuE: optional boolean (1|true|yes|on to enable)"
  echo "  fixedL: optional boolean (1|true|yes|on to enable)"
  exit 1
}

# Help
[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && usage

# Expect 5 or 6 args
if [ "$#" -lt 5 ] || [ "$#" -gt 7 ]; then
  usage
fi

UE42=$1; UM42=$2; UT42=$3; DM2=$4; NPOINTS=$5
NO_NUE_INPUT="${6:-}"
FIXED_L_INPUT="${7:-}"

# Normalize boolean
shopt -s nocasematch
case "$NO_NUE_INPUT" in
  1|true|yes|on) NO_NUE_FLAG="--noNuE" ;;
  ""|0|false|no|off) NO_NUE_FLAG="" ;;
  *) echo "Error: invalid noNuE value '$NO_NUE_INPUT'"; usage ;;
esac
shopt -u nocasematch

shopt -s nocasematch
case "$FIXED_L_INPUT" in
  1|true|yes|on) FIXED_L_FLAG="--fixedL" ;;
  ""|0|false|no|off) FIXED_L_FLAG="" ;;
  *) echo "Error: invalid fixedL value '$FIXED_L_INPUT'"; usage ;;
esac
shopt -u nocasematch



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

# Logs
mkdir -p logs
ts=$(date +%Y%m%d_%H%M%S)

# Launch in background
./DoTemplateFit --scanMode ue42_um42 --tgtPar "$UE42" "$UM42" "$UT42" "$DM2" \
  --nPoints "$NPOINTS" ${NO_NUE_FLAG} ${FIXED_L_FLAG} --job 0 --njobs 1 > "logs/log_ue_um_${ts}.txt" 2>&1 &

./DoTemplateFit --scanMode ue42_dm2  --tgtPar "$UE42" "$UM42" "$UT42" "$DM2" \
  --nPoints "$NPOINTS" ${NO_NUE_FLAG} ${FIXED_L_FLAG} --job 0 --njobs 1 > "logs/log_ue_dm_${ts}.txt" 2>&1 &

./DoTemplateFit --scanMode um42_dm2  --tgtPar "$UE42" "$UM42" "$UT42" "$DM2" \
  --nPoints "$NPOINTS" ${NO_NUE_FLAG} ${FIXED_L_FLAG} --job 0 --njobs 1 > "logs/log_um_dm_${ts}.txt" 2>&1 &

./DoTemplateFit --scanMode ut42_dm2  --tgtPar "$UE42" "$UM42" "$UT42" "$DM2" \
  --nPoints "$NPOINTS" ${NO_NUE_FLAG} ${FIXED_L_FLAG} --job 0 --njobs 1 > "logs/log_ut_dm_${ts}.txt" 2>&1 &

wait
echo "All scans done: tgtPar=($UE42 $UM42 $UT42 $DM2), nPoints=$NPOINTS, noNuE=${NO_NUE_FLAG:+on}${NO_NUE_FLAG:-off}, fixedL=${FIXED_L_FLAG:+on}${FIXED_L_FLAG:-off}"
