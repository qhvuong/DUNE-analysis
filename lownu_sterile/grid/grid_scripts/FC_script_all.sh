#! /usr/bin/env bash

export IFDH_CP_UNLINK_ON_ERROR=1
export IFDH_CP_MAXRETRIES=1
export IFDH_DEBUG=0

# mkdir -p workalike for ifdh cp
ifdh_mkdir_p() {
  local dir=$1
  local force=$2
  if [ `ifdh ls $dir 0 $force | wc -l` -gt 0 ] 
  then
      : # we're done
  else
      ifdh_mkdir_p `dirname $dir` $force
      ifdh mkdir $dir $force
  fi
}

#N=$1
FITSPERJOB=$1

# Set up software
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup ifdhc
setup root v6_22_08d -q e20:p392:prof
setup jobsub_client
setup cigetcert

# Copy the tarball with all of the inputs to the grid node
ifdh cp /pnfs/dune/persistent/users/qvuong/NuCut/NuCut.tar.gz NuCut.tar.gz
tar -xzf NuCut.tar.gz
mv NuCut/* .

#FIRSTN=$((${FIRST}+${PROCESS}))

# for loop
for RUN in $(seq 0 $(($FITSPERJOB-1)) )
do
  #POINT=$(($FIRSTN*$FITSPERJOB + $RUN))
  POINT=$(($PROCESS*$FITSPERJOB + $RUN))

  # Determine the oscillation parameter bin from one integer PROCESS which will be incremented in each job
  n=$(python GetN.py $POINT$ 2>/dev/null | tail -1)

  # run it
  echo "./DoTemplateFit --n ${n}"
  ./DoTemplateFit --n ${n}

  mv output.txt output_${POINT}.txt
  ifdh cp -D output_${POINT}.txt /pnfs/dune/scratch/users/qvuong/output/FC

done

# Make the output a text file with the fit results
#ifdh cp -D output_file_${PROCESS}.txt /pnfs/dune/scratch/users/qvuong/output/LvsE/nCov/dm2_1/

