#!/bin/bash
# properties = {properties}
#$ -terse
#$ -S /bin/bash
#$ -v CONDA_EXE,CONDA_PREFIX
#$ -cwd
#$ -j yes
set -euo pipefail

export PATH=$CONDA_PREFIX/bin:$(dirname $CONDA_EXE):$PATH

{exec_job}

echo $?
