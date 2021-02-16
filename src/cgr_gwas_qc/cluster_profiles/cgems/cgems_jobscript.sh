#!/bin/bash
# properties = {properties}
#$ -terse
#$ -S /bin/bash
#$ -v CONDA_EXE,CONDA_PREFIX
#$ -cwd
#$ -j yes
export PATH=${CONDA_PREFIX}/bin$(dirname ${CONDA_EXE}):${PATH}

set -euo pipefail

{exec_job}

echo $?
