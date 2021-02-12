#!/bin/bash
# properties = {properties}
#$ -terse
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

set -e

{exec_job}

echo $?
