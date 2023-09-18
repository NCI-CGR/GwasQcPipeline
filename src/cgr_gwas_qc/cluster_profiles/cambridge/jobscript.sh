#!/bin/bash
# properties = {properties}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --parsable
set -euo pipefail

{exec_job}

echo $?
