#!/bin/bash
{% if cgems %}
#$ -S /bin/bash
#$ -N GwasQcPipeline
#$ -v CONDA_EXE,CONDA_PREFIX
#$ -cwd
#$ -j yes
#$ -q {{ queue }}
#$ -l h_rt={{ time_hr }}:00:00
export PATH=$CONDA_PREFIX/bin:$(dirname $CONDA_EXE):$PATH
source /etc/profile.d/modules.sh; module load sge; unset module
{% endif %}
{% if biowulf %}
#SBATCH --job-name="GwasQcPipeline"
#SBATCH --partition="{{ queue }}"
#SBATCH --time={{ time_hr }}:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2gb
{% endif %}

set -euo pipefail

cd {{ working_dir }}
[[ -d logs ]] || mkdir -p logs

{{ python_executable }} -m cgr_gwas_qc snakemake --profile {{ profile }} {{ group_options }}
