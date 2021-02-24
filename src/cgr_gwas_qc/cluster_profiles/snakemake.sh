#!/bin/bash
{% if cgems %}
#$ -S /bin/bash
#$ -N GwasQcPipeline
#$ -v CONDA_EXE,CONDA_PREFIX
#$ -cwd
#$ -j yes
#$ -o gwas_qc_log.$JOB_ID
#$ -q {{ queue }}
#$ -l h_rt={{ time_hr }}:00:00
#$ -l mem={{ local_mem_mb }}
{% if local_tasks > 1 %}#$ -pe by_node {{ local_tasks }}{% endif %}
export PATH=$CONDA_PREFIX/bin:$(dirname $CONDA_EXE):$PATH
source /etc/profile.d/modules.sh; module load sge; unset module
{% endif %}
{% if biowulf %}
#SBATCH --job-name="GwasQcPipeline"
#SBATCH --partition="{{ queue }}"
#SBATCH --output=gwas_qc_log.%j
#SBATCH --time={{ time_hr }}:00:00
#SBATCH --nodes=1
#SBATCH --ntasks={{ local_tasks }}
#SBATCH --cpus-per-task=1
#SBATCH --mem={{ local_mem_mb }}
{% endif %}

set -euo pipefail

cd {{ working_dir }}
[[ -d logs ]] || mkdir -p logs

{{ python_executable }} -m cgr_gwas_qc snakemake --local-cores {{ local_tasks }} --profile {{ profile }} {{ group_options }}
