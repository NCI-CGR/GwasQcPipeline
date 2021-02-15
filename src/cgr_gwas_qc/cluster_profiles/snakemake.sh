#!/bin/bash
{% if cgems %}
#$ -S /bin/bash
#$ -N GwasQcPipeline
#$ -V
#$ -cwd
#$ -j yes
#$ -q {{ queue }}
#$ -l h_rt={{ time_hr }}:00:00
source /etc/profile.d/modules.sh;
module load sge
unset module

{% endif %}
{% if biowulf %}
#SBATCH --job-name="GwasQcPipeline"
#SBATCH --partition="{{ queue }}"
#SBATCH --time={{ time_hr }}:00:00
{% endif %}
set -euo pipefail

cd {{ working_dir }}
[[ -d logs ]] || mkdir -p logs

{{ python_executable }} -m cgr_gwas_qc snakemake --profile {{ profile }} {{ group_options }}
