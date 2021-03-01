#!/bin/bash
{% if cgems %}
#$ -S /bin/bash
#$ -N GwasQcPipeline
#$ -v CONDA_EXE,CONDA_PREFIX
#$ -cwd
#$ -terse
#$ -j yes
#$ -o gwas_qc_log.$JOB_ID
#$ -q {{ queue }}
#$ -l h_rt={{ time_hr }}:00:00
#$ -l mem_free={{ local_mem_mb }}
{% if local_tasks > 1 %}
#$ -pe by_node {{ local_tasks }}
{% endif %}

export PATH=$CONDA_PREFIX/bin:$(dirname $CONDA_EXE):$PATH
source /etc/profile.d/modules.sh
module load sge
unset module

CLUSTER_JOB_ID=${JOB_ID}
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

CLUSTER_JOB_ID=${SLURM_JOB_ID}
{% endif %}

shopt -s -o errexit pipefail nounset

cgr_get_time() {
    date +'%H:%M:%S %Y-%m-%d'
}

cgr_start_message() {
    printf "################################################################################\n"
    printf "# CGR SUBMIT: Starting GWAS QC Workflow\t" && cgr_get_time
    printf "################################################################################\n"
}

cgr_restart_message() {
    printf "\n\n"
    printf "################################################################################\n"
    printf "# CGR SUBMIT: Re-Starting workflow to finish incomplete tasks\t" && cgr_get_time
    printf "################################################################################\n"
    printf "\n\n"
}

cgr_exit_message() {
    printf "\n\n"
    printf "################################################################################\n"
    if [[ $? != 0 ]]; then
        printf "# CGR SUBMIT: There was an error running the workflow\t" && cgr_get_time
    else
        printf "# CGR SUBMIT: Workflow complete\t" && cgr_get_time
    fi
    printf "################################################################################\n"
}
trap cgr_exit_message EXIT  # capture exit signal and print exit message

# Make sure logs dir exists
cd {{ working_dir }}
[[ -d logs ]] || mkdir -p logs

# Run the workflow
cgr_start_message
run_workflow() {
    {{ python_executable }} -m cgr_gwas_qc snakemake \
        --local-cores {{ local_tasks }} \
        --profile {{ profile }} \
        {{ group_options }}
}
run_workflow

# Check the log and make sure everything completed as expected i.e. "(100%) done"
sleep 10  # in case of filesystem latency
PCT_DONE=$(tail -n 5 gwas_qc_log.$CLUSTER_JOB_ID | sed -nr "s/.*[[:digit:]]+ of [[:digit:]]+ steps \((.*)\%\) done.*/\1/p")

if [[ $PCT_DONE != "" && $PCT_DONE != 100 ]]; then
    # The pipeline ran successfully but ended early so restart
    cgr_restart_message
    run_workflow
fi

exit 0
