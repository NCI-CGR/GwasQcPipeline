#!/bin/bash
{% if cgems %}
#$ -S /bin/bash
#$ -N GwasQcPipeline
#$ -v CONDA_EXE,CONDA_PREFIX
#$ -cwd
#$ -terse
#$ -notify
#$ -j yes
#$ -o gwas_qc_log.$JOB_ID
#$ -q {{ queue }}
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

MAX_ATTEMPTS=5
ATTEMPT=0
CLUSTER_KILLED=0

UNICORN="\360\237\246\204"
PENGUIN="\xF0\x9F\x90\xA7"
DNA="\xF0\x9F\xA7\xAC"

cgr_bar(){
    printf "################################################################################\n"
}

cgr_get_time() {
    date +'%H:%M:%S %Y-%m-%d'
}

cgr_start_message() {
    cgr_bar
    printf "# CGR SUBMIT: Starting GWAS QC Workflow\n"
    printf "# User: %s\n" $USER
    printf "# Version: {{ version }}\n"
    printf "# Profile: {{ profile }}\n"
    printf "# Snakefile: {{ snakefile }}\n"
    printf "# Start Time: %s\n# Date: %s\n" $(cgr_get_time)
    cgr_bar
}

cgr_restart_message() {
    printf "\n"
    cgr_bar
    printf "# CGR SUBMIT: Re-Starting workflow to finish incomplete tasks\n"
    printf "# ATTEMPT: %d\n" $1
    printf "# Start Time: %s\n# Date: %s\n" $(cgr_get_time)
    cgr_bar
    printf "\n"
}

cgr_exit_message() {
    exit_code=$?
    printf "\n"
    cgr_bar
    if [[ $exit_code == 0 ]]; then
        printf "# CGR SUBMIT: Workflow complete %b%b%b\n" $UNICORN $PENGUIN $DNA
    elif [[ $CLUSTER_KILLED == 1 ]]; then
        printf "# CGR SUBMIT: The workflow was KILLED by the cluster. Check resource limits and re-run.\n"
        printf "#             For example, try increasing the walltime using the '--time-hr' option\n"
        printf "#             'cgr submit --time-hr 300' \n"
        printf "# Exit Code: %d\n" ${exit_code}
    else
        printf "# CGR SUBMIT: There was a workflow error, check logs and re-run.\n"
        printf "# Exit Code: %d\n" ${exit_code}
    fi
    printf "# End Time: %s\n# Date: %s\n" $(cgr_get_time)
    cgr_bar
    exit $exit_code
}
trap cgr_exit_message EXIT # capture exit signal and print exit message

cgr_cluster_killed() {
    exit_code=$?
    CLUSTER_KILLED=1
    exit $exit_code
}
# capture cluster exit signals
{% if cgems %}
trap cgr_cluster_killed USR1
trap cgr_cluster_killed USR2
{% endif %}
{% if biowulf %}
trap cgr_cluster_killed TERM
{% endif %}

run_workflow() {
    {{ python_executable }} -m cgr_gwas_qc snakemake \
        --local-cores {{ local_tasks }} \
        --profile {{ profile }} \
        {{ added_options }}
}

# Make sure logs dir exists
cd {{ working_dir }}
[[ -d logs ]] || mkdir -p logs

# Run the workflow
cgr_start_message
run_workflow

until [ $ATTEMPT -gt $MAX_ATTEMPTS ]; do
    sleep 10 # in case of filesystem latency

    if [ -e "GwasQcPipeline.complete" ]; then
        # NOTE: This will never happen if submitting subworkflows
        exit 0
    elif grep -q "subworkflow" <<< "{{ added_options }}"; then  # Subworkflow
        # For subworkflows we will check the logs for snakemake to tell us that everything is complete
        final_stage=$(grep -e "^[[:digit:]]\+ of [[:digit:]]\+ steps ([[:digit:]]\+%) done$" gwas_qc_log.${CLUSTER_JOB_ID} | tail -n1)
        if [[ ! -z $final_stage ]]; then
            num_done=$(echo $final_stage | sed -r "s/^([[:digit:]]+) of [[:digit:]]+ steps \([[:digit:]]+%\) done/\1/")
            num_steps=$(echo $final_stage | sed -r "s/^[[:digit:]]+ of ([[:digit:]]+) steps \([[:digit:]]+%\) done/\1/")
            if [[ $num_done == $num_steps ]]; then
                # Comparing number done vs number of steps b/c for large
                # workflow snakemake will report 100% even when not all the
                # steps are complete.
                exit 0
            fi
        fi
    fi

    ((ATTEMPT++))
    cgr_restart_message $ATTEMPT
    run_workflow
done

# Workflow never completed, probably needs intervention
exit 1
