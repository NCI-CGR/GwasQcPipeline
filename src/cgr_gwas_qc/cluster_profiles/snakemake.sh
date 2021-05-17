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

DATE="$(date +'%Y-%m-%d')"
TIME="$(date +'%H:%M:%S')"
MAX_ATTEMPTS=5
ATTEMPT=0
CLUSTER_KILLED=0
LOG="gwas_qc_log.${CLUSTER_JOB_ID}"

################################################################################
# Signal Traps
################################################################################
cgr_start_message() {
  printf "################################################################################\n"
  printf "# CGR SUBMIT: Starting GWAS QC Workflow\n"
  printf "# User: ${USER}\n"
  printf "# Version: {{ version }}\n"
  printf "# Profile: {{ profile }}\n"
  printf "# Snakefile: {{ snakefile }}\n"
  printf "# Cluster Job ID: ${CLUSTER_JOB_ID}\n"
  printf "# Start Date: ${DATE}\n"
  printf "# Start Time: ${TIME}\n"
  printf "################################################################################\n"
}

cgr_restart_message() {
  printf "\n"
  printf "################################################################################\n"
  printf "# CGR SUBMIT: Re-Starting workflow to finish incomplete tasks\n"
  printf "# ATTEMPT: %d\n" $1
  printf "# Start Time: %s\n# Date: %s\n" $(cgr_get_time)
  printf "# End Date: ${DATE}\n"
  printf "# End Time: ${TIME}\n"
  printf "################################################################################\n"
  printf "\n"
}

cgr_exit_success() {
  local unicorn="\360\237\246\204"
  local penguin="\xF0\x9F\x90\xA7"
  local dna="\xF0\x9F\xA7\xAC"

  printf "\n"
  printf "################################################################################\n"
  printf "# CGR SUBMIT: Workflow complete ${unicorn}${penguin}${dna}\n"
  printf "# End Date: ${DATE}\n"
  printf "# End Time: ${TIME}\n"
  printf "################################################################################\n"
}

cgr_exit_locked() {
  printf "\n"
  printf "################################################################################\n"
  printf "# CGR SUBMIT: Your working directory appears to be locked. You probably want to\n"
  printf "#             run: 'cgr snakemake -n --unlock'\n"
  printf "# End Date: ${DATE}\n"
  printf "# End Time: ${TIME}\n"
  printf "################################################################################\n"
}

cgr_exit_failed() {
  local exit_code="$1"
  printf "\n"
  printf "################################################################################\n"
  printf "# CGR SUBMIT: There was a workflow error, check logs and re-run.\n"
  printf "# Exit Code: ${exit_code}\n"
  printf "# End Date: ${DATE}\n"
  printf "# End Time: ${TIME}\n"
  printf "################################################################################\n"
}

cgr_exit_killed() {
  local exit_code="$1"
  printf "\n"
  printf "################################################################################\n"
  printf "# CGR SUBMIT: The workflow was KILLED by the cluster. Check resource limits and re-run.\n"
  printf "#             For example, try increasing the walltime using the '--time-hr' option\n"
  printf "#             'cgr submit --time-hr 300' \n"
  printf "# Exit Code: ${exit_code}\n"
  printf "# End Date: ${DATE}\n"
  printf "# End Time: ${TIME}\n"
  printf "################################################################################\n"
}

cgr_exit_message() {
  exit_code=$?

  if (( exit_code == 0 )); then
    cgr_exit_success
  elif (( exit_code == 42 )); then
    cgr_exit_locked
  elif (( CLUSTER_KILLED == 1 )); then
    cgr_exit_killed "${exit_code}"
  else
    cgr_exit_failed "${exit_code}"
  fi

  exit ${exit_code}
}
trap cgr_exit_message EXIT # capture exit signal and print exit message

cgr_cluster_killed() {
  exit_code=$?
  CLUSTER_KILLED=1
  exit ${exit_code}
}
# capture cluster exit signals
{% if cgems %}
trap cgr_cluster_killed USR1
trap cgr_cluster_killed USR2
{% endif %}
{% if biowulf %}
trap cgr_cluster_killed TERM
{% endif %}

################################################################################
# Running Functions
################################################################################
run_workflow() {
  local exit_status
  {{ python_executable }} -m cgr_gwas_qc snakemake \
    --local-cores {{ local_tasks }} \
    --profile {{ profile }} \
    {{ added_options }} || exit_status=$?  # Don't exit on failure, see why we exited

  if (( exit_status != 0 )) && log_says_locked; then
    exit 42
  fi

}

log_says_complete() {
  local final_stage="$(grep -e "^[[:digit:]]\+ of [[:digit:]]\+ steps ([[:digit:]]\+%) done$" "${LOG}" | tail -n1)"

  if [[ -n "${final_stage}" ]]; then
    local num_done=$(echo "${final_stage}" | sed -r "s/^([[:digit:]]+) of [[:digit:]]+ steps \([[:digit:]]+%\) done/\1/")
    local num_steps=$(echo "${final_stage}" | sed -r "s/^[[:digit:]]+ of ([[:digit:]]+) steps \([[:digit:]]+%\) done/\1/")
    # Comparing number done vs number of steps b/c for large
    # workflow snakemake will report 100% even when not all the
    # steps are complete.
    (( num_done == num_steps ))
    return
  fi

  false
}

log_says_locked() {
  $(grep -q "Error: Directory cannot be locked." "${LOG}")
}

running_subworkflow() {
    $(grep -q "subworkflow" <<< "{{ added_options }}")
}
################################################################################
# Main
################################################################################
main() {
  # Make sure logs dir exists
  cd {{ working_dir }}
  [[ -d logs ]] || mkdir -p logs

  # Run the workflow
  cgr_start_message
  run_workflow

  until (( ATTEMPT > MAX_ATTEMPTS )); do
    sleep 10 # in case of filesystem latency

    if [[ -e "GwasQcPipeline.complete" ]]; then
      # When running the entire workflow look for final output file
      exit 0
    elif running_subworkflow && log_says_complete; then  # Subworkflow
      # For subworkflows check the logs and make sure snakemake says everything is complete
      exit 0
    fi

    ATTEMPT=$(( ATTEMPT + 1 ))
    cgr_restart_message $ATTEMPT
    run_workflow
  done

  # Workflow never completed, probably needs intervention
  exit 1
}

main "$@"
