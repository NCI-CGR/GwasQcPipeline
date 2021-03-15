QC Report
=========

README
------

.. todo:: Describe the README

- **Inputs**:
- **Output**:
- **Script**:

Report
------

.. todo:: Describe the text report

- **Inputs**:
- **Output**:
- **Script**:

Summary Table
-------------

.. todo:: Describe the summary Excel file

- **Inputs**:
- **Output**:
- **Script**:

Supplemental Tables and Parts
-----------------------------

We create a number of potentially useful tables in order to generate the above reports. These tables and their locations are briefly described below.

Sample QC Summary
+++++++++++++++++

.. todo:: Describe the sample QC table.

- **Inputs**:
- **Output**:
- **Script**:

Population QC Summary
+++++++++++++++++++++

- **Inputs**:
    - ``population_level/{population}/subjects.het``
    - ``population_level/{population}/subjects_unrelated{pi}_maf{maf}_ld{ld}_pruned.eigenvec``
- **Output**: ``population_level/population_qc.csv``
- **Script**: ``cgr_gwas_qc/workflow/scripts/population_qc_table.py``

.. automodule:: cgr_gwas_qc.workflow.scripts.population_qc_table
