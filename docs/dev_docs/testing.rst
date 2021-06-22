
Working with Test Suite
=======================

Software testing is a technique where you isolate different parts of a system and test specific functionality.
Testing bioinformatics workflows is challenging because large amounts of data are often needed for certain steps to run.
Here we try to take a test driven approach and use both unit tests (tests targeting 1 specific piece of functionality) and integration tests (test targeting many steps at once).
Where possible we used very small "fake" datasets to ensure the test suite runs quickly.
However, much of the workflow could only be tested using real data.
We do not currently distribute our real test data but we have made the test suite module so that all developers can at least run the fake data tests.
We are using the powerful pytest_ to orchestrate running of our tests.

.. code-block:: bash

   .
   └── tests
      ├── __init__.py
      ├── cli/
      ├── cluster_profiles/
      ├── conftest.py
      ├── data/
      ├── models/
      ├── parsers/
      ├── reporting/
      ├── test_config.py
      ├── test_snakemake_job_grouping.py
      ├── testing/
      ├── validators/
      └── workflow/


Fake Data
---------

To run fake data tests you can simply run::

   $ pytest

Real Data
---------

To run all tests, including those with real data you can simply run::

   $ pytest --real-data
