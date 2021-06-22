Working with Test Suite
=======================

Software testing is a technique where you isolate different parts of a system and test specific functionality.
Testing in bioinformatics is challenging because we often need large amounts of data.
We used a test driven approach and use both unit tests (tests targeting 1 specific piece of functionality) and integration tests (test targeting many steps at once).
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
      ├── data/                                # This is a git submodule pointing to our fake test data at https://github.com/NCI-CGR/GwasQcPipeline-test-data
      ├── models/
      ├── parsers/
      ├── reporting/
      ├── test_config.py
      ├── test_snakemake_job_grouping.py
      ├── testing/
      ├── validators/
      └── workflow/

Notice that the ``test/`` structure closely mimics the ``src/cgr_gwas_qc/`` structure.
pytest_ requires a few naming conventions.

1. Any file that contains tests must be named ``test_*.py``
2. Any test function must be named ``def test_*()``
3. There is a special file ``conftest.py`` that is shared by all of the test files.
   The file ``tests/conftest.py`` is available to all tests.
   While a file ``tests/cli/conftest.py`` would be available to all tests in the ``tests/cli`` folder and subdirectories.

pytest_ also uses a concept called ``fixtures``, where you can define a function that generates some data/state and then use the return value in your tests.
For example,

.. code-block:: python

   import pytest

   @pytest.fixture
   def one():
      return 1

   def test_addition(one):
      assert 2 == one + 1


In this simple example, we create a fixture ``one`` that returns the number ``1``.
Our test function ``test_addition`` has the fixture name as one of its parameters.
pytest_ will replace the fixture name with the return value so you can use it directly as a variable in your test.
I use ``conftest.py`` as a way to store and share fixtures among tests.


Fake Data
---------

Our "fake" test data is stored at https://github.com/NCI-CGR/GwasQcPipeline-test-data.
We have this repository as a git submodule in the test directory.
So if you followed :ref:`dev-setup` where we ``git clone --recursive`` you already have this data at ``tests/data/``.

To run all tests that use these data simply::

   $ pytest


You will notice the use of ``FakeData()`` in our test cases.
This is a helper class for providing the correct paths to our test data.
See the ``cgr_gwas_qc.testing.data`` module for more information.

.. autoclass:: cgr_gwas_qc.testing.data.FakeData


Real Data
---------

We also have small set real data stored on CGEMs/CCAD at ``/DCEG/CGF/Bioinformatics/Production/data/cgr_gwas_qc``.
We do not distribute this data publicly but CGR users can sync this data to other NIH computers/servers using:

.. code-block::

   from cgr_gwas_qc.testing.data import RealData

   RealData(sync=True)

This will rsync the real data to ``<path to cloned version of repository>/.cache/cgr_gwas_qc/test_data``.
If down the road CGR shares this data with other users then they can save the data on their own internal server.
Then users would need to set the environmental variables ``TEST_DATA_USER``, ``TEST_DATA_SERVER``, and ``TEST_DATA_PATH`` to replicate the rsync capabilities.

To run all tests, including those with real data, simply::

   $ pytest --real-data

.. note::
   You will see the use of ``@pytest.mark.real_data`` as a decorator on test functions.
   For example,

   .. code-block:: python

      import pytest

      @pytest.mark.real_data
      def test_addition():
         pass

   This tells pytest_ that this test uses real data and should be skipped unless the ``--real-data`` flag is given.

.. autoclass:: cgr_gwas_qc.testing.data.RealData
