Documentation
=============

|cgr| uses the Sphinx_ python documentation generator.
Documentation is written in reStructuredText_, a markup language similar to Markdown but with more advanced features.

Writing Documentation
---------------------

All documentation is found in the ``docs/`` folder.

.. code-block::

   # GwasQcPipeline folder cloned from https://github.com/NCI-CGR/GwasQcPipeline
   .github/
   docs/ # This folder contains all documentation.
   src/
   tests/
   README.md  # This is displayed on the main GitHub page.
   ...


.. todo::

   Add note about docs folder structure.

.. note::

   If you add a new ``rst`` page, you must reference this page from ``docs/index.rst`` or from another page that is already linked.
   Otherwise, your new page will never be built.

Building Documentation
----------------------

Locally
^^^^^^^

To builds documentation locally you need to have your :doc:`development environment setup <setup>` correctly.
Then from the root folder run

.. code-block::

   make -C docs html

This will create a new folder ``docs/_build/html`` which contains the generated html pages.

My preferred way is to use VS Code with the `reStructuredText plugin`_.
Add the following to your VS Code `settings.json`:

.. code-block::

   {
      "restructuredtext.sphinxBuildPath": "/home/.../miniconda3/envs/cgr-dev/bin/sphinx-build",
      "restructuredtext.linter.executablePath": "/home/.../miniconda3/envs/cgr-dev/bin/rstcheck",
      "restructuredtext.confPath": "${workspaceFolder}/docs",
      "restructuredtext.languageServer.disabled": true,
   }

.. note::

   You must adjust the above paths ``...`` to match your installation.

Then if you preview the ``rst`` file in VS Code it will auto build the docs for you.

Publishing Documentation
^^^^^^^^^^^^^^^^^^^^^^^^

Publishing documentation is as easy as pushing a commit to GitHub prefixed by ``docs:``.
For example,

.. code-block::

   git add -A
   git commit -m "docs: updated the documentation"
   git push

This will automatically triggered a GitHub action to build and publish the documentation.
This action is also triggered when creating a package release.


.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _reStructuredText: https://www.sphinx-doc.org/en/master/
.. _`reStructuredText plugin`: https://github.com/vscode-restructuredtext/vscode-restructuredtext
