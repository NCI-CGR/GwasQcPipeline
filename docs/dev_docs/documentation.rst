.. _documentation:

Working with Documentation
==========================

|cgr| uses the Sphinx_ python documentation generator.
Documentation is written in reStructuredText_, a markup language similar to Markdown but with more advanced features.
The sphinx documentation is stored in the ``docs/`` folder.

.. code-block:: bash

   .
   ├── docs/
   │  ├── api/              # Implementation details for parsers and other helpers
   │  ├── common_urls.rst   # URLs used throughout documentation
   │  ├── conf.py           # Sphinx configuration file
   │  ├── dev_docs/         # Development documentation
   │  ├── getting_started/  # Getting started documentation
   │  ├── index.rst         # Main page, if you add a new page link to it here
   │  ├── Makefile          # Sphinx make file for building docs ``make html``
   │  ├── reference/        # User reference material, i.e., file-type reference
   │  ├── static/           # Images and styling
   │  └── sub_workflows/    # Workflow documentation

However, I also embed documentation directly from the source code using keywords like ``.. automodule``, ``.. pydantic``, and ``.. click``.
This allows me to have details in a single location as close to the code as possible.

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
