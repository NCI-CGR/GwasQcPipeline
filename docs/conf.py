# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
from cgr_gwas_qc.version import __version__

# -- Project information -----------------------------------------------------

project = "GwasQcPipeline"
copyright = "2020, Leidos Biomedical Research, Inc."
authors = ["Justin Fear", "Eric Karlins", "Jiahui Wang", "Cameron Palmer", "Bari Ballew", "Bin Zhu"]

# The full version, including alpha/beta/rc tags
release = (
    __version__.lstrip("v").replace("-alpha.", "a").replace("-beta.", "b").replace("-rc.", "rc")
)
pkg = f"cgr_gwas_qc-{release}-py3-none-any.whl"


# -- General configuration ---------------------------------------------------
# Add any paths that contain templates here, relative to this directory.
templates_path = ["templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build",
    "api/cgr_gwas_qc.parsers.illumina.IlluminaBeadArrayFiles.rst",
    "common_urls.rst",  # used in rst_epilog not
]

rst_prolog = f"""
.. |pkgurl| replace:: https://github.com/NCI-CGR/GwasQcPipeline/releases/download/{__version__}/{pkg}
.. |cgr| replace:: *CGR GwasQcPipeline*
"""

rst_epilog = ".. include:: /common_urls.rst"

# -- Extension Settings -------------------------------------------------
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "sphinx-pydantic",
    "sphinxcontrib.mermaid",
    "sphinx-prompt",
    "sphinx_substitution_extensions",
    "sphinx_click",
]

todo_include_todos = True


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
numfig = True

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["static"]

html_css_files = ["css/custom.css"]
