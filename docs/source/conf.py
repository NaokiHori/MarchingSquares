import os
import sys

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MarchingSquares'
copyright = '2022, Naoki Hori'
author = 'Naoki Hori'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

sys.path.append(os.path.abspath("./ext"))
extensions = [
        "myliteralinclude",
        "sphinx.ext.mathjax",
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_theme_options = {
    "fixed_sidebar": "true",
    "github_user": "NaokiHori",
    "github_repo": "MarchingSquares",
    "github_type": "true",
}

