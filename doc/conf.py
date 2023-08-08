# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os, sys
sys.path.insert(0, os.path.abspath('GRATIOSA'))
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


project = 'GRATIOSA'
copyright = '2023, Forquet, Pineau, Meyer'
author = 'Forquet, Pineau, Meyer'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
source_dir = './doc'
html_dir='./build'
extensions = [
    'sphinx.ext.githubpages',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon'
]

napoleon_use_param = True
napoleon_use_rtype = True
autoclass_content = 'both'
templates_path = ['_templates']
autodoc_member_order = 'bysource'
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
