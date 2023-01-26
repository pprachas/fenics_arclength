# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'fenics-arclength'
copyright = '2023, Peerasait Prachaseree, Saeed Mohammadzadeh, Emma Lejeune'
author = 'Peerasait Prachaseree, Saeed Mohammadzadeh, Emma Lejeune'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))
print(sys.path)

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'myst_nb'
]

myst_enable_extensions = [
    'amsmath',
    'dollarmath',
    'html_image',
    'substitution'
]
source_suffix = {
    '.rst': 'restructuredtext',
    '.ipynb': 'myst-nb',
    '.myst': 'myst-nb',
}

html_theme_options = {
    'repository_url': 'https://github.com/pprachas/fenics_arclength',
    'use_repository_button': True 
}
nb_execution_mode = 'off'
nb_number_source_lines = True
nb_render_markdown_format = 'myst'
templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'

# html_sidebars = {
#    '**': ['globaltoc.html', 'sourcelink.html', 'searchbox.html']
# }

html_static_path = ['_static']
