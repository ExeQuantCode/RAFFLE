# Configuration file for the Sphinx documentation builder.

# -- Project information
import datetime
import os
import sys

from unittest.mock import Mock

MOCK_MODULES = ["raffle._raffle"]  # List any other modules if needed
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

# sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'src', 'raffle')))  # Sets the base path to find your modules
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'src')))  # Sets the base path to find your modules

project = 'RAFFLE'
copyright = f'{datetime.date.today().year}, RAFFLE-developers'
# release = '1.0'
# version = '1.0.0'

# -- General configuration
master_doc = 'index'

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinxcontrib.bibtex',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

exclude_patterns = ['_build', '.DS_Store', 'build']


# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

html_logo = "RAFFLE_logo_no_background.png"
# html_favicon = 'favicon.ico'
html_theme_options = {
    'logo_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    # 'style_nav_header_background': 'white',
    'flyout_display': 'hidden',
    'version_selector': True,
    'language_selector': True,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}


html_context = {
    "display_github": True,
    "github_repo": "RAFFLE",
    "github_user": "ExeQuantCode",
    "github_version": "main",
    "conf_py_path": "/docs/source/",
}

autoclass_content="both"

bibtex_bibfiles = ['references.bib']
