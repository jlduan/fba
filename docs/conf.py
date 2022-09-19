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
import os
import sys

sys.path.insert(0, os.path.abspath("../fba"))
from datetime import datetime

from __init__ import __version__

# -- Project information -----------------------------------------------------

project = "fba"
author = "Jialei Duan"
copyright = f"2020-{datetime.today().year}, {author}"

# The full version, including alpha/beta/rc tags
release = __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["myst_parser", "sphinxext.opengraph", "sphinx_copybutton"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The name of the Pygments (syntax highlighting) style to use.
# from pygments.styles import STYLE_MAP; print(sorted(STYLE_MAP.keys()))
pygments_style = "perldoc"
pygments_dark_style = "monokai"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = "furo"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {}

# The name for this set of Sphinx documents.
# "<project> v<release> documentation" by default.
#
# html_title = u'test vtest'
html_title = "FBA"

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#
# html_logo = None

# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#
html_favicon = "logo.ico"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_theme_options = {
    "light_logo": "logo.svg",
    "dark_logo": "logo_dark.svg",
    #
    "sidebar_hide_name": True,
    "navigation_with_keys": True,
    #
}

# If not None, a 'Last updated on:' timestamp is inserted at every page
# bottom, using the given strftime format.
# The empty string is equivalent to '%b %d, %Y'.
#
# html_last_updated_fmt = None

# If true, links to the reST sources are added to the pages.
#
html_show_sourcelink = False

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#
html_show_sphinx = False

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#
html_show_copyright = True

# Strip and configure input prompts for code cells
# https://sphinx-copybutton.readthedocs.io/en/latest/
copybutton_prompt_text = "$"
copybutton_prompt_text = r"\$ |>>> |\.\.\. |In \[\d*\]: | {2,5}\.\.\.: "
copybutton_prompt_is_regexp = True
copybutton_line_continuation_character = "\\"

# Generate OpenGraph metadata
# https://github.com/wpilibsuite/sphinxext-opengraph
ogp_site_url = "https://jlduan.github.io/fba/_build/html/index.html"
ogp_type = "article"
ogp_custom_meta_tags = [
    '<meta property="og:ignore_canonical" content="true" />',
]
