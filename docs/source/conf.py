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


# -- Project information -----------------------------------------------------
from datetime import datetime
from sage.env import SAGE_ENV
from sage.docs.conf import *

try:
    import sage.all
except ImportError:
    raise RuntimeError("to build the documentation you need to be inside a Sage shell"
                       "(run first the command 'sage -sh' in a shell")

SAGE_DOC_SRC = SAGE_ENV["SAGE_DOC_SRC"]
SAGE_DOC = SAGE_ENV["SAGE_DOC"]
SAGE_SRC = SAGE_ENV["SAGE_SRC"]
MATHJAX_DIR = SAGE_ENV["MATHJAX_DIR"]

project = 'MQ Estimator'
year = datetime.now().year
copyright = f'{year}, Technology Innovation Institute LLC'
author = 'Emanuele Bellini, Rusydi H. Makarim and Javier Verbel'
package_name = 'mpkc'
package_folder = f'../../{package_name}'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = [] # ['_static']

# -- Options copied from Sagemath conf.py file -------------------------------

# We use MathJax to build the documentation unless the environment
# variable SAGE_DOC_MATHJAX is set to "no" or "False".  (Note that if
# the user does not set this variable, then the script sage-env sets
# it to "True".)

if os.environ.get('SAGE_DOC_MATHJAX', 'no') != 'no' and os.environ.get('SAGE_DOC_MATHJAX', 'no') != 'False':

    extensions.append('sphinx.ext.mathjax')
    mathjax_path = 'MathJax.js?config=TeX-AMS_HTML-full,../mathjax_sage.js'

    # from sage.misc.latex_macros import sage_mathjax_macros
    # html_theme_options['mathjax_macros'] = sage_mathjax_macros()

    # mathjax_relative = os.path.basename(MATHJAX_DIR)

    # It would be really nice if sphinx would copy the entire mathjax directory,
    # (so we could have a _static/mathjax directory), rather than the contents of the directory

    # html_static_path.append(MATHJAX_DIR)
    # exclude_patterns += ['**/'+os.path.join(mathjax_relative, i) for i in ('docs', 'README*', 'test',
    #                                                                       'unpacked', 'LICENSE')]
    from sage.env import SAGE_LOCAL, SAGE_SHARE
    html_static_path.append(SAGE_LOCAL + "/lib/mathjax")  # conda
    html_static_path.append(SAGE_SHARE + "/mathjax")  # sage distribution
else:
    extensions.append('sphinx.ext.imgmath')

# This is to make the verbatim font smaller;
# Verbatim environment is not breaking long lines
from sphinx.highlighting import PygmentsBridge
from pygments.formatters.latex import LatexFormatter


class CustomLatexFormatter(LatexFormatter):
    def __init__(self, **options):
        super(CustomLatexFormatter, self).__init__(**options)
        self.verboptions = r"formatcom=\footnotesize"


PygmentsBridge.latex_formatter = CustomLatexFormatter

latex_elements['preamble'] += r'''
% One-column index
\makeatletter
\renewenvironment{theindex}{
  \chapter*{\indexname}
  \markboth{\MakeUppercase\indexname}{\MakeUppercase\indexname}
  \setlength{\parskip}{0.1em}
  \relax
  \let\item\@idxitem
}{}
\makeatother
\renewcommand{\ttdefault}{txtt}
'''


def setup(app):
    app.connect('autodoc-process-docstring', process_docstring_cython)
    app.connect('autodoc-process-docstring', process_directives)
    app.connect('autodoc-process-docstring', process_docstring_module_title)
    app.connect('autodoc-process-docstring', process_dollars)
    #app.connect('autodoc-process-docstring', process_inherited)
    app.connect('autodoc-process-docstring', skip_TESTS_block)
    app.connect('autodoc-skip-member', skip_member)
    app.add_transform(SagemathTransform)

    # When building the standard docs, app.srcdir is set to SAGE_DOC_SRC +
    # 'LANGUAGE/DOCNAME', but when doing introspection, app.srcdir is
    # set to a temporary directory.  We don't want to use intersphinx,
    # etc., when doing introspection.
    if app.srcdir.startswith(SAGE_DOC_SRC):
        app.add_config_value('intersphinx_mapping', {}, False)
        app.add_config_value('intersphinx_cache_limit', 5, False)
        # We do *not* fully initialize intersphinx since we call it by hand
        # in find_sage_dangling_links.
        #   app.connect('missing-reference', missing_reference)
        app.connect('missing-reference', find_sage_dangling_links)
        import sphinx.ext.intersphinx
        app.connect('builder-inited', set_intersphinx_mappings)
        app.connect('builder-inited', sphinx.ext.intersphinx.load_mappings)
        app.connect('builder-inited', nitpick_patch_config)
