# ------------------------------------------------------------------------------
# Author(s): David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Configuration file for the Sphinx documentation builder.
# -----------------------------------------------------------------------------
# For a full list of options see the Sphinx documentation:
#   https://www.sphinx-doc.org/en/master/usage/configuration.html
# -----------------------------------------------------------------------------


# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'SUNDIALS Developers Guide'
copyright = '2020, The SUNDIALS Team'
author = 'David J. Gardner, Cody J. Balos, Daniel R. Reynolds, Carol W. Woodward'
version = 'v5.3.0'

# The master toctree document (needed with Sphinx 1.6.7 but not 3.1.1)
master_doc = 'index'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Set the date format (full-month-name day, full-year)
today_fmt = '%B %d, %Y'


# -- Options for LaTeX output ------------------------------------------------

# 1. the rst file name used to generate the LaTeX file
# 2. the name of the LaTeX file to generate (and the resulting PDF file name)
# 3. the document title
# 4. text for \author
# 5. the LaTeX theme
# 6. include the file from 1. in the output
tex_author = """
    David J. Gardner, Cody J. Balos, and Carol S. Woodward\\\\
    {\em Center for Applied Scientific Computing} \\\\
    {\em Lawrence Livermore National Laboratory} \\\\
    \\\\
    Daniel R. Reynolds\\\\
    {\em Department of Mathematics} \\\\
    {\em Southern Methodist University}
    """

latex_documents = [('index', 'sundials_developers_guide.tex', project,
                    tex_author, 'manual', False)]

# LaTeX customizations
latex_elements = {
# paper size option of the document class
'papersize': 'letterpaper',
# font size option of the document class
'pointsize': '10pt',
# set the version number/release name
'releasename' : version,
# arguments to the sphinxsetup macro
'sphinxsetup':
    # the color for titles
    'TitleColor={RGB}{0,0,0},'+
    # disable frames around code-blocks
    'verbatimwithframe=false,'+
    # do not wrap long lines in code-blocks
    'verbatimwrapslines=false,'+
    # background color for code-blocks
    'VerbatimColor={RGB}{240.0,240.0,240.0},'+
    # font used by heading
    'HeaderFamily=\\rmfamily\\bfseries',
# disable the fncychap package
'fncychap':'',
# figure alignment options
'figure_align': 'htbp',
# additional preamble content
'preamble': r'''
% =====================================================
% Start custom preamble (see latex_elements in conf.py)
% =====================================================

% Use ragged-right for the whole document
\usepackage[document]{ragged2e}

% Specify depths for section numbering and table of contents
\setcounter{tocdepth}{2}
\setcounter{secnumdepth}{3}

% Link a footnote to its location in the text
\usepackage{footnotebackref}

% ===================================================
% End custom preamble (see latex_elements in conf.py)
% ===================================================
''',
# custom maketitle
'maketitle': r'''
% =======================================================
% Start custom cover page (see latex_elements in conf.py)
% =======================================================

\makeatletter

% Start roman numbering
\pagenumbering{Roman}

% Title page
\begin{titlepage}
  \newpage
  \null
  \vskip 2em%
  \begin{center}%
    \let \footnote \thanks
    {\huge \rmfamily \@title \par}
    {\Large \rmfamily SUNDIALS \releasename \par}
    \vskip 3.0em%
    {\large \lineskip .5em%
     \begin{tabular}[t]{c}%
       \@author
     \end{tabular}\par}%
    \vskip 1em%
    {\large \@date \par}%
    \vfill
    {\includegraphics[width=0.5\textwidth]{../../../sundials/figures/doc_logo_blue}}
    \vfill
  \end{center}
  \par
  \vskip 1.5em
\end{titlepage}

\makeatother

\clearpage

% Disclaimer
\thispagestyle{empty}% no number of this page
\vglue5\baselineskip
\begin{center}
  {\bf DISCLAIMER}
\end{center}
\noindent
This document was prepared as an account of work sponsored by an agency of
the United States government. Neither the United States government nor
Lawrence Livermore National Security, LLC, nor any of their employees makes
any warranty, expressed or implied, or assumes any legal liability or responsibility
for the accuracy, completeness, or usefulness of any information, apparatus, product,
or process disclosed, or represents that its use would not infringe privately owned rights.
Reference herein to any specific commercial product, process, or service by trade name,
trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement,
recommendation, or favoring by the United States government or Lawrence Livermore National
Security, LLC. The views and opinions of authors expressed herein do not necessarily state
or reflect those of the United States government or Lawrence Livermore National Security, LLC,
and shall not be used for advertising or product endorsement purposes.

\vskip2\baselineskip
\noindent
This work was performed under the auspices of the U.S. Department of Energy by Lawrence
Livermore National Laboratory under Contract DE-AC52-07NA27344.
\vfill
\begin{center}
  Approved for public release; further dissemination unlimited
\end{center}

\clearpage

% Start arabic numbering
\pagenumbering{arabic}

% =====================================================
% End custom cover page (see latex_elements in conf.py)
% =====================================================
'''
}


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
