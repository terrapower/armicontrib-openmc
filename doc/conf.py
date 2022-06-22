# Copyright 2021 TerraPower, LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Configuration file for the Sphinx documentation builder.

This file only contains a selection of the most common options. For a full
list see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html

-- Path setup --------------------------------------------------------------

If extensions (or modules to document with autodoc) are in another directory,
add these directories to sys.path here. If the directory is relative to the
documentation root, use os.path.abspath to make it absolute, like shown here.
"""

import os
import sys


sys.path.insert(0, os.path.abspath(os.path.join("..")))


import armi
from armi.utils.dochelpers import PyReverse, ExecDirective


# Configure the baseline framework "App" for framework doc building
armi.configure(armi.apps.App())


def setup(app):
    """Add custom directives"""
    app.add_directive("pyreverse", PyReverse)
    app.add_directive("exec", ExecDirective)


# -- Project information -----------------------------------------------------

project = "ARMI-OpenMC Plugin"
copyright = "2022, TerraPower, LLC"
author = "Nick Touran"

# The full version, including alpha/beta/rc tags
release = "1.0"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.extlinks",
    "sphinx.ext.viewcode",
    "sphinxcontrib.apidoc",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.needs",
    "sphinxcontrib.plantuml",
]

todo_include_todos = True

APIDOC_REL = ".apidocs"
SOURCE_DIR = os.path.join("..", "armicontrib")
APIDOC_DIR = APIDOC_REL

# apidoc config
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
}  # , 'private-members']
autodoc_member_order = "bysource"
autoclass_content = "both"

apidoc_module_dir = os.path.join("..", "armicontrib")
apidoc_output_dir = ".apidocs"
# apidoc_excluded_paths = ["tests", "*/test*"]
apidoc_separate_modules = True
apidoc_module_first = True
apidoc_extra_args = ["--implicit-namespaces"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

numfig = True

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

modindex_common_prefix = ["armicontrib.armiopenmc."]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_last_updated_fmt = "%Y-%m-%d"

# -- Options for LaTeX output ------------------------------------------------
latex_show_urls = "footnote"
latex_engine = "xelatex"

TITLE_PAGE = r"""
\makeatletter
\renewcommand{\sphinxmaketitle}{%
  \let\sphinxrestorepageanchorsetting\relax
  \ifHy@pageanchor\def\sphinxrestorepageanchorsetting{\Hy@pageanchortrue}\fi
  \hypersetup{pageanchor=false}% avoid duplicate destination warnings
  \begin{titlepage}%
    \let\footnotesize\small
    \let\footnoterule\relax
    \noindent\rule{\textwidth}{1pt}\par
      \begingroup % for PDF information dictionary
       \def\endgraf{ }\def\and{\& }%
       \pdfstringdefDisableCommands{\def\\{, }}% overwrite hyperref setup
       \hypersetup{pdfauthor={\@author}, pdftitle={\@title}}%
      \endgroup
    \begin{flushright}%
      \sphinxlogo
      \py@HeaderFamily
      {\Huge \@title \par}
      {\itshape\LARGE \py@release\releaseinfo \par}
      \vfill
      {\LARGE
        \begin{tabular}[t]{c}
          \@author 
        \end{tabular}\kern-\tabcolsep
        \par}
        \vfill
      {\large
       \@date \par
       \vfill
      }%
      \includegraphics[scale=1.0]{../../_static/tplogo.png} 
      \vfill
      {\Large Copyright \textcopyright\ 2021, TerraPower, LLC \\
      \py@authoraddress  \par }
      \vfill
    \end{flushright}%\par
  \end{titlepage}%
  \setcounter{footnote}{0}%
  \let\thanks\relax\let\maketitle\relax
  %\gdef\@thanks{}\gdef\@author{}\gdef\@title{}
  \clearpage
  \sphinxrestorepageanchorsetting
}
\sphinxmaketitle
"""

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # Try really hard to get TOC on one page
    "preamble": r"""\authoraddress{15800 Northup Way\\Bellevue, WA 98008}
    \usepackage[notlot,nottoc,notlof]{}
    \addtocontents{toc}{\vskip -2.2cm}

    \makeatletter
    \fancypagestyle{normal}{
% this is the stuff in sphinx.sty
    \fancyhf{}
    \fancyfoot[LE,RO]{{\py@HeaderFamily\thepage}}
% we comment this out and
    %\fancyfoot[LO]{{\py@HeaderFamily\nouppercase{\rightmark}}}
    %\fancyfoot[RE]{{\py@HeaderFamily\nouppercase{\leftmark}}}
% add copyright stuff
    \fancyfoot[LO,RE]{{Copyright \textcopyright\ 2021, TerraPower, LLC. All Rights Reserved.}}
% again original stuff
    \fancyhead[LE,RO]{{\py@HeaderFamily \@title\sphinxheadercomma\py@release}}
    \renewcommand{\headrulewidth}{0.4pt}
    \renewcommand{\footrulewidth}{0.4pt}
    }
% this is applied to each opening page of a chapter
    \fancypagestyle{plain}{
    \fancyhf{}
    \fancyfoot[LE,RO]{{\py@HeaderFamily\thepage}}
    \renewcommand{\headrulewidth}{0pt}
    \renewcommand{\footrulewidth}{0.4pt}
% add copyright stuff for example at left of footer on odd pages,
% which is the case for chapter opening page by default
    \fancyfoot[LO,RE]{{Copyright \textcopyright\ 2021, TerraPower, LLC. All Rights Reserved.}}
    }
\makeatother
    """,
    # Latex figure (float) alignment
    #
    "figure_align": "H",
    "extraclassoptions": ",openany,oneside",
    "babel": "\\usepackage{babel}",
    "maketitle": TITLE_PAGE,
    "sphinxsetup": "vmargin={1in,1in}",
    # "printindex": "\\footnotesize\\raggedright\\printindex",
    # "atendofbody": r"\renewenvironment{theindex}{}",
    "printindex": "",
}

latex_domain_indices = False


intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "armi": ("https://terrapower.github.io/armi/", None),
}

# requirements options
needs_id_regex = "^R-[a-z_-]+"
