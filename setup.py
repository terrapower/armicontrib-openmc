# Copyright 2024 TerraPower, LLC
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
"""Setup.py script for the Open Source OpenMC ARMI plugin"""

from setuptools import setup, find_namespace_packages

with open("README.md") as f:
    README = f.read()

setup(
    name="armicontrib-openmc",
    version="1.1.0",
    description=("ARMI plugin for neutronics analysis with OpenMC."),
    author="TerraPower LLC",
    author_email="armi-devs@terrapower.com",
    packages=find_namespace_packages(),
    long_description=README,
    install_requires=[
        "armi",
    ],
    extras_require={
        "dev": [
            "pytest",
            "sphinx",
            "sphinx_rtd_theme",
            "sphinxcontrib-apidoc",
            "sphinxcontrib-needs",
            "sphinxcontrib-plantuml",
        ]
    },
    keywords=["ARMI, OpenMC", "Neutronics"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Information Analysis",
    ],
    test_suite="tests",
    include_package_data=True,
)
