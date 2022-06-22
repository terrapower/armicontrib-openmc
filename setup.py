"""Setup.py script for the Open Source OpenMC ARMI plugin"""

from setuptools import setup, find_namespace_packages

with open("README.md") as f:
    README = f.read()

setup(
    name="armicontrib-openmc",
    version="1.0",
    description=("ARMI plugin for neutronics analysis with OpenMC."),
    author="TerraPower LLC",
    author_email="armi-devs@terrapower.com",
    packages=find_namespace_packages(),
    package_data={
        "armicontrib.armiopenmc": ["templates/*", "templates/**"],
    },
    entry_points={"console_scripts": ["armi-openmc-demo=openmcdemo:main"]},
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
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Information Analysis",
    ],
    test_suite="tests",
    include_package_data=True,
)
