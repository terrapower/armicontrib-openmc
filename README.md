# Open-Source ARMI-OpenMC Plugin

This code connects TerraPower's ARMI® engineering automation code with ANL's OpenMC monte
carlo radiation transport code.


## Features

* Build OpenMC models based on ARMI representations of a reactor in Hex
  geometry

* Execute OpenMC

* Read OpenMC output information and store it back on the ARMI reactor model for
  coupling with other codes (e.g. for thermal/hydraulics, fuel
  performance etc.)

## Limitations

* Only supports a subset of OpenMC geometries

* Does not support all OpenMC features

## Installation

The ARMI-OpenMC plugin and demonstration application, along with ARMI and all other
required dependencies can be installed by running

    pip install -r requirements.txt

from within a Python environment. It is recommended to do this inside of a
Python [Virtual Environment](https://docs.python.org/3/tutorial/venv.html) to
gain better control of installed dependencies and their versions.

