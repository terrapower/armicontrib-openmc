# Open-Source ARMI-OpenMC Plugin

This code connects TerraPower's ARMI® engineering automation code with ANL's OpenMC monte
carlo radiation transport code.


## Features

* Build OpenMC models based on ARMI representations of a reactor in Hex or Cartesian
  geometry

* Execute OpenMC

* Read OpenMC output information and store it back on the ARMI reactor model for
  coupling with other codes (e.g. for thermal/hydraulics, fuel
  performance etc.)

## Limitations

* Only supports a subset of OpenMC geometries

* Does not support all OpenMC features

## Installation

The ARMI-OpenMC plugin and demonstration application, along with ARMI and all
non-OpenMC required dependencies can be installed by running

    pip install -r requirements.txt

from within a Python environment. It is recommended to do this inside of a
Python [Virtual Environment](https://docs.python.org/3/tutorial/venv.html) to
gain better control of installed dependencies and their versions.

The OpenMC dcoumentation has [brief](https://docs.openmc.org/en/stable/quickinstall.html)
and [more detailed](https://docs.openmc.org/en/stable/usersguide/install.html)
instructions for installing OpenMC.

## Running the Example Case

The plugin contains an example case that simulates the Fast Flux Test Facility (FFTF)
and generates two flux spectrum plots. To run the example case,

1. Switch to the "develop" branch of the ARMI-OpenMC repository.

2. Clone the [FFTF ARMI model](https://github.com/AidanMcDonald/fftf-isothermal-model)
and switch to the "openmcSettings" branch. Note that this fork differs from 
the [main FFTF ARMI model](https://github.com/terrapower/fftf-isothermal-model) only in
some execution settings in FFTF.yaml. It is recommended to check FFTF.yaml and choose
execution settings appropriate for your system.

3. From the FFTF case directory, execute the "run-fftf" entry point with
```
python -m openmcdemo run-fftf
```
This will instatiate the ARMI model, create the OpenMC model, and run the OpenMC model.
Two flux spectrum plots will be displayed after OpenMC finishes executing.
Instructions for exploring the OpenMC output in more detail can be found in the
[OpenMC documentation](https://docs.openmc.org/en/stable/usersguide/processing.html).
