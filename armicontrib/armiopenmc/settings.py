# Copyright 2022 TerraPower, LLC
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
Define OpenMC-specific ARMI user-configurable settings.

This module implements the :py:meth:`ArmiPlugin.defineSettings()
<armi:armi.plugins.ArmiPlugin.defineSettings()>` and
:py:meth:`ArmiPlugin.defineSettingsValidators()
<armi:armi.plugins.ArmiPlugin.defineSettingsValidators()>` Plugin APIs. Aside from the
settings that control OpenMC's behavior specifically, this provides new options to the
ARMI built-in ``neutronicsKernel`` setting:

    * "OpenMC": Enable OpenMC.
"""
import os
import shutil

from armi.physics.neutronics import settings as neutronicsSettings
from armi.settings import setting
from armi.operators import settingsValidation
from armi.operators.settingsValidation import Query
from armi.physics import neutronics


CONF_EPS_BURN_TIME = "epsBurnTime"
CONF_EPS_CYCLIC = "epsCyclic"
CONF_EPS_NDENS = "epsNdens"
CONF_NEUTRONICS_OUTPUTS_TO_SAVE = "neutronicsOutputsToSave"
CONF_OPENMC_DB = "writeOpenMCDb"
CONF_OPENMC_PATH = "OpenMCExePath"
CONF_N_PARTICLES = "nParticles"
CONF_N_BATCHES = "nBatches"
CONF_N_INACTIVE = "nInactiveBatches"
CONF_TALLY_MESH_DIMENSION = "tallyMeshDimension"
CONF_ENTROPY_MESH_DIMENSION = "entropyMeshDimension"
#CONF_ENERGY_GROUP_STRUCTURE = "energyGroupStructure"
CONF_OPENMC_VERBOSITY = "openmcVerbosity"
CONF_N_OMP_THREADS = "nOMPThreads"
CONF_N_MPI_PROCESSES = "nMPIProcesses"


CONF_OPT_OPENMC = "OpenMC"


KERNELS = {CONF_OPT_OPENMC}


def defineSettings():
    settings = [
        setting.Setting(
            CONF_OPENMC_DB,
            default=False,
            label="Output database with OpenMC neutronics mesh",
            description="If enabled, a database will be created containing the results "
            "of the most recent OpenMC invocation before converting back to the input "
            "mesh. This is useful for visualizing/analyzing the direct results of the "
            "OpenMC run without any mesh conversions taking place.",
        ),
        setting.Option(CONF_OPT_OPENMC, neutronicsSettings.CONF_NEUTRONICS_KERNEL),
        setting.Setting(
            CONF_EPS_BURN_TIME,
            default=1.0,
            label="Burn time eps",
            description=(
                "Burn time eps (Cycle length convergence.  "
                "Set to 1.0 if the cycle length is known.)"
            ),
        ),
        setting.Setting(
            CONF_NEUTRONICS_OUTPUTS_TO_SAVE,
            default="Input/Output",
            label="Save OpenMC Files",
            description=(
                "Defines outputs from OpenMC-based neutronics kernel to be copied from "
                "the fast path to the network drive for inspection, restarts, debugging, "
                "etc."
            ),
            options=["", "Input/Output", "Flux files", "Restart files", "All"],
        ),
        setting.Setting(
            CONF_OPENMC_PATH,
            default="openmc",
            label="OpenMC path",
            description="The path to the OpenMC executable",
            options=[],
        ),
        setting.Setting(
            CONF_N_BATCHES,
            default = 100,
            label = "Number of batches",
            description=("Defines number of batches to be run in an OpenMC simulation.")
        ),
        setting.Setting(
            CONF_N_INACTIVE,
            default = 10,
            label = "Number of inactive batches",
            description=("Defines number of inactive batches to run before recording results"
            "in an OpenMC simulation."
            )
        ),
        setting.Setting(
            CONF_N_PARTICLES,
            default = 1000,
            label = "Number of particles per generation",
            description=("Defines number of particles to be simulated per generation"
            "in OpenMC."
            )
        ),
        setting.Setting(
            CONF_TALLY_MESH_DIMENSION,
            default = [10, 10, 10],
            label = "Tally mesh dimension in OpenMC",
            description=("Defines number of mesh cells in each dimension for fission and"
            "heating tallies."
            )
        ),
        setting.Setting(
            CONF_ENTROPY_MESH_DIMENSION,
            default = [10, 10, 10],
            label = "Entropy mesh dimension in OpenMC",
            description=("Defines number of mesh cells in each dimension for Shannon entropy"
            "mesh. Used to measure source distribution convergence."
            )
        ),
        setting.Setting(
            CONF_OPENMC_VERBOSITY,
            default = 7,
            label = "Verbosity of OpenMC output",
            description=("OpenMC verbosity. From OpenMC documentation:"
                         "1. don’t display any output"
                         "2. only show OpenMC logo"
                         "3. all of the above + headers"
                         "4. all of the above + results"
                         "5. all of the above + file I/O"
                         "6. all of the above + timing statistics and initialization messages"
                         "7. all of the above + k by generation"
                         "9. all of the above + indicate when each particle starts"
                         "10. all of the above + event information"
            )
        ),
        setting.Setting(
            CONF_N_OMP_THREADS,
            default = os.cpu_count(),
            label = "Number of OpenMP threads",
            description=("Number of OpenMP threads to use in OpenMC."
                         "Defaults to number of hardware threads available."
            )
        ),
        setting.Setting(
            CONF_N_MPI_PROCESSES,
            default = 1,
            label = "Number of MPI processes",
            description=("Number of MPI processes to use in OpenMC."
                         "Each process will load its own copy of the problem geometry,"
                         "So high numbers of MPI processes can use a lot of memory for"
                         "complex geometry problems. Using shared-memory OpenMP threads"
                         "can help mitigate this."
            )
        )
    ]
    return settings


def defineSettingValidators(inspector):
    """Define OpenMC-related setting validations."""
    queries = [
        Query(
            lambda: inspector.cs[CONF_N_PARTICLES] < 0,
            "The number of particles is less than 0.",
            "",
            inspector.NO_ACTION,
        )
    ]
    return queries
