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

    * "OpenMC": Enable OpenM .
"""
import shutil

from armi.physics.neutronics import settings as neutronicsSettings
from armi.settings import setting
from armi.operators import settingsValidation
from armi.physics import neutronics


CONF_EPS_BURN_TIME = "epsBurnTime"
CONF_EPS_CYCLIC = "epsCyclic"
CONF_EPS_NDENS = "epsNdens"
CONF_NEUTRONICS_OUTPUTS_TO_SAVE = "neutronicsOutputsToSave"
CONF_OPENMC_DB = "writeOpenMCDb"


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
    ]
    return settings


def defineSettingValidators(inspector):
    """Define OpenMC-related setting validations."""
    queries = []
    return queries
