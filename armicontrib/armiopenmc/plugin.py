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
Registers elements of the OpenMC plugin with ARMI.

.. pyreverse:: armicontrib.armiopenmc.plugin -A
    :align: center
    :width: 90%
"""

from armi import plugins
from armi import interfaces
from armi.physics.neutronics.globalFlux import globalFluxInterface

from armicontrib.armiopenmc import settings
from armicontrib.armiopenmc import schedulers

ORDER = interfaces.STACK_ORDER.FLUX


class OpenMCPlugin(plugins.ArmiPlugin):
    """Plugin for OpenMC."""

    @staticmethod
    @plugins.HOOKIMPL
    def exposeInterfaces(cs):
        """Function for exposing interface(s) to other code"""
        nk = cs["neutronicsKernel"]
        kwargs = {"enabled": bool(cs["globalFluxActive"]), "bolForce": True}

        if nk in settings.KERNELS:
            return [
                interfaces.InterfaceInfo(
                    globalFluxInterface.ORDER, schedulers.OpenMCInterface, kwargs
                )
            ]

    @staticmethod
    @plugins.HOOKIMPL
    def defineSettings():
        """Define settings."""
        return settings.defineSettings()

    @staticmethod
    @plugins.HOOKIMPL
    def defineSettingsValidators(inspector):
        """Define settings inspections."""
        return settings.defineSettingValidators(inspector)
