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

"""App for demo."""
import armi
from armi import plugins


class ARMIOpenMCApp(armi.apps.App):
    """App that adds only the OpenMC plugin for testing purposes."""

    name = "armi-openmc-demo"

    def __init__(self):

        armi.apps.App.__init__(self)

        from armicontrib.armiopenmc.plugin import OpenMCPlugin

        self._pm.register(OpenMCPlugin)

    @property
    def splashText(self):
        return """
     =================================
     == ARMI-OpenMC Demo Application ==
     =================================
     ARMI Version: {}
""".format(
            armi.__version__
        )
