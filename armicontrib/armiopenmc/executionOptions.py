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

import os
import shutil

from armi import runLog
from armi.physics.neutronics.globalFlux import globalFluxInterface

from armi.settings import caseSettings
from armi.physics import neutronics
from armi.physics.neutronics import settings as gsettings

from . import settings


class OpenMCOptions(globalFluxInterface.GlobalFluxOptions):
    """Define options for one particular OpenMC execution."""

    def __init__(self, label=None):
        globalFluxInterface.GlobalFluxOptions.__init__(self, label)
        self.templatePath = None
        self.kernelName = "OpenMC"
        self.libDataFile = neutronics.ISOTXS
        self.label = label if label else "openmc"
        self.executablePath = None
        self.runDir = None
        self.numberMeshPerEdge = 1
        self.neutronicsOutputsToSave = None
        self.xsLibraryName = None
        self.existingFixedSource = None
        self.bcCoefficient = None
        self.detailedDb: Optional[str] = None
        # When writing the neutronics-mesh database, we need the actual case Settings
        # object in order for that database to be complete and capable of being read.
        # Storing the full cs object here seems contrary to the goal of Options classes
        # to decouple Executers from raw Settings. However, here we need the **object**,
        # rather than the settings that it contains. Future work to the ARMI framework
        # may relax the need for case settings in a database file to support loading, in
        # which case this would no longer be required.
        self.csObject: Optional[caseSettings.Settings] = None

    def fromUserSettings(self, cs: caseSettings.Settings):
        """Set options from user settings"""
        globalFluxInterface.GlobalFluxOptions.fromUserSettings(self, cs)
        self.setRunDirFromCaseTitle(cs.caseTitle)
        self.neutronicsOutputsToSave = cs[settings.CONF_NEUTRONICS_OUTPUTS_TO_SAVE]
        self.existingFixedSource = cs[gsettings.CONF_EXISTING_FIXED_SOURCE]
        self.epsFissionSourceAvg = cs[gsettings.CONF_EPS_FSAVG]
        self.epsFissionSourcePoint = cs[gsettings.CONF_EPS_FSPOINT]
        self.epsEigenvalue = cs[gsettings.CONF_EPS_EIG]
        self.numberMeshPerEdge = cs[gsettings.CONF_NUMBER_MESH_PER_EDGE]

    def fromReactor(self, reactor):
        """Set options from an ARMI composite to be modeled (often a ``Core``)"""
        globalFluxInterface.GlobalFluxOptions.fromReactor(self, reactor)
        self.inputFile = f"{self.label}.inp"
        self.outputFile = f"{self.label}.out"

    def resolveDerivedOptions(self):
        """
        Set other options that are dependent on previous phases of loading options.
        """
        globalFluxInterface.GlobalFluxOptions.resolveDerivedOptions(self)
        if self.existingFixedSource:
            self.extraInputFiles.append((self.existingFixedSource, self.existingFixedSource))
        if self.isRestart:
            for _label, fnames in fileSetsHandler.specifyRestartFiles(self).items():
                self.extraInputFiles.extend([(f, f) for f in fnames])
