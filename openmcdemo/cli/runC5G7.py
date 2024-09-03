# Copyright 2023 TerraPower, LLC
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
import sys
import numpy as np
import openmc
import matplotlib.pyplot as plt

from armi import runLog
from armi import operators
from armi.reactor import reactors
from armi.cli.entryPoint import EntryPoint
from armi.settings import caseSettings
from armi.physics.neutronics import energyGroups
from armi.reactor.flags import Flags
from armi.physics.neutronics.settings import CONF_GROUP_STRUCTURE
from armicontrib.armiopenmc.executers import OpenMCExecuter
from armicontrib.armiopenmc.inputWriters import parseEnergyGroupStructure
from armicontrib.armiopenmc import executionOptions
from armi.reactor.converters.blockConverters import MultipleComponentMerger

from armicontrib.armiopenmc.settings import (
    CONF_OPENMC_PATH,
    CONF_N_PARTICLES,
    CONF_N_BATCHES,
    CONF_N_INACTIVE,
    CONF_TALLY_MESH_DIMENSION,
    CONF_ENTROPY_MESH_DIMENSION,
    CONF_OPENMC_VERBOSITY,
    CONF_N_OMP_THREADS,
    CONF_N_MPI_PROCESSES,
    CONF_VERTICAL_SYMMETRY,
)


class RunC5G7(EntryPoint):
    """
    Use armiopenmc to simulate the C5G7 thermal reactor benchmark.
    """

    name = "run-c5g7"
    settingsArgument = None

    def addOptions(self):
        self.createOptionFromSetting(CONF_OPENMC_PATH)
        self.createOptionFromSetting(CONF_N_PARTICLES)
        self.createOptionFromSetting(CONF_N_BATCHES)
        self.createOptionFromSetting(CONF_N_INACTIVE)
        self.createOptionFromSetting(CONF_TALLY_MESH_DIMENSION)
        self.createOptionFromSetting(CONF_ENTROPY_MESH_DIMENSION)
        self.createOptionFromSetting(CONF_GROUP_STRUCTURE)
        self.createOptionFromSetting(CONF_OPENMC_VERBOSITY)
        self.createOptionFromSetting(CONF_N_OMP_THREADS)
        self.createOptionFromSetting(CONF_N_MPI_PROCESSES)
        self.createOptionFromSetting(CONF_VERTICAL_SYMMETRY)
        self.parser.add_argument(
            "--inputs-only",
            action="store_true",
            default=False,
            help="Just make input files for the suite; don't run",
        )

    @staticmethod
    def _initSettings():
        """
        Provide hard-coded settings rather than user-provided ones.

        Notes
        -----
        We need to override this method to set our own settings object.

        These settings are modified because the base case is used for many other
        tests across the ARMI ecosystem. We start from them but modify as necessary
        to get the fixture run set up. By re-using this input file, we will
        benefit from getting automatically-updated test inputs as the ARMI
        framework progresses in the future.
        """
        cs = caseSettings.Settings(os.path.join("c5g7-settings.yaml"))
        print("setting from _initSettings", cs[CONF_OPENMC_PATH])

        return cs

    def invoke(self):
        if shutil.which(self.cs[CONF_OPENMC_PATH]) is None:
            runLog.error(
                "The requested OpenMC executable, `{}` cannot be found".format(
                    self.cs[CONF_OPENMC_PATH]
                )
            )
            if not self.args.inputs_only:
                sys.exit(1)

        o = operators.factory(self.cs)
        o.r = reactors.loadFromCs(self.cs)

        opts = executionOptions.OpenMCOptions()
        opts.fromReactor(o.r)
        opts.fromUserSettings(o.cs)

        # Add custom tallies
        tallies = openmc.Tallies()

        # Example custom tally
        mesh = openmc.RegularMesh()
        mesh.lower_left = [0.0, 0.0, 0.0]
        mesh.upper_right = [
            64.26,
            64.26,
            max([assembly.getHeight() for assembly in o.r.core]),
        ]
        mesh.dimension = [1, 1500, 1500]
        meshFilter = openmc.MeshFilter(mesh=mesh)
        groups = parseEnergyGroupStructure(energyGroups.getGroupStructure(opts.groupStructure))
        energyFilter = openmc.EnergyFilter(groups, filter_id=35)
        meshFluxTally = openmc.Tally(1, name="custom tally")
        meshFluxTally.scores = ["flux"]
        meshFluxTally.filters = [meshFilter, energyFilter]
        tallies.append(meshFluxTally)

        opts.addTallies(tallies)

        e = OpenMCExecuter(reactor=o.r, options=opts)

        e.run()
