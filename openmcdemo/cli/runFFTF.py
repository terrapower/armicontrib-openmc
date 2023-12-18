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


from armiopenmc.executers import OpenMCExecuter
from armiopenmc import executionOptions

#from openmcdemo import ARMIOpenMCApp

from armi import runLog
from armi import operators
from armi.reactor import reactors
from armi.cli.entryPoint import EntryPoint
from armi.settings import caseSettings
from armi.reactor.tests import test_reactors
from armi.utils import directoryChangers

from armiopenmc.settings import (
    CONF_OPENMC_PATH,
    CONF_N_PARTICLES,
    CONF_N_BATCHES,
    CONF_N_INACTIVE,
    CONF_TALLY_MESH_DIMENSION,
    CONF_ENTROPY_MESH_DIMENSION,
    CONF_ENERGY_GROUP_STRUCTURE,
    CONF_OPENMC_VERBOSITY,
    CONF_N_OMP_THREADS,
    CONF_N_MPI_PROCESSES
)


class RunFFTF(EntryPoint):
    """
    Use armiopenmc to simulate the Fast Flux Test Facility.
    """

    name = "run-fftf"
    settingsArgument = None

    def addOptions(self):
        self.createOptionFromSetting(CONF_OPENMC_PATH)
        self.createOptionFromSetting(CONF_N_PARTICLES)
        self.createOptionFromSetting(CONF_N_BATCHES)
        self.createOptionFromSetting(CONF_N_INACTIVE)
        self.createOptionFromSetting(CONF_TALLY_MESH_DIMENSION)
        self.createOptionFromSetting(CONF_ENTROPY_MESH_DIMENSION)
        self.createOptionFromSetting(CONF_ENERGY_GROUP_STRUCTURE)
        self.createOptionFromSetting(CONF_OPENMC_VERBOSITY)
        self.createOptionFromSetting(CONF_N_OMP_THREADS)
        self.createOptionFromSetting(CONF_N_MPI_PROCESSES)
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
        cs = caseSettings.Settings(os.path.join("FFTF.yaml"))
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
        e = OpenMCExecuter(reactor=o.r, options=opts)

        e.run()

        """
        Analysis not yet implemented. For now, just print mgFlux from one block to verify it's in there.
        """
        print("blocks[460].p['mgFlux']: " + str(o.r.core.getBlocks()[460].p['mgFlux']))
        print("sum(blocks[460].p['mgFlux']): " + str(sum(o.r.core.getBlocks()[460].p['mgFlux'])))

