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
Components for reading OpenMC output files and updating ARMI state with results.

To apply OpenMC results to a ``reactor``, run something like::

    opts = OpenMCOptions()
    opts.fromUserSettings(cs)
    opts.fromReactor(reactor)
    reader = OpenMCReader(opts)
    reader.apply(reactor)
"""
import os

import numpy as np

from armi.reactor import reactors
from armi import runLog

from armi.bookkeeping.db import Database3

from .const import SolutionType
from .executionOptions import OpenMCOptions

# TODO: should be a GlobalFluxResultMapper subclass to get dpa, etc.
class OpenMCReader:
    """
    Read OpenMC output files and apply to ARMI.

    Uses the ARMI state plus a template.
    """

    def __init__(self, options: OpenMCOptions):
        """
        Initialize the reader and check the output file, but do not read it yet.

        See Also
        --------
        apply
            Applies the results to a reactor object.
        """
        self.opts = options
        self.r: reactors.Reactor = None
        self._check()

    def _check(self):
        """
        Check current directory for valid output files.

        Also check convergence flag and warn if not converged.
        """
        if not os.path.exists("OpenMC"):
            raise RuntimeError("No valid OpenMC output found. Check OpenMC stdout for errors.")
        self._openmcData = openmcFile.OpenMCStream.readBinary(openmcFile.OpenMC)
        runLog.info("Found OpenMC output with:\n\t" + "\n\t".join(self._openmcData.makeSummary()))

        if self._openmcData.convergence != openmcFile.Convergence.CONVERGED:
            runLog.warning(
                f"OpenMC run did not converge. Convergence state is {self._openmcData.convergence}."
            )

    def apply(self, reactor: reactors.Reactor):
        """
        Read data from binary interface files apply to the ARMI model.

        Generally, the armiObject is the Case's Reactor object.
        """
        self.r = reactor
        self._readGeometry()
        self._readKeff()
        self._readPower()
        self._readFluxes()
        self._readPeakFluxes()

        if self.opts.detailedDb is not None:
            with Database3(self.opts.detailedDb, "w") as dbo:
                dbo.writeInputsToDB(self.opts.csObject)
                dbo.writeToDB(self.r)

    def getKeff(self):
        """Return keff directly from output files without applying to reactor."""
        return self._openmcData.keff
