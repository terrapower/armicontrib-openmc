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

import openmc

from armi.reactor import reactors
from armi import runLog

from .executionOptions import OpenMCOptions


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
        self.nf = None
        self._sp: openmc.StatePoint = None
        self._geometry: openmc.Geometry = None
        self._check()

        if self.opts.energyMode == "multigroup":
            self.heatingTallyName = "fission"
        else:
            self.heatingTallyName = "heating-local"

    def _check(self):
        if not os.path.exists(self.opts.outputFile):
            raise RuntimeError("No valid OpenMC output found. Check OpenMC stdout for errors.")

    def apply(self, reactor: reactors.Reactor):
        """
        Read data from binary interface files apply to the ARMI model.

        Generally, the armiObject is the Case's Reactor object.
        """

        self.r = reactor
        self._readStatePoint()
        self.getNormalizationFactor()
        self._readGeometry()
        self._readKeff()
        self._readPower()
        self._readFluxes()

        # if self.opts.detailedDb is not None:
        #    with Database3(self.opts.detailedDb, "w") as dbo:
        #        dbo.writeInputsToDB(self.opts.csObject)
        #        dbo.writeToDB(self.r)

    def getKeff(self):
        """Return keff directly from output files without applying to reactor."""
        return self._sp.keff.nominal_value

    def getNormalizationFactor(self):
        """
        OpenMC returns tallies in units per source particle.
        Normalize to usable units with heating tally and known reactor power.
        """
        totalHeatingTally = (
            sum(self._sp.get_tally(scores=[self.heatingTallyName]).mean) * 1.602e-19
        )  # [J/(sourceParticle)]
        if self.heatingTallyName == "fission":
            totalHeatingTally = totalHeatingTally * 200e6
        if self.r is None:
            raise ValueError(
                "OpenMCReader.r must be set before normalization factor can be calculated."
            )
        power = self.opts.power  # [J/s]
        self.nf = power / totalHeatingTally  # [sourceParticle/s]

    def _readStatePoint(self):
        """Load statepoint output file into memory"""
        runLog.info("Reading OpenMC output file...")
        self._sp = openmc.StatePoint(self.opts.outputFile)

    def _readGeometry(self):
        """
        Read openmc geometry object into memory. This uses a lot of memory,
        but we need it to ensure blockFiltered tally results are written to the correct blocks.
        """
        openmc.reset_auto_ids()  # openmc will complain about reused ids if we don't reset
        self._geometry = openmc.Geometry.from_xml("geometry.xml")

    def _readKeff(self):
        """Store keff from the outputs onto the reactor model."""
        self.r.core.p.keff = self._sp.keff.nominal_value

    def _readPower(self):
        """Read power density"""
        powerTally = self._sp.get_tally(name="power")
        blockFilter = powerTally.find_filter(openmc.CellFilter)
        reshapedPowerTally = powerTally.get_reshaped_data()
        cells = self._geometry.get_all_cells()

        for i in range(len(reshapedPowerTally)):
            cellNumber = blockFilter.bins[i]
            blockName = cells[cellNumber].name
            b = self.r.core.getBlockByName(blockName)

            if self.heatingTallyName == "fission":
                blockPower = reshapedPowerTally[i] * self.nf * 1.602e-19 * 200e6  # [W]
            else:
                blockPower = reshapedPowerTally[i] * self.nf * 1.602e-19  # [W]

            b.p.power = blockPower  # [W]
            b.p.pdens = blockPower / b.getVolume()  # [W/cm^3]
            # Our best estimate for peak power is average power for now.
            # This can be improved by a finer meshed tally inside the block.
            b.p.ppdens = blockPower / b.getVolume()

    def _readFluxes(self):
        """
        Read fluxes from flux tally.
        """

        fluxTally = self._sp.get_tally(id=102)
        blockFilter = fluxTally.find_filter(openmc.CellFilter)
        reshapedFluxTally = fluxTally.get_reshaped_data()
        cells = self._geometry.get_all_cells()

        for i in range(len(reshapedFluxTally[:, 0])):
            cellNumber = blockFilter.bins[i]
            blockName = cells[cellNumber].name
            b = self.r.core.getBlockByName(blockName)

            # Energy groups in openmc are in ascending order, so we need to reverse first
            blockFluxData = list(reshapedFluxTally[i, :] * self.nf)
            blockFluxData.reverse()
            setattr(b.p, "mgFlux", blockFluxData)
            setattr(b.p, "flux", sum(blockFluxData) / b.getVolume())
            # Our best estimate for peak flux is average flux for now.
            # This can be improved by a finer meshed tally inside the block.
            setattr(b.p, "fluxPeak", sum(blockFluxData) / b.getVolume())
