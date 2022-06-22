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

"""
Components for writing OpenMC inputs from ARMI models.


.. pyreverse:: armicontrib.openmc.inputWriters -A
    :align: center
"""
from typing import NamedTuple

import jinja2

from armi.nucDirectory import nuclideBases
from armi import runLog
from armi.physics import neutronics
from armi.reactor import geometry
from armi.utils.units import ASCII_LETTER_A, ASCII_ZERO
from armi.utils import hexagon
from armi.reactor.flags import Flags


from . import const

import openmc
from openmc.model import hexagonal_prism, ZCylinder


class OpenMCWriter:
    """
    Write OpenMC data using the openmc python api.
    """

    def __init__(self, reactor, options):
        """Build the writer"""
        self._env = None
        self.r = reactor
        self.options = options
        self.geometryWriter = HexWriter(reactor, options)

    def write(self, stream):
        """Write the input file to a stream."""

        runLog.info(f"Writing OpenMC input based on: {self.r.core}")
        cells = []
        for a in self.r.core.getChildrenWithFlags(Flags.FUEL):
            for b in a.getChildrenWithFlags(Flags.FUEL)[0:1]:
                cells.extend(self.geometryWriter.addBlock(b))

        universe = openmc.Universe(cells=cells)
        import matplotlib.pyplot as plt

        #        universe.plot(width=(150, 150), pixels=(1200, 1200))
        universe.plot(width=(20, 20), pixels=(1200, 1200))
        plt.show()


class GeometryWriter:
    """Handles geometry-specific information on A.NIP3."""

    def __init__(self, reactor, options):
        self.r = reactor
        self.options = options

    def makeMat(self, component):
        """Make an OpenMC material from an ARMI component."""
        mat = openmc.Material(name=component.material.name)
        for nuc, dens in component.getNumberDensities().items():
            if "C" in nuc:
                continue
            mat.add_nuclide(nuc.capitalize(), dens)

        mat.set_density("g/cc", component.material.density(Tc=component.temperatureInC))

        return mat


class HexWriter(GeometryWriter):
    def addBlock(self, block):
        """
        Add a hexblock to the model.

        This just calls various OpenMC api calls
        """
        xc, yc, zc = block.spatialLocator.getGlobalCoordinates()
        side = hexagon.side(block.getPitch())

        intercoolant = block.getChildrenWithFlags(Flags.INTERCOOLANT)[0]
        coolant = block.getChildrenWithFlags(Flags.COOLANT)[0]
        duct = block.getChildrenWithFlags(Flags.DUCT)[0]
        clad = block.getChildrenWithFlags(Flags.CLAD)[0]
        wire = block.getChildrenWithFlags(Flags.WIRE)[0]
        fuel = block.getChildrenWithFlags(Flags.FUEL)[0]

        # surfaces
        interstitial = hexagonal_prism(
            edge_length=side,
            origin=(xc, yc),
            boundary_type="reflective",
            orientation="x",
        )
        duct_outer = hexagonal_prism(
            edge_length=hexagon.side(duct.p.op), origin=(xc, yc), orientation="x"
        )
        duct_inner = hexagonal_prism(
            edge_length=hexagon.side(duct.p.ip), origin=(xc, yc), orientation="x"
        )

        # cells
        cells = []

        # interstitial coolant
        cells.append(
            openmc.Cell(fill=self.makeMat(intercoolant), region=(interstitial & ~duct_outer))
        )
        # duct
        cells.append(openmc.Cell(fill=self.makeMat(duct), region=(~duct_inner & duct_outer)))

        fuelMat = self.makeMat(fuel)
        cladMat = self.makeMat(clad)

        clad_surfs = []

        for pinLoc in fuel.spatialLocator:
            xc, yc, zc = pinLoc.getGlobalCoordinates()

            fuel_od = ZCylinder(r=fuel.p.od / 2.0, x0=xc, y0=yc)
            clad_od = ZCylinder(r=clad.p.od / 2.0, x0=xc, y0=yc)
            clad_id = ZCylinder(r=clad.p.id / 2.0, x0=xc, y0=yc)

            cells.append(openmc.Cell(fill=cladMat, region=(-clad_od & +clad_id)))
            # fuel
            cells.append(openmc.Cell(fill=fuelMat, region=(-fuel_od)))

            clad_surfs.append(clad_od)

        cool_reg = duct_inner
        for clad_surf in clad_surfs:
            cool_reg = cool_reg & -clad_surf
        # cool_reg = duct_inner & +clad_od

        # coolant
        cells.append(openmc.Cell(fill=self.makeMat(coolant), region=cool_reg))

        return cells
