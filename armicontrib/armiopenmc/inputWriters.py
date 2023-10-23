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
import warnings

import armi
from armi.nucDirectory import nuclideBases
from armi import runLog
from armi.physics import neutronics
from armi.reactor import geometry
from armi.utils.units import ASCII_LETTER_A, ASCII_ZERO
from armi.utils import hexagon
from armi.reactor.flags import Flags
from armi.reactor.converters.blockConverters import MultipleComponentMerger


from . import const

import openmc
from openmc.model import hexagonal_prism, ZCylinder

class OpenMCWriter:
    """
    Write OpenMC data using the openmc python api.
    """

    def __init__(self, reactor, options):
        """Build the writer"""

        self.r = reactor
        self.options = options
        self.materialMap = None  # Written by writeMaterials

    '''
    def writeMaterials(self):
        """Write the openmc materials input file."""

        core = self.r.core
        materials = openmc.Materials()
        self.materialMap = dict()

        for assembly in core:
            for block in assembly:
                blockIndicesString = str(block.spatialLocator.getCompleteIndices())
                self.materialMap[blockIndicesString] = dict()
                for component in block:
                    componentMaterial = openmc.Material(name=component.material.name)
                    componentMaterial.set_density("g/cm3", component.material.getProperty("pseudoDensity", Tc=component.temperatureInC))

                    componentNuclides = component.getNuclides()
                    componentNuclideDensities = component.getNuclideNumberDensities(componentNuclides)
                    totalComponentNuclideDensity = sum(componentNuclideDensities)

                    for i, nuclideName in enumerate(componentNuclides):
                        nuclide = nuclideBases.byName[nuclideName]
                        if nuclide.a>0:  # ignore lumped fission products and dummy nuclides for now
                            nuclideGNDSName = openmc.data.gnds_name(Z=nuclide.z, A=nuclide.a)
                            componentMaterial.add_nuclide(nuclideGNDSName, componentNuclideDensities[i]/totalComponentNuclideDensity, 'ao')
                    self.materialMap[blockIndicesString][component.name] = componentMaterial
                    materials.append(componentMaterial)
        materials.export_to_xml()
    '''

    def writeGeometry(self):
        """Write the openmc geometry input file."""
        
        def blendHelixComponentsIntoCoolant(block, solventName='coolant'):
            """OpenMC doesn't support helixes, so blend them all into coolant"""
            helixComponentNames = []
            for component in block:
                if isinstance(component, armi.reactor.components.complexShapes.Helix):
                    helixComponentNames.append(component.getName()) 
            if len(helixComponentNames)>0:
                block = MultipleComponentMerger(sourceBlock=block, soluteNames=helixComponentNames, solventName=solventName).convert()
            return block
        
        def buildCellRegion(component, bottomPlane, topPlane):
            """Build region based on shape"""

            # Circle
            if isinstance(component, armi.reactor.components.basicShapes.Circle):
                innerCylinder = openmc.ZCylinder(r=component.getDimension("id")/2)
                outerCylinder = openmc.ZCylinder(r=component.getDimension("od")/2)
                region = +bottomPlane & -topPlane & +innerCylinder & -outerCylinder
                return region

            # Hexagon
            if isinstance(component, armi.reactor.components.basicShapes.Hexagon):
                innerHexPrism = openmc.model.hexagonal_prism(edge_length=component.getDimension("ip")/3**.5,
                                                             orientation='x') # Check on orientation
                outerHexPrism = openmc.model.hexagonal_prism(edge_length=component.getDimension("op")/3**.5,
                                                             orientation='x')
                region = +bottomPlane & -topPlane & ~innerHexPrism & outerHexPrism
                return region

            # Rectangle
            if isinstance(component, armi.reactor.components.basicShapes.Hexagon):
                innerRectPrism = openmc.model.rectangular_prism(width=component.getDimension("widthInner"),
                                                                height=component.getDimension("lengthInner")) # Check that width/height aren't flipped
                outerRectPrism = openmc.model.rectangular_prism(width=component.getDimension("widthOuter"),
                                                                height=component.getDimension("lengthOuter"))
                region = +bottomPlane & -topPlane & ~innerRectPrism & outerRectPrism
                return region

            # DerivedShape
            if isinstance(component, armi.reactor.components.DerivedShape):
                # DerivedShape is supported, but we need to set DerivedShape regions after all others in the block.
                # We simply return here and explicitly set region later
                return

            # Helix
            if isinstance(component, armi.reactor.components.complexShapes.Helix):
                # Note: Helix components are automatically blended into coolant.
                warnings.warn("Helix shape not supported by OpenMC. Ignoring helical component.")
                return

            # Others
            raise NotImplementedError("Shape type not supported yet")

        def cartesianToRing(cartesianIndices):
            """Convenience function for converting from Cartesian coordinate system to OpenMC's ring system"""
            x = cartesianIndices[0]
            y = cartesianIndices[1]

            if x*y>=0:
                ring = abs(x)+abs(y)
                if x>0 or y>0:
                    pos = 6*ring/6-x
                elif x<0 or y<0:
                    pos = 6*ring*3/6+abs(y)
                else:
                    pos = 0

            elif x*y<0:
                ring = max([abs(x), abs(y)])
                if abs(x)==ring:
                    if x<0:
                        pos = 6*ring*3/6-y
                    else:
                        pos = 6*ring-abs(y)
                else:
                    if y>0:
                        pos = 6*ring*1/6+abs(x)
                    else:
                        pos = 6*ring*4/6+x

            return [ring, int(pos)]

        def buildComponentMaterial(component):
            if component.material.name=="Void":
                return None
            componentMaterial = openmc.Material(name=component.material.name)
            componentMaterial.set_density("g/cm3", component.material.getProperty("pseudoDensity", Tc=component.temperatureInC))

            componentNuclides = component.getNuclides()
            componentNuclideDensities = component.getNuclideNumberDensities(componentNuclides)
            totalComponentNuclideDensity = sum(componentNuclideDensities)

            for i, nuclideName in enumerate(componentNuclides):
                nuclide = nuclideBases.byName[nuclideName]
                if nuclide.a>0:  # ignore lumped fission products and dummy nuclides for now
                    nuclideGNDSName = openmc.data.gnds_name(Z=nuclide.z, A=nuclide.a)
                    componentMaterial.add_nuclide(nuclideGNDSName, componentNuclideDensities[i]/totalComponentNuclideDensity, 'ao')
            return componentMaterial

        core = self.r.core

        materials = openmc.Materials()
        emptyMaterial = openmc.Material()

        if core.geomType == armi.reactor.geometry.GeomType.HEX:
            # Openmc uses awkward ring indexing system for hex lattices

            """There may be a cleaner way to do this with core.getAssembliesInRing(), etc."""
            numRings = 1
            #numRings = core.numRings

            boundingCylinderBottomPlane = openmc.ZPlane(z0=0.0, boundary_type='vacuum')
            boundingCylinderTopPlane = openmc.ZPlane(z0=max([assembly.getHeight() for assembly in core]), boundary_type='vacuum')
            boundingCylinder = openmc.ZCylinder(r=core.getBoundingCircleOuterDiameter()/2, boundary_type='vacuum')
            
            emptyCell = openmc.Cell(region= +boundingCylinderBottomPlane & -boundingCylinderTopPlane & -boundingCylinder)

            # Create blank list of assembly universes for the lattice - we will fill in universes individually later
            emptyUniverse = openmc.Universe(cells=[emptyCell])
            assemblyRings = [None]*numRings
            assemblyRings[0] = [emptyUniverse]
            ringLength = 6
            for i in range(1, numRings):
                assemblyRings[i] = [emptyUniverse]*ringLength
                ringLength += 6

        else:
            raise TypeError("Unsupported geometry type")

        for assembly in core[156:157]:
            # Write a universe for each assembly
            print("working on assembly " + str(assembly.getName()))

            # Create ZPlanes between blocks
            z0 = 0.0
            ZPlanes=[openmc.ZPlane(z0=z0, boundary_type='vacuum')]
            for block in assembly:
                z0 += block.getHeight()
                ZPlanes.append(openmc.ZPlane(z0=z0, boundary_type='transmission'))
            ZPlanes[-1].boundary_type = 'vacuum' # Reset top plane boundary condition to vacuum

            componentCellsInAssembly = []
            for i, block in enumerate(assembly):
                block = blendHelixComponentsIntoCoolant(block)
                componentCellsInBlock = [None]*len(block)
                for j, component in enumerate(block):
                    componentMaterial = buildComponentMaterial(component)
                    cell = openmc.Cell(name=component.getName(),
                                       fill=componentMaterial,
                                       region=buildCellRegion(component, bottomPlane=ZPlanes[i], topPlane=ZPlanes[i+1]))
                    if componentMaterial is not None:
                        materials.append(componentMaterial)
                    componentCellsInBlock[j] = cell

                # Set region for DerivedShape component - there should be a max of one per block
                for j, component in enumerate(block):
                    if isinstance(component, armi.reactor.components.DerivedShape):
                        componentCellsInBlock[j].region = ~openmc.Union([blockCell.region for blockCell in componentCellsInBlock if blockCell.region is not None])

                componentCellsInAssembly += componentCellsInBlock
            assemblyUniverse = openmc.Universe(name=assembly.name, cells=componentCellsInAssembly)

            # Place the assembly in the correct place in the core
            if core.geomType == armi.reactor.geometry.GeomType.HEX:
                ringIndices = cartesianToRing(assembly.spatialLocator.getCompleteIndices()[0:2])
                assemblyRings[ringIndices[0]][ringIndices[1]] = assemblyUniverse
            else:
                raise TypeError("Unsupported geometry type")

        # Create core lattice
        if core.geomType == armi.reactor.geometry.GeomType.HEX:
            lattice = openmc.HexLattice()
            lattice.pitch = [core.getAssemblyPitch()]
            lattice.center = (0,0)
            assemblyRings.reverse()
            lattice.universes = assemblyRings
            lattice.outer = emptyUniverse
        else:
            raise TypeError("Unsupported geometry type")

        # Create root universe
        rootCell = openmc.Cell(name="rootCell", fill=lattice, region= +boundingCylinderBottomPlane & -boundingCylinderTopPlane & -boundingCylinder)
        rootUniverse = openmc.Universe(cells=[rootCell])

        # Write geometry to xml
        geom = openmc.Geometry(rootUniverse)

        materials.export_to_xml()
        geom.export_to_xml()

    def writeSettings(self):
        """Write the openmc settings input file."""
        settings = openmc.Settings()
        settings.run_mode = 'eigenvalue'
        point = openmc.stats.Point(xyz=(0.0, 0.0, self.r.core[0].getHeight()/2))
        settings.src = openmc.Source(space=point)
        settings.batches = 10
        settings.inactive = 1
        settings.particles = 100
        settings.generations_per_batch = 1
        settings.temperature = {'method': 'interpolation', 'default': 350.0}
        settings.output = {'tallies': True, 'summary': False}
        entropyMesh = openmc.RegularMesh()
        bbWidth = self.r.core.getBoundingCircleOuterDiameter()/2
        bbHeight = max([assembly.getHeight() for assembly in self.r.core])
        entropyMesh.lower_left = [-bbWidth, -bbWidth, 0]
        entropyMesh.upper_right = [bbWidth, bbWidth, bbHeight]
        entropyMesh.dimension = (8,8,8)
        settings.entropy_mesh = entropyMesh
        settings.export_to_xml()

    def writeTallies(self):
        """Write the openmc tallies input file."""
        tallies = openmc.Tallies()
        dummyTally = openmc.Tally()
        dummyTally.scores = ['fission']
        tallies.append(dummyTally)
        tallies.export_to_xml()

    def writePlots(self):
        """Write the openmc plots input file."""
        plot = openmc.Plot()
        plot.basis = 'xy'
        plot.filename = self.r.getName()
        plot.width = (30, 30) #(self.r.core.getBoundingCircleOuterDiameter(), self.r.core.getBoundingCircleOuterDiameter())
        plot.pixels = (1000, 1000)
        plot.origin = (0.0, 0.0, 150.0)
        plot.color_by = 'material'
        #geometry = openmc.Geometry.from_xml('geometry.xml')
        plot.colorize(openmc.Geometry.from_xml('geometry.xml'))
        #plot.colors = {'dummy': 'red'}
        plots = openmc.Plots([plot])
        plots.export_to_xml()

