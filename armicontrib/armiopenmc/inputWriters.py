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
import math

import armi
from armi.nucDirectory import nuclideBases
from armi.utils import hexagon
from armi.reactor.converters.blockConverters import MultipleComponentMerger
from armi.reactor.components import basicShapes, complexShapes
from armi.physics.neutronics import energyGroups

import openmc

class OpenMCWriter:
    """
    Write OpenMC data using the openmc python api.
    """

    def __init__(self, reactor, options):
        """Build the writer"""

        self.r = reactor
        self.options = options

    def writeGeometry(self):
        """Write the openmc geometry, materials, and plots input file."""

        def blendHelixComponentsIntoCoolant(block, solventName='coolant'):
            """OpenMC doesn't support helixes, so blend them all into coolant"""
            helixComponentNames = []
            for component in block:
                if isinstance(component, complexShapes.Helix):
                    helixComponentNames.append(component.getName()) 
            if len(helixComponentNames)>0:
                block = MultipleComponentMerger(sourceBlock=block, soluteNames=helixComponentNames, solventName=solventName).convert()
            return block

        def buildCellRegion(component, bottomPlane, topPlane):
            """Build region based on shape"""

            # Circle
            if isinstance(component, basicShapes.Circle):
                innerCylinder = openmc.ZCylinder(r=component.getDimension("id")/2)
                outerCylinder = openmc.ZCylinder(r=component.getDimension("od")/2)
                region = +bottomPlane & -topPlane & +innerCylinder & -outerCylinder
                return region

            # Hexagon
            if isinstance(component, basicShapes.Hexagon):
                innerHexPrism = openmc.model.hexagonal_prism(edge_length=component.getDimension("ip")/3**.5,
                                                             orientation='x')
                outerHexPrism = openmc.model.hexagonal_prism(edge_length=component.getDimension("op")/3**.5,
                                                             orientation='x')
                region = +bottomPlane & -topPlane & ~innerHexPrism & outerHexPrism
                return region

            # Rectangle
            if isinstance(component, basicShapes.Rectangle):
                innerRectPrism = openmc.model.rectangular_prism(width=component.getDimension("widthInner"),
                                                                height=component.getDimension("lengthInner")) # Check that width/height aren't flipped
                outerRectPrism = openmc.model.rectangular_prism(width=component.getDimension("widthOuter"),
                                                                height=component.getDimension("lengthOuter"))
                region = +bottomPlane & -topPlane & ~innerRectPrism & outerRectPrism
                return region

            # DerivedShape
            if isinstance(component, armi.reactor.components.DerivedShape):
                # DerivedShape is supported, but we need to set DerivedShape regions after all others in the block.
                # We should never get here.
                warnings.warn("DerivedShape not supported in OpenMCWriter.buildCellRegion")
                return

            # Helix
            if isinstance(component, complexShapes.Helix):
                # Note: Helix components are automatically blended into coolant. We should never get here.
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

        def buildRings(numRings, universe):
            """Assemble rings for use in openmc.HexLattice"""
            # Note: Need to rings.reverse() before assigning to lattice.universes. Not doing it here keeps indexing easy.
            rings = [None]*numRings
            rings[0] = [universe]
            ringLength = 6
            for i in range(1, numRings):
                rings[i] = [universe]*ringLength
                ringLength += 6
            return rings

        def buildComponentMaterial(component):
            """Build OpenMC material for ARMI component"""
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
        plotColors = dict()
        colorLookup = {"HT9": "steelblue",
                       "UZr": "green",
                       "Sodium": "antiquewhite",
                       "B4C": "black"}

        if core.geomType == armi.reactor.geometry.GeomType.HEX:

            numRings = core.numRings
            boundingCylinderRadius = core.getCoreRadius()
            boundingCylinder = openmc.ZCylinder(r=boundingCylinderRadius, boundary_type='vacuum')

            if core.symmetry.domain==armi.reactor.geometry.DomainType.FULL_CORE:
                boundingCellBottomPlane = openmc.ZPlane(z0=0.0, boundary_type='vacuum')
                boundingCellTopPlane = openmc.ZPlane(z0=max([assembly.getHeight() for assembly in core]), boundary_type='vacuum')
                boundingCellRegion = -boundingCylinder & +boundingCellBottomPlane & -boundingCellTopPlane

            if core.symmetry.domain==armi.reactor.geometry.DomainType.THIRD_CORE:
                armi.reactor.converters.geometryConverters.EdgeAssemblyChanger().addEdgeAssemblies(core)
                boundingCellBottomPlane = openmc.ZPlane(z0=0.0, boundary_type='vacuum')
                boundingCellTopPlane = openmc.ZPlane(z0=max([assembly.getHeight() for assembly in core]), boundary_type='vacuum')
                
                periodicPlane0 = openmc.Plane(a=0.0,b=1.0,c=0.0,d=0.0)
                periodicPlane1 = openmc.Plane(a=3.0**.5,b=1.0,c=0.0,d=0.0)

                if core.symmetry.boundary==armi.reactor.geometry.BoundaryType.PERIODIC:
                    periodicPlane0.boundary_type = 'periodic'
                    periodicPlane1.boundary_type = 'periodic'

                if core.symmetry.boundary==armi.reactor.geometry.BoundaryType.REFLECTIVE:
                    periodicPlane0.boundary_type = 'reflective'
                    periodicPlane1.boundary_type = 'reflective'
                boundingCellRegion = -boundingCylinder & +periodicPlane0 & +periodicPlane1 & +boundingCellBottomPlane & -boundingCellTopPlane

            emptyCellRegion = -boundingCylinder & +boundingCellBottomPlane & -boundingCellTopPlane
            emptyCell = openmc.Cell(region = emptyCellRegion)

            # Create blank list of assembly universes for the lattice - we will fill in universes individually later
            emptyUniverse = openmc.Universe(cells=[emptyCell])
            assemblyRings = buildRings(numRings, emptyUniverse)

        else:
            raise TypeError("Unsupported geometry type")

        for assembly in core:
            # Write a universe for each assembly
            print("Working on assembly " + str(assembly.getName()))
            assemblyUniverse = openmc.Universe(name=assembly.name)

            # Create ZPlanes between blocks
            z0 = 0.0
            ZPlanes=[openmc.ZPlane(z0=z0, boundary_type='vacuum')]
            for block in assembly:
                z0 += block.getHeight()
                ZPlanes.append(openmc.ZPlane(z0=z0, boundary_type='transmission'))
            ZPlanes[-1].boundary_type = 'vacuum' # Reset top plane boundary condition to vacuum

            componentCellsInAssembly = []
            for i, block in enumerate(assembly):
                blockBottomPlane = ZPlanes[i]
                blockTopPlane = ZPlanes[i+1]

                blockWithoutHelices = blendHelixComponentsIntoCoolant(block)

                # Get DerivedShape component. We need to set its region last
                derivedShapeComponent = [component for component in blockWithoutHelices if isinstance(component, armi.reactor.components.DerivedShape)][0]
                derivedShapeComponentMaterial = buildComponentMaterial(derivedShapeComponent)
                if derivedShapeComponentMaterial is not None:
                    materials.append(derivedShapeComponentMaterial)
                    plotColors[derivedShapeComponentMaterial.id] = colorLookup[derivedShapeComponent.material.name]
                blockMinusDerivedShape = [component for component in blockWithoutHelices if not isinstance(component, armi.reactor.components.DerivedShape)]

                componentCellsInBlock = []
                # Divide all components into groups with same mult
                multGroups = dict()
                for component in blockMinusDerivedShape:
                    mult = int(component.getDimension("mult"))
                    if mult not in multGroups:
                        multGroups[mult] = []
                    multGroups[mult].append(component)

                # If any components have mult>1, we need a cell filled with a lattice of them
                if len(multGroups)>1:
                    for mult in multGroups:
                        if mult>1:
                            # Determine number of rings -> solve mult=1+6*(nRings-1)(nRings)/2
                            nRings = math.ceil(.5*(1+(1+4/3*(mult-1))**.5))
                            # NOTE: Currently supporting only 1 multGroup with mult>1
                            # NOTE: Need better latticePitch
                            if block.hasPinPitch():
                                latticePitch = block.getPinPitch()
                            else:
                               latticePitch = 1.2*max([component.getBoundingCircleOuterDiameter() for component in multGroups[mult]])
                            componentCellsInMultGroupUniverse = []

                            for component in multGroups[mult]:
                                componentMaterial = buildComponentMaterial(component)
                                cell = openmc.Cell(name=component.getName(),
                                                   fill=componentMaterial,
                                                   region=buildCellRegion(component, bottomPlane=blockBottomPlane, topPlane=blockTopPlane))
                                if componentMaterial is not None:
                                    materials.append(componentMaterial)
                                    plotColors[componentMaterial.id] = colorLookup[component.material.name]
                                componentCellsInMultGroupUniverse.append(cell)
                            # Fill unused space in lattice with derivedShapeComponentMaterial
                            derivedShapeComponentLatticeCell = openmc.Cell(name=derivedShapeComponent.getName(),
                                                                    fill=derivedShapeComponentMaterial,
                                                                    region=~openmc.Union([cell.region for cell in componentCellsInMultGroupUniverse]) & +blockBottomPlane & -blockTopPlane)
                            componentCellsInMultGroupUniverse.append(derivedShapeComponentLatticeCell)
                            multGroupUniverse = openmc.Universe(name="multGroup"+str(mult), cells=componentCellsInMultGroupUniverse)
                    # Set blockLattice
                    blockLattice = openmc.HexLattice()
                    blockLattice.pitch = [latticePitch]
                    blockLattice.orientation = 'x'
                    blockLattice.center = (0,0)
                    blockLatticeRings = buildRings(nRings, multGroupUniverse)
                    blockLatticeRings.reverse()
                    blockLattice.universes = blockLatticeRings
                    blockLatticeOuterCell = openmc.Cell(region=+blockBottomPlane & -blockTopPlane & -boundingCylinder,
                                                        fill=derivedShapeComponentMaterial)
                    blockLatticeOuterUniverse = openmc.Universe(cells=[blockLatticeOuterCell])
                    blockLattice.outer = blockLatticeOuterUniverse
                    # Need cell to fill with lattice
                    # Get smallest mult 1 component - blockLatticeCell will have region inside it
                    mult1ComponentInnerDiameters = [component.getCircleInnerDiameter() for component in multGroups[1]]
                    smallestMult1Component = multGroups[1][min(range(len(mult1ComponentInnerDiameters)), key=mult1ComponentInnerDiameters.__getitem__)]
                    if isinstance(smallestMult1Component, basicShapes.Circle):
                        innerCylinder = openmc.ZCylinder(r=smallestMult1Component.getDimension("id")/2-.05)
                        blockLatticeCellRegion = -innerCylinder & +blockBottomPlane & -blockTopPlane
                    elif isinstance(smallestMult1Component, basicShapes.Hexagon):
                        innerHexPrism = openmc.model.hexagonal_prism(edge_length=smallestMult1Component.getDimension("ip")/3**.5-.05, orientation='x')
                        blockLatticeCellRegion = innerHexPrism & +blockBottomPlane & -blockTopPlane
                    elif isinstance(smallestMult1Component, basicShapes.Rectangle):
                        innerRectPrism = openmc.model.rectangular_prism(width=smallestMult1Component.getDimension("widthInner")-.05,
                                                                        height=smallestMult1Component.getDimension("lengthInner")-.05)
                        blockLatticeCellRegion = innerRectPrism & +blockBottomPlane & -blockTopPlane
                    else:
                        raise NotImplementedError("Shape type not supported yet")
                    blockLatticeCell = openmc.Cell(name="blockLattice",
                                                   fill = blockLattice,
                                                   region = blockLatticeCellRegion)
                    componentCellsInBlock.append(blockLatticeCell)

                for component in multGroups[1]:
                    componentMaterial = buildComponentMaterial(component)
                    cell = openmc.Cell(name=component.getName(),
                                       fill=componentMaterial,
                                       region=buildCellRegion(component, bottomPlane=blockBottomPlane, topPlane=blockTopPlane))
                    if componentMaterial is not None:
                        materials.append(componentMaterial)
                        plotColors[componentMaterial.id] = colorLookup[component.material.name]
                    componentCellsInBlock.append(cell)

                # Set region for DerivedShape component - there should be a max of one per block
                derivedShapeComponentCell = openmc.Cell(name=derivedShapeComponent.getName(),
                                                        fill = derivedShapeComponentMaterial,
                                                        region = ~openmc.Union([blockCell.region for blockCell in componentCellsInBlock]) & +blockBottomPlane & -blockTopPlane)
                componentCellsInBlock.append(derivedShapeComponentCell)

                componentCellsInAssembly += componentCellsInBlock
            assemblyUniverse.add_cells(componentCellsInAssembly)

            # Place the assembly in the correct place in the core
            if core.geomType == armi.reactor.geometry.GeomType.HEX:
                ringIndices = cartesianToRing(assembly.spatialLocator.getCompleteIndices()[0:2])
                assemblyRings[ringIndices[0]][ringIndices[1]] = assemblyUniverse
            else:
                raise TypeError("Unsupported geometry type")

        # Create core lattice
        if core.geomType == armi.reactor.geometry.GeomType.HEX:
            lattice = openmc.HexLattice()
            lattice.pitch = [core.getAssemblyPitch()+0.0001]
            lattice.center = (0,0)
            assemblyRings.reverse()
            lattice.universes = assemblyRings
            lattice.outer = emptyUniverse
        else:
            raise TypeError("Unsupported geometry type")

        # Create root universe
        rootCell = openmc.Cell(name="rootCell", fill=lattice, region = boundingCellRegion)
        rootUniverse = openmc.Universe(cells=[rootCell])

        # Write geometry to xml
        geom = openmc.Geometry(rootUniverse)
        geom.merge_surfaces=True # merge redundant surfaces

        # Write plots xml file
        plot = openmc.Plot()
        plot.basis = 'xy'
        plot.filename = self.r.getName()
        plot.width = (boundingCylinderRadius, boundingCylinderRadius) #(self.r.core.getBoundingCircleOuterDiameter(), self.r.core.getBoundingCircleOuterDiameter())
        plot.pixels = (4000, 4000)
        plot.origin = (0.0, 0.0, 20.0)
        plot.color_by = 'material'
        plot.colors = plotColors
        plots = openmc.Plots([plot])

        print("Exporting to xml...")
        materials.export_to_xml()
        plots.export_to_xml()
        geom.export_to_xml()

    def writeSettings(self):
        """Write the OpenMC settings input file."""
        settings = openmc.Settings()
        settings.run_mode = 'eigenvalue'
        boundingCyliderRadius = self.r.core.getCoreRadius()
        point = openmc.stats.Box(lower_left=(-boundingCyliderRadius,-boundingCyliderRadius,0.0), upper_right=(boundingCyliderRadius,boundingCyliderRadius,100.0)) #xyz=(0.0, 0.0, 20.0))#self.r.core[0].getHeight()/2))
        settings.source = openmc.Source(space=point)
        settings.batches = 40
        settings.inactive = 10
        settings.particles = 100000
        settings.generations_per_batch = 1
        settings.temperature = {'method': 'interpolation', 'default': 350.0}
        settings.output = {'tallies': True, 'summary': True}
        settings.verbosity = 7
        entropyMesh = openmc.RegularMesh()
        bbWidth = boundingCyliderRadius
        bbHeight = max([assembly.getHeight() for assembly in self.r.core])
        entropyMesh.lower_left = [-bbWidth, -bbWidth, 0]
        entropyMesh.upper_right = [bbWidth, bbWidth, bbHeight]
        entropyMesh.dimension = (10, 10, 10)
        settings.entropy_mesh = entropyMesh
        settings.export_to_xml()

    def writeTallies(self):
        """Write the OpenMC tallies input file."""
        tallies = openmc.Tallies()
        
        # Fission tally
        fissionTally = openmc.Tally()
        fissionTally.scores = ['fission']
        fissionTally.nuclides = ['U235', 'U238']
        fissionTallyMesh = openmc.RegularMesh()
        bbWidth = self.r.core.getCoreRadius()
        bbHeight = max([assembly.getHeight() for assembly in self.r.core])
        fissionTallyMesh.lower_left = [-bbWidth, -bbWidth, 0]
        fissionTallyMesh.upper_right = [bbWidth, bbWidth, bbHeight]
        fissionTallyMesh.dimension = (1000, 1000, 1)
        fissionTally.filters = [openmc.MeshFilter(mesh=fissionTallyMesh)]
        tallies.append(fissionTally)
        
        # Multigroup Flux tally
        fluxTally = openmc.Tally()
        fluxTally.scores = ['flux']
        fluxTallyMesh = openmc.RegularMesh()
        bbWidth = self.r.core.getCoreRadius()
        bbHeight = max([assembly.getHeight() for assembly in self.r.core])
        fluxTallyMesh.lower_left = [-bbWidth, -bbWidth, 0]
        fluxTallyMesh.upper_right = [bbWidth, bbWidth, bbHeight]
        fluxTallyMesh.dimension = (1000, 1000, 1)
        energyGroupStructure = energyGroups.getGroupStructure("ANL33")
        energyGroupStructure.append(0.0)
        energyGroupStructure.reverse()
        fluxTallyEnergyFilter = openmc.EnergyFilter(energyGroupStructure)
        fluxTally.filters = [fluxTallyEnergyFilter, openmc.MeshFilter(mesh=fluxTallyMesh)]
        tallies.append(fluxTally)
        
        # Power tally
        powerTally = openmc.Tally()
        powerTally.scores = ['heating']
        powerTallyMesh = openmc.RegularMesh()
        bbWidth = self.r.core.getCoreRadius()
        bbHeight = max([assembly.getHeight() for assembly in self.r.core])
        powerTallyMesh.lower_left = [-bbWidth, -bbWidth, 0]
        powerTallyMesh.upper_right = [bbWidth, bbWidth, bbHeight]
        powerTallyMesh.dimension = (1000, 1000, 1)
        powerTally.filters = [openmc.MeshFilter(mesh=powerTallyMesh)]
        tallies.append(powerTally)
        
        tallies.export_to_xml()

