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
import os

import numpy as np

import armi
from armi import runLog
from armi.nucDirectory import nuclideBases
from armi.utils import hexagon
from armi.reactor.converters.blockConverters import MultipleComponentMerger
from armi.reactor.components import basicShapes, complexShapes
from armi.reactor import systemLayoutInput
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

        self.blockFilterCells = None

    def write(self):
        """Wrapper that writes all OpenMC input files"""
        self.writeGeometry()
        self.writeSettings()
        self.writeTallies()

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

        def buildCellRegion(component, bottomPlane, topPlane, origin=(0.0,0.0), outsideBuffer=0.0):
            """Build region based on shape"""

            # Circle
            if isinstance(component, basicShapes.Circle):
                innerCylinder = openmc.ZCylinder(x0=origin[0], y0=origin[1], r=component.getDimension("id")/2)
                outerCylinder = openmc.ZCylinder(x0=origin[0], y0=origin[1], r=component.getDimension("od")/2 + outsideBuffer)
                region = +bottomPlane & -topPlane & +innerCylinder & -outerCylinder
                return region

            # Hexagon
            if isinstance(component, basicShapes.Hexagon):
                innerHexPrism = openmc.model.hexagonal_prism(edge_length=component.getDimension("ip")/3**.5, orientation='x', origin=origin)
                outerHexPrism = openmc.model.hexagonal_prism(edge_length=component.getDimension("op")/3**.5 + outsideBuffer, orientation='x', origin=origin)
                region = +bottomPlane & -topPlane & ~innerHexPrism & outerHexPrism
                return region

            # Rectangle
            if isinstance(component, basicShapes.Rectangle):
                innerRectPrism = openmc.model.rectangular_prism(width=component.getDimension("widthInner"), height=component.getDimension("lengthInner"), origin=origin) # Check that width/height aren't flipped
                outerRectPrism = openmc.model.rectangular_prism(width=component.getDimension("widthOuter")+outsideBuffer, height=component.getDimension("lengthOuter")+outsideBuffer, origin=origin)
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

        def buildCell(component, material, block, bottomPlane, topPlane, origin=(0.0,0.0), outsideBuffer=0.0):
            if component.getDimension("mult")==1:
                        cell = openmc.Cell(name=component.getName(),
                                           fill=material,
                                           region=buildCellRegion(component, bottomPlane=blockBottomPlane, topPlane=blockTopPlane, origin=origin, outsideBuffer=outsideBuffer))
            else:
                if block.hasPinPitch():
                    latticePitch = block.getPinPitch()
                else:
                    latticePitch = 1.2*max([component.getBoundingCircleOuterDiameter() for component in block.getComponents() if component.getDimension("mult")!=1])

                cellRegions = []
                for location in component.spatialLocator:
                    if self.r.core.geomType == armi.reactor.geometry.GeomType.HEX:
                        origin = (location[0]*latticePitch+location[1]*.5*latticePitch, location[1]*3**.5*latticePitch)
                    elif self.r.core.geomType == armi.reactor.geometry.GeomType.CARTESIAN:
                        origin = (location[0]*latticePitch[0], location[1]*latticePitch[1])
                    else:
                        raise TypeError("Unsupported geometry type")

                    cellRegions.append(buildCellRegion(component, bottomPlane=blockBottomPlane, topPlane=blockTopPlane, origin=origin))
                cell = openmc.Cell(name=component.getName(),
                                   fill=material,
                                   region=openmc.Union(cellRegions))
            return cell

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

            # OpenMC's ring system starts at the top and goes clockwise
            lenRing = max([ring*6, 1])
            pos = (lenRing+ring-pos)%lenRing
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

        runLog.info("Writing geometry, materials, and plots...")
        core = self.r.core

        materials = openmc.Materials()
        emptyMaterial = openmc.Material()
        plotColors = dict()
        colorLookup = {"HT9": "steelblue",
                       "UZr": "green",
                       "UO2": "green",
                       "UraniumOxide": "green",
                       "Sodium": "antiquewhite",
                       "Custom": "purple",
                       "B4C": "black",
                       "SaturatedWater": "blue"}

        numRings = core.numRings
        boundingCellBottomPlane = openmc.ZPlane(z0=0.0, boundary_type='vacuum')
        boundingCellTopPlane = openmc.ZPlane(z0=max([assembly.getHeight() for assembly in core]), boundary_type='vacuum')

        if core.geomType == armi.reactor.geometry.GeomType.HEX:

            boundingCylinderRadius = core.getCoreRadius()
            boundingCylinder = openmc.ZCylinder(r=boundingCylinderRadius, boundary_type='vacuum')

            if core.symmetry.domain==armi.reactor.geometry.DomainType.FULL_CORE:
                boundingCellRegion = -boundingCylinder & +boundingCellBottomPlane & -boundingCellTopPlane

            if core.symmetry.domain==armi.reactor.geometry.DomainType.THIRD_CORE:
                armi.reactor.converters.geometryConverters.EdgeAssemblyChanger().addEdgeAssemblies(core)
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
            assemblyLatticeIndices = buildRings(numRings, emptyUniverse)

        elif core.geomType == armi.reactor.geometry.GeomType.CARTESIAN:
            boundingCylinderRadius = core.getAssemblyPitch()[0]*core.numRings*2**.5
            boundingCylinder = openmc.ZCylinder(r=boundingCylinderRadius, boundary_type='vacuum')

            if core.symmetry.domain==armi.reactor.geometry.DomainType.QUARTER_CORE:
                periodicPlane0 = openmc.Plane(a=0.0,b=1.0,c=0.0,d=0.0)
                periodicPlane1 = openmc.Plane(a=1.0,b=0.0,c=0.0,d=0.0)
                periodicPlane0.boundary_type = 'periodic'
                periodicPlane1.boundary_type = 'periodic'
                boundingCellRegion = -boundingCylinder & +periodicPlane0 & +periodicPlane1 & +boundingCellBottomPlane & -boundingCellTopPlane

            emptyCellRegion = -boundingCylinder & +boundingCellBottomPlane & -boundingCellTopPlane
            emptyCell = openmc.Cell(region = emptyCellRegion)

            emptyUniverse = openmc.Universe(cells=[emptyCell])
            assemblyLatticeIndices = []
            for ring in range(numRings):
                assemblyLatticeIndices.append([emptyUniverse]*numRings)
        else:
            raise TypeError("Unsupported geometry type")

        self.blockFilterCells = []
        for assembly in core:
            # Write a universe for each assembly
            assemblyUniverse = openmc.Universe(name=assembly.name)

            # Create ZPlanes between blocks
            z0 = 0.0
            ZPlanes=[openmc.ZPlane(z0=z0, boundary_type='vacuum')]
            for block in assembly:
                z0 += block.getHeight()
                ZPlanes.append(openmc.ZPlane(z0=z0, boundary_type='transmission'))
            ZPlanes[-1].boundary_type = 'vacuum' # Reset top plane boundary condition to vacuum

            blockCellsInAssembly = []
            for i, block in enumerate(assembly):
                blockBottomPlane = ZPlanes[i]
                blockTopPlane = ZPlanes[i+1]

                blockWithoutHelices = blendHelixComponentsIntoCoolant(block)

                blockHasDerivedShapeComponent = any([isinstance(component, armi.reactor.components.DerivedShape) for component in blockWithoutHelices])
                # Get DerivedShape component. We need to set its region last
                if blockHasDerivedShapeComponent:
                    derivedShapeComponent = [component for component in blockWithoutHelices if isinstance(component, armi.reactor.components.DerivedShape)][0]
                    derivedShapeComponentMaterial = buildComponentMaterial(derivedShapeComponent)
                    if derivedShapeComponentMaterial is not None:
                        materials.append(derivedShapeComponentMaterial)
                        plotColors[derivedShapeComponentMaterial.id] = colorLookup[derivedShapeComponent.material.name]
                    blockMinusDerivedShape = [component for component in blockWithoutHelices if not isinstance(component, armi.reactor.components.DerivedShape)]
                else:
                    blockMinusDerivedShape = blockWithoutHelices
                    derivedShapeComponentMaterial = None

                componentCellsInBlock = []
                # Divide all components into groups with same mult
                """
                Better way to get multGroups:
                mults = {c.getDimension("mult") for c in blockMinusDerivedShape.iterComponents()}
                """

                multGroups = dict()
                for component in blockMinusDerivedShape:
                    mult = int(component.getDimension("mult"))
                    if mult not in multGroups:
                        multGroups[mult] = []
                    multGroups[mult].append(component)

                # If any components have mult>1, we need a cell filled with a lattice of them
                if len(multGroups)==2:
                    if block.hasPinPitch():
                        latticePitch = block.getPinPitch()
                    else:
                        latticePitch = 1.2*max([[component.getBoundingCircleOuterDiameter() for component in multGroups[multGroup]] for multGroup in multGroups if multGroup!=1])[0]

                    if core.geomType == armi.reactor.geometry.GeomType.HEX:
                        totalSlots = [mult for mult in multGroups if mult>1][0]

                        # Determine number of rings -> solve mult=1+6*(nRings-1)(nRings)/2
                        nRings = math.ceil(.5*(1+(1+4/3*(totalSlots-1))**.5))
                        blockLattice = openmc.HexLattice()
                        blockLattice.pitch = [latticePitch]
                        blockLattice.orientation = 'x'
                        blockLattice.center = (0,0)
                        blockLatticeIndices = buildRings(nRings, emptyUniverse)

                    elif core.geomType == armi.reactor.geometry.GeomType.CARTESIAN:
                        boundingIndices = block.spatialGrid.getIndexBounds()
                        blockLatticeIndicesOffset = (boundingIndices[0][0],boundingIndices[1][0])
                        blockLattice = openmc.RectLattice()
                        blockLattice.pitch = latticePitch
                        blockLattice.lower_left = (latticePitch[0]*boundingIndices[0][0],latticePitch[1]*boundingIndices[1][0])
                        blockLatticeIndices = []
                        for ring in range(boundingIndices[1][1]-boundingIndices[1][0]-1):
                            blockLatticeIndices.append([emptyUniverse]*(boundingIndices[0][1]-boundingIndices[0][0]-1))

                    for mult in multGroups:
                        if mult>1:
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
                            if blockHasDerivedShapeComponent:
                                # Fill unused space in lattice with derivedShapeComponentMaterial
                                derivedShapeComponentLatticeCell = openmc.Cell(name=derivedShapeComponent.getName(),
                                                                               fill=derivedShapeComponentMaterial,
                                                                               region=~openmc.Union([cell.region for cell in componentCellsInMultGroupUniverse]) & +blockBottomPlane & -blockTopPlane)
                                componentCellsInMultGroupUniverse.append(derivedShapeComponentLatticeCell)

                            multGroupUniverse = openmc.Universe(name="multGroup"+str(mult), cells=componentCellsInMultGroupUniverse)

                            """if spatialLocator is a multiIndexLocation, use it. Otherwise, assume all blockLatticeUniverses are multGroupUniverse"""
                            multGroupIndices = multGroups[mult][0].spatialLocator
                            if not isinstance(multGroupIndices, armi.reactor.grids.MultiIndexLocation):
                                # Assume every universe in blockLattice is a multGroupUniverse
                                for rowIndex, row in enumerate(blockLatticeIndices):
                                    for colIndex in range(len(row)):
                                        blockLatticeIndices[rowIndex][colIndex] = multGroupUniverse
                            else:
                                if core.geomType == armi.reactor.geometry.GeomType.CARTESIAN:
                                    for location in multGroupIndices:
                                        blockLatticeIndices[location[0]][location[1]] = multGroupUniverse
                                if core.geomType == armi.reactor.geometry.GeomType.HEX:
                                    for location in multGroupIndices:
                                        ringIndices = cartesianToRing(location[0:2])
                                        blockLatticeIndices[ringIndices[0]][ringIndices[1]] = multGroupUniverse

                    if core.geomType == armi.reactor.geometry.GeomType.HEX:
                        blockLatticeIndices.reverse()
                    blockLattice.universes = blockLatticeIndices

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

                    remainingComponents = multGroups[1]
                else:
                    remainingComponents = [component for component in blockMinusDerivedShape]

                # If component cell surfaces are perfectly coincident with lattice boundaries, it's possible for a particle to get lost by hitting a corner just right.
                # Fix by adding a buffer to the outside of the largest component (e.g. intercoolant). Does not change effective geometry in the simulation.
                # This code should be migrated to proper plugin geometryTransformations because we do not undo this transformation anywhere.
                remainingComponentOuterDiameters = [component.getBoundingCircleOuterDiameter() for component in remainingComponents]
                largestRemainingComponent = remainingComponents[max(range(len(remainingComponentOuterDiameters)), key=remainingComponentOuterDiameters.__getitem__)]
                componentMaterial = buildComponentMaterial(largestRemainingComponent)
                if componentMaterial is not None:
                        materials.append(componentMaterial)
                        plotColors[componentMaterial.id] = colorLookup[component.material.name]
                cell = buildCell(largestRemainingComponent,
                                 material=componentMaterial,
                                 block=block,
                                 bottomPlane=blockBottomPlane,
                                 topPlane=blockTopPlane,
                                 origin=(0.0,0.0),
                                 outsideBuffer=0.01)
                componentCellsInBlock.append(cell)
                remainingComponents.remove(largestRemainingComponent)

                for component in remainingComponents:
                    componentMaterial = buildComponentMaterial(component)
                    if componentMaterial is not None:
                        materials.append(componentMaterial)
                        plotColors[componentMaterial.id] = colorLookup[component.material.name]
                    cell = buildCell(component,
                                     material=componentMaterial,
                                     block=block,
                                     bottomPlane=blockBottomPlane,
                                     topPlane=blockTopPlane,
                                     origin=(0.0,0.0),
                                     outsideBuffer=0.0)
                    componentCellsInBlock.append(cell)

                if blockHasDerivedShapeComponent:
                    # Set region for DerivedShape component - there should be a max of one per block
                    derivedShapeComponentCell = openmc.Cell(name=derivedShapeComponent.getName(),
                                                            fill = derivedShapeComponentMaterial,
                                                            region = ~openmc.Union([blockCell.region for blockCell in componentCellsInBlock]) & +blockBottomPlane & -blockTopPlane)
                    componentCellsInBlock.append(derivedShapeComponentCell)
                else:
                    # If there's no DerivedShape component, add an empty cell in unoccupied region in case ComponentCellsInBlock doesn't fill its assemblyLattice cell
                    emptyComponentCell = openmc.Cell(name='emptyComponent',
                                                          region = ~openmc.Union([blockCell.region for blockCell in componentCellsInBlock]) & +blockBottomPlane & -blockTopPlane)
                    componentCellsInBlock.append(emptyComponentCell)
                blockUniverse = openmc.Universe(cells = componentCellsInBlock)
                blockUniverseCell = openmc.Cell(name=block.getName(),
                                                fill=blockUniverse,
                                                region=openmc.Union([blockCell.region for blockCell in componentCellsInBlock]))

                blockCellsInAssembly.append(blockUniverseCell)
                self.blockFilterCells.append(blockUniverseCell)
            assemblyUniverse.add_cells(blockCellsInAssembly)

            # Place the assembly in the correct place in the core
            if core.geomType == armi.reactor.geometry.GeomType.HEX:
                ringIndices = cartesianToRing(assembly.spatialLocator.getCompleteIndices()[0:2])
                assemblyLatticeIndices[ringIndices[0]][ringIndices[1]] = assemblyUniverse
            elif core.geomType == armi.reactor.geometry.GeomType.CARTESIAN:
                spatialIndices = assembly.spatialLocator.getCompleteIndices()
                assemblyLatticeIndices[spatialIndices[0]][spatialIndices[1]] = assemblyUniverse
            else:
                raise TypeError("Unsupported geometry type")

        # Create core lattice
        if core.geomType == armi.reactor.geometry.GeomType.HEX:
            lattice = openmc.HexLattice()
            lattice.center = (0,0)

            lattice.pitch = [core.getAssemblyPitch()]
        elif core.geomType == armi.reactor.geometry.GeomType.CARTESIAN:
            lattice = openmc.RectLattice()
            if core.symmetry.domain==armi.reactor.geometry.DomainType.QUARTER_CORE:
                lattice.lower_left = (0.0, 0.0)
            else:
                lattice.lower_left = (-core.getAssemblyPitch()[0]*(core.numRings/2),-core.getAssemblyPitch()[1]*(core.numRings/2))
            lattice.pitch = [p for p in core.getAssemblyPitch()]
        else:
            raise TypeError("Unsupported geometry type")

        assemblyLatticeIndices.reverse()
        lattice.universes = assemblyLatticeIndices
        lattice.outer = emptyUniverse

        # Create root universe
        rootCell = openmc.Cell(name="rootCell", fill=lattice, region = boundingCellRegion)
        rootUniverse = openmc.Universe(cells=[rootCell])

        # Write geometry to xml
        geometry = openmc.Geometry(rootUniverse)
        geometry.merge_surfaces=True # merge redundant surfaces

        # Write plots xml file
        plot = openmc.Plot()
        plot.basis = 'xy'
        plot.filename = self.r.getName()
        plot.width = (boundingCylinderRadius, boundingCylinderRadius)
        plot.pixels = (100, 100)
        plot.origin = (0.0, 0.0, 20.0)
        plot.color_by = 'material'
        plot.colors = plotColors
        plots = openmc.Plots([plot])

        materials.export_to_xml()
        plots.export_to_xml()
        geometry.export_to_xml()

    def writeSettings(self):
        """Write the OpenMC settings input file."""
        runLog.info("Writing settings...")
        settings = openmc.Settings()
        settings.run_mode = 'eigenvalue'
        bbHeight = max([assembly.getHeight() for assembly in self.r.core])
        if self.r.core.geomType == armi.reactor.geometry.GeomType.HEX:
            boundingCylinderRadius = self.r.core.getCoreRadius()
        elif self.r.core.geomType == armi.reactor.geometry.GeomType.CARTESIAN:
            boundingCylinderRadius = self.r.core.getAssemblyPitch()[0]*self.r.core.numRings*2**.5
        point = openmc.stats.Box(lower_left=(-boundingCylinderRadius,-boundingCylinderRadius,0.0), 
                                 upper_right=(boundingCylinderRadius,boundingCylinderRadius,bbHeight),
                                 only_fissionable=True)
        settings.source = openmc.IndependentSource(space=point)
        settings.batches = self.options.nBatches
        settings.inactive = self.options.nInactiveBatches
        settings.particles = self.options.nParticles
        settings.generations_per_batch = 1
        settings.temperature = {'method': 'interpolation', 'default': 350.0}
        settings.output = {'tallies': False, 'summary': False}
        settings.verbosity = self.options.openmcVerbosity
        entropyMesh = openmc.RegularMesh()
        bbWidth = boundingCylinderRadius
        entropyMesh.lower_left = [-bbWidth, -bbWidth, 0]
        entropyMesh.upper_right = [bbWidth, bbWidth, bbHeight]
        entropyMesh.dimension = tuple(self.options.entropyMeshDimension)
        settings.entropy_mesh = entropyMesh
        settings.export_to_xml()

    def writeTallies(self):
        """Write the OpenMC tallies input file."""
        runLog.info("Writing tallies...")
        if self.options.Tallies is None:
            tallies = openmc.Tallies()
        else:
            tallies = self.options.Tallies

        if self.r.core.geomType == armi.reactor.geometry.GeomType.HEX:
            bbWidth = self.r.core.getCoreRadius()
        elif self.r.core.geomType == armi.reactor.geometry.GeomType.CARTESIAN:
            bbWidth = self.r.core.getAssemblyPitch()[0]*self.r.core.numRings

        tallyMesh = openmc.RegularMesh()
        bbHeight = max([assembly.getHeight() for assembly in self.r.core])
        tallyMesh.lower_left = [-bbWidth, -bbWidth, 0]
        tallyMesh.upper_right = [bbWidth, bbWidth, bbHeight]
        tallyMesh.dimension = tuple(self.options.tallyMeshDimension)
        meshFilter = openmc.MeshFilter(mesh=tallyMesh, filter_id=101)

        energyGroupStructure = parseEnergyGroupStructure(energyGroups.getGroupStructure(self.options.groupStructure))
        energyFilter = openmc.EnergyFilter(energyGroupStructure, filter_id=102)
        blockFilter = openmc.CellFilter(bins = self.blockFilterCells, filter_id=103)

        fissionTally = openmc.Tally(101, name="fission rate")
        fissionTally.scores = ['fission']
        fissionTally.nuclides = ['U235', 'U238']
        fissionTally.filters = [meshFilter]
        tallies.append(fissionTally)

        blockFluxTally = openmc.Tally(102, name="block filter multigroup flux")
        blockFluxTally.scores = ['flux']
        blockFluxTally.filters = [blockFilter, energyFilter]
        tallies.append(blockFluxTally)

        meshFluxTally = openmc.Tally(103, name="mesh filter multigroup flux")
        meshFluxTally.scores = ['flux']
        meshFluxTally.filters = [meshFilter, energyFilter]
        tallies.append(meshFluxTally)

        powerTally = openmc.Tally(104, name="power")
        powerTally.scores = ['heating-local']
        powerTally.filters = [blockFilter]
        tallies.append(powerTally)

        absorptionTally = openmc.Tally(105, name="absorption")
        absorptionTally.scores = ['absorption']
        absorptionTally.filters = [blockFilter]
        tallies.append(absorptionTally)

        tallies.export_to_xml()

def parseEnergyGroupStructure(energyGroupStructure):
    """Convert ARMI group structure to openmc group structure"""
    energyGroupStructure.append(0.0)
    energyGroupStructure.reverse()
    return energyGroupStructure

