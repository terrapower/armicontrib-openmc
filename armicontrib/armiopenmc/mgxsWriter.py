import numpy as np
import os
import openmc
from armi.nuclearDataIO.xsLibraries import CompxsLibrary, IsotxsLibrary
from armi.nuclearDataIO.cccc import compxs, isotxs

class MGXSWriter:
    def __init__(self, inputFile, inputFormat, regionNames=None):
        self.inputFile = inputFile
        self.inputType = inputFormat
        self.outputFilename = os.path.splitext(inputFile)[0] + ".h5"
        fileExtension = os.path.splitext(inputFile)[1]
        self.regionNames = {'0': 'fuel',
                            '1': 'mox low fuel',
                            '2': 'mox medium fuel',
                            '3': 'mox high fuel',
                            '4': 'moderator',
                            '5': 'guide tube',
                            '6': 'fission chamber'}
        if self.inputType == 'macro':
            if fileExtension == '.ascii':
                self.lib = compxs.readAscii(inputFile)
            elif fileExtension == '':
                self.lib = compxs.readBinary(inputFile)
            elif fileExtension == '.h5':
                self.lib = None
            else:
                raise ValueError("Unknown file extension")
        elif self.inputType == 'micro':
            if fileExtension == '.ascii':
                self.lib = isotxs.readAscii(inputFile)
            elif fileExtension == '':
                self.lib = isotxs.readBinary(inputFile)
            elif fileExtension == '.h5':
                self.lib = None
            else:
                raise ValueError("Unknown file extension")
        else:
            raise ValueError("Unknown cross section type")
    def write(self):
        '''
        Write the mgxs.h5 file used by openmc. Returns the output filename as a string
        If input file is in .h5 format, does nothing.
        '''
        # Compxs precomputed macroscopic cross sections
        if type(self.lib) is CompxsLibrary:
            return self.writeCompxs()
        # Isotxs microscopic cross sections
        if type(self.lib) is IsotxsLibrary:
            return self.writeIsotxs()
        # Openmc's hdf5 formatted cross sections
        if self.lib == None:
            return self.inputFile
    def writeCompxs(self):
        numGroups = self.lib.compxsMetadata["numGroups"]
        groups = openmc.mgxs.EnergyGroups(np.append(np.array([0.]), self.lib.neutronEnergyUpperBounds))
        mg_cross_sections_file = openmc.MGXSLibrary(groups)
        for region in self.lib.regions:
            region_xsdata = openmc.XSdata(self.regionNames[str(region.regionNumber)], groups)
            region_xsdata.order = 0
            for interaction in region.macros.__dict__:
                data = region.macros[interaction]
                if data is not None:
                    if interaction == 'total':
                        region_xsdata.set_total(data)
                    if interaction == 'absorption':
                        region_xsdata.set_absorption(data)
                    if interaction == 'totalScatter':
                        region_xsdata.set_scatter_matrix(np.reshape(data.toarray(),(numGroups,numGroups,1)))
                    if interaction == 'fission':
                        region_xsdata.set_fission(data)
                    if interaction == 'nuSigF':
                        region_xsdata.set_nu_fission(data)                    
                    if interaction == 'neutronsPerFission':
                        region_xsdata.set_nu_fission(data*region.macros[fission])
                    if interaction == 'chi':
                        region_xsdata.set_chi(np.reshape(data,(numGroups,)))
            mg_cross_sections_file.add_xsdatas([region_xsdata])
        mg_cross_sections_file.export_to_hdf5(filename=self.outputFilename)
        return self.outputFilename

    def writeIsotxs(self):
        raise NotImplementedError()
