import numpy as np
import openmc
from armi.nuclearDataIO.cccc import compxs

class MGXSWriter:
    def __init__(self, inputFile, regionNames=None):
        self.inputFile = inputFile
        self.inputType = 'compxs'
        self.inputFormat = 'ascii'
        self.regionNames = {'0': 'fuel',
                            '1': 'mox low fuel',
                            '2': 'mox medium fuel',
                            '3': 'mox high fuel',
                            '4': 'moderator',
                            '5': 'guide tube',
                            '6': 'fission chamber'}
        self.lib = None
    def readInput(self):
        if self.inputType == 'compxs':
            self.lib = compxs.readAscii(self.inputFile)
    def write(self):
        if self.inputType == 'compxs':
            self.writeCompxs()
    def writeCompxs(self):
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
                        region_xsdata.set_scatter_matrix(np.reshape(data.toarray(),(7,7,1)))
                    if interaction == 'fission':
                        region_xsdata.set_fission(data)
                    if interaction == 'nuSigF':
                        region_xsdata.set_nu_fission(data)                    
                    if interaction == 'neutronsPerFission':
                        region_xsdata.set_nu_fission(data*region.macros[fission])
                    if interaction == 'chi':
                        region_xsdata.set_chi(np.reshape(data,(7,)))
            mg_cross_sections_file.add_xsdatas([region_xsdata])
        mg_cross_sections_file.export_to_hdf5()
