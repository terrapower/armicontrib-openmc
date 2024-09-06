import numpy as np
from scipy.sparse import csc_matrix
from armi.nuclearDataIO.cccc import compxs
from armi.nuclearDataIO.xsLibraries import CompxsLibrary
from armi.nuclearDataIO.nuclearFileMetadata import RegionXSMetadata

c = CompxsLibrary()

c.compxsMetadata["numGroups"] = 7
c.compxsMetadata["fileWideChiFlag"] = 0
c.compxsMetadata["numFissComps"] = 5
c.compxsMetadata["maxUpScatterGroups"] = 7
c.compxsMetadata["maxDownScatterGroups"] = 7
c.compxsMetadata["numDelayedFam"] = 0
c.compxsMetadata["maxScatteringOrder"] = 0
c.compxsMetadata["reservedFlag1"] = 0
c.compxsMetadata["reservedFlag2"] = 0
c.compxsMetadata["minimumNeutronEnergy"] = 0.0
c.compxsMetadata["compFamiliesWithPrecursors"] = np.array([0,0,0,0,0,0,0])
c.compxsMetadata["fissionWattSeconds"] = np.array([-1.e+20,-1.e+20,-1.e+20,-1.e+20,-1.e+20,-1.e+20,-1.e+20])
c.compxsMetadata["captureWattSeconds"] = np.array([-1.e+20,-1.e+20,-1.e+20,-1.e+20,-1.e+20,-1.e+20,-1.e+20])
c.neutronVelocity = np.logspace(-5, 7, 8)[:-1]
c.neutronEnergyUpperBounds = np.logspace(-5, 7, 8)[1:]

uo2 = compxs.CompxsRegion(c, 0)
uo2.macros['numGroups'] = 7
uo2.macros['transport'] = np.array([1.779490E-01, 3.298050E-01, 4.803880E-01, 5.543670E-01, 3.118010E-01, 3.951680E-01, 5.644060E-01])
uo2.macros['absorption'] = np.array([8.02480E-03, 3.71740E-03, 2.67690E-02, 9.62360E-02, 3.00200E-02, 1.11260E-01, 2.82780E-01])
scatter_matrix = \
    [[[1.275370E-01, 4.237800E-02, 9.437400E-06, 5.516300E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 3.244560E-01, 1.631400E-03, 3.142700E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 4.509400E-01, 2.679200E-03, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 4.525650E-01, 5.566400E-03, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 1.252500E-04, 2.714010E-01, 1.025500E-02, 1.002100E-08],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 1.296800E-03, 2.658020E-01, 1.680900E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 8.545800E-03, 2.730800E-01]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.reshape(np.rollaxis(scatter_matrix, 0, 3), (7,7))
uo2.macros['totalScatter'] = csc_matrix(scatter_matrix)
uo2.macros['fission'] = np.array([7.212060E-03, 8.193010E-04, 6.453200E-03, 1.856480E-02, 1.780840E-02, 8.303480E-02, 2.160040E-01])
uo2.macros['nuSigF'] = np.array([2.005998E-02, 2.027303E-03, 1.570599E-02, 4.518301E-02, 4.334208E-02, 2.020901E-01, 5.257105E-01])
uo2.macros['neutronsPerFission'] = np.array([2.78145, 2.47443, 2.43383, 2.43380, 2.43380, 2.43380, 2.43380])
uo2.macros['chi'] = np.reshape(np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.000000E-00, 0.000000E-00, 0.000000E-00]), (7,1))

uo2.macros['total'] = uo2.macros['transport']
uo2.macros['removal'] = uo2.macros['absorption']
uo2.macros['n2n'] = np.array([0., 0., 0., 0., 0., 0., 0.])
uo2.macros['higherOrderScatter'] = uo2.macros['totalScatter']
uo2.metadata = RegionXSMetadata()
uo2.metadata["chiFlag"] = 1
uo2.metadata["numUpScatterGroups"] = [6,5,4,3,2,1,0]
uo2.metadata["numDownScatterGroups"] = [0,1,2,3,4,5,6]
uo2.metadata["powerConvMult"] = [1., 1., 1., 1., 1., 1., 1.]
uo2.metadata["d1Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
uo2.metadata["d1Additive"] = [0., 0., 0., 0., 0., 0., 0.]
uo2.metadata["d2Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
uo2.metadata["d2Additive"] = [0., 0., 0., 0., 0., 0., 0.]
uo2.metadata["d3Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
uo2.metadata["d3Additive"] = [0., 0., 0., 0., 0., 0., 0.]


mox43 = compxs.CompxsRegion(c, 1)
mox43.macros['numGroups'] = 7
mox43.macros['transport'] = np.array([1.787310E-01, 3.308490E-01, 4.837720E-01, 5.669220E-01, 4.262270E-01, 6.789970E-01, 6.828520E-01])
mox43.macros['absorption'] = np.array([8.43390E-03, 3.75770E-03, 2.79700E-02, 1.04210E-01, 1.39940E-01, 4.09180E-01, 4.09350E-01])
scatter_matrix = \
    [[[1.288760E-01, 4.141300E-02, 8.229000E-06, 5.040500E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 3.254520E-01, 1.639500E-03, 1.598200E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 4.531880E-01, 2.614200E-03, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 4.571730E-01, 5.539400E-03, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 1.604600E-04, 2.768140E-01, 9.312700E-03, 9.165600E-09],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 2.005100E-03, 2.529620E-01, 1.485000E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 8.494800E-03, 2.650070E-01]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.reshape(np.rollaxis(scatter_matrix, 0, 3), (7,7))
mox43.macros['totalScatter'] = csc_matrix(scatter_matrix)
mox43.macros['fission'] = np.array([7.62704E-03, 8.76898E-04, 5.69835E-03, 2.28872E-02, 1.07635E-02, 2.32757E-01, 2.48968E-01])
mox43.macros['nuSigF'] = np.array([2.175300E-02, 2.535103E-03, 1.626799E-02, 6.547410E-02, 3.072409E-02, 6.666510E-01, 7.139904E-01])
mox43.macros['neutronsPerFission'] = np.array([2.85209, 2.89099, 2.85486, 2.86073, 2.85447, 2.86415, 2.86780])
mox43.macros['chi'] = np.reshape(np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.000000E-00, 0.000000E-00, 0.000000E-00]), (7,1))

mox43.macros['total'] = mox43.macros['transport']
mox43.macros['removal'] = mox43.macros['absorption']
mox43.macros['n2n'] = np.array([0., 0., 0., 0., 0., 0., 0.])
mox43.macros['higherOrderScatter'] = mox43.macros['totalScatter']
mox43.metadata = RegionXSMetadata()
mox43.metadata["chiFlag"] = 1
mox43.metadata["numUpScatterGroups"] = [6,5,4,3,2,1,0]
mox43.metadata["numDownScatterGroups"] = [0,1,2,3,4,5,6]
mox43.metadata["powerConvMult"] = [1., 1., 1., 1., 1., 1., 1.]
mox43.metadata["d1Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
mox43.metadata["d1Additive"] = [0., 0., 0., 0., 0., 0., 0.]
mox43.metadata["d2Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
mox43.metadata["d2Additive"] = [0., 0., 0., 0., 0., 0., 0.]
mox43.metadata["d3Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
mox43.metadata["d3Additive"] = [0., 0., 0., 0., 0., 0., 0.]

mox70 = compxs.CompxsRegion(c, 2)
mox70.macros['numGroups'] = 7
mox70.macros['transport'] = np.array([1.813230E-01, 3.343680E-01, 4.937850E-01, 5.912160E-01, 4.741980E-01, 8.336010E-01, 8.536030E-01])
mox70.macros['absorption'] = np.array([9.06570E-03, 4.29670E-03, 3.28810E-02, 1.22030E-01, 1.82980E-01, 5.68460E-01, 5.85210E-01])
scatter_matrix = \
    [[[1.304570E-01, 4.179200E-02, 8.510500E-06, 5.132900E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 3.284280E-01, 1.643600E-03, 2.201700E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 4.583710E-01, 2.533100E-03, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 4.637090E-01, 5.476600E-03, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 1.761900E-04, 2.823130E-01, 8.728900E-03, 9.001600E-09],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 2.276000E-03, 2.497510E-01, 1.311400E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 8.864500E-03, 2.595290E-01]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.reshape(np.rollaxis(scatter_matrix, 0, 3), (7,7))
mox70.macros['totalScatter'] = csc_matrix(scatter_matrix)
mox70.macros['fission'] = np.array([8.25446E-03, 1.32565E-03, 8.42156E-03, 3.28730E-02, 1.59636E-02, 3.23794E-01, 3.62803E-01])
mox70.macros['nuSigF'] = np.array([2.381395E-02, 3.858689E-03, 2.413400E-02, 9.436622E-02, 4.576988E-02, 9.281814E-01, 1.043200E+00])
mox70.macros['neutronsPerFission'] = np.array([2.88498, 2.91079, 2.86574, 2.87063, 2.86714, 2.86658, 2.87539])
mox70.macros['chi'] = np.reshape(np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.000000E-00, 0.000000E-00, 0.000000E-00]), (7,1))

mox70.macros['total'] = mox70.macros['transport']
mox70.macros['removal'] = mox70.macros['absorption']
mox70.macros['n2n'] = np.array([0., 0., 0., 0., 0., 0., 0.])
mox70.macros['higherOrderScatter'] = mox43.macros['totalScatter']
mox70.metadata = RegionXSMetadata()
mox70.metadata["chiFlag"] = 1
mox70.metadata["numUpScatterGroups"] = [6,5,4,3,2,1,0]
mox70.metadata["numDownScatterGroups"] = [0,1,2,3,4,5,6]
mox70.metadata["powerConvMult"] = [1., 1., 1., 1., 1., 1., 1.]
mox70.metadata["d1Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
mox70.metadata["d1Additive"] = [0., 0., 0., 0., 0., 0., 0.]
mox70.metadata["d2Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
mox70.metadata["d2Additive"] = [0., 0., 0., 0., 0., 0., 0.]
mox70.metadata["d3Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
mox70.metadata["d3Additive"] = [0., 0., 0., 0., 0., 0., 0.]

mox87 = compxs.CompxsRegion(c, 3)
mox87.macros['numGroups'] = 7
mox87.macros['transport'] = np.array([1.830450E-01, 3.367050E-01, 5.005070E-01, 6.061740E-01, 5.027540E-01, 9.210280E-01, 9.552310E-01])
mox87.macros['absorption'] = np.array([9.48620E-03, 4.65560E-03, 3.62400E-02, 1.32720E-01, 2.08400E-01, 6.58700E-01, 6.90170E-01])
scatter_matrix = \
    [[[1.315040E-01, 4.204600E-02, 8.697200E-06, 5.193800E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 3.304030E-01, 1.646300E-03, 2.600600E-09, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 4.617920E-01, 2.474900E-03, 0.000000E-00, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 4.680210E-01, 5.433000E-03, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 1.859700E-04, 2.857710E-01, 8.397300E-03, 8.928000E-09],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 2.391600E-03, 2.476140E-01, 1.232200E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 8.968100E-03, 2.560930E-01]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.reshape(np.rollaxis(scatter_matrix, 0, 3), (7,7))
mox87.macros['totalScatter'] = csc_matrix(scatter_matrix)
mox87.macros['fission'] = np.array([8.67209E-03, 1.62426E-03, 1.02716E-02, 3.90447E-02, 1.92576E-02, 3.74888E-01, 4.30599E-01])
mox87.macros['nuSigF'] = np.array([2.518600E-02, 4.739509E-03, 2.947805E-02, 1.122500E-01, 5.530301E-02, 1.074999E+00, 1.239298E+00])
mox87.macros['neutronsPerFission'] = np.array([2.90426, 2.91795, 2.86986, 2.87491, 2.87175, 2.86752, 2.87808])
mox87.macros['chi'] = np.reshape(np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.000000E-00, 0.000000E-00, 0.000000E-00]), (7,1))

mox87.macros['total'] = mox87.macros['transport']
mox87.macros['removal'] = mox87.macros['absorption']
mox87.macros['n2n'] = np.array([0., 0., 0., 0., 0., 0., 0.])
mox87.macros['higherOrderScatter'] = mox43.macros['totalScatter']
mox87.metadata = RegionXSMetadata()
mox87.metadata["chiFlag"] = 1
mox87.metadata["numUpScatterGroups"] = [6,5,4,3,2,1,0]
mox87.metadata["numDownScatterGroups"] = [0,1,2,3,4,5,6]
mox87.metadata["powerConvMult"] = [1., 1., 1., 1., 1., 1., 1.]
mox87.metadata["d1Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
mox87.metadata["d1Additive"] = [0., 0., 0., 0., 0., 0., 0.]
mox87.metadata["d2Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
mox87.metadata["d2Additive"] = [0., 0., 0., 0., 0., 0., 0.]
mox87.metadata["d3Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
mox87.metadata["d3Additive"] = [0., 0., 0., 0., 0., 0., 0.]

moderator = compxs.CompxsRegion(c, 4)
moderator.macros['numGroups'] = 7
moderator.macros['transport'] = np.array([1.592060E-01, 4.129700E-01, 5.903100E-01, 5.843500E-01, 7.180000E-01, 1.254450E+00, 2.650380E+00])
moderator.macros['absorption'] = np.array([6.01050E-04, 1.57930E-05, 3.37160E-04, 1.94060E-03, 5.74160E-03, 1.50010E-02, 3.72390E-02])
scatter_matrix = \
    [[[4.447770E-02, 1.134000E-01, 7.234700E-04, 3.749900E-06, 5.318400E-08, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 2.823340E-01, 1.299400E-01, 6.234000E-04, 4.800200E-05, 7.448600E-06, 1.045500E-06],
      [0.000000E-00, 0.000000E-00, 3.452560E-01, 2.245700E-01, 1.699900E-02, 2.644300E-03, 5.034400E-04],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 9.102840E-02, 4.155100E-01, 6.373200E-02, 1.213900E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 7.143700E-05, 1.391380E-01, 5.118200E-01, 6.122900E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 2.215700E-03, 6.999130E-01, 5.373200E-01],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 1.324400E-01, 2.480700E+00]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.reshape(np.rollaxis(scatter_matrix, 0, 3), (7,7))
moderator.macros['totalScatter'] = csc_matrix(scatter_matrix)

moderator.macros['total'] = moderator.macros['transport']
moderator.macros['removal'] = moderator.macros['absorption']
moderator.macros['n2n'] = np.array([0., 0., 0., 0., 0., 0., 0.])
moderator.macros['higherOrderScatter'] = mox43.macros['totalScatter']
moderator.metadata = RegionXSMetadata()
moderator.metadata["chiFlag"] = 0
moderator.metadata["numUpScatterGroups"] = [6,5,4,3,2,1,0]
moderator.metadata["numDownScatterGroups"] = [0,1,2,3,4,5,6]
moderator.metadata["powerConvMult"] = [1., 1., 1., 1., 1., 1., 1.]
moderator.metadata["d1Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
moderator.metadata["d1Additive"] = [0., 0., 0., 0., 0., 0., 0.]
moderator.metadata["d2Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
moderator.metadata["d2Additive"] = [0., 0., 0., 0., 0., 0., 0.]
moderator.metadata["d3Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
moderator.metadata["d3Additive"] = [0., 0., 0., 0., 0., 0., 0.]

guideTube = compxs.CompxsRegion(c, 5)
guideTube.macros['numGroups'] = 7
guideTube.macros['transport'] = np.array([1.260320E-01, 2.931600E-01, 2.842400E-01, 2.809600E-01, 3.344400E-01, 5.656400E-01, 1.172150E+00])
guideTube.macros['absorption'] = np.array([5.11320E-04, 7.58010E-05, 3.15720E-04, 1.15820E-03, 3.39750E-03, 9.18780E-03, 2.32420E-02])
scatter_matrix = \
    [[[6.616590E-02, 5.907000E-02, 2.833400E-04, 1.462200E-06, 2.064200E-08, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 2.403770E-01, 5.243500E-02, 2.499000E-04, 1.923900E-05, 2.987500E-06, 4.214000E-07],
      [0.000000E-00, 0.000000E-00, 1.832970E-01, 9.239700E-02, 6.944600E-03, 1.080300E-03, 2.056700E-04],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 7.885110E-02, 1.701400E-01, 2.588100E-02, 4.929700E-03],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 3.733300E-05, 9.973720E-02, 2.067900E-01, 2.447800E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 9.172600E-04, 3.167650E-01, 2.387700E-01],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 4.979200E-02, 1.099120E+00]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.reshape(np.rollaxis(scatter_matrix, 0, 3), (7,7))
guideTube.macros['totalScatter'] = csc_matrix(scatter_matrix)

guideTube.macros['total'] = guideTube.macros['transport']
guideTube.macros['removal'] = guideTube.macros['absorption']
guideTube.macros['n2n'] = np.array([0., 0., 0., 0., 0., 0., 0.])
guideTube.macros['higherOrderScatter'] = uo2.macros['totalScatter']
guideTube.metadata = RegionXSMetadata()
guideTube.metadata["chiFlag"] = 0
guideTube.metadata["numUpScatterGroups"] = [6,5,4,3,2,1,0]
guideTube.metadata["numDownScatterGroups"] = [0,1,2,3,4,5,6]
guideTube.metadata["powerConvMult"] = [1., 1., 1., 1., 1., 1., 1.]
guideTube.metadata["d1Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
guideTube.metadata["d1Additive"] = [0., 0., 0., 0., 0., 0., 0.]
guideTube.metadata["d2Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
guideTube.metadata["d2Additive"] = [0., 0., 0., 0., 0., 0., 0.]
guideTube.metadata["d3Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
guideTube.metadata["d3Additive"] = [0., 0., 0., 0., 0., 0., 0.]

fissionChamber = compxs.CompxsRegion(c, 6)
fissionChamber.macros['numGroups'] = 7
fissionChamber.macros['transport'] = np.array([1.260320E-01, 2.931600E-01, 2.842500E-01, 2.810200E-01, 3.344600E-01, 5.656400E-01, 1.172140E+00])
fissionChamber.macros['absorption'] = np.array([5.11320E-04, 7.58130E-05, 3.16430E-04, 1.16750E-03, 3.39770E-03, 9.18860E-03, 2.32440E-02])
scatter_matrix = \
    [[[6.616590E-02, 5.907000E-02, 2.833400E-04, 1.462200E-06, 2.064200E-08, 0.000000E-00, 0.000000E-00],
      [0.000000E-00, 2.403770E-01, 5.243500E-02, 2.499000E-04, 1.923900E-05, 2.987500E-06, 4.214000E-07],
      [0.000000E-00, 0.000000E-00, 1.834250E-01, 9.228800E-02, 6.936500E-03, 1.079000E-03, 2.054300E-04],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 7.907690E-02, 1.699900E-01, 2.586000E-02, 4.925600E-03],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 3.734000E-05, 9.975700E-02, 2.067900E-01, 2.447800E-02],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 9.174200E-04, 3.167740E-01, 2.387600E-01],
      [0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 4.979300E-02, 1.09910E+00]]]
scatter_matrix = np.array(scatter_matrix)
scatter_matrix = np.reshape(np.rollaxis(scatter_matrix, 0, 3), (7,7))
fissionChamber.macros['totalScatter'] = csc_matrix(scatter_matrix)
fissionChamber.macros['fission'] = np.array([4.79002E-09, 5.82564E-09, 4.63719E-07, 5.24406E-06, 1.45390E-07, 7.14972E-07, 2.08041E-06])
fissionChamber.macros['nuSigF'] = np.array([1.323401E-08, 1.434500E-08, 1.128599E-06, 1.276299E-05, 3.538502E-07, 1.740099E-06, 5.063302E-06])
fissionChamber.macros['neutronsPerFission'] = np.array([2.76283, 2.46239, 2.43380, 2.43380, 2.43380, 2.43380, 2.43380])
fissionChamber.macros['chi'] = np.reshape(np.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0.000000E-00, 0.000000E-00, 0.000000E-00]), (7,1))

fissionChamber.macros['total'] = fissionChamber.macros['transport']
fissionChamber.macros['removal'] = fissionChamber.macros['absorption']
fissionChamber.macros['n2n'] = np.array([0., 0., 0., 0., 0., 0., 0.])
fissionChamber.macros['higherOrderScatter'] = uo2.macros['totalScatter']
fissionChamber.metadata = RegionXSMetadata()
fissionChamber.metadata["chiFlag"] = 1
fissionChamber.metadata["numUpScatterGroups"] = [6,5,4,3,2,1,0]
fissionChamber.metadata["numDownScatterGroups"] = [0,1,2,3,4,5,6]
fissionChamber.metadata["powerConvMult"] = [1., 1., 1., 1., 1., 1., 1.]
fissionChamber.metadata["d1Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
fissionChamber.metadata["d1Additive"] = [0., 0., 0., 0., 0., 0., 0.]
fissionChamber.metadata["d2Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
fissionChamber.metadata["d2Additive"] = [0., 0., 0., 0., 0., 0., 0.]
fissionChamber.metadata["d3Multiplier"] = [1., 1., 1., 1., 1., 1., 1.]
fissionChamber.metadata["d3Additive"] = [0., 0., 0., 0., 0., 0., 0.]

c.compxsMetadata["numComps"] = len(c.regions)

compxs.writeAscii(c,"C5G7_COMPXS.ascii")
