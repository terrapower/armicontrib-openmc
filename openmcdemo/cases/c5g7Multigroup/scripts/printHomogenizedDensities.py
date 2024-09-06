from math import pi

# Script to calculate number densities for use in homogenized-pin C5G7 case.
# These numbers are not used when the case is executed in multigroup mode.


def printNewDensity(label, oldDensity, oldVolume, newVolume):
    newDensity = oldDensity * oldVolume / newVolume
    print(label + ": " + str(f"{newDensity:.5}"))
    return newDensity


fuelOldVolume = pi * (0.8190 / 2) ** 2
fuelNewVolume = pi * (1.08 / 2) ** 2
print("UO2")
printNewDensity("U235", 8.65e-4, fuelOldVolume, fuelNewVolume)
printNewDensity("U238", 2.225e-2, fuelOldVolume, fuelNewVolume)
printNewDensity("O", 4.622e-2, fuelOldVolume, fuelNewVolume)
ZrOldVolume = pi * ((0.95 / 2) ** 2 - (0.836 / 2) ** 2)
printNewDensity("ZR", 4.3e-2, ZrOldVolume, fuelNewVolume)
AlOldVolume = pi * ((1.08 / 2) ** 2 - (0.97 / 2) ** 2)
printNewDensity("AL27", 6.00e-2, AlOldVolume, fuelNewVolume)
print("")
print("mox low")
printNewDensity("U235", 5.00e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("U238", 2.21e-2, fuelOldVolume, fuelNewVolume)
printNewDensity("PU238", 1.50e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("PU239", 5.80e-4, fuelOldVolume, fuelNewVolume)
printNewDensity("PU240", 2.40e-4, fuelOldVolume, fuelNewVolume)
printNewDensity("PU241", 9.80e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("PU242", 5.40e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("AM241", 1.30e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("O", 4.63e-2, fuelOldVolume, fuelNewVolume)
printNewDensity("ZR", 4.3e-2, ZrOldVolume, fuelNewVolume)
printNewDensity("AL27", 6.00e-2, AlOldVolume, fuelNewVolume)
print("")
print("mox medium")
printNewDensity("U235", 5.00e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("U238", 2.21e-2, fuelOldVolume, fuelNewVolume)
printNewDensity("PU238", 2.40e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("PU239", 9.30e-4, fuelOldVolume, fuelNewVolume)
printNewDensity("PU240", 3.90e-4, fuelOldVolume, fuelNewVolume)
printNewDensity("PU241", 1.52e-4, fuelOldVolume, fuelNewVolume)
printNewDensity("PU242", 8.40e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("AM241", 2.00e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("O", 4.63e-2, fuelOldVolume, fuelNewVolume)
printNewDensity("ZR", 4.3e-2, ZrOldVolume, fuelNewVolume)
printNewDensity("AL27", 6.00e-2, AlOldVolume, fuelNewVolume)
print("")
print("mox high")
printNewDensity("U235", 5.00e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("U238", 2.21e-2, fuelOldVolume, fuelNewVolume)
printNewDensity("PU238", 3.00e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("PU239", 1.16e-3, fuelOldVolume, fuelNewVolume)
printNewDensity("PU240", 4.90e-4, fuelOldVolume, fuelNewVolume)
printNewDensity("PU241", 1.90e-4, fuelOldVolume, fuelNewVolume)
printNewDensity("PU242", 1.05e-4, fuelOldVolume, fuelNewVolume)
printNewDensity("AM241", 2.50e-5, fuelOldVolume, fuelNewVolume)
printNewDensity("O", 4.63e-2, fuelOldVolume, fuelNewVolume)
printNewDensity("ZR", 4.3e-2, ZrOldVolume, fuelNewVolume)
printNewDensity("AL27", 6.00e-2, AlOldVolume, fuelNewVolume)
print("")
print("guide tube")
AlOldVolumeGT = pi * ((1.08 / 2) ** 2 - (0.68 / 2) ** 2)
printNewDensity("AL27", 6.00e-2, AlOldVolumeGT, fuelNewVolume)
modOldVolume = pi * (0.68 / 2) ** 2
printNewDensity("H", 6.70e-2, modOldVolume, fuelNewVolume)
printNewDensity("O", 3.35e-2, modOldVolume, fuelNewVolume)
printNewDensity("B", 2.78e-5, modOldVolume, fuelNewVolume)
printNewDensity("U235", 1.0e-8, modOldVolume, fuelNewVolume)
