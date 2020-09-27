import numpy as np
import pandas as pd
import csv,os
from matplotlib import pyplot as plt
from scipy import stats
from osgeo import gdal,osr,ogr
from scipy import signal
import ctypes
ctypes.WinDLL("kernel32.dll")
ctypes.WinDLL("msvcrt.dll")
ctypes.WinDLL("user32.dll")
ctypes.WinDLL(r"D:\msys64\mingw64\bin\libwinpthread-1.dll")
ctypes.WinDLL(r"D:\msys64\mingw64\bin\libgcc_s_seh-1.dll")
ctypes.WinDLL(r"D:\msys64\mingw64\bin\libgomp-1.dll")
ctypes.WinDLL(r"D:\msys64\home\zmhwh\gsl-2.6\cblas\.libs\libgslcblas-0.dll")
ctypes.WinDLL(r"D:\msys64\home\zmhwh\gsl-2.6\.libs\libgsl-25.dll")
ctypes.WinDLL(r"D:\msys64\home\zmhwh\libeemd\libeemd.so")
from pyeemd import ceemdan,eemd

def fit(y):
    return signal.detrend(y)

path = r"D:\Cornell\EthiopianDrought\0ExperimentData\ExtractPixelCSV\Agg_Mask50.csv"


def ExtractVIFromCSV(path,StartYear,EndYear):


    data = pd.read_csv(path)

    RFKS = ["RF{}Short".format(year) for year in range(StartYear,EndYear+1)]
    RFKL = ["RF{}Long".format(year) for year in range(StartYear,EndYear+1)]


    RFS = data[RFKS].to_numpy()
    RFL = data[RFKL].to_numpy()
    RFSA = (RFS - RFS.mean(axis=1)[:,np.newaxis])/RFS.std(axis=1)[:,np.newaxis]
    del RFS
    RFLA = (RFL - RFL.mean(axis=1)[:, np.newaxis]) / RFL.std(axis=1)[:, np.newaxis]
    del RFL

    EVIKS = ["EVI{}Short".format(year) for year in range(StartYear,EndYear+1)]
    EVIKL = ["EVI{}Long".format(year) for year in range(StartYear,EndYear+1)]
    EVIS = data[EVIKS].to_numpy()
    EVIL = data[EVIKL].to_numpy()
    EVISA = (EVIS - EVIS.mean(axis=1)[:, np.newaxis]) / EVIS.std(axis=1)[:, np.newaxis]
    del EVIS
    EVILA = (EVIL - EVIL.mean(axis=1)[:, np.newaxis]) / EVIL.std(axis=1)[:, np.newaxis]
    del EVIL

    SIFKS = ["SIF{}Short".format(year) for year in range(StartYear, EndYear + 1)]
    SIFKL = ["SIF{}Long".format(year) for year in range(StartYear, EndYear + 1)]
    SIFS = data[SIFKS].to_numpy()
    SIFL = data[SIFKL].to_numpy()

    results = map(fit, SIFS)

    deResults = np.array(list(results))
    SIFS = deResults
    del deResults
    results = map(fit, SIFL)
    deResults = np.array(list(results))
    SIFL = deResults
    del deResults

    SIFSA = (SIFS - SIFS.mean(axis=1)[:, np.newaxis]) / SIFS.std(axis=1)[:, np.newaxis]
    del SIFS
    SIFLA = (SIFL - SIFL.mean(axis=1)[:, np.newaxis]) / SIFL.std(axis=1)[:, np.newaxis]
    del SIFL

    PVIKS = ["PVI{}Short".format(year) for year in range(StartYear, EndYear + 1)]
    PVIKL = ["PVI{}Long".format(year) for year in range(StartYear, EndYear + 1)]
    PVIS = data[PVIKS].to_numpy()
    PVIL = data[PVIKL].to_numpy()
    PVISA = (PVIS - PVIS.mean(axis=1)[:, np.newaxis]) / PVIS.std(axis=1)[:, np.newaxis]
    del PVIS
    PVILA = (PVIL - PVIL.mean(axis=1)[:, np.newaxis]) / PVIL.std(axis=1)[:, np.newaxis]
    del PVIL
    COL = data["COL"].to_numpy()
    ROW = data["ROW"].to_numpy()
    DataDict = {}

    for index, year in enumerate(range(StartYear, EndYear + 1)):
        DataDict["RF" + str(year) + "Short"] = RFSA[:, index]
        DataDict["RF" + str(year) + "Long"] = RFLA[:, index]
        DataDict["EVI" + str(year) + "Short"] = EVISA[:, index]
        DataDict["EVI" + str(year) + "Long"] = EVILA[:, index]
        DataDict["SIF" + str(year) + "Short"] = SIFSA[:, index]
        DataDict["SIF" + str(year) + "Long"] = SIFLA[:, index]
        DataDict["PVI" + str(year) + "Short"] = PVISA[:, index]
        DataDict["PVI" + str(year) + "Long"] = PVILA[:, index]

    DataDict["ROW"] = ROW.flatten()
    DataDict["COL"] = COL.flatten()

    df = pd.DataFrame(DataDict)

    outpath = r"D:\Cornell\EthiopianDrought\0ExperimentData\ExtractPixelCSV\Agg_Mask50_Anom.csv"
    df.to_csv(outpath, index=False)
def ExtractVIFromCSV2(path,StartYear,EndYear):


    data = pd.read_csv(path)

    RFKS = ["RF{}Short".format(year) for year in range(StartYear,EndYear+1)]
    RFKL = ["RF{}Long".format(year) for year in range(StartYear,EndYear+1)]


    RFS = data[RFKS].to_numpy()
    RFL = data[RFKL].to_numpy()
    RFSA = (RFS - RFS.mean(axis=1)[:,np.newaxis])
    del RFS
    RFLA = (RFL - RFL.mean(axis=1)[:, np.newaxis])
    del RFL

    EVIKS = ["EVI{}Short".format(year) for year in range(StartYear,EndYear+1)]
    EVIKL = ["EVI{}Long".format(year) for year in range(StartYear,EndYear+1)]
    EVIS = data[EVIKS].to_numpy()
    EVIL = data[EVIKL].to_numpy()
    EVISA = (EVIS - EVIS.mean(axis=1)[:, np.newaxis])
    del EVIS
    EVILA = (EVIL - EVIL.mean(axis=1)[:, np.newaxis])
    del EVIL

    SIFKS = ["SIF{}Short".format(year) for year in range(StartYear, EndYear + 1)]
    SIFKL = ["SIF{}Long".format(year) for year in range(StartYear, EndYear + 1)]
    SIFS = data[SIFKS].to_numpy()
    SIFL = data[SIFKL].to_numpy()

    results = map(fit, SIFS)
    deResults = np.array(list(results))
    SIFS = deResults
    del deResults
    results = map(fit, SIFL)
    deResults = np.array(list(results))
    SIFL = deResults
    del deResults

    SIFSA = (SIFS - SIFS.mean(axis=1)[:, np.newaxis])
    del SIFS
    SIFLA = (SIFL - SIFL.mean(axis=1)[:, np.newaxis])
    del SIFL

    PVIKS = ["PVI{}Short".format(year) for year in range(StartYear, EndYear + 1)]
    PVIKL = ["PVI{}Long".format(year) for year in range(StartYear, EndYear + 1)]
    PVIS = data[PVIKS].to_numpy()
    PVIL = data[PVIKL].to_numpy()
    PVISA = (PVIS - PVIS.mean(axis=1)[:, np.newaxis])
    del PVIS
    PVILA = (PVIL - PVIL.mean(axis=1)[:, np.newaxis])
    del PVIL
    COL = data["COL"].to_numpy()
    ROW = data["ROW"].to_numpy()
    DataDict = {}

    for index, year in enumerate(range(StartYear, EndYear + 1)):
        DataDict["RF" + str(year) + "Short"] = RFSA[:, index]
        DataDict["RF" + str(year) + "Long"] = RFLA[:, index]
        DataDict["EVI" + str(year) + "Short"] = EVISA[:, index]
        DataDict["EVI" + str(year) + "Long"] = EVILA[:, index]
        DataDict["SIF" + str(year) + "Short"] = SIFSA[:, index]
        DataDict["SIF" + str(year) + "Long"] = SIFLA[:, index]
        DataDict["PVI" + str(year) + "Short"] = PVISA[:, index]
        DataDict["PVI" + str(year) + "Long"] = PVILA[:, index]

    DataDict["ROW"] = ROW.flatten()
    DataDict["COL"] = COL.flatten()

    df = pd.DataFrame(DataDict)

    outpath = r"D:\Cornell\EthiopianDrought\0ExperimentData\ExtractPixelCSV\Agg_Mask50_Anom.csv"
    df.to_csv(outpath, index=False)


ExtractVIFromCSV(path,2007,2018)
