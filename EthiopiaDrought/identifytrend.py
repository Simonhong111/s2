from osgeo import gdal,osr,ogr
import os
import glob
import numpy as np
import pandas as pd
import h5py
from netCDF4 import Dataset
from dateutil import rrule
from datetime import *
from matplotlib import cm
# from matplotlib import pyplot as plt
from pylab import plt
from PyEMD import EMD
x=[13,20]
y=[16,22]
x=[15,25]
y=[21,32]
def identifyTrend(chirpsdirectory):
    start = datetime.strptime("-".join(["2003", "01", "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2018-12-31", "%Y-%m-%d").date()
    value = []
    mT = []


    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):

        path = os.path.join(chirpsdirectory,
                                  "ET.v3.3a" + str(dt.year) + str(dt.month).zfill(
                                      2) + ".tif")

        if os.path.exists(path):
            data = gdal.Open(path).ReadAsArray(x[0], y[0], x[1] - x[0], y[1] - y[0])
            value.append(np.mean(data[data > -999]))
    return  value
            # mT.append(m)
if __name__ == "__main__":
    value = identifyTrend(r"D:\Cornell\EthiopianDrought\AData\ETClip")
    assert len(value) == 192, "missing data"
    mT = [i for i in range(192)]
    value = np.array(value)
    mT = np.array(mT)


    IMF = EMD(spline_kind="cubic").emd(value, mT,max_imf=5)
    sum = np.sum(IMF, axis=0)
    N = IMF.shape[0] + 1
    f = plt.figure(figsize=(12, 9))

    # Plot results
    plt.subplot(N, 1, 1)
    plt.title("Using Monthly ET data of all months from Jan to Dec")
    plt.plot(mT, value, 'r')
    plt.ylabel("ET")

    for n in range(N-2):
        plt.subplot(N, 1, n + 2)
        plt.plot(mT, IMF[n], 'g')
        plt.ylabel("eIMF %i" % (n + 1))
        plt.locator_params(axis='y', nbins=5)

    plt.tight_layout()

    plt.show()

    print(np.sum(IMF,axis=0).shape)