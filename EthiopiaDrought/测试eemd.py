from osgeo import gdal,osr,ogr
import os
import numpy as np
from dateutil import rrule
from datetime import *
from scipy import signal
from Experiment.Exp_Write_Image import write_Img
from matplotlib import pyplot as plt
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
from pyeemd import ceemdan
from pyeemd.utils import plot_imfs
if __name__ == "__main__":
    y = np.sin(2 * np.pi * np.linspace(0, 1, 1000)) + np.sin(20 * np.pi * np.linspace(0, 1, 1000))
    # plt.plot(y)
    # plt.show()
    imfs = ceemdan(y, num_imfs=4, S_number=4, num_siftings=50)
    # # plot_imfs(imfs)
    # # plt.show()
    # plt.plot(y-imfs[-1]-imfs[-2])
    # plt.show()
    print(imfs[-1])
