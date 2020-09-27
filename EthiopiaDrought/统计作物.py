from osgeo import gdal,osr,ogr
import numpy as np
import glob
import os
from dateutil import rrule
from datetime import *
import csv
import pandas as pd
from matplotlib import pyplot as plt


def ExtractCSV(path):
    data = pd.read_csv(path)
    CTp = ["WHEATOPH", "MAIZEOPH", "BARLEYOPH", "SORGHUMOPH", "TEFFOPH"]

    keys = []
    for year in range(2010,2017):
        for ct in CTp:
            Var = data[ct+str(year)].to_numpy()
            count = 0
            for v in Var:
               if np.isnan(v).any()==True or (v==0).any()==True:
                   continue
               else:
                   count +=1
            print(ct+str(year),count)





CTp = ["WHEATOPH","MAIZEOPH","BARLEYOPH","SORGHUMOPH","TEFFOPH"]
croppath = r"D:\Cornell\EthiopianDrought\CropCSV\AgSS_2010_2016_5_Crops.csv"

ExtractCSV(croppath)




