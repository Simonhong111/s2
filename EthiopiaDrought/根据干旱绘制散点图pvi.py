from osgeo import gdal,osr,ogr
import numpy as np
import glob
import os
from dateutil import rrule
from datetime import *
import csv
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
def monthAnalysis(path,monthtype,var1,var2,year):

    data = pd.read_csv(path)
    Var1List = data[monthtype+var1+str(year)].to_numpy()
    Var2List = data[monthtype+var2+str(year)].to_numpy()
    Freq = data[monthtype+"Fre"].to_numpy()
    return Var1List,Var2List,Freq


path = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\PolyGonAgg_anom.csv"
colors =['b','g','orange','brown','r','k']
DataType ="polygon"
for Year in range(2010,2017):

    fig = plt.figure(figsize=(15, 10))
    plt.title(str(Year) +" "+DataType+"\n",fontsize=16)

    for m in range(1, 7):

        ax = fig.add_subplot(3, 3, m)

        if m == 1:

            ylabel = "short PVI anom"
            xlabel = "short EVI anom"
            Var1List, Var2List, Freq = monthAnalysis(path, "Short", "EVI", "PVI", Year)
            droughtID = list(set(list(Freq)))
        if m == 2:
            ylabel = "short PVI anom"
            xlabel = "short GSIF anom"
            Var1List, Var2List, Freq = monthAnalysis(path, "Short", "GSIF", "PVI", Year)
            droughtID = list(set(list(Freq)))
        if m == 3:
            ylabel = "short PVI anom"
            xlabel = "short NSIF anom"
            Var1List, Var2List, Freq = monthAnalysis(path, "Short", "NSIF", "PVI", Year)
            droughtID = list(set(list(Freq)))

        if m == 4:
            ylabel = "long PVI anom:"
            xlabel = "long EVI anom"
            Var1List, Var2List, Freq = monthAnalysis(path, "Long", "EVI", "PVI", Year)
            droughtID = list(set(list(Freq)))
        if m == 5:
            ylabel = "long PVI anom"
            xlabel = "long GSIF Anomaly"
            Var1List, Var2List, Freq = monthAnalysis(path, "Long", "GSIF", "PVI", Year)
            droughtID = list(set(list(Freq)))
        if m == 6:
            print(m)
            ylabel = "long PVI anom"
            xlabel = "long NSIF anom"
            Var1List, Var2List, Freq =monthAnalysis(path, "Long", "NSIF", "PVI", Year)
            droughtID = list(set(list(Freq)))
        droughtID = sorted(droughtID)
        # print(droughtID)
        mcolor = []
        for i in range(Freq.shape[0]):

            mcolor.append(colors[Freq[i]])

        colors = np.array(colors)

        scatter = ax.scatter( Var1List,  Var2List, c=mcolor)
        titles = "R for freq "
        for freq in droughtID:
            titles = titles+str(int(freq))+" "
        titles=titles+":"
        for freq in droughtID:
            mask = np.where(Freq == freq)
            sample = mask[0].shape[0]
            slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[mask], Var2List[mask])
            print(m,mask[0].shape[0],r_value)
            titles = titles + str(round(r_value, 2)) + ","
        ax.tick_params(labelsize=8)
        ax.set_title(titles, fontsize=14)
        ax.set_xlabel(xlabel,fontsize=12)
        ax.set_ylabel(ylabel,fontsize=12)

    ax = fig.add_subplot(3, 3, 7)
    ax.scatter([0.96, 0.96, 0.96, 0.98, 0.98, 0.98], [1.5, 1, 0.5, 1.5, 1, 0.5],
               c=['b', 'g', 'orange', 'brown', 'r', 'k'])
    ax.text(0.965, 1.485, "Frequency 0")
    ax.text(0.965, 0.985, "Frequency 1")
    ax.text(0.965, 0.485, "Frequency 2")
    ax.text(0.985, 1.485, "Frequency 3")
    ax.text(0.985, 0.985, "Frequency 4")
    ax.text(0.985, 0.485, "Frequency 5")
    fig.tight_layout()  # 调整整体空白
    path2 = os.path.join(r"D:\Cornell\EthiopianDrought\CropCSV\Crop_image_drought_fre",
                         str(Year) + "polyPVI.jpg")
    print(path2)
    plt.savefig(path2)
    plt.close(fig)

# plt.show()
