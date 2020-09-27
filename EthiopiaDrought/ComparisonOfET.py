from osgeo import gdal,osr,ogr
from datetime import *
from dateutil import rrule
import os
from epdatapro import write_Img
# import datetime
import glob
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import calendar
from scipy.stats import gaussian_kde
from scipy import stats
MName =['Jan','Feb','Mar','Apr','May','Jun','Jul',"Aug",'Sep','Oct','Nov','Dec']

def scatter(year):
    fig = plt.figure(figsize=(12,10))

    for month in range(1,13):


        gleamET = gdal.Open(r'D:\Cornell\EthiopianDrought\AData\ETMonthV3.3a'+'\\'+"ET.v3.3a"+str(year)+
                        str(month).zfill(2)+'.tif').ReadAsArray()
        bessET = gdal.Open(r'D:\Cornell\BessAgg0.25V2'+'\\'+str(year)+str(month).zfill(2)+"01.tif").ReadAsArray()

        mask = (gleamET >-999.0) & (bessET > -999.0)
        gleamET = gleamET[mask]
        bessET = bessET[mask]

        slope, intercept, r_value, p_value, std_err = stats.linregress(gleamET, bessET)
        minx,maxx = gleamET.min(),gleamET.max()
        miny, maxy = bessET.min(), bessET.max()
        xy = np.vstack([gleamET, bessET])
        z = gaussian_kde(xy)(xy)
        ax = fig.add_subplot(4, 3, month)  # create 12-months subfigure

        ax.set_title(MName[month-1]+','+str(year))
        ax.scatter(gleamET, bessET, c=z, s=100, edgecolor='')
        ax.plot([minx,maxx],[slope*minx+intercept,slope*maxx+intercept])
        tex = r'${R}^2$='+str(round(r_value*r_value,3))
        ax.text(minx + (maxx-minx)*0.01,maxy - (maxy - miny)*0.1 , tex, fontsize=10)
        tex = r'$Y={}X+{}$'.format(round(slope,3),round(intercept,3))
        ax.text(minx + (maxx - minx) * 0.01, maxy - (maxy - miny) * 0.2, tex, fontsize=10)
        ax.set_xlabel('Monthly Gleam ET')
        ax.set_ylabel('Monthly Bess ET')
    fig.tight_layout()  # 调整整体空白
    plt.show()

# scatter(2013)

def Map(year,month):
    fig = plt.figure(figsize=(5,8))


    gleamET = gdal.Open(r'D:\Cornell\EthiopianDrought\AData\ETMonthV3.3a'+'\\'+"ET.v3.3a"+str(year)+
                    str(month).zfill(2)+'.tif').ReadAsArray()
    bessET = gdal.Open(r'D:\Cornell\BessAgg0.25V2'+'\\'+str(year)+str(month).zfill(2)+"01.tif").ReadAsArray()

    rainFall = gdal.Open(r'D:\Cornell\EthiopianDrought\Chirps2'+'\\'+'chirps-v2.0.{}.{}.tif'.
                         format(str(year),str(month).zfill(2))).ReadAsArray()
    mask = (gleamET >-999.0) & (bessET > -999.0)
    gleamET2 = gleamET[mask]
    bessET2 = bessET[mask]
    minx, maxx = gleamET2.min(), gleamET2.max()
    miny, maxy = bessET2.min(), bessET2.max()
    minX = min(minx,miny)
    maxX = max(maxx,maxy)

    # rainFall[rainFall < -999] == np.nan
    rainFall2 = rainFall[rainFall > -999]
    rmin,rmax = rainFall2.min(),rainFall2.max()


    from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
    cmap = plt.get_cmap("RdBu")

    ax1 = fig.add_subplot(3, 1, 1)
    ax1.set_title('Monthly Gleam ET {},{}'.format(MName[month-1],year))

    cax1 = ax1.imshow(gleamET, cmap=cmap, vmin=minX, vmax=maxX)
    cbar1 = plt.colorbar(cax1, ax=ax1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax2.set_title('Monthly Bess ET {},{}'.format(MName[month-1],year))

    cax2 = ax2.imshow(bessET, cmap=cmap, vmin=minX, vmax=maxX)
    cbar = plt.colorbar(cax2, ax=ax2)
    cmap.set_under(color='white')

    ax3 = fig.add_subplot(3, 1, 3)
    ax3.set_title('Monthly Precipitation {},{}'.format(MName[month - 1], year))

    cax3 = ax3.imshow(rainFall, cmap=cmap, vmin=rmin, vmax=rmax)
    cbar = plt.colorbar(cax3, ax=ax3)
    cmap.set_under(color='white')

    fig.tight_layout()  # 调整整体空白
    # plt.savefig(r"D:\Cornell\Slides\\"+str(year)+MName[month-1]+'_Spatial.png')
    plt.show()

for i in range(1,13):
    Map(2015, i)

def TSeries(GleamDir,BessDir):
    start = datetime.strptime("-".join(["2002", '03', "01"]), "%Y-%m-%d").date()
    stop = datetime.strptime("2015-12-31", "%Y-%m-%d").date()
    gMean = []
    bMean = []
    gStd = []
    bStd = []
    axisTime = []

    for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
        bessET = os.path.join(BessDir,
                                   "{}{}01.tif".format(str(dt.year),str(dt.month).zfill(2)))
        gleamET = os.path.join(GleamDir,
                              "ET.v3.3a{}{}.tif".format(str(dt.year), str(dt.month).zfill(2)))

        assert os.path.exists(bessET) and os.path.exists(gleamET),'no file {}{}'.format(bessET,gleamET)

        gleamETArr = gdal.Open(gleamET).ReadAsArray()
        bessETArr = gdal.Open(bessET).ReadAsArray()

        mask = (gleamETArr > -999.0) & (bessETArr > -999.0)
        gleamETArr = gleamETArr[mask]
        bessETArr = bessETArr[mask]

        gleamETArrStd = gleamETArr.std()
        bessETArrstd = bessETArr.std()
        gleamETArrMean = gleamETArr.mean()
        bessETArrMean = bessETArr.mean()
        gMean.append(gleamETArrMean)
        gStd.append(gleamETArrStd)
        bMean.append(bessETArrMean)
        bStd.append(bessETArrstd)
        axisTime.append(str(dt.year)+str(dt.month).zfill(2))

    gStd = np.array(gStd)
    gMean = np.array(gMean)
    bStd = np.array(bStd)
    bMean = np.array(bMean)
    TiLabel =[]
    TiValue = []
    T = [i for i in range(len(bStd))]
    for idx,ti in enumerate(axisTime):
        if ti.endswith('06'):
            TiLabel.append(ti)
            TiValue.append(T[idx])
    print(TiLabel,TiValue)

    plt.title("ET Time Series 200203-201512")

    plt.plot(T, np.array(gMean),label='Gleam ET',c='green')
    plt.fill_between(T, gMean + gStd, gMean - gStd, color='lightgreen')
    plt.plot(T, bMean, label='Bess ET', c='blue')
    plt.fill_between(T, bMean + bStd, bMean - bStd, color='lightblue')
    plt.xticks(TiValue,TiLabel)
    for xtick in plt.gca().get_xticklabels():
        xtick.set_rotation(45)
    plt.legend()
    plt.show()




# TSeries(r'D:\Cornell\EthiopianDrought\AData\ETMonthV3.3a',r'D:\Cornell\BessAgg0.25V2')


