import numpy as np
import pandas as pd
from osgeo import gdal,osr,ogr
from matplotlib import pyplot as plt
from scipy import stats
import matplotlib.gridspec as gridspec
path = r"D:\Cornell\EthiopianDrought\0ExperimentData\ExtractPixelCSV\Agg_Mask50.csv"

def ExtractVIFromCSV(path,SeasonType,var1,var2,year):

    assert SeasonType in ["Short","Long"], "SeasonType must be Short/Long"
    data = pd.read_csv(path)
    key1 = var1 + str(year) + SeasonType
    key2 = var2 + str(year) + SeasonType
    key3 = "PVI" + str(year) + SeasonType
    Var1List = data[key1].to_numpy()
    Var2List = data[key2].to_numpy()
    PVIList =  data[key3].to_numpy()
    COL = data["COL"].to_numpy()
    ROW = data["ROW"].to_numpy()

    return Var1List,Var2List,PVIList,COL,ROW
def GenValidRC():
    ref_path = r"D:\Cornell\EthiopianDrought\0ExperimentData\CropMask\agg_clip50.tif"
    ref_raster = gdal.Open(ref_path)
    geo_t = ref_raster.GetGeoTransform()

    daShapefile = r"D:\Cornell\EthiopianDrought\0ExperimentData\CropSHP\CropLandR.shp"

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(daShapefile, 0)
    layer = dataSource.GetLayer()
    feature = layer.GetFeature(0)
    geo = feature.GetGeometryRef()
    geo = str(geo).split("((")[1].split("))")[0].split(",")
    x2 = []
    y2 = []
    for term in geo:
        x2.append(float(term.split(" ")[0]))
        y2.append(float(term.split(" ")[1]))

    x2 = np.array(x2)
    y2 = np.array(y2)
    x2 = (x2 - geo_t[0]) / geo_t[1]
    y2 = (y2 - geo_t[3]) / geo_t[5]
    x2 = x2.tolist()
    y2 = y2.tolist()

    ring = ogr.Geometry(ogr.wkbLinearRing)
    for id, p in enumerate(x2):
        ring.AddPoint(x2[id], y2[id])
    ring.AddPoint(x2[0], y2[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    path = r"D:\Cornell\EthiopianDrought\0ExperimentData\ExtractPixelCSV\Agg_Mask50.csv"
    Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Short", "SIF", "RF", 2009)

    ROW2 = ROW.tolist()
    COL2 = COL.tolist()
    valid = []

    for id, r in enumerate(ROW2):

        point2 = ogr.Geometry(ogr.wkbPoint)
        point2.AddPoint(COL2[id], ROW2[id])

        if poly.Contains(point2):

            valid.append(id)
    return valid


valid = GenValidRC()
Year = 2010
fig2 = plt.figure(constrained_layout=True,figsize=(10,12))
spec2 = gridspec.GridSpec(ncols=2, nrows=3, figure=fig2)
f2_ax1 = fig2.add_subplot(spec2[0, 0])
f2_ax2 = fig2.add_subplot(spec2[0, 1])
f2_ax3 = fig2.add_subplot(spec2[1, 0])
f2_ax4 = fig2.add_subplot(spec2[1, 1])
f2_ax5 = fig2.add_subplot(spec2[2, 0])
f2_ax6 = fig2.add_subplot(spec2[2, 1])

# 1 **********************
Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Short", "SIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])

f2_ax1.set_ylabel("Seasonly P ($mm$)",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax1.set_xlabel("Seasonly SIF ($mWm^{-2}nm^{-1}sr^{-1}$)",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax1.scatter(Var1List[valid],Var2List[valid])
f2_ax1.text(0.5, 0.95,'(a)'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = f2_ax1.transAxes)
f2_ax1.text(0.05, 0.8,'R={}\np={}'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':14},
     horizontalalignment='left',
     verticalalignment='center',
     transform = f2_ax1.transAxes)

#2 ****************************
Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Short", "SIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], PVIList[valid])
f2_ax2.set_ylabel("Seasonly PVI",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax2.set_xlabel("Seasonly SIF ($mWm^{-2}nm^{-1}sr^{-1}$)",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax2.text(0.5, 0.95,'(b)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = f2_ax2.transAxes)

f2_ax2.scatter(Var1List[valid],PVIList[valid])
f2_ax2.text(0.95, 0.8,'R={}\nP={}'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':14},
     horizontalalignment='right',
     verticalalignment='center',
     transform = f2_ax2.transAxes)

# 3 **********************
Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Long", "SIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])

f2_ax3.set_ylabel("Seasonly P ($mm$)",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax3.set_xlabel("Seasonly SIF ($mWm^{-2}nm^{-1}sr^{-1}$)",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax3.scatter(Var1List[valid],Var2List[valid])
f2_ax3.text(0.5, 0.95,'(c)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = f2_ax3.transAxes)
f2_ax3.text(0.05, 0.8,'R={}\nP={}'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':14},
     horizontalalignment='left',
     verticalalignment='center',
     transform = f2_ax3.transAxes)

#4 ****************************
Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Long", "SIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], PVIList[valid])
f2_ax4.set_ylabel("Seasonly PVI",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax4.set_xlabel("Seasonly SIF ($mWm^{-2}nm^{-1}sr^{-1}$)",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax4.text(0.5, 0.95,'(d)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='left',
     verticalalignment='center',
     transform = f2_ax4.transAxes)

f2_ax4.scatter(Var1List[valid],PVIList[valid])
f2_ax4.text(0.95, 0.8,'R={}\np={}'.format(round(r_value,3),round(p_value,3)),fontdict={'fontname':'Times New Roman','fontsize':14},
     horizontalalignment='right',
     verticalalignment='center',
     transform = f2_ax4.transAxes)

# 5

SP = []
SPV = []

SPVI = []
SPVIV = []

LP = []
LPV =[]

LPVI = []
LPVIV = []

for year in range(2007,2019):

    Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Short", "SIF", "RF", year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
    SP.append(r_value)
    SPV.append(p_value)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], PVIList[valid])
    SPVI.append(r_value)
    SPVIV.append(p_value)

    # long
    Var1List,Var2List,PVIList,COL,ROW = ExtractVIFromCSV(path, "Long", "SIF", "RF", year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], Var2List[valid])
    LP.append(r_value)
    LPV.append(p_value)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List[valid], PVIList[valid])
    LPVI.append(r_value)
    LPVIV.append(p_value)

print("SP")
[print(y+2007,round(SP[y],3)," p-value ",round(SPV[y],3),round(SPVI[y],3),' p_value ',round(SPVIV[y],3)) for y in range(12)]
print("LP")
[print(y+2007,round(LP[y],3)," p-value ",round(LPV[y],3),round(LPVI[y],3),' p_value ',round(LPVIV[y],3)) for y in range(12)]
x = np.arange(12)
f2_ax5.bar(x-0.2,SP,width=0.4,label="P")
f2_ax5.bar(x+0.2,SPVI,width=0.4,label="PVI")
f2_ax5.legend()
f2_ax5.text(0.5, 0.95,'(e)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = f2_ax5.transAxes)
f2_ax5.set_xlabel('Year(2007-2018)',fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax5.set_ylabel("R",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax5.set_xticks(range(0,12,3))
f2_ax5.set_xticklabels([str(2007+y)[2:] for y in range(0,12,3)])

# 6 **************
x = np.arange(12)
f2_ax6.bar(x-0.2,LP,width=0.4,label="P")
f2_ax6.bar(x+0.2,LPVI,width=0.4,label="PVI")
f2_ax6.legend()
f2_ax6.text(0.5, 0.95,'(f)',fontdict={'fontname':'Times New Roman','fontsize':16},
     horizontalalignment='center',
     verticalalignment='center',
     transform = f2_ax6.transAxes)
f2_ax6.set_xlabel('Year(2007-2018)',fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax6.set_ylabel("R",fontdict={'fontname':'Times New Roman','fontsize':16})
f2_ax6.set_xticks(range(0,12,3))
f2_ax6.set_xticklabels([str(2007+y)[2:] for y in range(0,12,3)])

plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.2,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_2.png',bbox_inches='tight',dpi=fig2.dpi,pad_inches=0.05)

plt.show()