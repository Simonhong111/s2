import numpy as np
import pandas as pd
import csv,os
from matplotlib import pyplot as plt
from scipy import stats

#
# def monthAnalysis(path):
#
#     data = pd.read_csv(path)
#
#     short_evi_name = ["Short"+"EVI"+str(year) for year in range(2003,2019)]
#     long_evi_name = ["Long" + "EVI" + str(year) for year in range(2003, 2019)]
#     short_rf_name = ["Short" + "RF" + str(year) for year in range(2003, 2019)]
#     long_rf_name = ["Long" + "RF" + str(year) for year in range(2003, 2019)]
#     short_pvi_name = ["Short" + "PVI" + str(year) for year in range(2003, 2019)]
#     long_pvi_name = ["Long" + "PVI" + str(year) for year in range(2003, 2019)]
#     short_gsif_name = ["Short" + "GSIF" + str(year) for year in range(2003, 2019)]
#     long_gsif_name = ["Long" + "GSIF" + str(year) for year in range(2003, 2019)]
#     short_nsif_name = ["Short" + "NSIF" + str(year) for year in range(2003, 2019)]
#     long_nsif_name = ["Long" + "NSIF" + str(year) for year in range(2003, 2019)]
#     short_Fre_name = ["Short" + "Fre"]
#     long_Fre_name = ["Long" + "Fre" ]
#
#     short_evi = data[short_evi_name].to_numpy()
#     long_evi =data[long_evi_name].to_numpy()
#     short_rf =data[short_rf_name].to_numpy()
#     long_rf =data[long_rf_name].to_numpy()
#     short_pvi =data[short_pvi_name].to_numpy()
#     long_pvi =data[long_pvi_name].to_numpy()
#     short_gsif =data[short_gsif_name].to_numpy()
#     long_gsif =data[long_gsif_name].to_numpy()
#     short_nsif =data[short_nsif_name].to_numpy()
#     long_nsif =data[long_nsif_name].to_numpy()
#     short_Fre =data[short_Fre_name].to_numpy()
#     long_Fre = data[long_Fre_name].to_numpy()
#
#
#     return short_evi ,long_evi ,short_rf ,long_rf,short_pvi,long_pvi ,short_gsif ,long_gsif,short_nsif ,long_nsif
def monthAnalysis(path,monthtype,var1,var2,year):
    if monthtype == "Short":
        months = [2,3,4,5]
    if monthtype == "Long":
        months =[6,7,8,9]
    data = pd.read_csv(path)
    key1 = [var1+str(year)+str(months[0]).zfill(2),var1+str(year)+str(months[1]).zfill(2),
            var1+str(year)+str(months[2]).zfill(2),var1+str(year)+str(months[3]).zfill(2)]
    key2 = [var2 + str(year) + str(months[0]).zfill(2), var2 + str(year) + str(months[1]).zfill(2),
            var2 + str(year) + str(months[2]).zfill(2), var2 + str(year) + str(months[3]).zfill(2)]

    Var1List = data[key1].to_numpy().mean(axis=1)
    Var2List = data[key2].to_numpy().mean(axis=1)
    PVIList = data[monthtype+"PVI"+str(year)].to_numpy()
    return Var1List,Var2List,PVIList

Years = [str(year) for year in range(2003,2019)]

short_evi_rf = []
short_nsif_rf = []
short_gsif_rf = []
short_evi_pvi = []
short_nsif_pvi = []
short_gsif_pvi = []

long_evi_rf = []
long_nsif_rf = []
long_gsif_rf = []
long_evi_pvi = []
long_nsif_pvi = []
long_gsif_pvi = []


path = r"D:\Cornell\EthiopianDrought\CropCSV\DailyCrop\PolyGonAgg_Mask70.csv"
for Year in Years:
    Var1List1, Var2List1, Var3List1 = monthAnalysis(path, "Short", "EVI", "RF", Year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List1, Var2List1)
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(Var1List1, Var3List1)
    short_evi_rf.append(r_value)
    short_evi_pvi.append(r_value2*(-1))

    Var1List2, Var2List2, Var3List2 = monthAnalysis(path, "Short", "GSIF", "RF", Year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List2, Var2List2)
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(Var1List2, Var3List2)
    short_gsif_rf.append(r_value)
    short_gsif_pvi.append(r_value2*(-1))

    Var1List3, Var2List3, Var3List3 = monthAnalysis(path, "Short", "NSIF", "RF", Year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List3, Var2List3)
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(Var1List3, Var3List3)
    short_nsif_rf.append(r_value)
    short_nsif_pvi.append(r_value2*(-1))

    Var1List4, Var2List4, Var3List4 = monthAnalysis(path, "Long", "EVI", "RF", Year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List4, Var2List4)
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(Var1List4, Var3List4)
    long_evi_rf.append(r_value)
    long_evi_pvi.append(r_value2*(-1))

    Var1List5, Var2List5, Var3List5 = monthAnalysis(path, "Long", "GSIF", "RF", Year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List5, Var2List5)
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(Var1List5, Var3List5)
    long_gsif_rf.append(r_value)
    long_gsif_pvi.append(r_value2*(-1))

    Var1List6, Var2List6, Var3List6 = monthAnalysis(path, "Long", "NSIF", "RF", Year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List6, Var2List6)
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(Var1List6, Var3List6)
    long_nsif_rf.append(r_value)
    long_nsif_pvi.append(r_value2*(-1))
#
mask = "threshold 70%"
fig = plt.figure(1)
Years = range(2003,2019)
plt.title("2003-2018 short rains " + "Correlation Coefficient for {}\n".format(mask), fontsize=16)

print('short evi rf', short_evi_rf)
print('short evi pvi',short_evi_pvi)
print("short gsif rf",short_gsif_rf)
print("short gsif pvi",short_gsif_pvi)
print("short nsif rf",short_nsif_rf)
print("short nsif pvi",short_nsif_pvi)
plt.ylabel("R",fontsize=16)
plt.grid(b=True)
plt.plot(Years,short_evi_rf,"g",label="evi vs rainfall")
plt.plot(Years,short_evi_pvi,"g--",label="evi vs pvi")
plt.plot(Years,short_gsif_rf,"k",label="gsif vs rainfall")
plt.plot(Years,short_gsif_pvi,"k--",label="gsif vs pvi")
plt.plot(Years,short_nsif_rf,"r",label="nsif vs rainfall")
plt.plot(Years,short_nsif_pvi,"r--",label="nsif vs pvi")
plt.legend()


fig2 = plt.figure(2)

Years = range(2003,2019)
plt.title("2003-2018 long rains " + "Correlation Coefficient for {}\n".format(mask), fontsize=16)
print('long evi rf', long_evi_rf)
print('long evi pvi',long_evi_pvi)
print("long gsif rf",long_gsif_rf)
print("long gsif pvi",long_gsif_pvi)
print("long nsif rf",long_nsif_rf)
print("long nsif pvi",long_nsif_pvi)
plt.grid(b=True)
plt.plot(Years,long_evi_rf,"g",label="evi vs rainfall")
plt.plot(Years,long_evi_pvi,"g--",label="evi vs pvi")
plt.plot(Years,long_gsif_rf,"k",label="gsif vs rainfall")
plt.plot(Years,long_gsif_pvi,"k--",label="gsif vs pvi")
plt.plot(Years,long_nsif_rf,"r",label="nsif vs rainfall")
plt.plot(Years,long_nsif_pvi,"r--",label="nsif vs pvi")
plt.ylabel("R",fontsize=16)
plt.legend()

plt.show()






