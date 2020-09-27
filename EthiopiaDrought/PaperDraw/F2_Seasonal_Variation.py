import numpy as np
import pandas as pd
import csv,os
from matplotlib import pyplot as plt
from scipy import stats
path = r"D:\Cornell\EthiopianDrought\0ExperimentData\ExtractPixelCSV\Agg_Mask50.csv"

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
    PVIList = np.log(1+PVIList)
    return Var1List,Var2List,PVIList

fig = plt.figure(figsize=(12, 10))
plt.xticks([])
plt.yticks([])
plt.axis('off')
Year = 2009

ax = fig.add_subplot(3, 2, 1)
Var1List, Var2List, PVIList = monthAnalysis(path, "Short", "NSIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
ax.scatter(Var1List,Var2List)

ax2 = fig.add_subplot(3, 2, 2)
Var1List, Var2List, PVIList = monthAnalysis(path, "Short", "NSIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, PVIList)
ax2.scatter(Var1List,PVIList)

ax3 = fig.add_subplot(3, 2, 3)
Var1List, Var2List, PVIList = monthAnalysis(path, "Long", "NSIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
ax3.scatter(Var1List,Var2List)

ax4 = fig.add_subplot(3, 2, 4)
Var1List, Var2List, PVIList = monthAnalysis(path, "Long", "NSIF", "RF", Year)
slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, PVIList)
ax4.scatter(Var1List,PVIList)

SP = []
SPV = []

SPVI = []
SPVIV = []

LP = []
LPV =[]

LPVI = []
LPVIV = []

for year in range(2003,2019):

    Var1List, Var2List, PVIList = monthAnalysis(path, "Short", "NSIF", "RF", year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
    SP.append(r_value)
    SPV.append(p_value)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, PVIList)
    SPVI.append(r_value*-1)
    SPVIV.append(p_value)

    # long
    Var1List, Var2List, PVIList = monthAnalysis(path, "Long", "NSIF", "RF", year)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, Var2List)
    LP.append(r_value)
    LPV.append(p_value)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Var1List, PVIList)
    LPVI.append(r_value*-1)
    LPVIV.append(p_value)

ax5 = fig.add_subplot(3, 2, 5)
x = np.arange(16)
ax5.bar(x-0.2,SP,width=0.4,label="P")
ax5.bar(x+0.2,SPVI,width=0.4,label="PVI")
ax5.legend()

ax6 = fig.add_subplot(3, 2, 6)
x = np.arange(16)
ax6.bar(x-0.2,LP,width=0.4,label="P")
ax6.bar(x+0.2,LPVI,width=0.4,label="PVI")
ax6.legend()



plt.subplots_adjust(left=0.1,right=0.9,bottom=0,top=0.9,wspace=0.1,hspace=0.2)
plt.savefig(r'D:\Cornell\EthiopianDrought\0ExperimentData\Fig\Fig_2.png',bbox_inches='tight',dpi=fig.dpi,pad_inches=0.05)
plt.show()
