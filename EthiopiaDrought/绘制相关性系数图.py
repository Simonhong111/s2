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

# path = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\PolyGonAgg_Mask40_Mean.csv"
#
# short_evi40 ,long_evi40 ,short_rf40 ,long_rf40,short_pvi40,long_pvi40 ,short_gsif40 ,long_gsif40,short_nsif40 ,long_nsif40=monthAnalysis(path)
#
# path = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\PolyGonAgg_Mask50_Mean.csv"
#
# short_evi50 ,long_evi50 ,short_rf50 ,long_rf50,short_pvi50,long_pvi50 ,short_gsif50 ,long_gsif50,short_nsif50 ,long_nsif50=monthAnalysis(path)
#
# path = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\PolyGonAgg_Mask60_Mean.csv"
#
# short_evi60 ,long_evi60 ,short_rf60 ,long_rf60,short_pvi60,long_pvi60 ,short_gsif60 ,long_gsif60,short_nsif60 ,long_nsif60=monthAnalysis(path)
#
# path = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\PolyGonAgg_Mask70_Mean.csv"
#
# short_evi70 ,long_evi70 ,short_rf70 ,long_rf70,short_pvi70,long_pvi70 ,short_gsif70 ,long_gsif70,short_nsif70 ,long_nsif70=monthAnalysis(path)
#
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
#
#
# short_evi_rf5 = []
# short_nsif_rf5 = []
# short_gsif_rf5 = []
# short_evi_pvi5 = []
# short_nsif_pvi5 = []
# short_gsif_pvi5 = []
#
# long_evi_rf5 = []
# long_nsif_rf5 = []
# long_gsif_rf5 = []
# long_evi_pvi5 = []
# long_nsif_pvi5 = []
# long_gsif_pvi5 = []
#
# short_evi_rf6 = []
# short_nsif_rf6 = []
# short_gsif_rf6 = []
# short_evi_pvi6 = []
# short_nsif_pvi6 = []
# short_gsif_pvi6 = []
#
# long_evi_rf6 = []
# long_nsif_rf6 = []
# long_gsif_rf6 = []
# long_evi_pvi6 = []
# long_nsif_pvi6 = []
# long_gsif_pvi6 = []
#
# short_evi_rf7 = []
# short_nsif_rf7 = []
# short_gsif_rf7 = []
# short_evi_pvi7 = []
# short_nsif_pvi7 = []
# short_gsif_pvi7 = []
#
# long_evi_rf7 = []
# long_nsif_rf7 = []
# long_gsif_rf7 = []
# long_evi_pvi7 = []
# long_nsif_pvi7 = []
# long_gsif_pvi7 = []
#
# for idx,year in enumerate(Years):
#
#     slope, intercept, r_value, p_value, std_err = stats.linregress(short_evi40[:,idx], short_rf40[:,idx])
#     slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(long_evi40[:,idx], long_rf40[:,idx])
#     slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(short_evi50[:,idx], short_rf50[:,idx])
#     slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(long_evi50[:,idx], long_rf50[:,idx])
#     slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(short_evi60[:,idx], short_rf60[:,idx])
#     slope6, intercept6, r_value6, p_value6, std_err6 = stats.linregress(long_evi60[:,idx], long_rf60[:,idx])
#     slope7, intercept7, r_value7, p_value7, std_err7 = stats.linregress(short_evi70[:,idx], short_rf70[:,idx])
#     slope8, intercept8, r_value8, p_value8, std_err8 = stats.linregress(long_evi70[:,idx], long_rf70[:,idx])
#     short_evi_rf.append(r_value)
#     long_evi_rf.append(r_value2)
#     short_evi_rf5.append(r_value3)
#     long_evi_rf5.append(r_value4)
#     short_evi_rf6.append(r_value5)
#     long_evi_rf6.append(r_value6)
#     short_evi_rf7.append(r_value7)
#     long_evi_rf7.append(r_value8)
#
#     slope, intercept, r_value, p_value, std_err = stats.linregress(short_evi40[:,idx], short_pvi40[:,idx])
#     slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(long_evi40[:,idx], long_pvi40[:,idx])
#     slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(short_evi50[:,idx], short_pvi50[:,idx])
#     slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(long_evi50[:,idx], long_pvi50[:,idx])
#     slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(short_evi60[:,idx], short_pvi60[:,idx])
#     slope6, intercept6, r_value6, p_value6, std_err6 = stats.linregress(long_evi60[:,idx], long_pvi60[:,idx])
#     slope7, intercept7, r_value7, p_value7, std_err7 = stats.linregress(short_evi70[:,idx], short_pvi70[:,idx])
#     slope8, intercept8, r_value8, p_value8, std_err8 = stats.linregress(long_evi70[:,idx], long_pvi70[:,idx])
#     short_evi_pvi.append(r_value)
#     long_evi_pvi.append(r_value2)
#     short_evi_pvi5.append(r_value3)
#     long_evi_pvi5.append(r_value4)
#     short_evi_pvi6.append(r_value5)
#     long_evi_pvi6.append(r_value6)
#     short_evi_pvi7.append(r_value7)
#     long_evi_pvi7.append(r_value8)
#
#     slope, intercept, r_value, p_value, std_err = stats.linregress(short_gsif40[:,idx], short_rf40[:,idx])
#     slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(long_gsif40[:,idx], long_rf40[:,idx])
#     slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(short_gsif50[:,idx], short_rf50[:,idx])
#     slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(long_gsif50[:,idx], long_rf50[:,idx])
#     slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(short_gsif60[:,idx], short_rf60[:,idx])
#     slope6, intercept6, r_value6, p_value6, std_err6 = stats.linregress(long_gsif60[:,idx], long_rf60[:,idx])
#     slope7, intercept7, r_value7, p_value7, std_err7 = stats.linregress(short_gsif70[:,idx], short_rf70[:,idx])
#     slope8, intercept8, r_value8, p_value8, std_err8 = stats.linregress(long_gsif70[:,idx], long_rf70[:,idx])
#     short_gsif_rf.append(r_value)
#     long_gsif_rf.append(r_value2)
#     short_gsif_rf5.append(r_value3)
#     long_gsif_rf5.append(r_value4)
#     short_gsif_rf6.append(r_value5)
#     long_gsif_rf6.append(r_value6)
#     short_gsif_rf7.append(r_value7)
#     long_gsif_rf7.append(r_value8)
#
#     slope, intercept, r_value, p_value, std_err = stats.linregress(short_gsif40[:,idx], short_pvi40[:,idx])
#     slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(long_gsif40[:,idx], long_pvi40[:,idx])
#     slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(short_gsif50[:,idx], short_pvi50[:,idx])
#     slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(long_gsif50[:,idx], long_pvi50[:,idx])
#     slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(short_gsif60[:,idx], short_pvi60[:,idx])
#     slope6, intercept6, r_value6, p_value6, std_err6 = stats.linregress(long_gsif60[:,idx], long_pvi60[:,idx])
#     slope7, intercept7, r_value7, p_value7, std_err7 = stats.linregress(short_gsif70[:,idx], short_pvi70[:,idx])
#     slope8, intercept8, r_value8, p_value8, std_err8 = stats.linregress(long_gsif70[:,idx], long_pvi70[:,idx])
#     short_gsif_pvi.append(r_value)
#     long_gsif_pvi.append(r_value2)
#     short_gsif_pvi5.append(r_value3)
#     long_gsif_pvi5.append(r_value4)
#     short_gsif_pvi6.append(r_value5)
#     long_gsif_pvi6.append(r_value6)
#     short_gsif_pvi7.append(r_value7)
#     long_gsif_pvi7.append(r_value8)
#
#     slope, intercept, r_value, p_value, std_err = stats.linregress(short_nsif40[:,idx], short_rf40[:,idx])
#     slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(long_nsif40[:,idx], long_rf40[:,idx])
#     slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(short_nsif50[:,idx], short_rf50[:,idx])
#     slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(long_nsif50[:,idx], long_rf50[:,idx])
#     slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(short_nsif60[:,idx], short_rf60[:,idx])
#     slope6, intercept6, r_value6, p_value6, std_err6 = stats.linregress(long_nsif60[:,idx], long_rf60[:,idx])
#     slope7, intercept7, r_value7, p_value7, std_err7 = stats.linregress(short_nsif70[:,idx], short_rf70[:,idx])
#     slope8, intercept8, r_value8, p_value8, std_err8 = stats.linregress(long_nsif70[:,idx], long_rf70[:,idx])
#     short_nsif_rf.append(r_value)
#     long_nsif_rf.append(r_value2)
#     short_nsif_rf5.append(r_value3)
#     long_nsif_rf5.append(r_value4)
#     short_nsif_rf6.append(r_value5)
#     long_nsif_rf6.append(r_value6)
#     short_nsif_rf7.append(r_value7)
#     long_nsif_rf7.append(r_value8)
#
#     slope, intercept, r_value, p_value, std_err = stats.linregress(short_nsif40[:,idx], short_pvi40[:,idx])
#     slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(long_nsif40[:,idx], long_pvi40[:,idx])
#     slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(short_nsif50[:,idx], short_pvi50[:,idx])
#     slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(long_nsif50[:,idx], long_pvi50[:,idx])
#     slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(short_nsif60[:,idx], short_pvi60[:,idx])
#     slope6, intercept6, r_value6, p_value6, std_err6 = stats.linregress(long_nsif60[:,idx], long_pvi60[:,idx])
#     slope7, intercept7, r_value7, p_value7, std_err7 = stats.linregress(short_nsif70[:,idx], short_pvi70[:,idx])
#     slope8, intercept8, r_value8, p_value8, std_err8 = stats.linregress(long_nsif70[:,idx], long_pvi70[:,idx])
#     short_nsif_pvi.append(r_value)
#     long_nsif_pvi.append(r_value2)
#     short_nsif_pvi5.append(r_value3)
#     long_nsif_pvi5.append(r_value4)
#     short_nsif_pvi6.append(r_value5)
#     long_nsif_pvi6.append(r_value6)
#     short_nsif_pvi7.append(r_value7)
#     long_nsif_pvi7.append(r_value8)
#
# fig = plt.figure(figsize=(15, 10))
# plt.title("correlation coefficient for short rains, 2003 - 2018 \n")
# plt.xticks([])
# plt.yticks([])
#
# ax1 = fig.add_subplot(2, 3, 1)
# ax1.set_title("evi vs rainfall")
# ax1.grid(b=True)
# ax1.plot(range(2003,2019),short_evi_rf,label="threshold 40")
# ax1.plot(range(2003, 2019), short_evi_rf5, label="threshold 50")
# ax1.plot(range(2003, 2019), short_evi_rf6, label="threshold 60")
# ax1.plot(range(2003, 2019), short_evi_rf7, label="threshold 70")
# ax1.set_xticks(range(2003,2019))
# ax1.set_xticklabels(Years)
# # ax1.legend()
#
# ax2 = fig.add_subplot(2, 3, 2)
# ax2.set_title("evi vs pvi")
# ax2.grid(b=True)
# ax2.plot(range(2003, 2019), short_evi_pvi, label="threshold 40")
# ax2.plot(range(2003, 2019), short_evi_pvi5, label="threshold 50")
# ax2.plot(range(2003, 2019), short_evi_pvi6, label="threshold 60")
# ax2.plot(range(2003, 2019), short_evi_pvi7, label="threshold 70")
# ax2.set_xticks(range(2003, 2019))
# ax2.set_xticklabels(Years)
# # ax2.legend()
#
# ax3 = fig.add_subplot(2, 3, 3)
# ax3.set_title("gsif vs rainfall")
# ax3.grid(b=True)
# ax3.plot(range(2003, 2019), short_gsif_rf, label="threshold 40")
# ax3.plot(range(2003, 2019), short_gsif_rf5, label="threshold 50")
# ax3.plot(range(2003, 2019), short_gsif_rf6, label="threshold 60")
# ax3.plot(range(2003, 2019), short_gsif_rf7, label="threshold 70")
# ax3.set_xticks(range(2003, 2019))
# ax3.set_xticklabels(Years)
# # ax3.legend()
#
#
# ax4 = fig.add_subplot(2, 3, 4)
# ax4.set_title("gsif vs pvi")
# ax4.grid(b=True)
# ax4.plot(range(2003, 2019), short_gsif_pvi, label="threshold 40")
# ax4.plot(range(2003, 2019), short_gsif_pvi5, label="threshold 50")
# ax4.plot(range(2003, 2019), short_gsif_pvi6, label="threshold 60")
# ax4.plot(range(2003, 2019), short_gsif_pvi7, label="threshold 70")
# ax4.set_xticks(range(2003, 2019))
# ax4.set_xticklabels(Years)
# # ax4.legend()
#
# ax5 = fig.add_subplot(2, 3, 5)
# ax5.set_title("nsif vs rainfall")
# ax5.grid(b=True)
# ax5.plot(range(2003, 2019), short_nsif_rf, label="threshold 40")
# ax5.plot(range(2003, 2019), short_nsif_rf5, label="threshold 50")
# ax5.plot(range(2003, 2019), short_nsif_rf6, label="threshold 60")
# ax5.plot(range(2003, 2019), short_nsif_rf7, label="threshold 70")
# ax5.set_xticks(range(2003, 2019))
# ax5.set_xticklabels(Years)
# # ax5.legend()
#
# ax6 = fig.add_subplot(2, 3, 6)
# ax6.set_title("nsif vs pvi")
# ax6.grid(b=True)
# ax6.plot(range(2003, 2019), short_nsif_pvi, label="threshold 40")
# ax6.plot(range(2003, 2019), short_nsif_pvi5, label="threshold 50")
# ax6.plot(range(2003, 2019), short_nsif_pvi6, label="threshold 60")
# ax6.plot(range(2003, 2019), short_nsif_pvi7, label="threshold 70")
# ax6.set_xticks(range(2003, 2019))
# ax6.set_xticklabels(Years)
# ax6.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0.)
#
# # fig.tight_layout()#调整整体空白
# # plt.subplots_adjust(wspace =0, hspace =0.3)#调整子图间距
# plt.show()

path = r"D:\Cornell\EthiopianDrought\CropCSV\Crop\new\PolyGonAgg_Mask70.csv"
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






