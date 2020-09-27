import numpy as np
import xlrd,os
from scipy.signal import argrelextrema
from matplotlib import pyplot as plt
from  matplotlib import  cm
import  matplotlib
from pylab import mpl

mpl.rcParams['font.sans-serif'] = ['FangSong']
mpl.rcParams['axes.unicode_minus'] = False
global pos
def Resultload(resultFilePath):

    dataset = []
    workbook = xlrd.open_workbook(resultFilePath)
    table = workbook.sheets()[0]
    rows = table.nrows
    HeadTitle = table.row_values(0)
    for row in range(table.nrows)[1:]:
        dataset.append(table.row_values(row))
    mData  = np.array([d for d in dataset]).astype(np.float)
    return mData
def ALLday(Dataset,field):
    slope2 = []
    intercept2 = []
    rvalue2 = []
    pvalue2 = []
    time2 = []


    for Data in Dataset:
        slope = []
        intercept = []
        rvalue = []
        pvalue = []
        time = []
        if Data[0] == field:
            for s in range(12):
                slope.append(Data[s * 5 + 7])
            for inter in range(12):
                intercept.append(Data[inter * 5 + 8])
            for r in range(12):
                rvalue.append(Data[r * 5 + 9])
            for p in range(12):
                pvalue.append(Data[p * 5 + 10])
            time.append(Data[3])

            slope2.append(slope)
            intercept2.append(intercept)
            rvalue2.append(rvalue)
            pvalue2.append(pvalue)
            time2.append(time)


    slope3 = []
    intercept3 = []
    rvalue3 = []
    pvalue3 = []
    time3 = []
    for s in slope2:
        if not s in slope3 and len(s):
            slope3.append(s)
    for intr in intercept2:
        if not intr in intercept3:
            intercept3.append(intr)
    for r in rvalue2:
        if not r in rvalue3:
            rvalue3.append(r)
    for p in pvalue2:
        if not p in pvalue3:
            pvalue3.append(p)
    for t in time2:
        if not t in time3:
            time3.append(t)

    slope3 = np.array(slope3)
    intercept3 = np.array(intercept3)
    rvalue3 = np.array(rvalue3)
    pvalue3 = np.array(pvalue3)
    time3 = np.array(time3)
    print(slope3.shape)



    timestr = ""

    for tt in time3:

        if len(str(tt[0])) == 7:
            timestr = timestr + str(tt[0])[0] + ":" + str(tt[0])[1]+str(tt[0])[2] + " "

        else:
            timestr = timestr + str(tt[0])[0]+str(tt[0])[1] + ":" + str(tt[0])[2] + str(tt[0])[3]+ " "

    timeX = np.array(range(len(time3)))
    # timeX = range(len(time3))
    Label = ["FLD 斜率","FLD2 斜率","3FLD 斜率","3FLD2 斜率","NDVI 斜率","EVI 斜率","GRVI 斜率",
             "GNDVI 斜率","MSAVI 斜率","OSAVI 斜率","SAVI 斜率","WDRVI 斜率"]
    Label2 = ["FLD 截距", "FLD2 截距", "3FLD 截距", "3FLD2 截距", "NDVI 截距", "EVI 截距", "GRVI 截距",
             "GNDVI 截距", "MSAVI 截距", "OSAVI 截距", "SAVI 截距", "WDRVI 截距"]
    Label3 = ["FLD R", "FLD2 R", "3FLD R", "3FLD2 R", "NDVI R", "EVI R", "GRVI R",
              "GNDVI R", "MSAVI R", "OSAVI R", "SAVI R", "WDRVI R"]

    Label4 = ["FLD P", "FLD2 P", "3FLD P", "3FLD2 P", "NDVI P", "EVI P", "GPVI P",
              "GNDVI P", "MSAVI P", "OSAVI P", "SAVI P", "WDPVI P"]
    Label5 = ["FLD 比例", "FLD2 比例", "3FLD 比例", "3FLD2 比例", "NDVI 比例", "EVI 比例", "GRVI 比例",
              "GNDVI 比例", "MSAVI 比例", "OSAVI 比例", "SAVI 比例", "WDRVI 比例"]

    fig1 = plt.figure(1)
    # plt.title("时间点是{}".format(timestr), fontsize=20)
    # print(time3)
    # print(slope3)
    for idx,label in enumerate(Label):
        fig1.add_subplot(3, 4, idx+1)
        plt.plot(timeX,slope3[:,idx],label= label)
        plt.legend()
    # fig1.savefig(os.path.join(r"C:\Users\zmhwh\Desktop\Satellive",str(month) + str(day) + "slope" + ".png"))
    # plt.close(fig1)
    fig2 = plt.figure(2)
    # plt.title("时间点是{}".format(timestr), fontsize=20)
    for idx,label in enumerate(Label2):
        fig2.add_subplot(3, 4, idx+1)
        # print(timeX)
        # print(intercept3[:,idx])
        plt.plot(timeX,intercept3[:,idx],label= label)
        plt.legend()

    # fig2.savefig(os.path.join(r"C:\Users\zmhwh\Desktop\Satellive", str(month) + str(day) + "intercept" + ".png"))
    # plt.close(fig2)
    fig3 = plt.figure(3)
    # plt.title("时间点是{}".format(timestr), fontsize=20)
    for idx,label in enumerate(Label3):
        fig3.add_subplot(3, 4, idx+1)
        plt.plot(timeX,rvalue3[:,idx],label= label)
        plt.legend()
    # fig3.savefig(os.path.join(r"C:\Users\zmhwh\Desktop\Satellive", str(month) + str(day) + "R" + ".png"))
    # plt.close(fig3)
    fig4 = plt.figure(4)
    # plt.title("时间点是{}".format(timestr), fontsize=20)
    for idx,label in enumerate(Label4):
        fig4.add_subplot(3, 4, idx+1)
        plt.plot(timeX,pvalue3[:,idx],label= label)
        plt.plot([timeX.min(),timeX.max()],[0.05,0.05])
        plt.legend()
    # fig4.savefig(os.path.join(r"C:\Users\zmhwh\Desktop\Satellive", str(month) + str(day) + "Pvalue" + ".png"))
    # plt.close(fig4)
    fig5 = plt.figure(5)

    # plt.title("时间点是{}".format(timestr), fontsize=20)
    for idx, label in enumerate(Label5):
        fig5.add_subplot(3, 4, idx + 1)
        plt.ylim(0, 100)
        L = len(pvalue3[:, idx])
        Tg = len(np.where(pvalue3[:, idx] < 0.05)[0])/L*100
        NTg = len(np.where(pvalue3[:, idx] >= 0.05)[0])/L*100

        plt.bar([1,2], [Tg,NTg], label=label,color =['g','b'])
        plt.xticks([1,2],["通过","未通过"])

        plt.legend()
    plt.show()
def plotD(resultFilePath,field):
    Dataset = Resultload(resultFilePath)
    ALLday(Dataset,field)

resultFilePath = r"D:\Satellive\Estimation\maifei.xlsx"

# Dataset = Resultload(resultFilePath)
# MD = [[d[1],d[2]] for d in Dataset]
# MD2 = []
# for md in MD:
#     if not md in MD2:
#         MD2.append(md)
# print(MD2)

# for md in MD2:
#     plotD(resultFilePath,md[0],md[1])

plotD(resultFilePath,2)