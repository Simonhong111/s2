import numpy as np
import xlrd,os
from scipy.signal import argrelextrema
from matplotlib import pyplot as plt
from  matplotlib import  cm
import  matplotlib
import pandas as pd
from pylab import mpl
from scipy import stats
from sklearn.linear_model import LinearRegression
mpl.rcParams['font.sans-serif'] = ['FangSong']
mpl.rcParams['axes.unicode_minus'] = False
global pos
def Resultload(resultFilePath):
    """
    :param resultFilePath:download the input data to cal or val the model
    :return:
    """
    dataset = []
    workbook = xlrd.open_workbook(resultFilePath)
    table = workbook.sheets()[0]
    rows = table.nrows
    HeadTitle = table.row_values(0)
    for row in range(table.nrows)[1:]:
        dataset.append(table.row_values(row))
    mData  = np.array([d for d in dataset]).astype(np.float)
    return mData


def getPara(Dataset,Field,Month,Day,Time):
    """
    :param Dataset: get the input of the model from VIS
    :param Field:  the Conditions that the input data needs to meet
    :param Month:   the conditions
    :param Day:     the conditios
    :param Time:    the conditions
    :return:
    """

    VI = []
    Area = []
    for data in Dataset:

        temp = []
        # print("***",data[0]  == Field and data[1] == Month and data[2] == Day and data[3] == Time)
        if data[1] == Month and data[2] == Day and data[3] == Time:

            temp.append(data[6])

            for i in range(12):
                temp.append(data[7+5*i])
            VI.append(temp)
            Area.append(data[5])
    return np.array(VI).astype(np.float),np.array(Area).astype(np.float)

def getPara2019(Dataset,Month,Day):
    """
    :param Dataset: get the input VIS in 2019 to val the model
    :param Month:
    :param Day:
    :return:
    """
    VI = []
    Area = []
    for data in Dataset:
        temp = []


        if data[1] == Month and data[2] == Day:
            temp.append(data[6])

            for i in range(12):
                temp.append(data[7+i])
            VI.append(temp)
            Area.append(data[5])
    return np.array(VI).astype(np.float),np.array(Area).astype(np.float)



FitFilePath = r"D:\Satellive\Estimation\maifei.xlsx"
FitFilePath2 = r"D:\Satellive\Estimation\maifeiIncludingA31A32-1120A5.xlsx"
# TestPath = r"D:\Satellive\Estimation\maifeiIncludingA31A32.xlsx"
TestPath = r"D:\Satellive\Estimation\maifeiIncludingA31A32-1120.xlsx"
# TestPath = r"D:\Satellive\Estimation\maifeiIncludingA31A32-1120A5.xlsx"
FitDataset = Resultload(FitFilePath)
FitDataset2 = Resultload(FitFilePath2)
TestDataset = Resultload(TestPath)


FitData,FitArea  =  getPara(FitDataset,2,9,5,93858)
# FitData,FitArea2  =  getPara2019(FitDataset2,7,18)
# print(FitArea)
# FitData = np.vstack((FitData,FitData2))

# FitArea = np.vstack((FitArea,FitArea2))
# print("fitdata",FitData)
# print("fitdata2",FitData2)
# FitData  =  getPara2019(FitDataset,6,30)
TestData,TestArea =  getPara2019(TestDataset,7,22)

model = LinearRegression()
model.fit(FitData[:,[3,5,9]] , FitData[:,0] )
predictions = model.predict(TestData[:,[3,5,9]])

# #
from sklearn import linear_model
# model = linear_model.Lasso(alpha=0.1)
# model.fit(FitData[:,[3,5,9]] , FitData[:,0] )
# predictions = model.predict(TestData[:,[3,5,9]])
# print(model.coef_)


# predictions = model.predict(TestData[:,[4,5,9]])

for i, prediction in enumerate(predictions):
    print('Predicted: %s, Target: %s and Difference %s' % (prediction , TestData[:,0][i],TestData[:,0][i]-prediction))

area = TestArea
ObservedCropYield = TestData[:,0]*TestArea
predictedCropYield = predictions * area

for idx,_ in enumerate(predictedCropYield):
    print(predictedCropYield[idx], ObservedCropYield[idx],(ObservedCropYield[idx] -predictedCropYield[idx]))


s, inter, r, p, std_err = stats.linregress(predictedCropYield,ObservedCropYield)
#
print("**",s,r,p)
plt.scatter(predictedCropYield,ObservedCropYield)
xmin = predictedCropYield.min()
xmax = predictedCropYield.max()
plt.plot([xmin,xmax],[s*xmin+inter,s*xmax+inter],"r")
plt.title("预测产量与实际产量关系图",fontsize = 20)
plt.xlabel("预测产量值 (kg)",fontsize =18)
plt.ylabel("实际产量值 (kg)",fontsize = 18)
plt.text(2610,2640,"Y = {} * X {}".format(round(s,3),round(inter,3)),fontsize =18)
plt.text(1000,3800,"P < 0.05",fontsize =18)
plt.text(1000,3900,"R = {}".format(round(r,3)),fontsize =18)
plt.show()
