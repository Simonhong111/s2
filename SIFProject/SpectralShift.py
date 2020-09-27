import numpy as np
import xlrd
from scipy.signal import argrelextrema
from matplotlib import pyplot as plt
import pandas as pd
from  matplotlib import  cm
import  matplotlib
import os
def writecsv(filename,White_Reflectance = 0.73):
    """
    :param filename: xlsx 数据记录文件路径
    :return: 返回float类型的第31列之后的数据，这些数据都是数字
    """
    fwhm = []
    workbook = xlrd.open_workbook(filename)
    table = workbook.sheets()[0]
    rows = table.nrows
    HeadTitle = table.row_values(0)[13:]
    Dark = table.row_values(1)[13:]
    White = table.row_values(2)[13:]

    White = np.array(White).astype(np.float)
    Dark = np.array(Dark).astype(np.float)
    HeadTitle = np.array(HeadTitle).astype(np.float)
    White = (White-Dark)/White_Reflectance

    for idx, item in enumerate(White):
        if idx ==0:
            fwhm.append(1.5)
        else:
            fwhm.append("")


    # for row in range(table.nrows)[3:]:
    #     # print(table.row_values(row))
    #     dataset.append(table.row_values(row))


    dataframe = pd.DataFrame(
        {'WL': HeadTitle, 'Counts': White,"fwhm":fwhm})

    # 将DataFrame存储为csv,index表示是否显示行名，default=True
    dataframe.to_csv(os.path.join(r"D:\SpecCal-master\EXAMPLE_INPUT_DATA", "maifei.csv"), index=False, sep=';')
writecsv(r"D:\HUBEIFIELDDATA\20180711airline\49535149514848535049955750955057_20180711143813.77.xlsx")