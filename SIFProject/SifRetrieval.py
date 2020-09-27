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
def fileload(filename,White_Reflectance = 0.73):
    """
    :param filename: xlsx 数据记录文件路径
    :return: 返回float类型的第31列之后的数据，这些数据都是数字
    """
    dataset = []
    workbook = xlrd.open_workbook(filename)
    table = workbook.sheets()[0]
    rows = table.nrows
    HeadTitle = table.row_values(0)
    Dark = table.row_values(1)
    White = table.row_values(2)
    for row in range(table.nrows)[3:]:
        # print(table.row_values(row))
        dataset.append(table.row_values(row))
    mHT = np.array(HeadTitle[30:]).astype(np.float)
    mDark = np.array(Dark[30:]).astype(np.float)
    mWhite = (np.array(White[30:]).astype(np.float) - mDark)/White_Reflectance
    mData  = np.array([d[30:] for d in dataset]).astype(np.float) - mDark

    return mData,mHT,mWhite
def fileload2(filename):
    """
    :param filename: xlsx 数据记录文件路径
    :return: 返回float类型的第31列之后的数据，这些数据都是数字
    """
    dataset = []
    workbook = xlrd.open_workbook(filename)
    table = workbook.sheets()[0]
    rows = table.nrows
    HeadTitle = table.row_values(0)
    Dark = table.row_values(1)
    White = table.row_values(2)
    for row in range(table.nrows)[3:]:
        # print(table.row_values(row))
        dataset.append(table.row_values(row))
    mHT = np.array(HeadTitle)
    mData  = np.array(dataset)

    return mData,mHT


def RtrsFLD(Data,HeadTitle,White,LeftStart =755,LeftEnd = 759,InStart=759,InEnd = 764):

   """
   the retrieval of SIF at O-A absorption band at 760 nm
   :param Ein:Incident Radiance inside the absorption window 
   :param Lin: Reflected Radiance inside the absorption window
   :param Eout: Incident Radiance outside the absorption window
   :param Lout: Reflected Radiance outside the absorption window
   :return: sif
   """

   LeftExtent,LeftIndex = getMaxValue(HeadTitle,White,Start=LeftStart,End=LeftEnd)
   Eout, Lout = White[LeftExtent][LeftIndex],Data[LeftExtent][LeftIndex]

   InExtent, InIndex = getMinValue(HeadTitle, White,Start=InStart,End=InEnd)
   Ein, Lin  =  White[InExtent][InIndex],Data[InExtent][InIndex]

   R = (Lout - Lin)/(Eout - Ein)

   F = (Eout*Lin - Ein*Lout)/(Eout - Ein)

   # print("FLD Retrieval Wave ：Out {} and In {}".format(HeadTitle[LeftExtent][LeftIndex], HeadTitle[InExtent][InIndex]))
   # print("FLD Retrieval Eout {},Lout {},Ein {},Lin {}".format(Eout,Lout,Ein,Lin))
   # print("FLD Retrieval Reflectance：{} and SIF {}".format(R,F))

   return F


def RtrsFLD_AVE(Data, HeadTitle, White,LeftStart= 755,LeftEnd=759,InStart=759,InEnd = 764):
    """
    the retrieval of SIF at O-A absorption band at 760 nm
    :param Ein:Incident Radiance inside the absorption window
    :param Lin: Reflected Radiance inside the absorption window
    :param Eout: Incident Radiance outside the absorption window
    :param Lout: Reflected Radiance outside the absorption window
    :return: sif
    """

    InExtent, InIndex = getMinValue(HeadTitle, White, Start=InStart, End=InEnd)
    Ein, Lin = White[InExtent][InIndex], Data[InExtent][InIndex]

    smooth_White = smooth(White)
    smooth_Data = smooth(Data)

    LeftExtent, LeftIndex = getMaxValue(HeadTitle, smooth_White, Start=LeftStart, End=LeftEnd)
    Eout, Lout = smooth_White[LeftExtent][LeftIndex], smooth_Data[LeftExtent][LeftIndex]

    # print("FLD_AVE Retr wave center at {},Eout {} and Lout {}".format(HeadTitle[LeftExtent][LeftIndex], Eout,Lout))


    R = (Lout - Lin) / (Eout - Ein)

    F = (Eout * Lin - Ein * Lout) / (Eout - Ein)

    # print("FLD Retrieval Wave ：Out {} and In {}".format(HeadTitle[LeftExtent][LeftIndex], HeadTitle[InExtent][InIndex]))
    # print("FLD Retrieval Eout {},Lout {},Ein {},Lin {}".format(Eout, Lout, Ein, Lin))
    # print("FLD Retrieval Reflectance：{} and SIF {}".format(R, F))

    return F

def Rtr3FLD(Data, HeadTitle, White,LeftStart= 755,LeftEnd=759,InStart=759,InEnd = 764,RightStart= 755,RightEnd=759):

    LeftExtent, LeftIndex = getMaxValue(HeadTitle, White, Start=LeftStart, End=LeftEnd)
    EL, LL = White[LeftExtent][LeftIndex], Data[LeftExtent][LeftIndex]

    InExtent, InIndex = getMinValue(HeadTitle, White, Start=InStart, End=InEnd)
    Ein, Lin = White[InExtent][InIndex], Data[InExtent][InIndex]

    RightExtent, RightIndex = getMinValue(HeadTitle, White, Start=RightStart, End=RightEnd)
    ER, LR = White[RightExtent][RightIndex], Data[RightExtent][RightIndex]

    LamdaL = HeadTitle[LeftExtent][LeftIndex]

    LamdaIn = HeadTitle[InExtent][InIndex]

    LamdaR = HeadTitle[RightExtent][RightIndex]

    wleft = (LamdaR - LamdaIn)/(LamdaR - LamdaL)
    wright = (LamdaIn - LamdaL)/(LamdaR - LamdaL)
    Lout = wleft*LL + wright*LR
    Eout = wleft*EL + wright*ER
    F = (Eout*Lin - Ein*Lout)/(Eout - Ein)
    #
    # print("3FLD Retrieval Wave ：Out {} and In {}".format(HeadTitle[LeftExtent][LeftIndex], HeadTitle[InExtent][InIndex]))
    # print("3FLD Retrieval EL {},LL {},Ein {},Lin {},ER {} LR {}".format(EL, LL, Ein, Lin,ER,LR))
    # print("3FLD Retrieval Eout {} Lout {} and SIF {}".format(Eout,Lout,F))
    return F


def Rtr3FLD_AVE(Data, HeadTitle, White, LeftStart=755, LeftEnd=759, InStart=759, InEnd=764, RightStart=755, RightEnd=759):

    smooth_White = smooth(White)
    smooth_Data = smooth(Data)

    LeftExtent, LeftIndex = getMaxValue(HeadTitle, smooth_White, Start=LeftStart, End=LeftEnd)
    EL, LL = smooth_White[LeftExtent][LeftIndex], smooth_Data[LeftExtent][LeftIndex]

    InExtent, InIndex = getMinValue(HeadTitle, White, Start=InStart, End=InEnd)
    Ein, Lin = White[InExtent][InIndex], Data[InExtent][InIndex]

    RightExtent, RightIndex = getMinValue(HeadTitle, smooth_White, Start=RightStart, End=RightEnd)
    ER, LR = smooth_White[RightExtent][RightIndex], smooth_Data[RightExtent][RightIndex]

    LamdaL = HeadTitle[LeftExtent][LeftIndex]

    LamdaIn = HeadTitle[InExtent][InIndex]

    LamdaR = HeadTitle[RightExtent][RightIndex]

    wleft = (LamdaR - LamdaIn) / (LamdaR - LamdaL)
    wright = (LamdaIn - LamdaL) / (LamdaR - LamdaL)
    Lout = wleft * LL + wright * LR
    Eout = wleft * EL + wright * ER
    F = (Eout * Lin - Ein * Lout) / (Eout - Ein)

    # print(
    #     "3FLD Retrieval Wave ：Out {} and In {}".format(HeadTitle[LeftExtent][LeftIndex], HeadTitle[InExtent][InIndex]))
    # print("3FLD Retrieval EL {},LL {},Ein {},Lin {},ER {} LR {}".format(EL, LL, Ein, Lin, ER, LR))
    # print("3FLD Retrieval Eout {} Lout {} and SIF {}".format(Eout, Lout, F))
    return F


def NDVI(Data, HeadTitle):

    Red,Nir = getRedNir(HeadTitle)

    r = np.where(HeadTitle == float(Red))
    nir = np.where(HeadTitle == float(Nir))

    RedV = float(Data[r])
    NirV = float(Data[nir])

    return (NirV-RedV)/(NirV+RedV)

def EVI(Data, HeadTitle):

    Red,Nir = getRedNir(HeadTitle)
    Green,Blue,Redege = getGBRedege(HeadTitle)

    r = np.where(HeadTitle == float(Red))
    nir = np.where(HeadTitle == float(Nir))
    b = np.where(HeadTitle == float(Blue))

    RedV = float(Data[r])
    NirV = float(Data[nir])
    BlueV = float(Data[b])

    return (NirV - RedV)/(NirV + 6 * RedV - 7.5 * BlueV + 1)

def GRVI(Data, HeadTitle):

    Red,Nir = getRedNir(HeadTitle)
    Green,Blue,Redege = getGBRedege(HeadTitle)

    r = np.where(HeadTitle == float(Red))
    nir = np.where(HeadTitle == float(Nir))
    b = np.where(HeadTitle == float(Blue))
    g = np.where(HeadTitle == float(Green))

    RedV = float(Data[r])
    NirV = float(Data[nir])
    BlueV = float(Data[b])
    GreenV = float(Data[g])

    return (GreenV - RedV)/(GreenV + RedV)


def GNDVI(Data, HeadTitle):

    Red,Nir = getRedNir(HeadTitle)
    Green,Blue,Redege = getGBRedege(HeadTitle)

    r = np.where(HeadTitle == float(Red))
    nir = np.where(HeadTitle == float(Nir))
    b = np.where(HeadTitle == float(Blue))
    g = np.where(HeadTitle == float(Green))

    RedV = float(Data[r])
    NirV = float(Data[nir])
    BlueV = float(Data[b])
    GreenV = float(Data[g])

    return (NirV - GreenV)/(NirV + GreenV)

def MSAVI(Data, HeadTitle):

    Red,Nir = getRedNir(HeadTitle)
    Green,Blue,Redege = getGBRedege(HeadTitle)

    r = np.where(HeadTitle == float(Red))
    nir = np.where(HeadTitle == float(Nir))
    b = np.where(HeadTitle == float(Blue))
    g = np.where(HeadTitle == float(Green))

    RedV = float(Data[r])
    NirV = float(Data[nir])
    BlueV = float(Data[b])
    GreenV = float(Data[g])

    return (2*NirV +1 - np.sqrt((2*NirV +1)*(2*NirV +1) - 8*(NirV - RedV)))/2

def OSAVI(Data, HeadTitle):

    Red,Nir = getRedNir(HeadTitle)
    Green,Blue,Redege = getGBRedege(HeadTitle)

    r = np.where(HeadTitle == float(Red))
    nir = np.where(HeadTitle == float(Nir))
    b = np.where(HeadTitle == float(Blue))
    g = np.where(HeadTitle == float(Green))

    RedV = float(Data[r])
    NirV = float(Data[nir])
    BlueV = float(Data[b])
    GreenV = float(Data[g])

    return (NirV - RedV)/(NirV + RedV + 0.16)

def SAVI(Data, HeadTitle):

    Red,Nir = getRedNir(HeadTitle)
    Green,Blue,Redege = getGBRedege(HeadTitle)

    r = np.where(HeadTitle == float(Red))
    nir = np.where(HeadTitle == float(Nir))
    b = np.where(HeadTitle == float(Blue))
    g = np.where(HeadTitle == float(Green))

    RedV = float(Data[r])
    NirV = float(Data[nir])
    BlueV = float(Data[b])
    GreenV = float(Data[g])

    return 1.5*(NirV - RedV)/(NirV + RedV + 0.5)


def WDRVI(Data, HeadTitle):

    Red,Nir = getRedNir(HeadTitle)
    Green,Blue,Redege = getGBRedege(HeadTitle)

    r = np.where(HeadTitle == float(Red))
    nir = np.where(HeadTitle == float(Nir))
    b = np.where(HeadTitle == float(Blue))
    g = np.where(HeadTitle == float(Green))

    RedV = float(Data[r])
    NirV = float(Data[nir])
    BlueV = float(Data[b])
    GreenV = float(Data[g])

    return (0.1*NirV - RedV)/(0.1*NirV + RedV)

def VIS(Data, HeadTitle):
    # Data = smooth(Data,7)
    Red, Nir = getRedNir(HeadTitle)
    Green, Blue, Redege = getGBRedege(HeadTitle)

    r = np.where(HeadTitle == float(Red))
    nir = np.where(HeadTitle == float(Nir))
    b = np.where(HeadTitle == float(Blue))
    g = np.where(HeadTitle == float(Green))

    RedV = float(Data[r])
    NirV = float(Data[nir])
    BlueV = float(Data[b])
    GreenV = float(Data[g])

    NDVI = (NirV - RedV) / (NirV + RedV)
    # print("(NirV + 6 * RedV - 7.5 * BlueV + 1)", NirV, RedV, BlueV)
    EVI =  (NirV - RedV) / (NirV + 6 * RedV - 7.5 * BlueV + 1)

    GRVI =   (GreenV - RedV) / (GreenV + RedV)
    GNDVI = (NirV - GreenV)/(NirV + GreenV)
    MSAVI = (2 * NirV + 1 - np.sqrt((2 * NirV + 1) * (2 * NirV + 1) - 8 * (NirV - RedV))) / 2
    OSAVI = (NirV - RedV)/(NirV + RedV + 0.16)
    SAVI = 1.5*(NirV - RedV)/(NirV + RedV + 0.5)
    WDRVI = (0.1 * NirV - RedV) / (0.1 * NirV + RedV)
    return NDVI,EVI,GRVI,GNDVI,MSAVI,OSAVI,SAVI,WDRVI



# def RtriFLD(Title,White,Dark,Data):
#     Lmdleft, Lmdin, Lmdright = getFLDParam([str(d) for d in Title])
#     print("IFLD left inside right",Lmdleft, Lmdin, Lmdright)
#     E_hat,R_hat = getIFLDParam(Title, White, Dark, Data)
#
#     HT = np.array([float(d) for d in Title[30:]])
#     mWhite = np.array([float(d) for d in White[30:]])
#     mDark = np.array([float(d) for d in Dark[30:]])
#     mData = np.array([float(d) for d in Data[30:]])
#
#
#     mData = mData - mDark
#     mWhite = (mWhite - mDark)/0.6
#
#     wout = np.where(HT == float(Lmdleft))
#
#     win = np.where(HT == float(Lmdin))
#     print("E_hat",E_hat,R_hat)
#
#     Rapp = mData[wout]/mWhite[wout]
#     alphaR = Rapp/R_hat
#
#     alphaF = (mWhite[wout]/E_hat)*alphaR
#
#     F = (alphaR*mWhite[wout]*mData[win] -mData[wout]*mWhite[win])/(alphaR*mWhite[wout]-alphaF*mWhite[win])
#     R = (mData[win]-F)/mWhite[win]
#     print("ifld Refl",F)
#     return F




def getMaxValue(HeadTitle,Data,Start=755,End=759):

    Extent = np.where((HeadTitle > Start) & (HeadTitle < End))
    Index = np.argmax(Data[Extent])

    return Extent,Index


def getMinValue(HeadTitle, Data, Start=759, End=764):
    Extent = np.where((HeadTitle > Start) & (HeadTitle < End))
    Index = np.argmin(Data[Extent])

    return Extent, Index


# def getFLDParam(Title):
#     HeadTitle = Title
#
#     Lmdin = "761.9"
#     Lmdin_idx = 0
#     Lmdin_Flag = False
#
#     Lmdleft = "755.5"
#     Lmdleft_idx = 0
#     Lmdleft_Flag = False
#
#     Lmdright = "769.7"
#     Lmdright_idx = 0
#     Lmdright_Flag = False
#
#     for idx,Lmd in enumerate(HeadTitle):
#
#         if any([wave in str(Lmd) for wave in ["757.5","757.6","757.7","757.8","757.9"]]):
#             Lmdleft = str(Lmd)
#             Lmdleft_idx = idx
#             Lmdleft_Flag = True
#
#
#         if any([wave in str(Lmd) for wave in ["760.1","760.2","760.3","760.4"]]):
#             Lmdin = str(Lmd)
#             Lmdin_idx = idx
#             Lmdin_Flag = True
#
#         if any([wave in str(Lmd) for wave in ["769.7","769.8","769.9"]]):
#             Lmdright = str(Lmd)
#             Lmdright_idx = idx
#             Lmdright_Flag = True
#
#     assert Lmdleft_Flag, "please see the xlsx file to guarentee that the wavelengt at {} exsits".format(Lmdleft)
#     assert Lmdin_Flag, "please see the xlsx file to guarentee that the wavelengt at {} exsits".format(Lmdin)
#     assert Lmdright_Flag, "please see the xlsx file to guarentee that the wavelengt at {} exsits".format(Lmdright)
#
#     return Lmdleft,Lmdin,Lmdright


def getRedNir(Title):
    HeadTitle = Title

    Red = "761.9"
    Red_idx = 0
    Red_Flag = False

    Nir = "755.5"
    Nir_idx = 0
    Nir_Flag = False


    for idx, Lmd in enumerate(HeadTitle):

        # if any([wave in str(Lmd) for wave in ["645.0","645.1","645.2","645.3","645.4"]]):
        if any([wave in str(Lmd) for wave in ["668.0", "668.1", "668.2", "668.3", "668.4"]]):
            Red = str(Lmd)
            Red_idx = idx
            Red_Flag = True

        if any([wave in str(Lmd) for wave in ["820.0","820.1","820.2","820.3","820.4"]]):
            Nir = str(Lmd)
            Nir_idx = idx
            Nir_Flag = True



    assert Red_Flag, "please see the xlsx file to guarentee that the wavelengt at {} exsits".format(Red)
    assert Nir_Flag, "please see the xlsx file to guarentee that the wavelengt at {} exsits".format(Nir)



    return Red,Nir


def getGBRedege(Title):
    HeadTitle = Title

    Green = "560"
    Green_idx = 0
    Green_Flag = False

    Blue = "475"
    Blue_idx = 0
    Blue_Flag = False

    Redege = "717"
    Redege_idx = 0
    Redege_Flag = False

    for idx, Lmd in enumerate(HeadTitle):

        if any([wave in str(Lmd) for wave in ["560.0","560.1", "560.2", "560.3", "560.4"]]):
            Green = str(Lmd)
            Green_idx = idx
            Green_Flag = True

        if any([wave in str(Lmd) for wave in ["475.0","475.1", "475.2", "475.3", "475.4"]]):
            Blue = str(Lmd)
            Blue_idx = idx
            Blue_Flag = True
        if any([wave in str(Lmd) for wave in ["717.0","717.1", "717.2", "717.3", "717.4"]]):

            Redege = str(Lmd)
            Redege_idx = idx
            Redege_Flag = True

    assert Green_Flag, "please see the xlsx file to guarentee that the wavelengt at {} exsits".format(Green)
    assert Blue_Flag, "please see the xlsx file to guarentee that the wavelengt at {} exsits".format(Blue)
    assert Redege_Flag, "please see the xlsx file to guarentee that the wavelengt at {} exsits".format(Redege)
    return Green, Blue,Redege




# def getIFLDParam(Title,White,Dark,Data):
#
#     HeadTitle = Title
#     mWhite = np.array([float(d) for d in White[30:]])
#     mHT = np.array([float(d) for d in HeadTitle[30:]])
#     mDark = np.array([float(d) for d in Dark[30:]])
#     mData = smooth(np.array([float(d) for d in Data[30:]])- mDark,1)
#     mWhite = smooth((mWhite - mDark)/0.73,1)
#
#
#
#     ext_left = (mHT > 750) & (mHT < 759)
#     ext_left = np.where(ext_left)
#
#     ext_inside = (mHT > 760) & (mHT < 763)
#     ext_inside = np.where(ext_inside)
#     lmd_inside = mHT[ext_inside][np.argmin(mWhite[ext_inside])]
#     print("lmd_inside",lmd_inside)
#
#     ext_right = (mHT > 770) & (mHT < 780)
#     ext_right = np.where(ext_right)
#
#
#
#
#     E_l_maxs = mWhite[ext_left][argrelextrema(mWhite[ext_left],np.greater)]
#     E_l_lmds = mHT[ext_left][argrelextrema(mWhite[ext_left],np.greater)]
#
#     E_r_maxs = mWhite[ext_right][argrelextrema(mWhite[ext_right],np.greater)]
#     E_r_lmds = mHT[ext_right][argrelextrema(mWhite[ext_right],np.greater)]
#
#
#     E_lmds = np.hstack((E_l_lmds, E_r_lmds))
#     E_maxs = np.hstack((E_l_maxs, E_r_maxs))
#     z1 = np.polyfit(E_lmds, E_maxs, 2)  # 用2次多项式拟合
#     p1 = np.poly1d(z1)
#     Ein_hat = p1(lmd_inside)
#
#
#
#
#
#     ext_left2 = (mHT > 755) & (mHT < 759)
#     ext_left2 = np.where(ext_left2)
#
#     ext_right2 = (mHT > 770) & (mHT < 775)
#     ext_right2 = np.where(ext_right2)
#
#     E_l_maxs2 = mWhite[ext_left2]
#     E_l_lmds2 = mHT[ext_left2]
#
#     E_r_maxs2 = mWhite[ext_right2]
#     E_r_lmds2 = mHT[ext_right2]
#
#     E_lmds2 = np.hstack((E_l_lmds2, E_r_lmds2))
#     E_maxs2 = np.hstack((E_l_maxs2, E_r_maxs2))
#
#
#     L_l_maxs2 = mData[ext_left2]
#     L_r_maxs2 = mData[ext_right2]
#     L_maxs2 = np.hstack((L_l_maxs2,L_r_maxs2))
#
#     R_maxs2 = L_maxs2/E_maxs2
#
#     z2 = np.polyfit(E_lmds2, R_maxs2, 3)  # 用2次多项式拟合
#     p2 = np.poly1d(z2)
#     R_hat = p2(lmd_inside)
#
#
#     return Ein_hat,R_hat


def smooth(arr,kernel= 3):
    smooth_data = []
    arr_len = len(arr)
    for idx,data in enumerate(arr):
        if (idx < int(kernel/2)) or (idx >= arr_len - int(kernel/2)):

            smooth_data.append(data)
        else:

            mean = 0
            start = idx - int(kernel/2)
            for i in range(kernel):
                mean  += arr[start + i]

            mean = mean/kernel

            smooth_data.append(mean)

    return np.array(smooth_data)

# fig = plt.figure()


def readATCParam(path_txt=r"D:\Satellive\SimonHong_Code\SIFProject\ATCParams.txt"):
    waveLen = []
    Trans = []
    with open(path_txt,"r") as f:
        lines = f.readlines()
        for line in lines:



            waveLen.append(10**7/float(' '.join(line.split()).split(" ")[0]))
            Trans.append(' '.join(line.split()).split(" ")[2])



    return np.array(waveLen).astype(np.float),np.array(Trans).astype(np.float)



def plotRef(Data,HeadTitle,White,LeftStart =755,LeftEnd = 759,InStart=759,InEnd = 764,weight=1):

    extent = (HeadTitle > 450) & (HeadTitle < 830)
    ext = np.where(extent)
    wavelen,trans = readATCParam()
    waveext = (wavelen> 450) & (wavelen < 830)
    waveext = np.where(waveext)

    mwavelen = wavelen[waveext]
    mtrans = trans[waveext]


    Extent, Index =getMinValue(HeadTitle,White)


    mWhite = White[ext]
    mData = Data[ext]
    mHT = HeadTitle[ext]
    # cmap = plt.cm.get_cmap('RdYlBu')
    # norm = matplotlib.colors.Normalize(vmin=0,vmax=287)

    mmmx = argrelextrema(mWhite,np.greater)

    trs = []
    for idx,item in enumerate(mHT-1.328):
        trs.append(mtrans[np.argmin(np.abs(mwavelen-item))])
        print(item,mwavelen[np.argmin(np.abs(mwavelen-item))])

    trs = np.array(trs).astype(np.float)

    plt.plot(mHT,mData/mWhite/trs,"--",label=u"大气校正后表观反射率")

    plt.plot(mHT,  mData/mWhite,label=u"实际表观反射率")
    # plt.plot(mHT, mWhite, label=u"实际参考板光谱 DN")
    # plt.plot(mHT, mData, "--",label=u"农田光谱 DN")
    # plt.plot(mHT,  mData/mWhite/trs,"--",label=u"大气校正后表观反射率")


    # plt.plot(mHT, mWhite,label=u"参考板 DN：反射率 0.73")
    # plt.scatter(mHT[mmmx],mData[mmmx]/mWhite[mmmx])


    # cc = [norm(weight) for i in range(287)]
    # plt.scatter(mHT,mData/mWhite,c=cc, marker="_",alpha=30,vmin=0, vmax=1, s=10, cmap=cmap)
    # plt.plot(mHT,mData/mWhite,c =cmap(norm(weight)))


    # plt.plot(mHT, mWhite)
    # plt.scatter([HeadTitle[Extent][Index]],[Data[Extent][Index]/White[Extent][Index]],c="r",marker="s")
    # plt.scatter([HeadTitle[Extent][Index],HeadTitle[Extent][Index]], [0,0.2])
    # plt.scatter([HeadTitle[Extent][Index] ], [White[Extent][Index]])
    # plt.plot(mHT, mWhite*0.0001)

    # ax1 =fig.add_add_subplot(111)
    # ax1.plot(mHT,mData)
    # ax2 =  ax1.twinx()
    # ax2.plot(wavelen,)







# Data,HT,White =fileload(r"D:\Satellive\SimonHong_Code\HUBEIFIELDDATA\20180830airline1-5\15350964458442982835611128394267_20180830113249.06.xlsx",White_Reflectance = 0.73)
# print(len(Data))
# intr = 0
# endd = intr+1
# for i in range(len(Data))[intr:endd]:
#
#     plotRef(Data[i], HT, White, LeftStart=755, LeftEnd=759, InStart=759, InEnd=764, weight=i)
# # fig.colorbar(plt.gci())
# plt.legend()
# plt.show()

# F = RtrsFLD(mData[0],mHT,mWhite,LeftStart =755,LeftEnd = 759,InStart=759,InEnd = 764)
#
# Ext,Index = getMaxValue(mHT,mWhite,Start=755,End=759)
# print(Ext)
# print(Index)
# print("ee",mHT[Ext][Index])
#
# Ext,Index = getMinValue(mHT,mWhite)
# print(Ext)
# print(Index)
# print(mWhite[Ext][Index])


# plt.plot()
