from osgeo import gdal,osr,ogr
from GetPointFromShp import *
from SifRetrieval2019 import *
import xlrd
from scipy import stats
import pandas as pd
path = r"D:\Satellive\cropyieldforecast\boundarypolygon.shp"



def retrieval(spec_path,cropyield,cpid,cparea):

    Dataset,HeadTitle,White = fileload(spec_path)

    Dataset2, HeadTitle2 = fileload2(spec_path) # 包含了所有的列数据



    # fig = plt.figure()
    Crop_Sif = []
    assert len(Dataset) == len(Dataset2), "mispatch occurs"


    obser_num = 0
    Sum_FLD, Sum_FLD_AVE,Sum_3FLD,Sum_3FLD_AVE = 0.0,0.0,0.0,0.0
    mSum_FLD, mSum_FLD_AVE, mSum_3FLD, mSum_3FLD_AVE = -99, -99, -99, -99

    Sum_NDVI,Sum_EVI,Sum_GRVI,Sum_GNDVI,Sum_MSAVI,Sum_OSAVI,Sum_SAVI,Sum_WDRVI = 0,0,0,0,0,0,0,0
    mSum_NDVI,mSum_EVI,mSum_GRVI,mSum_GNDVI,mSum_MSAVI,mSum_OSAVI,mSum_SAVI,mSum_WDRVI = \
        -99,-99,-99,-99,-99,-99,-99,-99   # 这些是为了求均值的呢

    for idx in range(len(Dataset)):

        Sum_FLD += RtrsFLD(Dataset[idx],HeadTitle,White,LeftStart=757) # FLD反演的荧光

        Sum_FLD_AVE += RtrsFLD_AVE(Dataset[idx],HeadTitle,White) # FLD平滑后反演的荧光

        Sum_3FLD += Rtr3FLD(Dataset[idx],HeadTitle,White) # 3FLD反演的荧光

        Sum_3FLD_AVE += Rtr3FLD_AVE(Dataset[idx], HeadTitle, White) # 3FLD平滑后反演的荧光
        NDVI2, EVI2, GRVI2, GNDVI2, MSAVI2, OSAVI2, SAVI2, WDRVI2 = VIS(Dataset[idx], HeadTitle)


        Sum_NDVI += NDVI2
        Sum_EVI      += EVI2
        Sum_GRVI     += GRVI2
        Sum_GNDVI    += GNDVI2
        Sum_MSAVI    += MSAVI2
        Sum_OSAVI    += OSAVI2
        Sum_SAVI     += SAVI2
        Sum_WDRVI    += WDRVI2

        obser_num  += 1  #观测数量+1



    if obser_num < 5:
        mSum_FLD = -99
        mSum_FLD_AVE = -99
        mSum_3FLD = -99
        mSum_3FLD_AVE = -99
        mSum_NDVI = -99
        mSum_EVI  = -99
        mSum_GRVI = -99
        mSum_GNDVI = -99
        mSum_MSAVI = -99
        mSum_OSAVI = -99
        mSum_SAVI = -99
        mSum_WDRVI = -99

    else:

        mSum_FLD = Sum_FLD/obser_num   # mSum_FLD means the average of the sum
        mSum_FLD_AVE = Sum_FLD_AVE/obser_num
        mSum_3FLD = Sum_3FLD/obser_num
        mSum_3FLD_AVE = Sum_3FLD_AVE/obser_num
        mSum_NDVI = Sum_NDVI/obser_num
        mSum_EVI = Sum_EVI/obser_num
        mSum_GRVI = Sum_GRVI /obser_num
        mSum_GNDVI = Sum_GNDVI/obser_num
        mSum_MSAVI = Sum_MSAVI/obser_num
        mSum_OSAVI = Sum_OSAVI/obser_num
        mSum_SAVI = Sum_SAVI/obser_num
        mSum_WDRVI = Sum_WDRVI/obser_num

        info_dict ={}
        info_dict["cropId"] = cpid
        info_dict["yield"] = cropyield
        info_dict['area'] = cparea
        info_dict["FLD"] = mSum_FLD
        info_dict["FLD_AVE"] = mSum_FLD_AVE
        info_dict["3FLD"] = mSum_3FLD
        info_dict["3FLD_AVE"] = mSum_3FLD_AVE
        info_dict["NDVI"] = mSum_NDVI
        info_dict["EVI"] = mSum_EVI
        info_dict["GRVI"] = mSum_GRVI
        info_dict["GNDVI"] = mSum_GNDVI
        info_dict["MSAVI"] = mSum_MSAVI
        info_dict["OSAVI"] = mSum_OSAVI
        info_dict["SAVI"] = mSum_SAVI
        info_dict["WDRVI"] = mSum_WDRVI

        Crop_Sif.append(info_dict)

    print(Crop_Sif)
    return Crop_Sif





def plotD(measurements,cropyield,cpid,cparea):


    cp_sif = retrieval(measurements, cropyield,cpid,cparea)
    Area = np.array([float(cps['area']) for cps in cp_sif])
    FLD = np.array([float(cps["FLD"]) for cps in cp_sif])
    FLD_AVE = np.array([float(cps["FLD_AVE"]) for cps in cp_sif])
    TrFLD = np.array([float(cps["3FLD"]) for cps in cp_sif])
    TrFLD_AVE = np.array([float(cps["3FLD_AVE"]) for cps in cp_sif])
    NDVI = np.array([float(cps["NDVI"]) for cps in cp_sif])
    EVI = np.array([float(cps["EVI"]) for cps in cp_sif])
    GRVI = np.array([float(cps["GRVI"]) for cps in cp_sif])
    GNDVI = np.array([float(cps["GNDVI"]) for cps in cp_sif])
    MSAVI = np.array([float(cps["MSAVI"]) for cps in cp_sif])
    OSAVI = np.array([float(cps["OSAVI"]) for cps in cp_sif])
    SAVI = np.array([float(cps["SAVI"]) for cps in cp_sif])
    WDRVI = np.array([float(cps["WDRVI"]) for cps in cp_sif])

    YLD = np.array([float(cps["yield"]) for cps in cp_sif])

    masks = np.where(FLD != -99)
    mFLD, mFLD_AVE, m3FLD, m3FLD_AVE,mNDVI,mEVI,mGRVI ,mGNDVI ,mMSAVI ,mOSAVI ,mSAVI ,mWDRVI = FLD[masks], FLD_AVE[masks], TrFLD[masks], TrFLD_AVE[masks],NDVI[masks], \
                                             EVI[masks],GRVI[masks] ,GNDVI[masks] ,MSAVI[masks] ,OSAVI[masks] ,SAVI[masks] ,WDRVI[masks]
    Area, YLD = Area[masks], YLD[masks]

    # mFLD = mFLD * Area * 0.001
    # mFLD_AVE = mFLD_AVE * Area * 0.001
    # m3FLD = m3FLD * Area * 0.001
    # m3FLD_AVE = m3FLD_AVE * Area * 0.001
    # mNDVI = mNDVI*Area*0.001
    # mEVI  = mEVI*Area*0.001
    # mGRVI = mGRVI*Area*0.001
    # mGNDVI = mGNDVI*Area*0.001
    # mMSAVI = mMSAVI*Area*0.001
    # mOSAVI = mOSAVI*Area*0.001
    # mSAVI = mSAVI*Area*0.001
    # mWDRVI = mWDRVI*Area*0.001


    # print("NDVI",mNDVI)
    #
    # slope1, intercept1, r_value_fld, p_value_fld, std_err = stats.linregress(mFLD*Area,YLD)
    # print("FLD cor", r_value_fld ** 2,slope1)
    #
    # slope2, intercept2, r_value_fld2, p_value_fld2, std_err = stats.linregress(mFLD_AVE*Area,YLD)
    # print("FLD_AVE cor", r_value_fld2 ** 2,slope2)
    #
    # slope3, intercept3, r_value_3fld, p_value_3fld, std_err = stats.linregress(m3FLD*Area,YLD, )
    # print("3FLD cor", r_value_3fld ** 2,slope3)
    #
    # slope4, intercept4, r_value_3fld2, p_value_3fld2, std_err = stats.linregress(m3FLD_AVE*Area,YLD)
    # print("3FLD_AVE cor", r_value_3fld2 ** 2,slope4)
    #
    # slope5, intercept5, r_value_ndvi, p_value_ndvi, std_err = stats.linregress(mNDVI*Area,YLD)
    # print("NDVI cor", r_value_ndvi ** 2,slope5)
    #
    # slope6, intercept6, r_value_evi, p_value_evi, std_err = stats.linregress(mEVI*Area,YLD)
    # print("EVI cor", r_value_evi ** 2, slope6)
    #
    # slope7, intercept7, r_value_grvi, p_value_grvi, std_err = stats.linregress(mGRVI*Area,YLD)
    # print("GRVI cor", r_value_grvi ** 2, slope7)
    #
    # slope8, intercept8, r_value_gndvi, p_value_gndvi, std_err = stats.linregress(mGNDVI*Area,YLD)
    # print("GNDVI cor", r_value_gndvi ** 2, slope8)
    #
    # slope9, intercept9, r_value_msavi, p_value_msavi, std_err = stats.linregress(mMSAVI*Area,YLD)
    # print("MSAVI cor", r_value_msavi ** 2, slope9)
    #
    # slope10, intercept10, r_value_osavi, p_value_osavi, std_err = stats.linregress(mOSAVI*Area,YLD)
    # print("OSAVI cor", r_value_osavi ** 2, slope10)
    #
    # slope11, intercept11, r_value_savi, p_value_savi, std_err = stats.linregress(mSAVI*Area,YLD)
    # print("SAVI cor", r_value_savi ** 2, slope11)
    #
    # slope12, intercept12, r_value_wdrvi, p_value_wdrvi, std_err = stats.linregress(mWDRVI*Area,YLD)
    # print("WDRVI cor", r_value_wdrvi ** 2, slope12)
    #
    # xmin,xmax = np.min(YLD) ,np.max(YLD)
    print("dat",os.path.basename(measurements).split("_")[1][:-5])
    month = os.path.basename(measurements).split("_")[1][:-5][5:6]
    day = os.path.basename(measurements).split("_")[1][:-5][6:8]
    time = os.path.basename(measurements).split("_")[1][:-5][8:14]
    areaname = os.path.dirname(measurements)
    # if areaname.endswith("1-5"):
    #     Field =1
    # if areaname.endswith("6-10"):
    #     Field = 2
    dataframe = pd.DataFrame(
        {"Field":cpid,"month":month,"day":day,"time":time,
         'Yield': YLD, "CropArea": Area,
         "FLD": mFLD,
         "FLD_AVE": mFLD_AVE,
         "3FLD": m3FLD,
         "3FLD_AVE": m3FLD_AVE,
         "NDVI":mNDVI,
         "mEVI":mEVI,
         "mGRVI":mGRVI ,
         "mGNDVI":mGNDVI,
         "mMSAVI":mMSAVI,
         "mOSAVI":mOSAVI,
         "mSAVI":mSAVI,
         "mWDRVI":mWDRVI
         })

    # 将DataFrame存储为csv,index表示是否显示行名，default=True
    dataframe.to_csv(os.path.join(r"D:\Satellive\Estimation", "maifei2019.csv"), index=False, sep=',',mode="a")



import glob



Datapath = r"D:\Satellive\20190517-20190827\Data"
subDatapath = os.listdir(Datapath)
secondSubDatapath =["A1","A2","A3","A4","A6"]
FNS =[]
for subpath in subDatapath:
    for sec in secondSubDatapath:
        fns = glob.glob(os.path.join(Datapath,subpath,sec,"*.xlsx"))
        for f in fns:
            FNS.append(f)
print(FNS)


for md in FNS:
    if "A1" in md:
        cropyield, cpid, cparea = 2245,1,2075
    if "A2" in md:
        cropyield, cpid, cparea = 2530,2,2242
    if "A3" in md:
        cropyield, cpid, cparea = 2720,3,2458
    if "A4" in md:
        cropyield, cpid, cparea = 3205,4,2620
    if "A6" in md:
        cropyield, cpid, cparea = 3830,6,3460

    plotD(md, cropyield, cpid, cparea)


