import csv,os
import numpy as np
import glob
from dateutil import rrule
from datetime import *
import time
from pylab import mpl

mpl.rcParams['font.sans-serif'] = ['FangSong']
mpl.rcParams['axes.unicode_minus'] = False
from matplotlib import pyplot as plt

global interval
interval = 10
def fileload(filename,key="ndvi"):
    """
    :param filename: xlsx 数据记录文件路径
    :return: 返回float类型的第31列之后的数据，这些数据都是数字
    """
    Timeset = []
    Dataset = []
    with open(filename, encoding='UTF-8') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            Timeset.append(row["date"])
            Dataset.append(row[key])

    return Timeset,Dataset


def getWater(path):

    mTime, mWater = fileload(path,"water")
    assert len(mTime) == len(mWater), "注意，文件{}列有空缺值".format(path)

    Tset = np.array(mTime)
    Wset = np.array(mWater).astype(np.float)

    mask = np.where(Wset != -1)

    Tset = Tset[mask]
    Wset = Wset[mask]

    return Tset, Wset


def getBatchData(dir,key):

    Files = glob.glob(os.path.join(dir,"*csv"))

    NDVISet = []
    TimeSet = []
    for mfile in Files:
        print("**",mfile)
        mTime,mNDVI = fileload(mfile,key)
        assert len(mTime) == len(mNDVI),"注意，文件{}列有空缺值".format(mfile)
        TimeSet.extend(mTime)
        NDVISet.extend(mNDVI)
    Tset = np.array(TimeSet)
    Nset = np.array(NDVISet).astype(np.float)
    if key == "ndvi":
        print("ndvi")
        mask = (Nset <1) & (Nset > -1) & (Nset != 0)
    if key == "preciption":
        print("precipation")
        mask = (Nset > 0) & (Nset < 999)
    if key == "avg":
        print("avg")
        mask = (Nset < 999)
    if key == "realavg":
        print("realavg")
        mask = (Nset > 0)
    if key == "streamflow":
        print("streamflow")
        mask = (Nset > 0)
    maskNdvi = np.where(mask)


    Tset = Tset[maskNdvi]
    Nset = Nset[maskNdvi]


    return Tset,Nset

def converTime(Tset):
    start = datetime.strptime("2006-1-1", "%Y-%m-%d").date()
    Timeset = []
    for dt in Tset:
        # print(dt)
        mdt = datetime.strptime(dt, "%Y-%m-%d").date() - start
        mdt = int(mdt.days)
        Timeset.append(mdt)
    return np.array(Timeset)
def converTimeWater(Tset):
    start = datetime.strptime("2006-1-1", "%Y-%m-%d").date()
    Timeset = []
    for dt in Tset:
        mdt = datetime.strptime(dt, "%Y/%m/%d").date() - start
        mdt = int(mdt.days)
        Timeset.append(mdt)
    return np.array(Timeset)


def getDate():

    start = datetime.strptime("2006-1-1", "%Y-%m-%d").date()
    stop = datetime.strptime("2015-12-31", "%Y-%m-%d").date()
    TM      = []
    TMLabel = []
    for dt in rrule.rrule(rrule.YEARLY, interval= 1, dtstart=start, until=stop):
        mdt = (dt.date() - start)
        mdt = int(mdt.days)
        TM.append(mdt)
        TMLabel.append(str(dt.date()))

    return TM,TMLabel



def drawPlot(T1,D1,Name1,T2,D2,Name2,extent):
    TM, TML = getDate()
    fig = plt.figure(figsize=(20,8))

    ax1 = fig.add_subplot(111)
    ax1.plot(T1, D1, label=u"{}".format(Name1))
    ax1.set_ylabel(u"{}天{}值".format(interval,Name1))
    ax1.set_title(u"{}天{}与{}时序图".format(interval,Name1,Name2))
    ax1.legend(loc=1)
    ax1.set_xticks(TM)
    ax1.set_xticklabels(TML)
    ax2 = ax1.twinx()  # this is the important function
    ax2.plot(T2, D2, c="r", label=u"{}".format(Name2))
    # ax2.set_xlim([0, np.e])
    ax2.set_ylabel(u'{}天{}'.format(interval,Name2))
    ax2.set_xlabel(u'时间')
    ax2.legend(loc=2)
    for xtick in ax1.get_xticklabels():
        xtick.set_rotation(50)
    fig.savefig(os.path.join(r"C:\Users\Administrator\Desktop\finalexcel_hzm\Tre",str(interval)+Name1+"_"+Name2+extent+".png"))
    plt.close()
    # plt.show()

def TimeAnalysis(dir,waterlevel,extent):

    WTset,Wset = getWater(waterlevel)
    Tset,Nset = getBatchData(dir,"ndvi")
    LTset,LSTset = getBatchData(dir,"avg")
    pTset, pset = getBatchData(dir, "preciption")
    # STset, Sset = getBatchData(dir, "streamflow")

    start = datetime.strptime("2006-1-1", "%Y-%m-%d").date()
    stop = datetime.strptime("2015-12-31", "%Y-%m-%d").date()

    NT = converTime(Tset)
    WT = converTimeWater(WTset)
    LT = converTime(LTset)
    PT = converTime(pTset)
    # ST = converTimeWater(STset)

    Tm =   []
    Ndvi = []
    LTm = []
    LST = []
    Water = []
    WTm = []
    PTm = []
    Prec =[]
    STm = []
    Stream = []

    for dt in rrule.rrule(rrule.DAILY, interval=interval, dtstart=start, until=stop):

        mdt = (dt.date() - start)
        mdt = int(mdt.days)

        mask = (NT >= mdt) & (NT < mdt+interval)
        mask = np.where(mask)

        if len(mask[0]) > 0:
            Tm.append(mdt)
            Ndvi.append(np.mean(Nset[mask]))

        wmask = (WT >= mdt) & (WT < mdt+interval)
        wmask = np.where(wmask)

        if len(wmask[0]) > 0:
            WTm.append(mdt)
            Water.append(np.mean(Wset[wmask]))

        Lmask = (LT >= mdt) & (LT < mdt + interval)
        Lmask = np.where(Lmask)

        if len(Lmask[0]) > 0:
            LTm.append(mdt)
            LST.append(np.mean(LSTset[Lmask]))

        Pmask = (PT >= mdt) & (PT < mdt + interval)
        Pmask = np.where(Pmask)

        if len(Pmask[0]) > 0:
           PTm.append(mdt)
           Prec.append(np.mean(pset[Pmask]))

        # Smask = (ST >= mdt) & (ST < mdt + interval)
        # Smask = np.where(Smask)

        # if len(Smask[0]) > 0:
        #     STm.append(mdt)
        #     Stream.append(np.mean(Sset[Smask]))

        # drawPlot(STm, Stream, "径流", Tm,Ndvi,"ndvi")
        # drawPlot(STm, Stream, "径流", LTm, LST, "温度")
        # drawPlot(STm, Stream, "径流", PTm, Prec, "降雨")

        drawPlot(WTm, Water, "水位", PTm, Prec, "降雨",extent)
        drawPlot(WTm, Water, "水位", LTm, LST, "温度",extent)
        drawPlot(PTm, Prec, "降雨", Tm, Ndvi, "ndvi",extent)
        drawPlot( LTm, LST, "温度",Tm, Ndvi, "ndvi",extent)















TimeAnalysis(r"C:\Users\Administrator\Desktop\finalexcel_hzm\threegeorge",\
             r"C:\Users\Administrator\Desktop\三峡日平均水位.csv","")
TimeAnalysis(r"C:\Users\Administrator\Desktop\finalexcel_hzm\4个buffer\buffer25",\
             r"C:\Users\Administrator\Desktop\三峡日平均水位.csv","25")
TimeAnalysis(r"C:\Users\Administrator\Desktop\finalexcel_hzm\4个buffer\buffer50",\
             r"C:\Users\Administrator\Desktop\三峡日平均水位.csv","50")
TimeAnalysis(r"C:\Users\Administrator\Desktop\finalexcel_hzm\4个buffer\buffer75",\
             r"C:\Users\Administrator\Desktop\三峡日平均水位.csv","75")
TimeAnalysis(r"C:\Users\Administrator\Desktop\finalexcel_hzm\4个buffer\buffer100",\
             r"C:\Users\Administrator\Desktop\三峡日平均水位.csv","100")