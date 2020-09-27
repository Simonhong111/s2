import numpy as np
from GridGen import GridDefination
def Pairs(s2point_path,sifpoint_path,hb_shp_path,save_path):
    """
    :param s2point_path:
    :param sifpoint_path: 2018-5.npz 某一月的sif 数据
    :param hb_shp_path:湖北省
    :return:
    """
    # hb_shp_grid = GridDefination()
    # hb_shp_grid.getGrid(hb_shp_path)
    # cols,rows = hb_shp_grid.GridToalNum()
    cols,rows = 7200,3600
    print(cols,rows)
    s2point = np.load(s2point_path)
    sifpoint = np.load(sifpoint_path)
    # kwargs = {"gridId": GridID, "gridEtRing": GridEtRing, "Date": Date, "mean": BandArr}
    s2grid = s2point["gridId"]
    s2date = s2point["Date"]
    s2Refl = s2point["mean"]
    s2EtRing=s2point["gridEtRing"]
    sifdate =sifpoint["date"]
    sifid = sifpoint["sifgridId"]
    sif757 = sifpoint["sif757"]
    sif771 = sifpoint["sif771"]



    pairDate =[]
    pairID =[]
    pairS2 =[]
    pairSIF757 =[]
    pairSIF771 =[]
    pairEtRing =[]


    for idx,dt in enumerate(sifdate):

        for sidx,sdt in enumerate(s2date):
            print("***", dt, sdt)
            print(sifid[idx], s2grid[sidx])
            if dt == sdt and sifid[idx] == s2grid[sidx]:

                pairDate.append(dt)
                pairID.append(sifid[idx])
                pairS2.append(s2Refl[sidx])
                pairSIF757.append(sif757[idx])
                pairSIF771.append(sif771[idx])
                pairEtRing.append(s2EtRing[sdt])

    kwargs = {"date":pairDate,"gridId": pairID, "S2Refl": pairS2, "sif757": pairSIF757,"sif771":pairSIF771,"EtRing":pairEtRing}
    # print(ValTile)
    np.savez(save_path, **kwargs)
