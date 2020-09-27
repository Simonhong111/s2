import os
import glob



def validateSentinel2(SAFEFile):
    L1_sub = glob.glob(os.path.join(SAFEFile,"*"))
    L1_sub_names = [os.path.basename(sub) for sub in L1_sub]
    L1_sub_all = ["AUX_DATA","DATASTRIP","GRANULE","HTML","rep_info","INSPIRE.xml","manifest.safe","MTD_MSIL2A.xml"]
    for valid in L1_sub_all:
        if valid not in L1_sub_names:
            print("File {} is not complete lossing the subfile {}".format(SAFEFile,valid))
            return SAFEFile
    L2_DATASTRIP = glob.glob(os.path.join(SAFEFile,"DATASTRIP","*"))
    if len(L2_DATASTRIP) !=1:
        print("File {} is not complete lossing the subfile in {}".format(SAFEFile, "DATASTRIP"))
        return SAFEFile
    else:
        L2_DATASTRIP_QI_DATA = glob.glob(os.path.join(SAFEFile, "DATASTRIP", "*","QI_DATA"))
        if len(L2_DATASTRIP_QI_DATA) !=1:
            print("File {} is not complete lossing the subfile in {}".format(SAFEFile, os.path.basename(L2_DATASTRIP)))
            return SAFEFile
        else:
            QI = glob.glob(os.path.join(SAFEFile, "DATASTRIP", "*","QI_DATA","*"))
            QI_names = [os.path.basename(name) for name in QI]
            QI_name_all =["FORMAT_CORRECTNESS.xml","GENERAL_QUALITY.xml","GEOMETRIC_QUALITY.xml",
                          "RADIOMETRIC_QUALITY.xml","SENSOR_QUALITY.xml"]
            for name in QI_name_all:
                if name not in QI_names:
                    print("File {} is not complete lossing the {} in DATASTRIP/*/QI_DATA".format(SAFEFile,name))
                    return SAFEFile
        L2_DATASTRIP_MTD_DS = glob.glob(os.path.join(SAFEFile, "DATASTRIP", "*", "MTD_DS.xml"))
        if len(L2_DATASTRIP_MTD_DS) != 1:
            print("File {} is not complete lossing MTD_DS.xml in DATASTRIP".format(SAFEFile))
            return SAFEFile
    L2_HTML = glob.glob(os.path.join(SAFEFile,"HTML","*"))
    HTML_name = [os.path.basename(name) for name in L2_HTML]
    HTML_name_all = ["banner_1.png","banner_2.png","banner_3.png","star_bg.jpg","UserProduct_index.html","UserProduct_index.xsl"]
    for name in HTML_name_all:
        if name not in HTML_name:
            print("File {} is not complete lossing the subfile {} in HTML".format(SAFEFile,name))
            return SAFEFile
    L2_rep_info = glob.glob(os.path.join(SAFEFile,"rep_info","*"))
    rep_info_name = [os.path.basename(name) for name in L2_rep_info]
    rep_info_name_all = ["S2_PDI_Level-2A_Datastrip_Metadata.xsd","S2_PDI_Level-2A_Tile_Metadata.xsd","S2_User_Product_Level-2A_Metadata.xsd"]
    for name in rep_info_name_all:
        if name not in rep_info_name:
            print("File {} is not complete lossing the subfile {} in repo_info".format(SAFEFile, name))
            return SAFEFile
    L2_GRANULE = glob.glob(os.path.join(SAFEFile,"GRANULE","*"))
    if len(L2_GRANULE) !=1:
        print("File {} is not complete lossing the subfile in GRANULE".format(SAFEFile))
        return SAFEFile
    else:
        if len(glob.glob(os.path.join(SAFEFile,"GRANULE","*","AUX_DATA","AUX_ECMWFT"))) !=1:
            print("File {} is not complete lossing the subfile in GRANULE/*/AUX_DATA".format(SAFEFile))
            return SAFEFile
        if len(glob.glob(os.path.join(SAFEFile,"GRANULE","*","MTD_TL.xml"))) !=1:
            print("File {} is not complete lossing the MTD_TL.xml in GRANULE".format(SAFEFile))
            return SAFEFile
        if len(glob.glob(os.path.join(SAFEFile,"GRANULE","*","QI_DATA"))) !=1:
            print("File {} is not complete lossing the QI_DATA in GRANULE".format(SAFEFile))
            return SAFEFile
        if len(glob.glob(os.path.join(SAFEFile,"GRANULE","*","IMG_DATA"))) !=1:
            print("File {} is not complete lossing the IMG_DATA in GRANULE".format(SAFEFile))
            return SAFEFile
        else:
            L2_GRANULE_10m = glob.glob(os.path.join(SAFEFile,"GRANULE","*","IMG_DATA","R10m"))
            L2_GRANULE_20m = glob.glob(os.path.join(SAFEFile, "GRANULE", "*", "IMG_DATA", "R20m"))
            L2_GRANULE_60m = glob.glob(os.path.join(SAFEFile, "GRANULE", "*", "IMG_DATA", "R60m"))

            if len(L2_GRANULE_10m) !=1:
                print("File {} is not complete lossing the R10m in GRANULE/*/IMG_DATA".format(SAFEFile))
                return SAFEFile
            else:
                L2_GRANULE_10m_bands = glob.glob(os.path.join(SAFEFile, "GRANULE", "*", "IMG_DATA", "R10m","*jp2"))
                bands_name = [os.path.basename(name).split("_")[-2]+"_"+os.path.basename(name).split("_")[-1] for name in L2_GRANULE_10m_bands]
                bands_name_all = ["AOT_10m.jp2","B02_10m.jp2","B03_10m.jp2","B04_10m.jp2","B08_10m.jp2","TCI_10m.jp2",
                                  "WVP_10m.jp2"]
                for name in bands_name_all:
                    if name not in bands_name:
                        print("File {} is not complete lossing the {} in GRANULE/*/IMG_DATA/R10m".format(SAFEFile,name))
                        return SAFEFile
                # 20m
            if len(L2_GRANULE_20m) != 1:
                print("File {} is not complete lossing the R20m in GRANULE/*/IMG_DATA".format(SAFEFile))
                return SAFEFile
            else:
                L2_GRANULE_20m_bands = glob.glob(os.path.join(SAFEFile, "GRANULE", "*", "IMG_DATA", "R20m", "*jp2"))
                bands_name = [os.path.basename(name).split("_")[-2]+"_"+os.path.basename(name).split("_")[-1] for name in L2_GRANULE_20m_bands]
                bands_name_all = ["AOT_20m.jp2", "B02_20m.jp2", "B03_20m.jp2", "B04_20m.jp2", "B05_20m.jp2",
                                  "B06_20m.jp2","B07_20m.jp2","B8A_20m.jp2",
                                  "B11_20m.jp2","B12_20m.jp2","SCL_20m.jp2","TCI_20m.jp2",
                                  "WVP_20m.jp2"]
                for name in bands_name_all:
                    if name not in bands_name:
                        print(
                            "File {} is not complete lossing the {} in GRANULE/*/IMG_DATA/R20m".format(SAFEFile, name))
                        return SAFEFile

                        # 20m
            if len(L2_GRANULE_60m) != 1:
                print("File {} is not complete lossing the R60m in GRANULE/*/IMG_DATA".format(SAFEFile))
                return SAFEFile
            else:
                L2_GRANULE_60m_bands = glob.glob(
                    os.path.join(SAFEFile, "GRANULE", "*", "IMG_DATA", "R60m", "*jp2"))
                bands_name = [
                    os.path.basename(name).split("_")[-2] + "_" + os.path.basename(name).split("_")[-1] for
                    name in L2_GRANULE_60m_bands]
                bands_name_all = ["AOT_60m.jp2", "B01_60m.jp2","B02_60m.jp2", "B03_60m.jp2", "B04_60m.jp2", "B05_60m.jp2",
                                  "B06_60m.jp2", "B07_60m.jp2", "B8A_60m.jp2","B09_60m.jp2",
                                  "B11_60m.jp2", "B12_60m.jp2", "SCL_60m.jp2", "TCI_60m.jp2",
                                  "WVP_60m.jp2"]
                for name in bands_name_all:
                    if name not in bands_name:
                        print(
                            "File {} is not complete lossing the {} in GRANULE/*/IMG_DATA/R60m".format(
                                SAFEFile, name))
                        return SAFEFile
            Granule_QI_DATA = glob.glob(os.path.join(SAFEFile,"GRANULE","*","QI_DATA","*"))
            QI_name = [os.path.basename(name) for name in Granule_QI_DATA]
            QI_name.append("PVI.jp2")
            QI_name_all =["FORMAT_CORRECTNESS.xml","GENERAL_QUALITY.xml","GEOMETRIC_QUALITY.xml","RADIOMETRIC_QUALITY.xml",
                          "SENSOR_QUALITY.xml","PVI.jp2","MSK_CLDPRB_20m.jp2","MSK_CLDPRB_60m.jp2",
                          "MSK_CLOUDS_B00.gml",
                          "MSK_DEFECT_B01.gml","MSK_DEFECT_B02.gml","MSK_DEFECT_B03.gml","MSK_DEFECT_B04.gml",
                          "MSK_DEFECT_B05.gml","MSK_DEFECT_B06.gml","MSK_DEFECT_B07.gml","MSK_DEFECT_B08.gml",
                          "MSK_DEFECT_B8A.gml","MSK_DEFECT_B09.gml","MSK_DEFECT_B10.gml","MSK_DEFECT_B11.gml",
                          "MSK_DEFECT_B12.gml",
                          "MSK_DETFOO_B01.gml", "MSK_DETFOO_B02.gml", "MSK_DETFOO_B03.gml", "MSK_DETFOO_B04.gml",
                          "MSK_DETFOO_B05.gml", "MSK_DETFOO_B06.gml", "MSK_DETFOO_B07.gml", "MSK_DETFOO_B08.gml",
                          "MSK_DETFOO_B8A.gml", "MSK_DETFOO_B09.gml", "MSK_DETFOO_B10.gml", "MSK_DETFOO_B11.gml",
                          "MSK_DETFOO_B12.gml",
                          "MSK_NODATA_B01.gml", "MSK_NODATA_B02.gml", "MSK_NODATA_B03.gml", "MSK_NODATA_B04.gml",
                          "MSK_NODATA_B05.gml", "MSK_NODATA_B06.gml", "MSK_NODATA_B07.gml", "MSK_NODATA_B08.gml",
                          "MSK_NODATA_B8A.gml", "MSK_NODATA_B09.gml", "MSK_NODATA_B10.gml", "MSK_NODATA_B11.gml",
                          "MSK_NODATA_B12.gml",
                          "MSK_SATURA_B01.gml", "MSK_SATURA_B02.gml", "MSK_SATURA_B03.gml", "MSK_SATURA_B04.gml",
                          "MSK_SATURA_B05.gml", "MSK_SATURA_B06.gml", "MSK_SATURA_B07.gml", "MSK_SATURA_B08.gml",
                          "MSK_SATURA_B8A.gml", "MSK_SATURA_B09.gml", "MSK_SATURA_B10.gml", "MSK_SATURA_B11.gml",
                          "MSK_SATURA_B12.gml","MSK_SNWPRB_20m.jp2","MSK_SNWPRB_60m.jp2",
                          "MSK_TECQUA_B01.gml", "MSK_TECQUA_B02.gml", "MSK_TECQUA_B03.gml", "MSK_TECQUA_B04.gml",
                          "MSK_TECQUA_B05.gml", "MSK_TECQUA_B06.gml", "MSK_TECQUA_B07.gml", "MSK_TECQUA_B08.gml",
                          "MSK_TECQUA_B8A.gml", "MSK_TECQUA_B09.gml", "MSK_TECQUA_B10.gml", "MSK_TECQUA_B11.gml",
                          "MSK_TECQUA_B12.gml"
                          ]
            for name in QI_name_all:
                if name not in QI_name:
                    print(
                        "File {} is not complete lossing the {} in GRANULE/*/QI_DATA".format(
                            SAFEFile, name))
                    return SAFEFile
    return "Complete"


DataDirectory =r'D:\aria2\hongboshi'

TileDirs = glob.glob(os.path.join(DataDirectory,"T*"))
print(TileDirs)
InvalidTiles = []
for TileD in TileDirs:
    Tiles = glob.glob(os.path.join(DataDirectory, TileD,"*"))
    for name in Tiles:
        # print(name)
        Tile = validateSentinel2(name)
        # print(Tile)
        if Tile != "Complete":
            InvalidTiles.append(Tile)



import shutil
path =r'D:\aria2\hongboshi\S2A_MSIL2A_20190105T030111_N0211_R032_T49RGP_20190105T064753'

outpath = os.path.join("D:\\",os.path.basename(path))
shutil.move(path,outpath)
# print(outpath)