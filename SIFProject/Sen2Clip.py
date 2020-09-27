from osgeo import gdal
from geo2mapxy import lonlat2geo
def GetRoi(infile,outfile,scope):
        """
        :param outfile:需要裁剪的文件
        :param scope:（ulx,uly,lrx,lry）,左上右下，地理坐标，需要转成投影坐标能用
        :return:
        """

        pszSrcFile = gdal.Open(infile,gdal.GA_ReadOnly)

        if pszSrcFile is None:
            raise Exception("The source file is empty, please check the self.RasterData ")
        else:
            uper = lonlat2geo(pszSrcFile, scope[0], scope[1])
            bttom = lonlat2geo(pszSrcFile, scope[2], scope[3])
            pszDstFile = gdal.Translate(outfile, pszSrcFile, projWin=(uper[0],uper[1],bttom[0],bttom[1]))
            pszDstFile = None





