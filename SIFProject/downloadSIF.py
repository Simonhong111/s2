# encoding:utf-8
# 下载普特英语文音频件的脚本
import urllib.request as urllib2
import os
from ftplib import FTP


def searchUrl():
    source = r"ftp://fluo.gps.caltech.edu/data/OCO2/sif_lite_B8100/"
    ftp = FTP('fluo.gps.caltech.edu')
    ftp.login()
    ftp.cwd("data/OCO2/sif_lite_B8100")
    yList = ["2015","2016","2017","2018",]
    mList = ["01","02","03","04","05","06","07","08","09","10","11","12"]
    ymList = []
    urlList = []
    for y in yList:
        if y =="2015":
            for m in mList[8:]:
                ymList.append(y +"/"+m)
        else:
            if y == "2017":
                for m in ["01","02","03","04","05","06","07","09","10","11","12"]:
                    ymList.append(y +"/"+m)
            else:
                for m in mList:
                    ymList.append(y +"/"+m)

    for ym in ymList:
        ftp.cwd(ym)
        temp = None
        temp = ftp.nlst()
        for u in temp:
            urlList.append(source+ym+"/"+u)
        print(temp)
        ftp.cwd("../../")
        print(urlList)



    return urlList


def download(urldict,output):
    for u in urlList:
        if not os.path.exists(os.path.join(output,u.split("/")[6])):
            print("**year",os.path.join(output,u.split("/")[6]))
            os.makedirs(os.path.join(output,u.split("/")[6]))
        if not os.path.exists(os.path.join(output,u.split("/")[6],u.split("/")[7])):
            os.makedirs(os.path.join(output, u.split("/")[6],u.split("/")[7]))
        f=urllib2.urlopen(u)
        data=f.read()
        with open(os.path.join(output,u.split("/")[6],u.split("/")[7],u.split("/")[-1]),"wb")\
                as file:
            file.write(data)


        savepath  = os.path.join(os.path.join(output,u.split("/")[6],u.split("/")[7],u.split("/")[-1]))
        print(savepath)


urlList = searchUrl()
print(urlList)
download(urlList,"D:\Satellive\SIFdata")

# url="ftp://fluo.gps.caltech.edu/data/OCO2/sif_lite_B8100/2015/08/oco2_LtSIF_150801_B8100r_171011082205s.nc4"
# f=urllib2.urlopen(url)
# data=f.read()
#
# with open(r"D:\Satellive\SIFdata\oco2_LtSIF_150801_B8100r_171011082205s.nc4","wb") as file:
#     file.write(data)