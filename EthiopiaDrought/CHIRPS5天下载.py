import urllib.request
import os
import glob
from dateutil import rrule
from datetime import *
import time
DownloadDir = r"D:\Cornell\EthiopianDrought\CHIRPS10Day"
#/workdir/zh55/Drought/CHIRPS5
urlpath =r"https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_pentad/tifs"



start = datetime.strptime("-".join(["2003", "01", "01"]), "%Y-%m-%d").date()
stop = datetime.strptime("-".join(["2018", "12", "31"]), "%Y-%m-%d").date()

baseyear = 2003
for dt in (rrule.rrule(rrule.MONTHLY, interval=1, dtstart=start, until=stop)):
    if dt.year != baseyear:
        time.sleep(20)
        baseyear = dt.year
        print(baseyear,"sleep")


    for i in range(1,7):

        chirps_file = "/".join([urlpath,
                               "chirps-v2.0.{}.{}.{}.tif.gz"
                               .format(str(dt.year),str(dt.month).zfill(2),str(i))])
        print(chirps_file)


        chirps_out = os.path.join(DownloadDir,os.path.basename(chirps_file))


        if not os.path.exists(chirps_out):
            urllib.request.urlretrieve(chirps_file, chirps_out)
            print(os.path.basename(chirps_out))

# chirps_file = "/".join([urlpath,str(2004),
#                                "chirps-v2.0.2004.01.22.tif.gz"])
# chirps_out = os.path.join(DownloadDir,os.path.basename(chirps_file))
# urllib.request.urlretrieve(chirps_file, chirps_out)
