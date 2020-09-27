import os
import glob
path = r"D:\MOS\T49RGQ"
names = glob.glob(os.path.join(path,"*"))
for name in names:
    nm = os.path.basename(name).split("_")

    newname = "_".join([nm[0],nm[1],nm[2],nm[4],nm[5]])
    newname = os.path.join(path,newname)
    os.rename(name,newname)


