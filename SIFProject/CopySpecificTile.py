import os
import glob
import shutil


directory = r"G:\2018S2docker"

sen2files = glob.glob(os.path.join(directory,"T49RCN","*201808*"+".SAFE"))



print(sen2files.__len__())


for file in sen2files:
    if os.path.exists(os.path.join("D:\S2Composite",os.path.basename(file))):
        print("The {} already exists, skip".format(os.path.basename(file)))
        continue
    shutil.copytree(file,os.path.join(r"D:\S2Composite",os.path.basename(file)))
    print("done")
print("finished")

