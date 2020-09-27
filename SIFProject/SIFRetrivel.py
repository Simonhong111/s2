import numpy as np

import csv



def readSpectrum(sp_path,title1,title2,threshold=720):
    files = open(sp_path)

    spectrum = csv.DictReader(files)

    bandcenter = []
    reflectance =[]

    for read in spectrum:

        bandcenter.append(float(read[title1]))
        reflectance.append(float(read[title2]))

    bandcenter = np.array(bandcenter)
    reflectance = np.array(reflectance)

    ind = np.where(bandcenter > threshold)
    bandcenter = bandcenter[ind]
    reflectance = reflectance[ind]
    return bandcenter,reflectance



bandcenter,reflectance = readSpectrum(r"D:\Sen2Projecton\spectrum.csv","Title","White2",650)




import matplotlib.pyplot as plt

fig = plt.figure()
plt.plot(bandcenter,reflectance)
# for id in ["1","2","3","4"]:
for id in ["White1", "White2", "White3", "4"]:
    bandcenter2, reflectance2 = readSpectrum(r"D:\Sen2Projecton\spectrum.csv", "Title", id, 650)


    plt.plot(bandcenter2,reflectance2)
plt.show()


