from osgeo import gdal
import numpy as np
from matplotlib import  pyplot as plt
path = path  = r'C:\Users\zmhwh\Documents\Tencent Files\1475598891\FileRecv\WuHanCity201901.tif'
bands = gdal.Open(path)
W ,H = bands.RasterXSize,bands.RasterYSize

w = 0
h = 0
stepx = int(W/60)
stepy = int(H/80)
print(stepx,stepy,W,H)
i = 0
while h < H:
    while w < W:
        print(w,h)
        w += stepx
        i +=1
    print("n",i)
    i =0
    w = 0
    h += stepy

# band1 = bands.GetRasterBand(1).ReadAsArray()
# band2 = bands.GetRasterBand(2).ReadAsArray()
# band3 = bands.GetRasterBand(3).ReadAsArray()
#
# band1[band1 >10000] = 10000
# band2[band2 >10000] = 10000
# band3[band3 >10000] = 10000
#
# max1,min1 = band1.max(),band1.min()-1
# max2,min2 = band2.max()+100,band2.min()-2
# max3,min3 = band3.max()+10,band3.min()
#
#
# max1,min1 = 10000,band1.min()
# max2,min2 = 10000,band2.min()
# max3,min3 = 10000,band3.min()
#
# print(max1,min1)
# print(max2,min2)
# print(max3,min3)
# print(W,H)
#
# # band1 = np.array((band1 -min1)*255/(max1-min1)).astype(np.bytes)
# # band2 = np.array((band2 -min2)*255/(max2-min2)).astype(np.bytes)
# # band3 = np.array((band3 -min3)*255/(max3-min3)).astype(np.bytes)
#
#
#
#
# step =20
# NW,NH = int(W/step)+1,int(H/step)+1
# print(NW,NH)
# arr1 = np.zeros(shape=(NH,NW),dtype=np.int)
# arr2 = np.zeros(shape=(NH,NW),dtype=np.int)
# arr3 = np.zeros(shape=(NH,NW),dtype=np.int)
# h = 0
# w = 0
#
# while h < NH-1:
#     # print("**",h)
#     while w < NW -1:
#
#         arr1[h][w] = int((band1[step * h][step * w] - min1) *255/ (max1 - min1))
#
#         arr2[h][w] = int((band2[step * h][step * w] - min2) / (max2 - min2))
#         arr3[h][w] = int((band3[step * h][step * w] - min3) / (max3 - min3))
#         # print("hang lie",step*h,step*w)
#         w += 1
#     w = 0
#     h += 1
#
#
# a = np.zeros(shape=(NH,NW,3),dtype=np.int)
# a[:,:,0] = arr1
# a[:,:,1] = arr2
# a[:,:,2] = arr3
#
# print(a[:,:,0].max(),a[:,:,0].min())
# print("d")
# plt.imshow(a)
# plt.show()