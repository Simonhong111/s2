import os
import numpy as np
import glob
from tqdm import tqdm
# from geo2mapxy import  *
import numba as nb
import matplotlib.pyplot as plt
#impor
# def Delta(x0,y0,x1,y1,n,m,m1,n1,m2,n2,m3,n3):
#     Dx1 = x1 - x0
#     Dy1 = y1 - y0
#     x2 = x0 - n * Dx1 + m*Dy1 -m1*(Dx1*Dx1 - Dy1*Dy1) -2*n1*Dx1*Dy1 -\
#         m2*Dy1*(3*Dx1*Dx1 - Dy1*Dy1) + n2*Dx1*(Dx1*Dx1 - 3*Dy1*Dy1)+\
#         m3*(Dx1*Dx1*Dx1*Dx1 - 6*Dx1*Dx1*Dy1*Dy1) + 4*n3*Dx1*Dy1*(Dx1*Dx1 - Dy1*Dy1)
#
#     y2 = y0 + m * Dx1 + n*Dy1 -n1*(Dx1*Dx1 - Dy1*Dy1) + 2*m1*Dx1*Dy1 -\
#         n2*Dy1*(3*Dx1*Dx1 - Dy1*Dy1) - m2*Dx1*(Dx1*Dx1 - 3*Dy1*Dy1)+\
#         n3*(Dx1*Dx1*Dx1*Dx1 + Dy1*Dy1*Dy1*Dy1- 6*Dx1*Dx1*Dy1*Dy1) - 4*m3*Dx1*Dy1*(Dx1*Dx1 - Dy1*Dy1)
import  time
from osgeo import osr,gdal,ogr
from scipy.optimize import fsolve

def geo2lonlat(srs, x, y):
    '''
    将投影坐标转为经纬度坐标（具体的投影坐标系由给定数据确定）
    :param dataset: GDAL地理数据
    :param x: 投影坐标x
    :param y: 投影坐标y
    :return: 投影坐标(x, y)对应的经纬度坐标(lon, lat)
    '''
    prosrs, geosrs = srs,srs.CloneGeogCS()
    ct = osr.CoordinateTransformation(prosrs, geosrs)
    coords = ct.TransformPoint(x, y)
    return coords[:2]
def lonlat2geo(srs, lon, lat):
    '''
    将经纬度坐标转为投影坐标（具体的投影坐标系由给定数据确定）
    :param dataset: GDAL地理数据
    :param lon: 地理坐标lon经度
    :param lat: 地理坐标lat纬度
    :return: 经纬度坐标(lon, lat)对应的投影坐标
    '''
    prosrs, geosrs = srs,srs.CloneGeogCS()
    ct = osr.CoordinateTransformation(geosrs, prosrs)
    coords = ct.TransformPoint(lon, lat)
    return coords[:2]


def Ret(X0, Y0, x0, y0, x, y, X, Y):
    P1 = x - x0
    Q1 = y - y0
    P2 = P1*P1 -Q1*Q1
    Q2 = Q1*P1 +P1*Q1
    P3 = P1*P2 -Q1*Q2
    Q3 = Q1 * P2 + P1 * Q2
    P4 = P1 * P3 - Q1 * Q3
    Q4 = Q1 * P3 + P1 * Q3
    P5 = P1 * P4 - Q1 * Q4
    Q5 = Q1 * P4 + P1 * Q4
    A = np.zeros((8, 8))
    b = np.zeros((8,))
    X = X - X0
    Y = Y - Y0
    for idx, _ in enumerate(P1):
        A[2 * idx, :] = [P1[idx], (-1) * Q1[idx], P2[idx], (-1) * Q2[idx], P3[idx], (-1) * Q3[idx], P4[idx],
                         (-1) * Q4[idx]]
        A[2 * idx + 1, :] = [Q1[idx], P1[idx], Q2[idx], P2[idx], Q3[idx], P3[idx], Q4[idx], P4[idx]]
        b[2 * idx] = X[idx]
        b[2 * idx+1] = Y[idx]
    return np.linalg.solve(A, b)
from numba import jit
# @jit(nopython=True)
def Trans(X0, Y0, x0, y0, x, y, X, Y,Par):
    P1 = x - x0
    Q1 = y - y0
    P2 = P1 * P1 - Q1 * Q1
    Q2 = Q1 * P1 + P1 * Q1
    # P3 = P1 * P2 - Q1 * Q2
    # Q3 = Q1 * P2 + P1 * Q2
    # P4 = P1 * P3 - Q1 * Q3
    # Q4 = Q1 * P3 + P1 * Q3
    # P5 = P1 * P4 - Q1 * Q4
    # Q5 = Q1 * P4 + P1 * Q4

    MX = X0 + Par[0]*P1 - Par[1]*Q1 + Par[2]*P2 - Par[3]*Q2# + Par[4]*P3 - Par[5]*Q3 + Par[6]*P4 - Par[7]*Q4
    MY = Y0 + Par[1] * P1 + Par[0] * Q1 + Par[3] * P2 + Par[2] * Q2 #+ Par[5] * P3 + Par[4] * Q3 #+ Par[7] * P4 + Par[6] * Q4

    return MX,MY




srs49  = osr.SpatialReference()
srs49.ImportFromEPSG(32649)
srs50  = osr.SpatialReference()
srs50.ImportFromEPSG(32650)


geo49 = [699960.0, 10.0, 0.0, 3500040.0, 0.0, -10.0]
geo50 = [199980.0, 10.0, 0.0, 3500040.0, 0.0, -10.0]
x0 = geo49[0] + 10000 * geo49[1]
y0 = geo49[3] + 5000 * geo49[5]
X0,Y0 = lonlat2geo(srs50,geo2lonlat(srs49,x0,y0)[0],geo2lonlat(srs49,x0,y0)[1])
Col = [9000,9500,10500,11000]
Row = [5100,5500,4500,4800]
X = []
Y = []
x = []
y = []

for i in range(4):
    mx = geo49[0] + Col[i] * geo49[1]
    my = geo49[3] + Row[i] * geo49[5]
    mX0, mY0 = lonlat2geo(srs50, geo2lonlat(srs49, mx, my)[0], geo2lonlat(srs49, mx, my)[1])
    x.append(mx)
    y.append(my)
    X.append(mX0)
    Y.append(mY0)
x = np.array(x)
y = np.array(y)
X = np.array(X)
Y = np.array(Y)
Par =Ret(X0,Y0,x0,y0,x,y,X,Y)




Col = [10000,11000,13000,15000,15000,20000]
Row = [5000,5200,5200,5030,5300,20000]

X = []
Y = []
x = []
y = []

for i in range(3):
    mx = geo49[0] + Col[i] * geo49[1]
    my = geo49[3] + Row[i] * geo49[5]
    mX0, mY0 = lonlat2geo(srs50, geo2lonlat(srs49, mx, my)[0], geo2lonlat(srs49, mx, my)[1])
    x.append(mx)
    y.append(my)
    X.append(mX0)
    Y.append(mY0)
x = np.array(x)
y = np.array(y)
X = np.array(X)
Y = np.array(Y)

s = time.time()
for i in range(100000):
    TrX,TrY = Trans(X0, Y0, x0, y0, x, y, X, Y,Par)
e = time.time()
print("the cost of time is ",e-s)
for i,_  in enumerate(TrX):
    print("W->E",TrX[i],X[i])
    print("N->S",TrY[i],Y[i])
X,Y = 100000,100000
H,W = 100001,100001

mx0 = geo49[0] + X * geo49[1]
my0 = geo49[3] + Y * geo49[5]

mx1 = geo49[0] + W * geo49[1]
my1 = geo49[3] + Y * geo49[5]

mx2 = geo49[0] + W * geo49[1]
my2 = geo49[3] + H * geo49[5]


mx3 = geo49[0] + X * geo49[1]
my3 = geo49[3] + H * geo49[5]

mxc = geo49[0] + (X+0.5) * geo49[1]
myc = geo49[3] + (Y+0.5) * geo49[5]

mX0, mY0 = lonlat2geo(srs50, geo2lonlat(srs49, mx0, my0)[0], geo2lonlat(srs49, mx0, my0)[1])
mX1, mY1 = lonlat2geo(srs50, geo2lonlat(srs49, mx1, my1)[0], geo2lonlat(srs49, mx1, my1)[1])
mX2, mY2 = lonlat2geo(srs50, geo2lonlat(srs49, mx2, my2)[0], geo2lonlat(srs49, mx2, my2)[1])
mX3, mY3 = lonlat2geo(srs50, geo2lonlat(srs49, mx3, my3)[0], geo2lonlat(srs49, mx3, my3)[1])
mXC, mYC = lonlat2geo(srs50, geo2lonlat(srs49, mxc, myc)[0], geo2lonlat(srs49, mxc, myc)[1])


#
# print(np.sqrt((mx0 - mx01)*(mx0 - mx01) + (my0 - my01)*(my0 - my01)))
# print(np.sqrt((mX0 - mX01)*(mX0 - mX01) + (mY0 - mY01)*(mY0 - mY01)))
# fig, ax1 = plt.subplots()
# ax1.plot([mx0,mx1,mx2,mx3,mx0],[my0,my1,my2,my3,my0])
# ax2 =  ax1.twinx()
# ax2.
plt.plot([mX0,mX1,mX2,mX3,mX0],[mY0,mY1,mY2,mY3,mY0])
plt.scatter([mXC],[mYC])
plt.show()
print(mx0-mx1,mX0-mX1)
print(mx1-mx2,mX1-mX2)
print(mx2-mx3,mX2-mX3)
print(mx3-mx0,mX3-mX0)
print((mx0+mx1)/2-mx0,(mx2+mx3)/2-mx0,mXC-mX3)

print("**",my0-my1,mY0-mY1)
print(my1-my2,mY1-mY2)
print(my2-my3,mY2-mY3)
print(my3-my0,mY3-mY0)
print((my0+my2)/2-my3,(my0+my3)/2-my3,mYC-mY3)


