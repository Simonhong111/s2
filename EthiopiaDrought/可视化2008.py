from epdatapro20082015 import *
from matplotlib import pyplot as plt
import os
import numpy as np
import glob
from matplotlib import cm

months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

yy = "2017"
rp = r"D:\Cornell\EthiopianDrought\Chirps2"
rpp = r"D:\Cornell\EthiopianDrought\AData\Chirps2Pars"
ep = r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia"
epp = r"D:\Cornell\EthiopianDrought\AData\MOD13C2.006EthiopiaAnomalyPars"
pvip = r"D:\Cornell\EthiopianDrought\AData\PVI2003-2018"
pvipp = r"D:\Cornell\EthiopianDrought\AData\PviParms"
fig = plt.figure(figsize=(14, 3))
plt.title("{} PVI Anomaly Map".format(yy) + '\n', fontsize=14)

# for m in range(1, 13):
#     ramn = chirpsAnomalyMap(rp, rpp, yy, str(m))
#     ax = fig.add_subplot(3, 4, m )
#     ax.set_title("Precipitaion_" + months[m - 1])
#     mask = np.where(ramn > -9999)
#     ramn[ramn == -9999] = np.nan
#     vmin = ramn[mask].min()
#     vmax = ramn[mask].max()
#     print("maxmin value", vmin, vmax)
#     cmap = plt.get_cmap("RdBu")
#
#     cax = ax.imshow(ramn, cmap=cmap, vmin=-2.5, vmax=2.5)
#     cbar = plt.colorbar(cax,ax=ax,fraction=0.036,pad=0.04)
#     ax.set_xticks([])
#     ax.set_yticks([])


# for m in range(1, 13):
#     eamn = eviAnomalyMap(ep, epp, yy, str(m))
#     ax = fig.add_subplot(3, 4, m)
#     ax.set_title("EVI_" + months[m - 1])
#     mask = np.where(eamn > -9999)
#     eamn[eamn == -9999] = np.nan
#     vmin = eamn[mask].min()
#     vmax = eamn[mask].max()
#     #     print("maxmin value",vmin,vmax)
#     cax = ax.imshow(eamn, cmap=plt.get_cmap("RdBu"), vmin=-2.5, vmax=2.5)
#     cbar = plt.colorbar(cax,ax=ax,fraction=0.036,pad=0.04)
#     ax.set_xticks([])
#     ax.set_yticks([])


eamn = PviAnomalyMap(pvip, pvipp, yy,pvitype="all")
ax = fig.add_subplot(1,3,1)
ax.set_title("PVI_Yearly")
mask = np.where(eamn > -9999)
eamn[eamn == -9999] = np.nan
vmin = eamn[mask].min()
vmax = eamn[mask].max()
print("maxmin value",vmin,vmax)
cax = ax.imshow(eamn, cmap=plt.get_cmap("RdBu"), vmin=-2.5, vmax=2.5)
cbar = plt.colorbar(cax,ax=ax,fraction=0.036,pad=0.04)
ax.set_xticks([])
ax.set_yticks([])

eamn = PviAnomalyMap(pvip, pvipp, yy,pvitype="short")
ax = fig.add_subplot(1,3,2)
ax.set_title("PVI_Short_rains")
mask = np.where(eamn > -9999)
eamn[eamn == -9999] = np.nan
vmin = eamn[mask].min()
vmax = eamn[mask].max()
print("maxmin value",vmin,vmax)
cax = ax.imshow(eamn, cmap=plt.get_cmap("RdBu"), vmin=-2.5, vmax=2.5)
cbar = plt.colorbar(cax,ax=ax,fraction=0.036,pad=0.04)
ax.set_xticks([])
ax.set_yticks([])

eamn = PviAnomalyMap(pvip, pvipp, yy,pvitype="long")
ax = fig.add_subplot(1,3,3)
ax.set_title("PVI_Long_rains")
mask = np.where(eamn > -9999)
eamn[eamn == -9999] = np.nan
vmin = eamn[mask].min()
vmax = eamn[mask].max()
print("maxmin value",vmin,vmax)
cax = ax.imshow(eamn, cmap=plt.get_cmap("RdBu"), vmin=-2.5, vmax=2.5)
cbar = plt.colorbar(cax,ax=ax,fraction=0.036,pad=0.04)
ax.set_xticks([])
ax.set_yticks([])

fig.tight_layout()  # 调整整体空白
plt.show()