from epdatapro20082015v2 import *
from matplotlib import pyplot as plt
from matplotlib import cm
months =["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]


rp = r"D:\Cornell\EthiopianDrought\Chirps2"
rpp = r"D:\Cornell\EthiopianDrought\AData\Chirps2ParsSeanson"

ep = r"D:\Cornell\EthiopianDrought\MOD13C2.006Ethiopia"
epp = r"D:\Cornell\EthiopianDrought\AData\MOD13C2.006EthiopiaAnomalyParsSeason"

sp = r"D:\Cornell\EthiopianDrought\AData\NewSIF0318AnomSplit"
spp = r"D:\Cornell\EthiopianDrought\AData\NewSIFMonthParsSeason"

gsif =r'D:\Cornell\EthiopianDrought\AData\GOSIF0317AnomSplit'
gsifp =r'D:\Cornell\EthiopianDrought\AData\GOSIFMonthParsSeason'

pvip = r"D:\Cornell\EthiopianDrought\AData\PVI2003-2018"
pvipp = r"D:\Cornell\EthiopianDrought\AData\PviParms"

yy = "2017"
mm = "long"
fig = plt.figure(figsize=(14, 8))
plt.title("{} Anomaly Map".format(yy) + '\n', fontsize=16)

ramn = chirpsAnomalyMap(rp, rpp, yy, str(mm))
ax = fig.add_subplot(2, 3, 1)
ax.set_title("Precipitaion Anomaly " + mm +" rains")
mask = np.where(ramn > -9999)
ramn[ramn == -9999] = np.nan
vmin = ramn[mask].min()
vmax = ramn[mask].max()
print("maxmin value", vmin, vmax)
cmap = plt.get_cmap("RdBu")

cax = ax.imshow(ramn, cmap=cmap, vmin=-2.5, vmax=2.5)
cbar = plt.colorbar(cax,ax=ax,fraction=0.036,pad=0.04)
ax.set_xticks([])
ax.set_yticks([])




samn = NewSIFAnomalyMap(sp, spp, yy, str(mm))
ax1 = fig.add_subplot(2, 3, 2)
ax1.set_title("New SIF Anomaly " + mm +" rains")
mask = np.where(samn > -9999)
samn[samn == -9999] = np.nan
vmin = samn[mask].min()
vmax = samn[mask].max()
print("maxmin value", vmin, vmax)
cax1 = ax1.imshow(samn, cmap=plt.get_cmap("RdBu"), vmin=-2.5, vmax=2.5)
cbar1 = plt.colorbar(cax1,ax=ax1,fraction=0.036,pad=0.04)
ax1.set_xticks([])
ax1.set_yticks([])



eamn = eviAnomalyMap(ep, epp, yy, str(mm))
ax2 = fig.add_subplot(2, 3, 3)
ax2.set_title("EVI Anomaly"  + mm +" rains")
mask = np.where(eamn > -9999)
eamn[eamn == -9999] = np.nan
vmin = eamn[mask].min()
vmax = eamn[mask].max()
print("maxmin value",vmin,vmax)
cax2 = ax2.imshow(eamn, cmap=plt.get_cmap("RdBu"), vmin=-2.5, vmax=2.5)
cbar2 = plt.colorbar(cax2,ax=ax2,fraction=0.036,pad=0.04)
ax2.set_xticks([])
ax2.set_yticks([])

eamn = GOSIFAnomalyMap(gsif, gsifp, yy, str(mm))
ax3 = fig.add_subplot(2, 3, 4)
ax3.set_title("GOSIF Anomaly " + mm +" rains")
mask = np.where(eamn > -9999)
eamn[eamn == -9999] = np.nan
vmin = eamn[mask].min()
vmax = eamn[mask].max()
print("maxmin value",vmin,vmax)
cax3 = ax3.imshow(eamn, cmap=plt.get_cmap("RdBu"), vmin=-2.5, vmax=2.5)
cbar3 = plt.colorbar(cax3,ax=ax3,fraction=0.036,pad=0.04)
ax3.set_xticks([])
ax3.set_yticks([])

eamn = PviAnomalyMap(pvip, pvipp, yy,pvitype=mm)
ax4 = fig.add_subplot(2,3,5)
ax4.set_title("PVI Anomaly "+mm +" rains")
mask = np.where(eamn > -9999)
eamn[eamn == -9999] = np.nan
vmin = eamn[mask].min()
vmax = eamn[mask].max()
print("maxmin value",vmin,vmax)
cax4 = ax4.imshow(eamn, cmap=plt.get_cmap("RdBu"), vmin=-2.5, vmax=2.5)
cbar4 = plt.colorbar(cax4,ax=ax4,fraction=0.036,pad=0.04)
ax4.set_xticks([])
ax4.set_yticks([])

fig.tight_layout()  # 调整整体空白
plt.show()