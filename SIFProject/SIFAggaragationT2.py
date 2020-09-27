from SIFAggaragation import *

shp_path = r"D:\Satellive\SHP\Hubei\Hubei.shp"
# grid(shp_path,"2015-08-01","2018-10-31",r"D:\Satellive\NPZ")
sif = np.load(r"D:\Satellive\NPZ\2015-2018-SifAgg.npz")
print(sif.files)
print(sif['vec_SIF_avg_norm'].shape)
print("************************")
# for idx,_ in enumerate(sif["sif757"]):
#
#     print("zijid",sif["sif757"][idx])
#     print(idx)
#     print(sif['date'][idx])
#     # print(sif["EtRing"][30])
# print("************************")
# fin = h5py.File(r"D:\Satellive\SIFdata\2015\08\oco2_LtSIF_150802_B8100r_171011082243s.nc4","r")
# lat_in = fin['latitude'][:]
#
# lon_in = fin['longitude'][:]
# daily_corr_in = fin["daily_correction_factor"][:]
# sif_771nm_in = fin["SIF_771nm"][:]
# sif_757nm_in = fin["SIF_757nm"][:]
# mode = fin['measurement_mode'][:]
#
#
# for idx,s757 in enumerate(sif["sif757"][0:1]):
#     print("zijid",sif["sif757"][idx])
#     print("date",sif['date'][idx])
#     print("extent",sif["EtRing"][idx])
#     s757 = sif["sif757"][idx]
#     d = sif['date'][idx]
#     et = sif["EtRing"][idx]
#     ind2 = (lat_in >= et[5]) & (lat_in < et[1]) & \
#            (lon_in < et[2]) & (lon_in >= et[0]) & (mode == 0)
#     print("shiji",sif_757nm_in[ind2])
#     s7 = np.array(sif_757nm_in[ind2])
#     print("jisuande mean",s7.mean())


















