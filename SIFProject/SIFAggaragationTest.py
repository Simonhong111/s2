from SIFAggaragation import *





shp_path = r"D:\Sen2Projecton\ShapeFile\Hubei\Hubei.shp"


grid(shp_path,"2018-05-01","2018-05-31",r"D:\Sen2Projecton\NPZ")



sif = np.load(r"D:\Sen2Projecton\NPZ\2018-5SifAgg.npz")
print(sif.files)
print(sif["sif757"])
print(sif['date'])
# from osgeo import gdal,ogr,osr
# def write_img(filename, im_data):
#
#
#
#     datatype = gdal.GDT_Float32
#
#     im_width = im_data.shape[1]
#     im_height = im_data.shape[0]
#
#
#     # 创建文件
#     driver = gdal.GetDriverByName("GTiff")  # 数据类型必须有，因为要计算需要多大内存空间
#     dataset = driver.Create(filename, im_width, im_height, 1, datatype)
#
#      # 写入投影
#     dataset.GetRasterBand(1).WriteArray(im_data)  # 写入数组数据
#
#
#
#
#
#
# # for idx,sif757 in enumerate(sif["sif757"]):
# #     write_img(os.path.join(r"D:\Sen2Projecton\SIF\SIFImg",sif["date"][idx]+".tif"),sif757)
#
# # a = np.array([1,2,3,4])
# # w = a>3
# # print("**",np.where(w)[0])
#
# import h5py
#
from netCDF4 import Dataset
# fin = h5py.File(r"D:\Sen2Projecton\SIF\2018\05\oco2_LtSIF_180501_B8100r_180703004855s.nc4","r")
# lat_in = fin['latitude'][:]
# lon_in = fin['longitude'][:]
# daily_corr_in = fin["daily_correction_factor"][:]
# sif_771nm_in = fin["SIF_771nm"][:]
# sif_757nm_in = fin["SIF_757nm"][:]
# mode = fin['measurement_mode'][:]
#
# MAX_LONGITUDE = 180
# MIN_LONGITUDE = -180
# MAX_LATITUDE = 90
# MIN_LATITUDE = -90
#
# # 经纬度间隔
# LON_INTERVAL = 0.05
# LAT_INTERVAL = 0.05
#
# driver = ogr.GetDriverByName('ESRI Shapefile')
# dataSource = driver.Open(shp_path, 0)

# 记住这里值获取一个layer和第一个feature，
# 所以正如上面所说，需要将湖北省矢量都聚合了，只要外部边界，只有一个特征，
# 在envi或者arcgis上面做

layer = dataSource.GetLayer()

# # 湖北省矢量边界的外接矩形
# extent = layer.GetExtent()
# # extent    minlon,maxlon,minlat,maxlat
# print("湖北省外界矩形", extent)
#
# USER_MAX_LONGITUDE = extent[1]
# USER_MIN_LONGIUDE = extent[0]
# USER_MAX_LATITUDE = extent[3]
# USER_MIN_LATITUDE = extent[2]
#
# user_min_Lon_id = int(np.floor((USER_MIN_LONGIUDE - MIN_LONGITUDE) / LON_INTERVAL))
# user_max_Lon_id = int(np.ceil((USER_MAX_LONGITUDE - MIN_LONGITUDE) / LON_INTERVAL))
# user_min_Lat_id = int(np.floor((USER_MIN_LATITUDE - MIN_LATITUDE) / LAT_INTERVAL))
# user_max_Lat_id = int(np.ceil((USER_MAX_LATITUDE - MIN_LATITUDE) / LAT_INTERVAL))
#
# user_min_Lon = user_min_Lon_id * LON_INTERVAL + MIN_LONGITUDE
# user_max_Lon = user_max_Lon_id * LON_INTERVAL + MIN_LONGITUDE
# user_min_Lat = user_min_Lat_id * LAT_INTERVAL + MIN_LATITUDE
# user_max_Lat = user_max_Lat_id * LAT_INTERVAL + MIN_LATITUDE
# lat = np.arange(user_min_Lat, user_max_Lat - 1e-8, LAT_INTERVAL)
# lon = np.arange(user_min_Lon, user_max_Lon - 1e-8, LON_INTERVAL)
#
# print("lat", lat, lat.shape)
# ind2 = (lat_in >= user_min_Lat) & (lat_in <= user_max_Lat) & (
#                 lon_in <= user_max_Lon) & (lon_in >= user_min_Lon) & (mode == 0)
# ind = np.array(np.where(ind2)[0])
# iLat = np.asarray(np.floor(
#                     ((np.asarray(lat_in[ind]) - user_min_Lat) / (user_max_Lat - user_min_Lat) * len(lat))),
#                                   dtype=int)
# iLon = np.asarray(np.floor(
#                     ((np.asarray(lon_in[ind]) - user_min_Lon) / (user_max_Lon - user_min_Lon) * len(lon))),
#                                   dtype=int)
# vec_n = np.zeros((len(lat),len(lon)))
# vec_sif757 = np.zeros((len(lat),len(lon)))
# vec_sif771 = np.zeros((len(lat),len(lon)))
# vec_sif_avg = np.zeros((len(lat),len(lon)))
# vec_sif_norm = np.zeros((len(lat),len(lon)))
# for i in range(len(ind)):
#     print(iLat[i],iLon[i])
#     vec_n[iLat[i],iLon[i]]  += 1
#     vec_sif757[iLat[i],iLon[i]]  += sif_757nm_in[ind[i]]
#     vec_sif771[iLat[i], iLon[i]] += sif_771nm_in[ind[i]]
#     vec_sif_avg[iLat[i], iLon[i]] += (sif_757nm_in[ind[i]]*0.5+1.5/2*sif_771nm_in[ind[i]])
#     vec_sif_norm[iLat[i], iLon[i]] += (sif_757nm_in[ind[i]] * 0.5 + 1.5 / 2 * sif_771nm_in[ind[i]]) * daily_corr_in[ind[i]]
#     print(vec_n[iLat[i],iLon[i]])
#
# vec_n[vec_n < 3] =1
# vec_sif757[vec_n == 1] = -99
# vec_sif757 /= vec_n
# vec_sif771[vec_n == 1] = -99
# vec_sif771 /= vec_n
# vec_sif_avg[vec_n == 1] = -99
# vec_sif_avg /= vec_n
#
# vec_sif_norm[vec_n == 1] = -99
# vec_sif_norm /= vec_n
#
#
# W = np.where(vec_sif757 !=-99)
# for i in range(len(W[0])):
#
#     row = 86-W[0][i]
#     col = W[1][i]
#     print(row,col)
#     id = row * len(lon) + col+1
#
#     print(id,vec_sif757[W[0][i],W[1][i]])

# print(vec_sif757[W])
# print(vec_sif771[W])
# print(vec_sif_avg[W])
# print(vec_sif_norm[W])

# print("*",vec_sif757)

# print("757",sif["sif757"][np.where(sif['date']== '20180501')])
# print("771",sif["sif771"][np.where(sif['date']== '20180501')])
# print("id",sif["sifgridId"][np.where(sif['date']== '20180501')])


# print("avg",sif["vec_SIF_avg"][np.where(sif['date']== '20180501')])
# print("norm",sif["vec_SIF_avg_norm"][np.where(sif['date']== '20180501')])