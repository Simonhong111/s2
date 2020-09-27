#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
import numpy as np
import h5py
import time
from netCDF4 import Dataset
from dateutil import rrule
from datetime import *
import time
import glob
from numba import jit
from osgeo import gdal,ogr,osr


dict_names = {
    'source': 'OCO-2 L1B',
    'references': '...',
    'Conventions': 'CF-1.6',
    'product_version': 'v1.0',
    'summmary': 'Fraunhofer-line based SIF retrievals',
    'keyword': 'satellite, OCO-2, Solar Induced Fluorescence, SIF',
    'keywords_vocabulary': 'NASA Global Change Master Directory (GCMD)',
    'cdm_data_type': 'grids',
    'comment': 'These data were produced at JPL/Caltech',
    'date_created': 'Created ' + time.ctime(time.time()),
    'creator_name': 'Caltech, Christian Frankenberg',
    'creator_email': 'cfranken@caltech.edu',
    'project': 'OCO-2 NASA/JPL',
    'geospatial_lat_min': '-90.0f; // float',
    'geospatial_lat_max': '90.0f; // float',
    'geospatial_lat_units': 'degrees_north',
    'geospatial_lon_min': '-180.0f; // float',
    'geospatial_lon_max': '180.0f; // float',
    'geospatial_lon_units': 'degrees_east',
    'geospatial_vertical_min': '0.0f; // float',
    'geospatial_vertical_max': '100000.0; // float',
    'standard_name_vocabulary': 'NetCDF Climate and Forecast (CF) Metadata Conventions Version 1.6',
    'platform': 'OCO-2',
    'sensor': 'OCO-2',
    'spatial_resolution': '2km x 1.3km at nadir (along-track x across-track)',
    '                   CoordSysBuilder': 'ucar.nc2.dataset.conv.CF1Convention',
}

dict_l2 = {
    'sif_757nm': 'SIF_757nm',
    'sif_771nm': 'SIF_771nm',
    'sif_757_1sigma': 'SIF_757nm_uncert',
    'sif_771_1sigma': 'SIF_771nm_uncert',
    #  'latitude': 'latitude',
    #  'longitude': 'longitude',
    'sif_757_rel': 'SIF_757nm_relative',
    'sif_771_rel': 'SIF_771nm_relative',
    'cont_radiance_757': 'continuum_radiance_757nm',
    'cont_radiance_771': 'continuum_radiance_771nm',
    'sza': 'solar_zenith_angle',
    'lza': 'sensor_zenith_angle',
    'co2_ratio': 'Cloud/co2_ratio',
    'o2_ratio': 'Cloud/o2_ratio',
    'daily_corr': 'daily_correction_factor',
    'T2m': 'Meteo/2m_temperature',
    'VPD_2m': 'Meteo/vapor_pressure_deficit',
    'Tskin': 'Meteo/skin_temperature',
    'altitude': 'surface_altitude',
}
from GridGen import GridDefination

def grid(shp_path,start,stop,save_sif_path,dTime=1):
    MAX_LONGITUDE = 180
    MIN_LONGITUDE = -180
    MAX_LATITUDE = 90
    MIN_LATITUDE = -90

    # 经纬度间隔
    LON_INTERVAL = 0.05
    LAT_INTERVAL = 0.05

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shp_path, 0)

    # 记住这里值获取一个layer和第一个feature，
    # 所以正如上面所说，需要将湖北省矢量都聚合了，只要外部边界，只有一个特征，
    # 在envi或者arcgis上面做
    layer = dataSource.GetLayer()
    # 湖北省矢量边界的外接矩形
    extent = layer.GetExtent()
    # extent    minlon,maxlon,minlat,maxlat
    print("湖北省外界矩形", extent)

    USER_MAX_LONGITUDE = extent[1]
    USER_MIN_LONGIUDE = extent[0]
    USER_MAX_LATITUDE = extent[3]
    USER_MIN_LATITUDE = extent[2]

    user_min_Lon_id = int(np.floor((USER_MIN_LONGIUDE - MIN_LONGITUDE) / LON_INTERVAL))
    user_max_Lon_id = int(np.ceil((USER_MAX_LONGITUDE - MIN_LONGITUDE) / LON_INTERVAL))
    user_min_Lat_id = int(np.floor((USER_MIN_LATITUDE - MIN_LATITUDE) / LAT_INTERVAL))
    user_max_Lat_id = int(np.ceil((USER_MAX_LATITUDE - MIN_LATITUDE) / LAT_INTERVAL))

    user_min_Lon = user_min_Lon_id * LON_INTERVAL + MIN_LONGITUDE
    user_max_Lon = user_max_Lon_id * LON_INTERVAL + MIN_LONGITUDE
    user_min_Lat = user_min_Lat_id * LAT_INTERVAL + MIN_LATITUDE
    user_max_Lat = user_max_Lat_id * LAT_INTERVAL + MIN_LATITUDE


    print("llll",user_min_Lon,user_max_Lon,user_min_Lat,user_max_Lat)
    # Define spatial grid:
    # 因为
    lat = np.arange(user_min_Lat, user_max_Lat + 1e-8, LAT_INTERVAL)
    lon = np.arange(user_min_Lon, user_max_Lon + 1e-8, LON_INTERVAL)

    print("lat",lat,lat.shape)
    print("lon",lon,lon.shape)
    start = datetime.strptime(start, "%Y-%m-%d").date()
    stop = datetime.strptime(stop, "%Y-%m-%d").date()


    vec_sif_757nm = np.zeros((len(lat)-1, len(lon)-1))
    vec_sif_771nm = np.zeros((len(lat)-1, len(lon)-1))
    vec_daily_corr = np.zeros((len(lat)-1, len(lon)-1))
    vec_n = np.zeros((len(lat)-1, len(lon)-1))
    vec_SIF_avg = np.zeros((len(lat)-1, len(lon)-1))
    vec_SIF_avg_norm = np.zeros((len(lat)-1, len(lon)-1))

    date = []
    sif757matrix = []
    sif771matrix = []
    sifgridId = []
    vec_SIF_avg_mat = []
    vec_SIF_avg_norm_mat =[]
    ETRing =[]
    samp_num = 0
    for dt in (rrule.rrule(rrule.DAILY, interval=dTime, dtstart=start, until=stop)):

        vec_n[:] = 0
        vec_sif_757nm[:] = 0
        vec_sif_771nm[:] = 0
        vec_daily_corr[:] = 0
        vec_SIF_avg[:] = 0
        vec_SIF_avg_norm[:] = 0

        files = glob.glob(os.path.join(r"D:\Satellive\SIFdata",dt.strftime('%Y'),dt.strftime('%m'),'oco2_*_??' + dt.strftime('%m%d') + '_*.nc4'))

        for file in files:

            fin = h5py.File(file,"r")
            lat_in = fin['latitude'][:]

            lon_in = fin['longitude'][:]
            daily_corr_in = fin["daily_correction_factor"][:]
            sif_771nm_in = fin["SIF_771nm"][:]
            sif_757nm_in = fin["SIF_757nm"][:]
            mode = fin['measurement_mode'][:]

            # 只要湖北省的数据
            ind2 = (lat_in >= user_min_Lat) & (lat_in <= user_max_Lat) & (
                lon_in <= user_max_Lon) & (lon_in >= user_min_Lon) & (mode == 0)

            ind = np.array(np.where(ind2)[0])

            if len(ind) > 0:

                iLat = np.asarray(np.floor(
                    ((np.asarray(lat_in[ind]) - user_min_Lat) / (user_max_Lat - user_min_Lat) * (len(lat)-1))),
                                  dtype=int)
                iLon = np.asarray(np.floor(
                    ((np.asarray(lon_in[ind]) - user_min_Lon) / (user_max_Lon - user_min_Lon) * (len(lon)-1))),
                                  dtype=int)

                favg(vec_n, iLat, iLon, ind, np.ones(len(sif_757nm_in), ), len(ind))
                # 第一个参数用来存储数据，以左下角为起点
                favg(vec_sif_757nm, iLat, iLon, ind, sif_757nm_in, len(ind))
                favg(vec_sif_771nm, iLat, iLon, ind, sif_771nm_in, len(ind))
                favg(vec_daily_corr, iLat, iLon, ind, daily_corr_in, len(ind))
                favg(vec_SIF_avg, iLat, iLon, ind, 0.5 * sif_757nm_in + 1.5 / 2. * sif_771nm_in, len(ind))
                favg(vec_SIF_avg_norm, iLat, iLon, ind, (0.5 * sif_757nm_in + 1.5 / 2. * sif_771nm_in) * daily_corr_in,
                     len(ind))
                # 上面的例子讲解的是计算了每一个矩阵里面的SIF数据的总和
                fin.close()

        if len(ind) >0:
            vec_n[vec_n < 3] = 1
            if len(np.where(vec_n != 1)[0]) > 0:
                print("valid file name is ", file)
                samp_num += len(np.where(vec_n != 1)[0])
                vec_sif_757nm[vec_n == 1] = -999.9
                vec_sif_757nm /= vec_n
                vec_sif_771nm[vec_n == 1] = -999.9
                vec_sif_771nm /= vec_n
                vec_daily_corr[vec_n == 1] = -999.9
                vec_daily_corr /= vec_n
                vec_SIF_avg[vec_n == 1] = -999.9
                vec_SIF_avg /= vec_n
                vec_SIF_avg_norm[vec_n == 1] = -999.9
                vec_SIF_avg_norm /= vec_n
                vec_sif_757nm_copy=vec_sif_757nm.copy()
                vec_sif_771nm_copy=vec_sif_771nm.copy()
                vec_SIF_avg_copy = vec_SIF_avg.copy()
                vec_SIF_avg_norm_copy = vec_SIF_avg_norm.copy()
                W = np.where(vec_sif_757nm_copy != -999.9)
                # print("W",W)

                for idx,row in enumerate(W[0]):
                    row1 = len(lat) - row
                    col = W[1][idx]
                    ulx = lon[0] + col*0.05
                    uly = lat[0] + row*0.05 + 0.05
                    lrx = ulx + 0.05
                    lry = uly - 0.05
                    temp =[ulx,uly,lrx,uly,lrx,lry,ulx,lry,ulx,uly]
                    ETRing.append(temp)
                    temp =[]
                    id = row*(len(lon)-1) + col +1
                    sifgridId.append(id)
                    month = str(dt.month)
                    if len(str(dt.month)) ==1:
                        month = "0" + month
                    day = str(dt.day)
                    if len(str(dt.day)) ==1:
                        day = "0" + day
                    date.append(str(dt.year)+month+day)
                    sif757matrix.append(vec_sif_757nm_copy[row,col])
                    sif771matrix.append(vec_sif_771nm_copy[row,col])
                    vec_SIF_avg_mat.append(vec_SIF_avg_copy[row,col])
                    vec_SIF_avg_norm_mat.append(vec_SIF_avg_norm_copy[row,col])
                    print("one day has benn handled")

    kwargs = {"date": date,"sif757":sif757matrix,"sif771":sif771matrix,"sifgridId":sifgridId,"vec_SIF_avg":vec_SIF_avg_mat,"vec_SIF_avg_norm":vec_SIF_avg_norm_mat,\
              "EtRing":ETRing}
    print(os.path.join(save_sif_path,str(dt.year)+"-"+str(dt.month))+"SifAgg.npz")
    # np.savez(os.path.join(save_sif_path,str(dt.year)+"-"+str(dt.month))+"SifAgg.npz", **kwargs)
    np.savez(os.path.join(save_sif_path, "2015-2018-SifAgg.npz"), **kwargs)
    print("num",samp_num)

# @jit(nopython=True, parallel=True)

def favg(arr, ix, iy, iz, inp, s):
    for i in range(s):
        arr[ix[i], iy[i]] += inp[iz[i]]

    return arr