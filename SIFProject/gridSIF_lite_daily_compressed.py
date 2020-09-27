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

                  
def grid(options, args):
    
    # Define spatial grid:
    dlat = np.arange(options.latMin, options.latMax+1e-8, options.dLat)
    dlon = np.arange(options.lonMin, options.lonMax+1e-8, options.dLon)
    lat = np.arange(options.latMin+options.dLat/2., options.latMax+1e-8-options.dLat/2., options.dLat)
    lon = np.arange(options.lonMin+options.dLon/2., options.lonMax+1e-8-options.dLat/2., options.dLon)
    #print lat,lon
    # how many time-slices for now
    start = datetime.strptime(options.start, "%Y-%m-%d").date()
    stop = datetime.strptime(options.stop, "%Y-%m-%d").date()
    print(start,stop)
    nT = 0
    # Check out times to be interated
    for dt in rrule.rrule(rrule.DAILY, interval=options.dTime, dtstart=start, until=stop):
      #print dt
      nT+=1

    print('Time dimension ' , nT)
    # Generate simple numpy arrays for averaging
    
    vec_sif_757nm=np.zeros((len(lat),len(lon)))
    vec_sif_771nm=np.zeros((len(lat),len(lon)))
    vec_sif_757_1sigma=np.zeros((len(lat),len(lon)))
    vec_sif_771_1sigma=np.zeros((len(lat),len(lon)))
    vec_sif_757_rel=np.zeros((len(lat),len(lon)))
    vec_sif_771_rel=np.zeros((len(lat),len(lon)))
    vec_cont_radiance_757=np.zeros((len(lat),len(lon)))
    vec_cont_radiance_771=np.zeros((len(lat),len(lon)))
    vec_sza=np.zeros((len(lat),len(lon)))
    vec_lza=np.zeros((len(lat),len(lon)))
    vec_co2_ratio=np.zeros((len(lat),len(lon)))
    vec_o2_ratio=np.zeros((len(lat),len(lon)))
    vec_daily_corr=np.zeros((len(lat),len(lon)))
    vec_T2m=np.zeros((len(lat),len(lon)))
    vec_VPD_2m=np.zeros((len(lat),len(lon)))
    vec_Tskin=np.zeros((len(lat),len(lon)))
    vec_altitude=np.zeros((len(lat),len(lon)))
    vec_n = np.zeros((len(lat),len(lon)))  
    vec_SIF_avg = np.zeros((len(lat),len(lon)))  
    vec_SIF_avg_norm = np.zeros((len(lat),len(lon)))    
    # create netCDF4 file:
    f = Dataset(options.outFile, 'w', format='NETCDF4')
    time = f.createDimension('time', None)
    lati = f.createDimension('lat', len(lat))
    loni = f.createDimension('lon', len(lon))
    times = f.createVariable('time','f8',('time',))
    latitudes = f.createVariable('lat','f8',('lat'))
    longitudes = f.createVariable('lon','f8',('lon'))
    latitudes[:] =  lat
    longitudes[:] = lon
    # Add units, long names, etc.
    times.units = 'days since 2014-1-1 0:0:0'
    latitudes.units="degrees_north"
    longitudes.units="degrees_east"
    latitudes.standard_name = "latitude"
    longitudes.standard_name="longitude"
    latitudes.axis = "Y"
    longitudes.axis="X"
    times.long_name = "time"
    latitudes.long_name="latitude"
    longitudes.long_name="longitude"

    Tskin=f.createVariable("Tskin","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    o2_ratio=f.createVariable("o2_ratio","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    lza=f.createVariable("lza","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    sif_757_rel=f.createVariable("sif_757_rel","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    cont_radiance_757=f.createVariable("cont_radiance_757","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    sif_771_1sigma=f.createVariable("sif_771_1sigma","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    sif_771_rel=f.createVariable("sif_771_rel","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    cont_radiance_771=f.createVariable("cont_radiance_771","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    VPD_2m=f.createVariable("VPD_2m","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    sif_771nm=f.createVariable("sif_771nm","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    altitude=f.createVariable("altitude","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    co2_ratio=f.createVariable("co2_ratio","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    sif_757_1sigma=f.createVariable("sif_757_1sigma","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    T2m=f.createVariable("T2m","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    sif_757nm=f.createVariable("sif_757nm","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    daily_corr=f.createVariable("daily_corr","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    sza=f.createVariable("sza","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)
    
    #for i in dict_names:
    #    cmd = 'f.'+ i + '="' +dict_names[i]+'"'
    #    try:
    #        exec(cmd)
    ##    except:
     #       print('Error executing ' + cmd)
     
    # Some basic stuff
    
    
    
    # Basic variables to be generated:
    #for i in dict_l2:
    #    cmd = i + '=f.createVariable("'+i+'","f4",("time","lat","lon",), zlib=True,least_significant_digit=4,fill_value=-999.9)'
    #    try:
    #        print(cmd)
    #    except:
    #        print('Error executing ' + cmd)
    n = f.createVariable('N','f4',('time','lat','lon',), zlib=True,least_significant_digit=4,fill_value=-999.9)
    SIF_a = f.createVariable('SIF_avg','f4',('time','lat','lon',), zlib=True,least_significant_digit=4,fill_value=-999.9)
    SIF_norm = f.createVariable('SIF_avg_daily','f4',('time','lat','lon',), zlib=True,least_significant_digit=4,fill_value=-999.9)
    
    # loop over files:
    counter = 0
    counter_time = 0

    # Start looping over files:
    counter = -1
    for dt in (rrule.rrule(rrule.DAILY, interval=options.dTime,  dtstart=start, until=stop)):
        # Set to 0 again
        vec_n[:]=0
        vec_sif_757nm[:]=0
        vec_sif_771nm[:]=0
        vec_sif_771_1sigma[:]=0
        vec_sif_757_1sigma[:]=0
        vec_sif_757_rel[:]=0
        vec_sif_771_rel[:]=0
        vec_cont_radiance_757[:]=0
        vec_cont_radiance_771[:]=0
        vec_sza[:]=0
        vec_lza[:]=0
        vec_co2_ratio[:]=0
        vec_daily_corr[:]=0
        vec_T2m[:]=0
        vec_VPD_2m[:]=0
        vec_Tskin[:]=0
        vec_altitude[:]=0
        vec_SIF_avg[:]=0
        vec_SIF_avg_norm[:]=0

        files = glob.glob(options.folder + dt.strftime('%Y/%m/')+ 'oco2_*_??'+dt.strftime('%m%d')+'_*.nc4')
        #print(options.folder + dt.strftime('%Y/%m/')+ 'oco2_*_??'+dt.strftime('%m%d')+'_*.nc4')
        counter +=1
        #date_center = dt + timedelta(days=np.floor(options.dTime/2.))
        date_start = datetime(2014,1,1)
        delta = dt - date_start
        times[counter] = delta.days
        for file in files:
            #print(file)
            fin = h5py.File(file,'r')

            lat_in = fin['latitude'][:]
            lon_in = fin['longitude'][:]
            biome = fin['IGBP_index'][:]
            mode = fin['measurement_mode'][:]
            eps =0.01
            ind2 = (lat_in>options.latMin+eps)&(lat_in<options.latMax-eps)&(lon_in<options.lonMax-0.01)&(lon_in>options.lonMin+0.01)& ((mode== options.mode) | (options.mode == -1))& ((biome== options.biome) | (options.biome == -1))
            ind = np.array(np.where(ind2)[0])
            print(file, counter , len(ind),)
            if len(ind)>0:
               # print len(ind)
                iLat = np.asarray(np.floor(((np.asarray(lat_in[ind])-options.latMin)/(options.latMax-options.latMin)*len(lat))),dtype=int)
                iLon = np.asarray(np.floor(((np.asarray(lon_in[ind])-options.lonMin)/(options.lonMax-options.lonMin)*len(lon))),dtype=int)
                #wo = np.where((iLon>119))[0]
               # print iLat
                #print iLat[wo], iLon[wo], lat_in[ind[wo]], lon_in[ind[wo]]
                #print np.max(iLon), np.max(iLat), len(lon), len(lat), np.max(lon_in)
                iT = counter
                #index_vector = np.asarray((iLon*len(lat)+iLat), dtype=int);
                
                o2_ratio_in=fin["Cloud/o2_ratio"][:]
                sif_771_1sigma_in=fin["SIF_771nm_uncert"][:]
                lza_in=fin["sensor_zenith_angle"][:]
                Tskin_in=fin["Meteo/skin_temperature"][:]
                cont_radiance_757_in=fin["continuum_radiance_757nm"][:]
                cont_radiance_771_in=fin["continuum_radiance_771nm"][:]
                daily_corr_in=fin["daily_correction_factor"][:]
                co2_ratio_in=fin["Cloud/co2_ratio"][:]
                longitude_in=fin["longitude"][:]
                latitude_in=fin["latitude"][:]
                altitude_in=fin["surface_altitude"][:]
                sif_771nm_in=fin["SIF_771nm"][:]
                VPD_2m_in=fin["Meteo/vapor_pressure_deficit"][:]
                sif_757_rel_in=fin["SIF_757nm_relative"][:]
                sif_757_1sigma_in=fin["SIF_757nm_uncert"][:]
                sif_771_rel_in=fin["SIF_771nm_relative"][:]
                sza_in=fin["solar_zenith_angle"][:]
                sif_757nm_in=fin["SIF_757nm"][:]
                T2m_in=fin["Meteo/2m_temperature"][:]
                
                
                favg(vec_n,iLat,iLon,ind,np.ones(len(ind),),len(ind))         
                favg(vec_sif_757nm,iLat,iLon,ind,sif_757nm_in,len(ind))
                favg(vec_sif_771nm,iLat,iLon,ind,sif_771nm_in,len(ind))
                favg(vec_sif_771_1sigma,iLat,iLon,ind,sif_771_1sigma_in,len(ind))
                favg(vec_sif_757_1sigma,iLat,iLon,ind,sif_757_1sigma_in,len(ind))
                favg(vec_sif_757_rel,iLat,iLon,ind,sif_757_rel_in,len(ind))
                favg(vec_sif_771_rel,iLat,iLon,ind,sif_771_rel_in,len(ind))
                favg(vec_cont_radiance_771,iLat,iLon,ind,cont_radiance_771_in,len(ind))
                favg(vec_sza,iLat,iLon,ind,sza_in,len(ind))
                favg(vec_lza,iLat,iLon,ind,lza_in,len(ind))
                favg(vec_co2_ratio,iLat,iLon,ind,co2_ratio_in,len(ind))
                favg(vec_daily_corr,iLat,iLon,ind,daily_corr_in,len(ind))
                favg(vec_T2m,iLat,iLon,ind,T2m_in,len(ind))
                favg(vec_VPD_2m,iLat,iLon,ind,VPD_2m_in,len(ind))
                favg(vec_Tskin,iLat,iLon,ind,Tskin_in,len(ind))
                favg(vec_altitude,iLat,iLon,ind,altitude_in,len(ind))
                favg(vec_SIF_avg,iLat,iLon,ind,0.5*sif_757nm_in+1.5/2.*sif_771nm_in,len(ind))
                favg(vec_SIF_avg_norm,iLat,iLon,ind,(0.5*sif_757nm_in+1.5/2.*sif_771nm_in)*daily_corr_in,len(ind))
                
                
                print('.. averaged')
                #print(np.max(vec_n))
                #for j in range(len(ind)):

                    
                    
                    #print iT
                    
                
                

            fin.close()
        # Demand a minimum of 5 points per grid cell
        if len(ind)>0:
            vec_n[vec_n<3]=1
            wo = np.where(vec_n>1)
            #print(wo)
            vec_sif_757_1sigma[vec_n==1]=-999.9
            vec_sif_757_1sigma/=vec_n
            vec_Tskin[vec_n==1]=-999.9
            vec_Tskin/=vec_n
            vec_sza[vec_n==1]=-999.9
            vec_sza/=vec_n
            vec_sif_757nm[vec_n==1]=-999.9
            vec_sif_757nm/=vec_n
            vec_VPD_2m[vec_n==1]=-999.9
            vec_VPD_2m/=vec_n
            vec_o2_ratio[vec_n==1]=-999.9
            vec_o2_ratio/=vec_n
            vec_sif_771nm[vec_n==1]=-999.9
            vec_sif_771nm/=vec_n
            vec_daily_corr[vec_n==1]=-999.9
            vec_daily_corr/=vec_n
            vec_cont_radiance_757[vec_n==1]=-999.9
            vec_cont_radiance_757/=vec_n
            vec_cont_radiance_771[vec_n==1]=-999.9
            vec_cont_radiance_771/=vec_n
            vec_sif_757_rel[vec_n==1]=-999.9
            vec_sif_757_rel/=vec_n
            vec_co2_ratio[vec_n==1]=-999.9
            vec_co2_ratio/=vec_n
            vec_lza[vec_n==1]=-999.9
            vec_lza/=vec_n
            vec_sif_771_1sigma[vec_n==1]=-999.9
            vec_sif_771_1sigma/=vec_n
            vec_altitude[vec_n==1]=-999.9
            vec_altitude/=vec_n
            vec_sif_771_rel[vec_n==1]=-999.9
            vec_sif_771_rel/=vec_n
            vec_T2m[vec_n==1]=-999.9
            vec_T2m/=vec_n
            #print(sif_757_1sigma[counter,wo[0],wo[1]])
            for j in range(len(wo[0])):
                sif_757_1sigma[counter,wo[0][j], wo[1][j]] = vec_sif_757_1sigma[wo[0][j], wo[1][j]]
                Tskin[counter,wo[0][j],wo[1][j]] = vec_Tskin[wo[0][j],wo[1][j]]
                sza[counter,wo[0][j],wo[1][j]] = vec_sza[wo[0][j],wo[1][j]]
                sif_757nm[counter,wo[0][j],wo[1][j]] = vec_sif_757nm[wo[0][j],wo[1][j]]
                VPD_2m[counter,wo[0][j],wo[1][j]] = vec_VPD_2m[wo[0][j],wo[1][j]]
                o2_ratio[counter,wo[0][j],wo[1][j]] = vec_o2_ratio[wo[0][j],wo[1][j]]
                sif_771nm[counter,wo[0][j],wo[1][j]] = vec_sif_771nm[wo[0][j],wo[1][j]]
                daily_corr[counter,wo[0][j],wo[1][j]] = vec_daily_corr[wo[0][j],wo[1][j]]
                cont_radiance_757[counter,wo[0][j],wo[1][j]] = vec_cont_radiance_757[wo[0][j],wo[1][j]]
                cont_radiance_771[counter,wo[0][j],wo[1][j]] = vec_cont_radiance_771[wo[0][j],wo[1][j]]
                sif_757_rel[counter,wo[0][j],wo[1][j]] = vec_sif_757_rel[wo[0][j],wo[1][j]]
                co2_ratio[counter,wo[0][j],wo[1][j]] = vec_co2_ratio[wo[0][j],wo[1][j]]
                lza[counter,wo[0][j],wo[1][j]] = vec_lza[wo[0][j],wo[1][j]]
                sif_771_1sigma[counter,wo[0][j],wo[1][j]] = vec_sif_771_1sigma[wo[0][j],wo[1][j]]
                altitude[counter,wo[0][j],wo[1][j]] = vec_altitude[wo[0][j],wo[1][j]]
                sif_771_rel[counter,wo[0][j],wo[1][j]] = vec_sif_771_rel[wo[0][j],wo[1][j]]
                T2m[counter,wo[0][j],wo[1][j]] = vec_T2m[wo[0][j],wo[1][j]]
            #print(np.max(Tskin[counter,:,:]))
            #for i in dict_l2:
            #    cmd1 = 'vec_'+i+'[vec_n==1]=-999.9'
            #    cmd2 = 'vec_' + i + '/=vec_n'
                #cmd3 = i + '[:,:,:] = vec_' + i+'[:,:,:]'
            #    print(cmd1)
            #    print(cmd2)
        #    vec_SIF_avg[vec_n==1]=-999.9
        #    vec_SIF_avg_norm[vec_n==1]=-999.9
            vec_SIF_avg/=vec_n
            vec_SIF_avg_norm/=vec_n
            
            SIF_a[counter,:,:]=-999.9
            SIF_norm[counter,:,:]=-999.9
            n[counter,:,:] = vec_n[:,:]
            #for i in dict_l2:
            #  cmd = i + '[counter,:,:] = vec_' + i+'[:,:]'
            #  #print(cmd)
            #  print(cmd)
            SIF_a[counter,:,:]=vec_SIF_avg[:,:]
            SIF_norm[counter,:,:]=vec_SIF_avg_norm[:,:]
            #print('###############')
    
    
    
    
    
    f.close()
      #  print ' done'

@jit(nopython=True, parallel=True)
def favg(arr,ix,iy,iz,inp,s):
    for i in range(s):
        arr[ix[i],iy[i]]+=inp[iz[i]]
    return arr        
        
def standalone_main():
    parser = OptionParser(usage="usage: %prog l2_file2")
    parser.add_option( "-o","--outFile", dest="outFile",
                       default='L3/OCO2_SIF_monthly_Allmodes_B8100.nc',
                       help="output filename (default OCO2_SIF_map.nc)")
    parser.add_option( "--latMin", dest="latMin",
                       type=float,
                       default=-90,
                       help="min latitude region")
    parser.add_option( "--dLat", dest="dLat",
                       type=float,
                       default=1,
                       help="latitude resolution (1 degree default)")
    parser.add_option( "--dLon", dest="dLon",
                       type=float,
                       default=1,
                       help="longitude resolution (1 degree default)")
    parser.add_option( "--startTime", dest="start",
                       default='2014-09-06',
                       help="default 2014-09-06")
    parser.add_option( "--stopTime", dest="stop",
                       default='2018-08-31',
                       help="default 2015-01-01")
    parser.add_option( "--dTime", dest="dTime",
                       default=1,
                       type=int,
                       help="default 1 month (determines time window size for each time step)")
    
    parser.add_option('--mode', dest='mode', type=int, default=-1,
                      help='mode (0=ND, 1=GL, 2=TG, etc, -1 for all)')

    parser.add_option('--biome', dest='biome', type=int, default=-1,
                      help='IGBP biome type (-1 default for ALL)')

    parser.add_option( "--latMax", dest="latMax",
                       type=float,
                       default=90,
                       help="max latitude region")
    parser.add_option( "--lonMin", dest="lonMin",
                       type=float,
                       default=-180,
                       help="min longitude region")
    parser.add_option( "--lonMax", dest="lonMax",
                       type=float,
                       default=180,
                       help="max longitude region")
    parser.add_option( "--folder", dest="folder",
                       default='/export/data1/ftp/data/OCO2/sif_lite_B8100/',
                       help="Default folder DIR root where data is stored in DIR/YYYY/MM/DD/LtSIF/oco2*.nc4")
   # /data/oco2/scf/product/B7???r/r0?/'
    # Parse command line arguments
    (options, args) = parser.parse_args()
    
    # start gridding
    grid(options, args)

if __name__ == "__main__":
    standalone_main()

