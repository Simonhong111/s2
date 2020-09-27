import argparse
import datetime
import glob
import numpy as np
import os
import re
import time
import sentinelsat
import zipfile
import time
import pdb
import pandas
from sentinelsat import SentinelAPI, read_geojson, geojson_to_wkt
from pandas import DataFrame
import pandas as pd
import hashlib

#用于获取下载哨兵二指定tile或者指定矢量的哨兵影像所需的URL地址、md5校验码及长度
#用于检验那些未下载完的链接

def getSentinelUrlByTiles():
    scihub_api = SentinelAPI('hust_hw_316', '12345678', 'https://scihub.copernicus.eu/dhus')
    startdate = sentinelsat.format_query_date('20190101')
    enddate = sentinelsat.format_query_date('20191231')
    # Search data, filtering by options.
    level = '2A'
    # tile = '50RLV'
    specify_tiles = ['49RGQ', '49RGP', '50RKU', '50RKV', '50RLV']  # 指定需要下载的tiles

    final_urls = DataFrame()
    url_list = []
    md5_list = []
    for tile in specify_tiles:
        products = scihub_api.query(beginposition=(startdate, enddate),
                                    platformname='Sentinel-2',
                                    producttype='S2MSI%s' % level,
                                    cloudcoverpercentage=(0, 100), filename='*T%s*' % tile)

        # Your IP has been temporarily locked
        products_df = scihub_api.to_dataframe(products)
        # for i in range(len(products_df)):
        # scihub_api.get_product_odata()
        # print(products_df['link'])
        # print(products_df['size'])

        # {'id': 'd9c11a08-9548-423a-9d5b-3790bdd1949b', 'title': 'S2B_MSIL1C_20180322T031539_N0206_R118_T49SCS_20180322T080731', 'size': 557305492, 'md5': 'A4DEE331BABF5D81FEEFED7E8F152E1A', 'date': datetime.datetime(2018, 3, 22, 3, 15, 39, 27000), 'footprint': 'POLYGON((109.12450894645511 32.434247217341294,109.15551946020656 32.55886040012398,109.19500574164893 32.7064511260414,109.23422411890954 32.85418162235375,109.27416742227145 33.00170389409306,109.31402177718599 33.149262433169554,109.35311404383268 33.297123960618194,109.38830946718173 33.42770696377061,110.02964494156487 33.435778798872015,110.04040193865697 32.44547285973725,109.12450894645511 32.434247217341294))', 'url': "https://scihub.copernicus.eu/dhus/odata/v1/Products('d9c11a08-9548-423a-9d5b-3790bdd1949b')/$value", 'Online': False, 'Creation Date': datetime.datetime(2018, 3, 22, 11, 36, 6, 645000), 'Ingestion Date': datetime.datetime(2018, 3, 22, 11, 25, 45, 119000)}
        # 第二种，通过矢量查询数据
        for i in range(len(products_df)):
            product = scihub_api.get_product_odata(products_df.iloc[i]['uuid'])
            print(product['url'])
            url_list.append(product['url'])
            print(product['md5'])
            md5_list.append((product['md5']))

    final_urls['url'] = url_list
    final_urls['md5'] = md5_list

    final_urls.to_csv("finale_urls_md5.csv", index=False)
    # 这里的md5一定要记录下来
    # 下面是以某个geojosn
    # footprint = geojson_to_wkt(read_geojson('test.geojson'))
    # products = scihub_api.query(footprint,
    #                     beginposition=(startdate, enddate),
    #                     platformname='Sentinel-2')

    # convert to Pandas DataFrame
    # products_df = scihub_api.to_dataframe(products)
    # print("............")
    # print(products_df)
    # pass


def printPath(level, path):
    allFileNum = 0
    ''''' 
    ´òÓ¡Ò»¸öÄ¿Â¼ÏÂµÄËùÓÐÎÄ¼þ¼ÐºÍÎÄ¼þ 
    '''
    # ËùÓÐÎÄ¼þ¼Ð£¬µÚÒ»¸ö×Ö¶ÎÊÇ´ÎÄ¿Â¼µÄ¼¶±ð
    dirList = []
    # ËùÓÐÎÄ¼þ
    fileList = []
    # ·µ»ØÒ»¸öÁÐ±í£¬ÆäÖÐ°üº¬ÔÚÄ¿Â¼ÌõÄ¿µÄÃû³Æ(google·­Òë)
    files = os.listdir(path)
    # ÏÈÌí¼ÓÄ¿Â¼¼¶±ð
    dirList.append(str(level))
    for f in files:
        if (os.path.isdir(path + '/' + f)):
            # ÅÅ³ýÒþ²ØÎÄ¼þ¼Ð¡£ÒòÎªÒþ²ØÎÄ¼þ¼Ð¹ý¶à
            if (f[0] == '.'):
                pass
            else:
                # Ìí¼Ó·ÇÒþ²ØÎÄ¼þ¼Ð
                dirList.append(f)
        if (os.path.isfile(path + '/' + f)):
            # Ìí¼ÓÎÄ¼þ
            fileList.append(f)
            # µ±Ò»¸ö±êÖ¾Ê¹ÓÃ£¬ÎÄ¼þ¼ÐÁÐ±íµÚÒ»¸ö¼¶±ð²»´òÓ¡
    i_dl = 0
    for dl in dirList:
        if (i_dl == 0):
            i_dl = i_dl + 1
        else:
            # ´òÓ¡ÖÁ¿ØÖÆÌ¨£¬²»ÊÇµÚÒ»¸öµÄÄ¿Â¼
            # ´òÓ¡Ä¿Â¼ÏÂµÄËùÓÐÎÄ¼þ¼ÐºÍÎÄ¼þ£¬Ä¿Â¼¼¶±ð+1
            printPath((int(dirList[0]) + 1), path + '/' + dl)
    final_list=[]
    for fl in fileList:
        # ´òÓ¡ÎÄ¼þ
        # Ëæ±ã¼ÆËãÒ»ÏÂÓÐ¶àÉÙ¸öÎÄ¼þ
        allFileNum = allFileNum + 1
        final_list.append(os.path.join(path, fl))
    return final_list



def checkwhichURL_not_done(tmp_filelist, urls_list_df):
    not_download_url=[] #最后未被下载的url

    df = urls_list_df.copy()
    df['flag']=0
    df = df.copy()

    for file in tmp_filelist:
        md5_hash = hashlib.md5()
        with open(file, "rb") as f:
            # Read and update hash in chunks of 4K
            for byte_block in iter(lambda: f.read(4096*1024*10), b""):
                md5_hash.update(byte_block)
            md5=str.upper(md5_hash.hexdigest())
            print(".................")
            #print(md5)
            mask = (df['md5'] == md5)
            df.loc[mask, 'flag'] = 1

    #df=df.fillna(0).copy()
    #最后找到所有那些flag=0.0的
    mask = (df['flag'] ==0)
    df_extra = df.loc[mask]
    df_extra = df_extra.copy()
    print(df_extra['url'])
    return df_extra







if __name__ == "__main__":
    strPath=r"D:\aria2\hongboshi"
    tmp_filelist=printPath(1,strPath)
    #print(tmp_filelist)
    urls_list_df = pd.read_csv(r"D:\aria2\finale_urls_md5.csv")
    df = checkwhichURL_not_done(tmp_filelist,urls_list_df )
     

    df.to_csv(r"D:\aria2\undone_usrl.csv")

    pass

