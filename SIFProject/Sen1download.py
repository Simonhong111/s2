from sentinelsat import SentinelAPI, read_geojson, geojson_to_wkt

api = SentinelAPI('', '',"https://scihub.copernicus.eu/dhus")
footprint = geojson_to_wkt(read_geojson(r'D:\Sentinel-2Data\Hubei.json'))
print(footprint)

products = api.query(footprint,producttype='SLC',date=("20180601","20180615"),orbitdirection='ASCENDING',platformname = 'Sentinel-1')
print(len(products))

for prod in products:
    print(api.get_product_odata(prod))




# api.download_all(products,directory_path="D:\Sentinel-2Data")
# api.download_all(products)