from epdatapro20082015v2 import *
from spatialmode import  *
from matplotlib import pyplot as plt
from matplotlib import cm
months =["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]


crop_path = r"D:\Cornell\EthiopianDrought\CropType2015\agg_clip60.tif"
crop = gdal.Open(crop_path)
geo_t = crop.GetGeoTransform()
Width,Height = crop.RasterXSize,crop.RasterYSize
crop_raster = crop.ReadAsArray()
mask = np.where(crop_raster == 255)
crop_raster[crop_raster != 255] = np.nan
# print(mask)
vRow = mask[0]
vCol = mask[1]
ind = vRow*Width + vCol

pvipath =r"D:\Cornell\EthiopianDrought\AData\PVIDaily"



# 计算矢量边界

daShapefile = r"D:\Cornell\EthiopianDrought\ETH_outline_SHP\ETH_outline.shp"

driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(daShapefile, 0)
layer = dataSource.GetLayer()
feature = layer.GetFeature(0)
geo = feature.GetGeometryRef()
geo = str(geo).split("((")[1].split("))")[0].split(",")
x = []
y = []
for term in geo:
    x.append(float(term.split(" ")[0]))
    y.append(float(term.split(" ")[1]))

x = np.array(x)
y = np.array(y)
x = (x - geo_t[0]) / geo_t[1]
y = (y - geo_t[3]) / geo_t[5]

plt.imshow(crop_raster)
# plt.colorbar()
plt.plot(x,y)
plt.show()

for year in range(2003,2019):
    yy = str(year)
    short_pvi = gdal.Open(os.path.join(pvipath, "short_pvi_{}.tif".format(yy))).ReadAsArray() * (1.0)
    long_pvi = gdal.Open(os.path.join(pvipath, "long_pvi_{}.tif".format(yy))).ReadAsArray() * (1.0)
    short_pvi_list = np.take(short_pvi,ind)
    long_pvi_list = np.take(long_pvi,ind)
    # print("short min max",short_pvi_list.min(),short_pvi_list.max())
    print("long min max", long_pvi_list.min(), long_pvi_list.max())


    fig = plt.figure(figsize=(8, 3))
    plt.title("{}  PVI Map ".format(yy) + '\n', fontsize=16)
    plt.xticks([])
    plt.yticks([])

    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_title("Short Rains PVI Map")
    mask1 = np.where(short_pvi > -9999)
    short_pvi[short_pvi == -9999] = np.nan
    vmin = short_pvi[mask1].min()
    vmax = short_pvi[mask1].max()
    # print("short maxmin value", vmin, vmax)

    cax1 = ax1.imshow(short_pvi, cmap=plt.get_cmap("rainbow"), vmin=0, vmax=1)
    cbar1 = plt.colorbar(cax1, ax=ax1, fraction=0.036, pad=0.04)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.imshow(crop_raster)
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_title("long Rains PVI Map")
    mask2 = np.where(long_pvi > -9999)
    long_pvi[long_pvi == -9999] = np.nan
    vmin = long_pvi[mask2].min()
    vmax = long_pvi[mask2].max()
    # print("maxmin value", vmin, vmax)

    cax2 = ax2.imshow(long_pvi, cmap=plt.get_cmap("rainbow"), vmin=0, vmax=1.5)
    cbar2 = plt.colorbar(cax2, ax=ax2, fraction=0.036, pad=0.04)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.imshow(crop_raster)
    fig.tight_layout()  # 调整整体空白
    path2 = os.path.join(r"D:\Cornell\EthiopianDrought\CropCSV\PVIMap", yy + "PVI60.jpg")
    plt.savefig(path2)
    plt.close()
    # plt.show()


