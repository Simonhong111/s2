from matplotlib import pyplot as plt
import numpy as np
import s2sphere as s2

points = [[116.9997801738,35.2430777118],[118.2065067255,35.2370638902],
          [118.192198189,34.2471243765],[116.9997827816,34.2529215885] ,
         [116.9997801738,35.2430777118]]
# points2 = [
# [116.999780173827, 35.2430777118033],
#     [116.999782781608, 34.2529215885203],[118.192198189045, 34.247124376463]
#     ,[118.206506725524, 35.2370638901883],[116.999780173827, 35.2430777118033]
# ]



import s2sphere

r = s2sphere.RegionCoverer()
p1 = s2sphere.LatLng.from_degrees(34.2471243765,118.192198189)
p2 = s2sphere.LatLng.from_degrees(35.2430777118,116.9997801738)
r.min_level = 8
r.max_level = 8
r.max_cells =100
print("ma",r.max_cells)
cell_ids = r.get_covering(s2sphere.LatLngRect.from_point_pair(p1, p2))
# print(cell_ids)
points2 = []
for i in cell_ids:
    print(i.face())
    p0,p1 = s2sphere.CellId(i.id()).to_lat_lng().lat().degrees,s2sphere.CellId(i.id()).to_lat_lng().lng().degrees
    print(p0,p1)
    points2.append([float(p1),float(p0)])
    # print("*",s2sphere.CellId(i.id()))


cellids =[
3874239171631513600,
3874274356003602432,
3874872490329112576,
3874942859073290240,
3874978043445379072,
3875013227817467904,
3875048412189556736,
3875083596561645568,
3875118780933734400,
3875259518422089728]
points3 = []

for i in cellids:

    p0,p1 = s2sphere.CellId(i).to_lat_lng().lat().degrees,s2sphere.CellId(i).to_lat_lng().lng().degrees
    print(p0,p1)
    points3.append([float(p1),float(p0)])
    # print("*",s2sphere.CellId(i.id()))



points = np.array(points)
x = points[:,0]
y = points[:,1]

points2 = np.array(points2)
x2 = points2[:,0]
y2 = points2[:,1]
points3 = np.array(points3)
x3 = points3[:,0]
y3 = points3[:,1]

plt.plot(x,y)
# plt.plot(x2,y2,'r*')
# plt.scatter([116.9997801738],[35.2430777118],c='g')
# plt.scatter([118.192198189],[34.2471243765],c='g')
# plt.scatter([117.36164757079143],[34.57137699422246],c='k',s=100)
plt.scatter(x3,y3,c='y')

plt.show()




