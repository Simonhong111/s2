from osgeo import gdal,osr,ogr
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def extractPointFromFeature(feature):
    feature = str(feature)
    polyNum = len(feature.split("))"))
    print("feature number is ",polyNum)
    feature = feature.split("((")[1].split("))")[0].split(",")
    points = [[float(f.split(" ")[1]),float(f.split(" ")[0])] for f in feature]
    return points

def is_collinear(points):
    for id in range(len(points)):
        id0 = id%len(points)
        id1 = (id+1)%len(points)
        id2 =(id+2)%len(points)

        v1 = [points[id1][0] -points[id0][0] , points[id1][1] -points[id0][1]]
        v2 = [points[id2][0] -points[id1][0] , points[id2][1] -points[id1][1]]
        if np.sqrt(v1[0]*v1[0] + v1[1]*v1[1]) < 0.01:
            print(id, points[id0], points[id1], points[id2])
            print("*1")
            return True
        if  np.sqrt(v2[0]*v2[0] + v2[1]*v2[1])< 0.01:
            print(id,points[id0],points[id1],points[id2])
            print("*2")
            return True
        if v1[0]*v2[1] == v1[1]*v2[0]:
            if points[id0][0] == points[id1][0] and points[id0][1] < points[id1][1]:
                if points[id2][1] < points[id1][1]:
                    print("*3")
                    return True
            if points[id0][0] == points[id1][0] and points[id0][1] > points[id1][1]:
                if points[id2][1] > points[id1][1]:
                    print("*4")
                    return True
            if points[id0][0] != points[id1][0] and points[id0][0] < points[id1][0]:
                if points[id2][0] < points[id1][0]:
                    print("*5")
                    return True
            if points[id0][0] != points[id1][0] and points[id0][0] > points[id1][0]:
                if points[id2][0] > points[id1][0]:
                    print("*6")
                    return True
    return False


if __name__ == '__main__':

    path = r'C:\Users\zmhwh\Desktop\Temp\s2latlon2.txt'

    # read txt method two
    f = open(path)
    coord = []
    for line2 in f:
        coord.append(line2)
    line = coord[0]
    line = line.strip().split(" ")
    # print(line)
    coord = []
    for term in line:
        term = term.split(",")
        coord.append([float(term[0]), float(term[1])])

    # coord.append(coord[0])
    points = coord
    # print(coord)


    # assert is_collinear(points[:-1]) == False, "is collinear"
    # print(points[223])
    # print(points[224])
    # print(points[0])
    Y,X = [float(p[0]) for p in points], [float(p[1]) for p in points]
    name_dict = {
        'X': X,
        'Y': Y
    }

    df = pd.DataFrame(name_dict)
    #
    # df.to_csv(r"C:\Users\zmhwh\Desktop\Temp\mm3.csv")

    plt.plot(X, Y)
    plt.plot(X[0:150], Y[0:150])
    for x in X:
        print(x)

    plt.show()


