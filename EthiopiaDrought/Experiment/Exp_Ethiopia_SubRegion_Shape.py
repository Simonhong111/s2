import shapefile as shp
def boundary(geo):
    sp = r"D:\Cornell\EthiopianDrought\ET_Region\et_region.shp"
    sf = shp.Reader(sp)
    X = []
    Y = []
    for shape in sf.shapeRecords():
        for i in range(len(shape.shape.parts)):
            i_start = shape.shape.parts[i]
            if i == len(shape.shape.parts) - 1:
                i_end = len(shape.shape.points)
            else:
                i_end = shape.shape.parts[i + 1]
            x = [(i[0]-geo[0])/geo[1] for i in shape.shape.points[i_start:i_end]]
            y = [(i[1]-geo[3])/geo[5] for i in shape.shape.points[i_start:i_end]]
            X.append(x)
            Y.append(y)
    return X,Y