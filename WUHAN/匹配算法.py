import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

MapX =np.random.randint(0,10,size=[50,50]).flatten()
MapY =np.random.randint(0,10,size=[50,50]).flatten()
T = np.array(range(MapY.shape[0]))


fig=plt.figure()
ax1 = Axes3D(fig)

# ax1.scatter3D(MapX,MapY,T, cmap='Blues')  #绘制散点图
ax1.plot3D(MapX,MapY,T,'gray')    #绘制空间曲线
plt.show()
