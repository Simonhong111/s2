import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

path  = r'C:\Users\zmhwh\Desktop\Temp\s2latlon2.txt'
f = open(path)
coord = []
for line2 in f:
    coord.append(line2)
line = coord[0]
line = line.strip().split(" ")
coord = []
for term in line:
    term = term.split(",")
    coord.append([float(term[0]),float(term[1])])

X = [p[1] for p in coord]
Y = [p[0] for p in coord]
# print(np.array(X).max(),np.array(X).min())
# print(np.array(Y).max(),np.array(Y).min())
# print(coord)
# plt.plot(X,Y)

def data_gen():
    for idx,item in enumerate(X):

        yield term, Y[idx]
fig, ax = plt.subplots()
ax.set_xlim(-44.9904, -4.28256)
ax.set_ylim(-31.8188, 20.3993)
xdata, ydata = [], []
ln, = plt.plot([], [], 'r')


def init():

    return ln,

def update(i):


    xdata.append(X[i])
    ydata.append(Y[i])

    ln.set_data(xdata, ydata)
    return ln,

ani = FuncAnimation(fig, update,
                    init_func=init,repeat=False, blit=True)

plt.show()