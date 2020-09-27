import numpy as np

N = 10000
a = np.random.normal(size=N)
b = np.random.normal(size=N)
c = np.random.normal(size=N)


V = a + 23
mask = np.where(V >=0)
valid = mask[0][0:10000]
a[valid]
b[valid]
c[valid]
