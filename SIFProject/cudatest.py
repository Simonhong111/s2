from pycuda.curandom import rand as curand
a_gpu = curand((50,))
b_gpu = curand((50,))
print(a_gpu)
