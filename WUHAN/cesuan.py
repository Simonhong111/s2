import numpy as np


a = np.array([144290,174730,187700,196710])
b = np.array([256420,260384,261192,262200])
c =np.array([0,172156,173123,174105])
h = a.sum()+b.sum()+c.sum()

# 高考
gk16 = 705
gk17 = 700
gk18 = 791
gk19 = 820
# 硕士
ss16 = 59
ss17 = 72
ss18 = 76.3


rate = (ss16/gk16+ss17/gk17+ss18/gk18)/3
bs16 = 7.7
bs17 = 7.7
bs18= 9.6
brate = (bs16/gk16+bs17/gk17+bs18/gk18)/3
rs = h + h/4*rate*3+h/4*brate*4

gzrate = (1 -rate)*0.6
print(h/4*gzrate*5+ rs)
zr = h/4*gzrate*5+ rs

print("个人",1670000000/zr,"人数",rs)

# 咨询费
zfzxf = 75
zfzxs = 50

# 系统集成
zfxtjc = 300
zfxtjcs = 50

# 地面战
zfdmz = 300
zfdmzs = 10

zffy = zfzxf*zfzxs + zfxtjc*zfxtjcs + zfdmz*zfdmzs
print("政府",zffy)
# 订阅
dy = 0.2
dys = 1000

# 合作收费
hz = 0.2
hzs = 100

# 地面战
qydmz = 300
qydmzs = 50
qyf = dy*dys + qydmz*qydmzs
print("企业",qyf)

print((1500000000-217500000-152000000)/zr,285*zr)
