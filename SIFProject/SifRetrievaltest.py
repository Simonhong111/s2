from SifRetrieval import *

dataset,HeadTitle,White,Dark = fileload(r"D:\HUBEIFIELDDATA\20180830airline6-10\15350963458812982835431128394280_20180830102641.08.xlsx")




mdataset = [da[30:] for da in dataset]
dataset = np.array(mdataset)

HT = np.array(HeadTitle[30:])
White = np.array(White[30:])
Dark = np.array(Dark[30:])
mDark = smooth(Dark,1)
mdataset = smooth(dataset[0],1) - mDark
mWhite = (smooth(White,1) - mDark)/0.73



w = (HT >750) & (HT <770)
w = np.where(w)


fig = plt.figure()
# plt.plot(HT,mWhite)
# plt.plot(HT,mdataset)
plt.subplot(121)
plt.plot(HT[w],mdataset[w]/mWhite[w])
plt.subplot(122)
plt.plot(HT[w],mWhite[w])

plt.show()


