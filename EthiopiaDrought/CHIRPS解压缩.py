import gzip
import os
import glob

def un_gz(file_name):
    # 获取文件的名称，去掉后缀名
    f_name = file_name.replace(".gz", "")
    # 开始解压
    g_file = gzip.GzipFile(file_name)
    # 读取解压后的文件，并写入去掉后缀名的同名文件（即得到解压后的文件）
    open(f_name, "wb+").write(g_file.read())
    g_file.close()


filenames = glob.glob(r"D:\Cornell\EthiopianDrought\CHIRPS5\*gz")
for file_name in filenames:
    print(file_name)
    un_gz(file_name)
