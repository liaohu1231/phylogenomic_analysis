#!/usr/bin/python
#coding:utf-8

import os

f1 = open('novel_order.csv','r')

name = f1.readlines()

path = "/hwfssz5/ST_CANCER/CGR/USER/liaohu/data/genome_4000/novel_order/"

files = os.listdir(path) #得到文件夹下的所有文件名称

for file in files:
    print(file)
    if file.split(".")[-1] == "fna":
        file_name = file.split('.')[0]
        print (file_name)
        for key in name:
            keyname = key.split(',')[5]
            keyname2=keyname.split('.')[0]
            print (keyname)
            taxon_id = key.split(',')[0]
            if file_name == keyname2:
                try:
                    os.rename(path+file,path+str(taxon_id)+'.fna')
                except OSError:
                    msg = " sorry, the file "+ file +" dose not exist."
                    #print(msg)
