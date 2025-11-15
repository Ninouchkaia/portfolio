#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

with open("colnames_annotated_2023.csv", 'r') as file_read :
    data = file_read.readlines()
    print(data)
    annot_list = []
    for line in data[1:]:
        line = line.split(";")
        annot = line[1].replace("\n","")
        annot_list.append(annot)
    print(annot_list)

    print(len(annot_list))

annot_list = list(set(annot_list))
print(len(annot_list))


stop = timeit.default_timer()
print(stop - start)  