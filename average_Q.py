#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 12:05:55 2020

@author: niibarkl
"""

import numpy as np
import tfs_io
import calc_vars
import os

out_file="../LES_outer_sword/Zf_hbCm3_dbs_hirenozz_v2_b10000/Qf_p10_zslice_av"
Q_FILE = "../LES_outer_sword/Zf_hbCm3_dbs_hirenozz_v2_b10000/Qf_p10_zslice_"#Zf_hbCm3_dbs_hirenozz_v2_b10000/Q_hbCm3_dbs_hirenozz_v2_b1_"
Q_FILE_ini = "../LES_outer_sword/Qf_ini_merged"
t_start=300
t_max=850

q_av=[]
count=0
for t in range(t_start,t_max+1):
    Q_tmp=Q_FILE+str(t)
    try:
        [q,time]=tfs_io.readqdata(Q_tmp)
        if t==t_start:
            q_av=q
        else:
            for num,b in enumerate(q_av):
                q_av[num]=q_av[num]+q[num]
        count+=1
    except OSError:
        print('Gibts nicht!')
    
for b in range(len(q_av)):
    q_av[b]=q_av[b]/count
#q_av=q_av/(t_max-t_start+1)

tfs_io.writeqdata(out_file,q_av)