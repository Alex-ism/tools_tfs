#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 12:54:55 2020

@author: abarklage
"""

import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
import calc_vars
import tfs_io

X_FILE = "../LES_outer_sword/Zf_hbCm3_dbs_hirenozz_v2_b10000/X_hbCm3_dbs_hirenozz_v2_b1_slice"
Q_FILE = "../LES_outer_sword_NPR3/Zf_hbCm3_dbs_hirenozz_v2_b10000/Qf_p10_zslice_"
t_start=1
t_max=1250
boxes=[17,18]
x_duesenende=1.36666667461395264


turbulent_outerflow={"NPR_low":3,"NPR_high":10, "delta_t":59.551900251890460,"delta_hold":2.247241518939263}



def progressBar(current, total, barLength = 20):
    percent = float(current) * 100 / total
    arrow   = '-' * int(percent/100 * barLength - 1) + '>'
    spaces  = ' ' * (barLength - len(arrow))

    print('Progress: [%s%s] %d %%' % (arrow, spaces, percent), end='\r')
    

'Reading X-file'

x=tfs_io.readgrid(X_FILE,boxes)
x_lip=tfs_io.readgrid(X_FILE,[13,14])


non_exist_file=[]
x_abl=np.zeros(t_max-t_start+1)
p_lip=np.zeros(t_max-t_start+1)
NPR=np.zeros(t_max-t_start+1)
for t in range(t_start,t_max+1):
    Q_tmp=Q_FILE+str(t)
    if os.path.exists(Q_tmp):
        [q,time]=tfs_io.readqdata(Q_tmp,boxes,qvars=[2])
        [q_lip,time]=tfs_io.readqdata(Q_tmp,[13,14])
        NPR[t-t_start]=calc_vars.getNPR(time,turbulent_outerflow)
        progressBar(t,t_max+1)
        
        'Duesenlippe'
        [ijk_lip,block_lip]=tfs_io.getPosition(x_lip,1.366667,0.280679792165756226,0)
        p_lip[t-t_start]=calc_vars.calcpressure(q_lip,block=block_lip,ijk=ijk_lip)
        for box in range(len(boxes)):
            
            kmax=x[box].shape[1]
            for k in range(kmax):
                i=0
                while x[box][0,k,1,i]<x_duesenende and q[box][1,k,1,i]>0:
                    i += 1
                x_abl[t-t_start] += x[box][0,0,1,i]
                
    else:
        x_abl[t-t_start]=x_abl[t-t_start-1]
        NPR[t-t_start]=NPR[t-t_start-1]
        p_lip[t-t_start]=p_lip[t-t_start-1]
        print('File of timestep %i does not exist' % t)
        non_exist_file.append(t)
                

x_abl = x_abl/2
plt.plot(NPR,x_abl)

# grouid=h5py.h5g.open(f, GROUP)


# dsetid = h5py.h5d.open(f, DATASET)

# newdata=0

# dsetid.read(h5py.h5s.ALL,h5py.h5s.ALL, newdata)
