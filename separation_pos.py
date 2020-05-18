#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 12:54:55 2020

@author: abarklage
"""

import numpy as np
import h5py
import os

X_FILE = "../LES_outer_sword/Zf_hbCm3_dbs_hirenozz_v2_b10000/X_hbCm3_dbs_hirenozz_v2_b1_slice"
Q_FILE = "../LES_outer_sword_NPR3/Zf_hbCm3_dbs_hirenozz_v2_b10000/Qf_p10_zslice_"
t_start=1
t_max=490
boxes=[17,18]
x_duesenende=1.36666667461395264

ONPR_low          = 3
ONPR_high         = 10
delta_t          = 59.551900251890460
delta_hold       = 2.247241518939263



def readgrid(X_FILE,boxes=None):
    f = h5py.File(X_FILE,'r')
    
    grp = f['/X/']
    
    box_names = list(grp.keys())
        
    if boxes is None:
        grid=[0]*len(box_names)
        boxes=list(range(1,len(box_names)+1))
    else:
        if max(boxes) > len(box_names):
            print('Requested box-range larger than in grid-data!')
        else:
            grid=[0]*len(boxes)
    
    index=0
    for i in boxes:
        name='box' + str(i)
        box=grp[name]
        x_var=box['x']
        grid_shape=[3]
        grid_shape.extend(x_var.shape)
        grid[index]=np.zeros(grid_shape,dtype=np.float)
        x_var.read_direct(grid[index][0])
        x_var=box['y']
        x_var.read_direct(grid[index][1])
        x_var=box['z']
        x_var.read_direct(grid[index][2])
        index += 1

    return grid

non_exist_file=0

def readqdata(Q_FILE,boxes=None,qvars=None):
    f = h5py.File(Q_FILE,'r')
    time=f.attrs.__getitem__('time')
    grp = f['/Q/']
    
    box_names = list(grp.keys())
        
    if boxes is None:
        q_data=[0]*len(box_names)
        boxes=list(range(1,len(box_names)+1))
    else:
        if max(boxes) > len(box_names):
            print('Requested box-range larger than in grid-data!')
        else:
            q_data=[0]*len(boxes)
    
    if qvars is None:
        index=0
        for i in boxes:
            name='box' + str(i)
            box=grp[name]
            q_names=list(box.keys())
            q_dspace=box[q_names[0]]
            q_shape=[len(q_names)]
            q_shape.extend(q_dspace.shape)
            q_data[index]=np.zeros(q_shape,dtype=np.float)
            qvars=list(range(1, len(q_names)+1))
            for q in qvars:
                q_dspace=box[q_names[q-1]]
                q_dspace.read_direct(q_data[index][q-1])
            index += 1
        
    else:
        index=0
        for i in boxes:
            name='box' + str(i)
            box=grp[name]
            q_names=list(box.keys())
            q_dspace=box[q_names[0]]
            q_shape=[len(qvars)]
            q_shape.extend(q_dspace.shape)
            q_data[index]=np.zeros(q_shape,dtype=np.float)
            j=0
            for q in qvars:
                q_dspace=box[q_names[q-1]]
                q_dspace.read_direct(q_data[index][j])
                j += 1
            index += 1
     
        
    return q_data,time

def progressBar(current, total, barLength = 20):
    percent = float(current) * 100 / total
    arrow   = '-' * int(percent/100 * barLength - 1) + '>'
    spaces  = ' ' * (barLength - len(arrow))

    print('Progress: [%s%s] %d %%' % (arrow, spaces, percent), end='\r')
    
#def getNPR(time,NPR_)

'Reading X-file'

x=readgrid(X_FILE,boxes)

'Reading grid'

x_abl=np.zeros(t_max-t_start+1)
for t in range(t_start,t_max+1):
    Q_tmp=Q_FILE+str(t)
    if os.path.exists(Q_tmp):
        [q,time]=readqdata(Q_tmp,boxes,qvars=[2])
        progressBar(t,t_max+1)
        
        'Ende der Düse'
        for box in range(len(boxes)):
            
            kmax=x[box].shape[1]
            for k in range(kmax):
                i=0
                while x[box][0,k,1,i]<x_duesenende and q[box][0,k,1,i]>0:
                    i += 1
                x_abl[t-t_start] += x[box][0,0,1,i]
                
    else:
        x_abl[t-t_start]=x_abl[t-t_start-1]
        print('File of timestep %i does not exist' % t)
        non_exist_file += 1
                

x_abl = x_abl/2

# grouid=h5py.h5g.open(f, GROUP)


# dsetid = h5py.h5d.open(f, DATASET)

# newdata=0

# dsetid.read(h5py.h5s.ALL,h5py.h5s.ALL, newdata)
