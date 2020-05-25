#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 14:08:44 2020

@author: niibarkl
"""

import numpy as np
import h5py
import os

def loaddat(DAT_DIR,file_spec=None):
    f=[]
    for root, dirs, files in os.walk(DAT_DIR):
        f.append(files)
    
#    files=files[1:3]
    if file_spec==None:
        files=files
    else:
        files=file_spec
    
    sensor_data=[0]*len(files)
    time=[0]*len(files)
    
    i=0
    for file in files:
        DAT_FILE=root+file
        sensor_data[i] = np.transpose(np.genfromtxt(DAT_FILE, dtype=float))
        time[i]=sensor_data[i][0,:]
        sensor_data[i]=sensor_data[i][1:,:]
        i += 1
    return sensor_data,time,files

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
    
#    if qvars is None:
    index=0
    for i in boxes:
        name='box'+str(i)
        box=grp[name]
        q_names=list(box.keys())
        q_dspace=box[q_names[0]]
        q_shape=[len(q_names)]
        q_shape.extend(q_dspace.shape)
        q_data[index]=np.zeros(q_shape,dtype=np.float)
        if qvars==None :
            qvars=list(range(1, len(q_names)+1))
        for q in qvars:
            q_dspace=box[q_names[q-1]]
            q_dspace.read_direct(q_data[index][q-1])
        index += 1

    return q_data,time

def getPosition(x,xq,yq,zq):
    ijk_q=[0]*3
    blocks=len(x)
    maxCoord=np.zeros((3,blocks));
    minCoord=np.zeros((3,blocks));
    index=0
    for block in x:
        maxCoord[0,index]=np.amax(block[0,:])
        maxCoord[1,index]=np.amax(block[1,:])
        maxCoord[2,index]=np.amax(block[2,:])
        minCoord[0,index]=np.amin(block[0,:])
        minCoord[1,index]=np.amin(block[1,:])
        minCoord[2,index]=np.amin(block[2,:])
        index += 1
        
    for block in range(len(x)):
        if (xq >= minCoord[0,block]) and (xq <= maxCoord[0,block]) \
        and (yq >= minCoord[1,block]) and (yq <= maxCoord[1,block]) \
        and (zq >= minCoord[2,block]) and (zq <= maxCoord[2,block]):
            blockq=block
    iiter=x[0].shape[3]
    jiter=x[0].shape[2]
    kiter=x[0].shape[1]
    
    dist_test=1
    for i in range(iiter):
        for j in range(jiter):
            for k in range(kiter):
                dist=np.abs((x[blockq][0,k,j,i]-xq))+np.abs((x[blockq][1,k,j,i]-yq)) \
                +np.abs((x[blockq][2,k,j,i]-zq))
                if dist<dist_test:
                    dist_test=dist
                    ijk_q[0]=k  
                    ijk_q[1]=j  
                    ijk_q[2]=i  
    
    return  ijk_q, blockq

def writeqdata(OUT_File,q_data,boxes=None,qvars=None,time=None):
    f = h5py.File(OUT_File,'w-')
    f.attrs.create('time',time)
    grp = f.create_group("/Q/")
        
    if boxes is None:
        boxes=list(range(1,len(q_data)+1))
        box_names=['box'+str(box) for box in boxes]
    else:
        if max(boxes) > len(box_names):
            print('Requested box-range larger than in grid-data!')
        else:
            box_names=['box'+str(box) for box in boxes]
    
    #    if qvars is None:
    for i in boxes:
        name=box_names[i-1]
        box=grp.create_group(name)
        if len(q_data[i-1].shape[:])>3:
            q_shape=q_data[i-1].shape[1:]
        else:
            q_shape=q_data[i-1].shape[:]
        if qvars==None :
            qvars=list(range(1, q_data[i-1].shape[0]+1))
            q_names=['Q'+str(q) for q in qvars]
        else:
            q_names=['Q'+str(q) for q in qvars]
        for q in qvars:
            q_dspace=box.create_dataset(q_names[q-1],q_shape,dtype='float64')
            q_dspace.write_direct(q_data[i-1][q-1])
    return None