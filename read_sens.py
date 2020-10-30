#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 14:54:12 2020

@author: niibarkl
"""
import numpy as np
import os

DAT_DIR="../LES_outer_sword_NPR3/sensors_dat/"


def loaddat(DAT_DIR):
    f=[]
    for root, dirs, files in os.walk(DAT_DIR):
        f.append(files)
    
    files=files[1:3]
    
    sensor_data=[0]*len(files)
    
    i=0
    for file in files:
        DAT_FILE=root+file
        sensor_data[i] = np.genfromtxt(DAT_FILE, dtype=float)
        i += 1
    
    return sensor_data,files

[data,files]=loaddat(DAT_DIR)

for d in range(len(data)):
    i=0
    while i < len(data[d]):
        if np.isnan(data[d][i,0])==True:
            data[d] = np.delete(data[d], i, axis=0)
        i+=1
        
time=data[0][:,0]   
nan_elements=np.isnan(time)
istart=[0]*len(data)
iend=[0]*len(data)
t_break=0

for d in range(len(data)):
    number_samples=data[d].shape[0]
    time=data[d][:,0]
    istart[d]=[]
    iend[d]=[-1]
    indicator=True
    for i in range(1,number_samples):
        if (time[-i-1] > time[-i] or (time[-i]-time[-i-1])>0.2) and indicator:
            t_break=time[-i+number_samples]
            istart[d].append(-i+number_samples)
            indicator=False
    
        elif time[-i-1] < t_break < time[-i+1] and np.abs(time[-i-1]-time[-i+1])<0.2 and indicator == False:
            iend[d].append(-i+number_samples-1)
            indicator=True
    if len(istart[d])<len(iend[d]):
        istart[d].append(0)
        
'Zuschneiden Signal'
new_data=[0]*len(data)

for d in range(len(data)):
    new_data[d]=np.zeros((1,6),dtype=np.float)
    for i in range(len(istart[d])):
        new_data[d]=np.append(new_data[d],data[d][istart[d][-i-1]:iend[d][-i-1],:],axis=0)

for d in range(len(data)):
    new_data[d] = np.delete(new_data[d], 0, axis=0)
'Schreiben der neuen Files'

#for d in range(len(data)):
    
#    filename='../LES_outer_sword/sens_cut/'+str(files[d])
#    np.savetxt(filename, new_data[d])
#        
        
        
        
        
