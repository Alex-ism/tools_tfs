#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 14:54:12 2020

@author: niibarkl
"""
import numpy as np
import os
import scipy
import math

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
indicator=True

for d in range(len(data)):
    number_samples=data[d].shape[0]
    time=data[d][:,0]
    istart[d]=[]
    iend[d]=[0]
    for i in range(1,number_samples):
        if time[i-1] > time[i] and indicator:
            t_break=time[i-1]
            istart[d].append(i-1)
            indicator=False
    
        elif time[i-1] < t_break < time[i+1] and indicator == False:
            iend[d].append(i-1)
            indicator=True
        
'Zuschneiden Signal'
new_data=[0]*len(data)

for d in range(len(data)):
    new_data[d]=np.zeros((1,6),dtype=np.float)
    for i in range(len(istart[d])):
        new_data[d]=np.append(new_data[d],data[d][iend[d][i]:istart[d][i],:],axis=0)


'Schreiben des neuen Files'

#for d in range(len(data)):
#    
#    filename='../LES_outer_sword/sens_cut/'+str(files[d])
#    np.savetxt(filename, new_data[d])
        
        
        
        
        
