#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 13:46:25 2020

@author: niibarkl
"""

import numpy as np

def calcpressure(q,block=None,ijk=None,Ma=3,kappa=1.4):
    if ijk==None:
        p=[0]*len(q)
        for block in range(len(q)):
            p[block]=q[block][4,:]-0.5*(q[block][1,:]**2+q[block][2,:]**2 \
                    +q[block][3,:]**2)/q[block][0,:]
        if Ma!=None:
            for block in range(len(q)):
                p[block][:]=p[block][:]*kappa*(kappa-1)*(1+(kappa-1)/2*Ma**2)**(kappa/(kappa-1))      
    else:
        p=q[block][4,ijk[0],ijk[1],ijk[2]]-0.5*(q[block][1,ijk[0],ijk[1],ijk[2]]**2\
           +q[block][2,ijk[0],ijk[1],ijk[2]]**2 \
           +q[block][3,ijk[0],ijk[1],ijk[2]]**2)/q[block][0,ijk[0],ijk[1],ijk[2]]
        if Ma!=None:
            p=p*kappa*(kappa-1)*(1+(kappa-1)/2*Ma**2)**(kappa/(kappa-1))  

    
    return p

def getArea(slice_ijk):
    dir1_max=slice_ijk.shape[1]
    dir2_max=slice_ijk.shape[2]
    A=np.zeros([dir1_max,dir2_max],dtype=np.float)
    n=np.zeros(slice_ijk.shape,dtype=np.float)
    for i in range(dir1_max):
        for j in range(dir2_max):
            ijk_pos=slice_ijk[:,i,j]
            if i!=dir1_max-1:
                vec1=ijk_pos+0.5*(slice_ijk[:,i+1,j]-ijk_pos)
            else:
                vec1=np.zeros(3)
                
            if i!=0:
                vec2=ijk_pos+0.5*(slice_ijk[:,i-1,j]-ijk_pos)
            else:
                vec2=np.zeros(3)
                
            if j!=dir2_max-1:
                vec3=ijk_pos+0.5*(slice_ijk[:,i,j+1]-ijk_pos)
            else:
                vec3=np.zeros(3)
                
            if j!=0:
                vec4=ijk_pos+0.5*(slice_ijk[:,i,j-1]-ijk_pos)
            else:
                vec4=np.zeros(3)
                    
            crossprod1=np.cross(vec1,vec3)
            crossprod2=np.cross(vec4,vec1)
            crossprod3=np.cross(vec3,vec2)
            crossprod4=np.cross(vec2,vec4)
            
            n[:,i,j]=crossprod1+crossprod2+crossprod3+crossprod4
            n[:,i,j]=n[:,i,j]/np.linalg.norm(n[:,i,j])
            A[i,j]=np.linalg.norm(crossprod1)+np.linalg.norm(crossprod2)+np.linalg.norm(crossprod3)+np.linalg.norm(crossprod4)
    return A,n