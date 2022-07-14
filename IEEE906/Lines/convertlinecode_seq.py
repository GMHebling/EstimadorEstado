#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 10:38:36 2021

@author: vitormelo
"""

import pandas as pd
import numpy as np


df_linecodes=pd.read_csv("LineCodes.csv",skiprows=(1))


j=complex(0,1)
a=1*np.exp(j*120*np.pi/180)
T = np.array([[1,1,1],[1,a**2, a],[1, a, a**2]]);
Tinv=np.linalg.inv(T)

R_abc=np.zeros(((len(df_linecodes["Name"]),3,3)))
X_abc=np.zeros(((len(df_linecodes["Name"]),3,3)))

for index, row in df_linecodes.iterrows():
    Zseq=[]
    Z_abc=[]
    Zseq=np.array([[row["R0"]+j*row["X0"],0,0],[0,row["R1"]+j*row["X1"],0],[0,0,row["R1"]+j*row["X1"]]])
    Z_abc=Tinv@Zseq@T
    R_abc[index]=np.real(Z_abc)
    X_abc[index]=np.imag(Z_abc)
    
    
df_linecodes["Raa"]=R_abc[:,0,0]*(abs(R_abc[:,0,0])>1e-15)
df_linecodes["Xaa"]=X_abc[:,0,0]*(abs(X_abc[:,0,0])>1e-15)
df_linecodes["Rab"]=R_abc[:,0,1]*(abs(R_abc[:,0,1])>1e-15)
df_linecodes["Xab"]=X_abc[:,0,1]*(abs(X_abc[:,0,1])>1e-15)
df_linecodes["Rac"]=R_abc[:,0,2]*(abs(R_abc[:,0,2])>1e-15)
df_linecodes["Xac"]=X_abc[:,0,2]*(abs(X_abc[:,0,2])>1e-15)
df_linecodes["Rbb"]=R_abc[:,1,1]*(abs(R_abc[:,1,1])>1e-15)
df_linecodes["Xbb"]=X_abc[:,1,1]*(abs(X_abc[:,1,1])>1e-15)
df_linecodes["Rbc"]=R_abc[:,1,2]*(abs(R_abc[:,1,2])>1e-15)
df_linecodes["Xbc"]=X_abc[:,1,2]*(abs(X_abc[:,1,2])>1e-15)
df_linecodes["Rcc"]=R_abc[:,2,2]*(abs(R_abc[:,2,2])>1e-15)
df_linecodes["Xcc"]=X_abc[:,2,2]*(abs(X_abc[:,2,2])>1e-15)


df_linecodes["Caa"]=np.zeros(len(df_linecodes["Name"]))
df_linecodes["Cab"]=np.zeros(len(df_linecodes["Name"]))
df_linecodes["Cac"]=np.zeros(len(df_linecodes["Name"]))
df_linecodes["Cbb"]=np.zeros(len(df_linecodes["Name"]))
df_linecodes["Cbc"]=np.zeros(len(df_linecodes["Name"]))
df_linecodes["Ccc"]=np.zeros(len(df_linecodes["Name"]))

columns=list(range(9,27))
columns.append(0)
columns.sort()
dfnew=df_linecodes.iloc[:,columns]

dfnew.to_csv("linecodes_estimador.csv",float_format='%.7f')