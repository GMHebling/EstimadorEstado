#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 13:48:56 2021

@author: vitormelo
"""
import pandas as pd
import numpy as np



dflines=pd.read_csv("Lines.csv",skiprows=(1))
dflinecodes=pd.read_csv("linecodes_estimador.csv")

parm=["Raa","Xaa","Rab","Xab","Rac","Xac","Rbb","Xbb","Rbc","Xbc","Rcc","Xcc","Caa","Cab","Cac","Cbb","Cbc","Ccc"]
dflines[parm]=0.000

for index,row in dflinecodes.iterrows():
    dflines[dflines["LineCode"].isin([row["Name"]])]
    for item in parm:
        dflines.loc[dflines["LineCode"].isin([row["Name"]]),item]=row[item]
        
        
        
dflines["Phases"]=dflines["Phases"].replace("ABC",7)
dflines["Phases"]=dflines["Phases"].replace("BC",6)
dflines["Phases"]=dflines["Phases"].replace("CA",5)
dflines["Phases"]=dflines["Phases"].replace("AB",4)
dflines["Phases"]=dflines["Phases"].replace("C",3)
dflines["Phases"]=dflines["Phases"].replace("B",2)
dflines["Phases"]=dflines["Phases"].replace("A",1)

dflines["Length"]=dflines["Length"]/1000
columns=["Bus1","Bus2","Phases","Length","Raa","Xaa","Rab","Xab","Rac","Xac","Rbb","Xbb","Rbc","Xbc","Rcc","Xcc","Caa","Cab","Cac","Cbb","Cbc","Ccc"]
df_DLIN=dflines[columns]
df_DLIN.to_csv("DLIN.csv",float_format='%.7f',header=None,index=None)
