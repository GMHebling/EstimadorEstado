#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 18:08:23 2022

@author: vitorhpmelo
"""


import pandas as pd
import numpy as np




df_DLIN=pd.read_csv("../DLIN.csv",header=None)
df_DBAR=pd.read_csv("../DBAR.csv",header=None)





df_DLIN.columns = ["DE", "PARA", "FASES", "COMPRIMENTO", "Raa","Xaa", "Rab","Xab", "Rac","Xac", "Rbb","Xbb",\
            "Rbc","Xbc", "Rcc","Xcc",\
            "Baa", "Bab", "Bac", "Bbb", "Bbc", "Bcc" ]
df_DBAR.columns = ["ID", "LIGACAO", "FASES", "TENSAO_NOM", "PNOM_A", "PNOM_B", "PNOM_C", "QNOM_A", "QNOM_B", "QNOM_C", \
            "ZIP", "VNOM_A", "VNOM_B", "VNOM_C", "ANG_VNOM_A", "ANG_VNOM_B", "ANG_VNOM_C"]
       
df_nadj=pd.DataFrame(df_DBAR["ID"])
df_nadj["nadj"]=0

for index, row in df_DLIN.iterrows():
    mask=df_nadj["ID"]==row["DE"]
    df_nadj.at[mask,"nadj"]=df_nadj[mask]["nadj"]+1
    mask=df_nadj["ID"]==row["PARA"]
    df_nadj.at[mask,"nadj"]=df_nadj[mask]["nadj"]+1
    
df_nadj=df_nadj.sort_values(by="nadj",ascending=(False))   
