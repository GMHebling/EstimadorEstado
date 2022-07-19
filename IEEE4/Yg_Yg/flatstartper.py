#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


dfVin=pd.read_csv("Vinicial1.csv",header=None)
vin=dfVin.values
x=np.random.rand(len(vin),len(vin[1]))
vout=vin
vout[:,1:]=vin[:,1:]-x[:,1:]*0.2
vout=pd.DataFrame(vout)
vout.to_csv("Vinicial.csv", sep = ',', index=False, header=False)