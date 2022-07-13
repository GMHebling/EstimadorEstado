import numpy as np
import pandas as pd
import math

df_bc = pd.read_csv('state_BC2_34.txt', header=None, delimiter='\t')
df_bc.columns = ['id', 'val']
df_pf = pd.read_csv('state_PF_34.txt', header=None, delimiter='\t')
df_pf.columns = ['id','val']


err = sum(pow(df_bc['val'] - df_pf['val'],2))/len(df_bc)
print(err)

rmse = math.sqrt(err)



err_mae = df_bc['val'] - df_pf['val']
mae = abs(np.nanmean(err_mae))
print('mae: ', mae)