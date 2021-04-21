import pandas as pd

med = pd.read_csv('DMED_col.csv')
med_not_zero = med[med['f'] != 0]
med_zero = med[med['f'] == 0]

med_not_zero.to_csv('DMED_py.csv', index=False)
med_zero.to_csv("DVMED_py.csv", index=False)