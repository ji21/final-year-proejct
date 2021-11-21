import pandas as pd
import numpy as np
#import sklearn as sk


df = pd.read_csv('data.csv')
df.drop('energy')
df.drop('homo')
df.drop('lumo')
print(df)
