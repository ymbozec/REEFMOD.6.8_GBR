#%%
import pandas as pd
import numpy as np
import scipy.io

#%%
mat = scipy.io.loadmat('R0_BCA_Moore_MIROC5_45.mat')

#%%
meta = mat['META']
vals = meta
keys = meta.dtype.descr
dic = {}
# Assemble the keys and values into variables with the same name as that used in MATLAB
for i in range(len(keys)):
    key = keys[i][0]
    val = np.squeeze(vals[key][0][0])  # squeeze is used to covert matlat (1,n) arrays into numpy (1,) arrays. 
    dic[key] = val

#%%
df_reef_id = pd.DataFrame({'Reef_ID':dic['reef_ID']})

# %%
