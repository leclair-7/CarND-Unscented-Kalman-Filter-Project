
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')


# In[2]:


data = pd.read_csv("../build/nis_array_.log",delimiter="\n")
t = [7.8 for i in range(498)]
ts = np.arange(0,498,1)


# In[3]:


plt.plot(ts, t, label='first plot')
plt.plot(ts, data, label='second plot')
plt.legend


# If the curve is way under, we're overestimating the uncertainty in the system; if half of the curve is over, we're underestimating the uncertainty
