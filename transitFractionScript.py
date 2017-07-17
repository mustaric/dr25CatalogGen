# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 10:12:41 2016

@author: smullall
"""

import pandas as p
import numpy as np
import rvIO as io
import matplotlib.pylab as plt

opsdata,dfinop,dfoutop,tcedfop=io.loadRVInOut(61710,type="OPS")
injdata,dfininj,dfoutinj,tcedfinj=io.loadRVInOut(61710,type="INJ-PlanetOn")
invdata,dfininv,dfoutinv,tcedfinv=io.loadRVInOut(61710,type="INV")



trfropsfile="/home/smullall/Kepler/RoboVetter/DR25/numberTransits/fracTransit-OPS-Dv-20160823.txt"

opstrfr=p.read_csv(trfropsfile,names=("tce","bkjd","frac"),delim_whitespace=True)

trfrinjfile="/home/smullall/Kepler/RoboVetter/DR25/numberTransits/fracTransit-INJ-Dv-20160823.txt"

injtrfr=p.read_csv(trfrinjfile,names=("tce","bkjd","frac"),delim_whitespace=True)

trfrinvfile="/home/smullall/Kepler/RoboVetter/DR25/numberTransits/fracTransit-INV-Dv-20160823.txt"

invtrfr=p.read_csv(trfrinvfile,names=("tce","bkjd","frac"),delim_whitespace=True)

#%%
plt.subplot(311)
plt.hist(opstrfr['frac'],bins=np.arange(.1,.95,.05))
plt.title('ops')
plt.subplot(312)
plt.hist(invtrfr['frac'],bins=np.arange(.1,.95,.05))
plt.title('inv')
plt.subplot(313)
plt.hist(injtrfr['frac'],bins=np.arange(.1,.95,.05))
plt.title('inj')

#%%
def numRemoved(tce,df):
    """
    Return the number removed given the dataframe
    """
    want=df['tce']==tce
    wantfrac=(df.loc[want,'frac']<0.75) & (df.loc[want,'frac']>0)
    
    return len(wantfrac[wantfrac])

#%%
#Plot the number of transits removed vs period.
L=len(opsdata.index)
ntransits=np.zeros(L)
periods=np.zeros(L)
duration=np.zeros(L)
pns=np.zeros(L)
mes=np.zeros(L)
df=opstrfr
data=opsdata;

for i,tce in enumerate(opsdata.index):
    ntransits[i]=numRemoved(tce,df)
    periods[i]=data.loc[tce,'Period']
    duration[i]=data.loc[tce,'duration']
    pns[i]=data.loc[tce,'pn']
    mes[i]=data.loc[tce,'mes']

plt.figure()
plt.plot(periods,1500/periods-ntransits,'b.')
plt.figure()
plt.plot(duration,1500/periods-ntransits,'r.')
