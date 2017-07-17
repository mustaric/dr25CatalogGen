# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 13:44:45 2017

Determine range of reliability in a certain box.

@author: smullall
"""


#Randomly choose 20% of the targets 100 times and calculate the effectiveness and reliability.

import pandas as p
import rvIO as io
import plotRobovetter as prv
import matplotlib.pyplot as plt
import numpy as np
import createRVPerfMetrics as crpm


#Change the revision number
id=62339
rtype="ALL"  #ALL, SS1 or INV

ops,a,b,c=io.loadRVInOut(id)
inj,a,b,c=io.loadRVInOut(id,type="INJ-PlanetOn")

if rtype == 'ALL':
    invdata,a,b,c=io.loadRVInOut(id,type='INV')
    invrvdata=crpm.trimINVData(invdata,"INV")
    ss1data,a,b,c=io.loadRVInOut(id,type='SS1')
    ss1rvdata=crpm.trimINVData(ss1data,"SS1")
    newinvdata=p.concat((invrvdata,ss1rvdata),ignore_index=True)
        
else:
   invdata,a,b,c=io.loadRVInOut(id,type=rtype)
   newinvdata=crpm.trimINVData(invdata,rtype)

inv=newinvdata


#%%
metric1='period'
metric2='mes'
range1=[100,350]
range2=[7,10]

niter=1000
frac=0.05

Ra=np.zeros(niter)
Ea=np.zeros(niter)

opswant=prv.inBox(ops,metric1,range1,metric2,range2)
invwant=prv.inBox(inv,metric1,range1,metric2,range2)

op=ops[opswant]
iv=inv[invwant]

for i in np.arange(1,niter,1):
    
    pick=np.random.rand(len(iv))<frac
    randinv=iv[pick]
    
    R,E=crpm.estimateReliability(randinv,op,s=(0,1.0))
    Ra[i]=R
    Ea[i]=E

#%%
bins=np.arange(.97,1,.002)
#plt.figure(figsize=(12,7))
plt.subplot(121)
plt.hist(Ra,20,histtype='step')
plt.xlabel('Reliability')
plt.subplot(122)
plt.hist(Ea,bins=bins,histtype='step')
plt.xlabel('Effectiveness')
plt.title(' 1000 iterations at 80 per cent of inversions')

    
    