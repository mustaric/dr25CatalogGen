# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 15:38:08 2017

@author: smullall
"""

import pandas as p
import numpy as np
import plotRobovetter as prv
import createRVPerfMetrics as crpm
import matplotlib.pyplot as plt
import createRVPerfMetricsPlots as pmp
import rvIO as io



def getData(id,rtype='ALL'):
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

    return ops,inj,inv


#Change the revision number
rtype="ALL"  #ALL, SS1 or INV    

ids=(62338,62339,62340)
fpscores=(.05,.1,.2,.5,.8,1.0)
pcscores=np.arange(0,1,.05)
colors=('k','b','r')
mark=('>','o','<')

#metric1='period'
#metric2='mes'
#range1=(200,500)
#range2=(7,10)

metric1='srad'
metric2='rplanet'
range1=(1.7,2.5)
range2=(1.7,2.5)

metric1b='logg'
range1b=(4.0,12.0)
metric2b='tstar'
range2b=(5300,6500)

#Initialize plot
fig1=plt.figure(figsize=(8,10))
ax1=fig1.add_subplot(211)
plt.xlabel('Reliability')
plt.ylabel('Completeness')

#fig2=plt.figure(figsize=(8,10))
ax2=fig1.add_subplot(212)
plt.xlabel('Reliability')
plt.ylabel('Npc * R /C')

for i,id in enumerate(ids):
    
    ops,inj,inv = getData(id,rtype=rtype)
    opsbox=prv.inBox(ops,metric1,range1,metric2,range2)
    invbox=prv.inBox(inv,metric1,range1,metric2,range2)
    injbox=prv.inBox(inj,metric1,range1,metric2,range2)
   
    opsbox2=prv.inBox(ops,metric1b,range1b,metric2b,range2b)
    invbox2=prv.inBox(inv,metric1b,range1b,metric2b,range2b)
    injbox2=prv.inBox(inj,metric1b,range1b,metric2b,range2b)
    
    for fps in fpscores:
        for pcs in pcscores:
  
            scores=(pcs,fps)
            injpc=prv.passes(inj,s=scores) 
            opspc=prv.passes(ops,s=scores)
            
            C=np.float(len(inj[injbox & injpc & injbox2])) / np.float(len(inj[injbox & injbox2]))
        
            R,E=pmp.estimateReliability(inv[invbox & invbox2],ops[opsbox & opsbox2],s=scores)
            if R<0:
                R=0
            
            adjPC=np.float(len(opspc[opspc & opsbox & opsbox2])) * R /C 
            
            ax1.scatter(R,C,color=colors[i],marker=mark[i])
            ax2.scatter(R,adjPC,color=colors[i],marker=mark[i])

plt.title('Adjust Scores for r62338(k), r62339(r), r62340(b)')

plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/compare/compare3Final-Npc-hz-feb02.png')


#%%

metric1='srad'
metric2='rplanet'
range1=(2.0,2.5)
range2=(1.0,1.8)

metric1b='logg'
range1b=(4.0,12.0)
metric2b='tstar'
range2b=(4000,7000)

score=(0,1.0)
for i,id in enumerate(ids):
    
    ops,inj,inv = getData(id,rtype=rtype)
    
    opsbox=prv.inBox(ops,metric1,range1,metric2,range2)  
    opsbox2=prv.inBox(ops,metric1b,range1b,metric2b,range2b)
    pcs=prv.passes(ops,s=score)
    
    ops[opsbox & opsbox2][['period','mes','score','disp','rplanet','srad','tstar','flags']].to_csv('/home/smullall/r%i-hz.txt' % id)