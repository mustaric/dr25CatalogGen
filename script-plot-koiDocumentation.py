# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 14:19:11 2017

Script to create figures for the KOI Documentation based on the final Robovetter

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
    rops,a,b,c=io.loadRVInOut(id)
    banned=p.read_csv('/soc/nfs/so-nfs/DR25/OPS/DATA/DR25-TCE-DirtyMulti-BanList.txt',comment='#')
    #rops=rops.drop(banned.tce)
    ops=io.getNewCandidates(rops)
    
    
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


#%%
id=62353
outdir='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/PlotsForDocuments'

ops,inj,inv=getData(id,rtype='ALL')


#banned=p.read_csv('/soc/nfs/so-nfs/DR25/OPS/DATA/DR25-TCE-DirtyMulti-BanList.txt',comment='#')
#ops=opskois.drop(banned.tce)

#%%
prv.plotHistResults(ops,nb=130)
plt.ylim((0,1400))
plt.savefig('%s/dr25ops-hist.png' % outdir)

#%%

nEBs= len(ops[ops.isdr25koi & (ops.S==1) & (ops.C==0) & (ops.E==0) & (ops.N==0)])
nPCs=len(ops[(ops.disp=='PC') & (ops.isdr25koi)])
nKOIs=len(ops[(ops.isdr25koi)])
NtlKOIs=len(ops[(ops.N==0) & ops.isdr25koi])
nCEN= len(ops[ops.isdr25koi & (ops.C==0) & (ops.E==0) & (ops.N==0)])


print "Number of EBs (NSCE=0100)=%u" % nEBs
print "Number of PCs = %u" % (nPCs)
print "number of KOIs = %u" % (nKOIs)
print "number of TL KOIs = %u" % (NtlKOIs)
print "number of onTarget KOIs (CEN=000) = %u" % (nCEN)

 #%%
metric1='period'
metric2='mes'
range1=(100,500)
range2=(7,10)
prv.plotRelCompScore(ops,inv,inj,metric1,metric2,range1,range2,scores=np.arange(0,1,.1),\
    fpscore=0.75,Rlim=(.2,1.02),Clim=(0.8,.2))
plt.annotate("Period:100-500 days\nMES:7--10",(0.4,0.988),xycoords="data",fontsize=13)
plt.title("Completeness and Reliablity for different Score Thresholds")

plt.savefig('%s/CR-adjustScore-DR25.png' % outdir)