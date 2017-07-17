# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 16:23:54 2017

@author: smullall
"""

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
fpsc=0.3
metric1='period'
metric2='mes'
range1=(200,500)
range2=(7,10)
scores=np.arange(0,1,.05)

prv.plotRelCompScore(ops,inv,inj,metric1,metric2,range1,range2,fpscore=fpsc,Rlim=(.2,1.02),Clim=(1.0,.2))
plt.title('Box2 Scores (fpscore=%.1f)  r%s' % (fpsc,str(id)))
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/adjScore/rvc-20170202-box2-%s.png' % str(id))

metric1='period'
metric2='mes'
range1=(10,200)
range2=(7,10)
scores=np.arange(0,1,.05)

prv.plotRelCompScore(ops,inv,inj,metric1,metric2,range1,range2,fpscore=fpsc)
plt.title('Box1 Scores (fpscore=%.1f) r%s' % (fpsc,str(id)))
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/adjScore/rvc-20170202-box1-%s.png' % str(id))

range1=(200,500)
range2=(10,20)
scores=np.arange(0,1,.05)

prv.plotRelCompScore(ops,inv,inj,metric1,metric2,range1,range2,fpscore=fpsc)
plt.title('Box5 Scores (fpscore=1.0) r%s' % str(id))
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/adjScore/rvc-20170202-box5-%s.png' % str(id))