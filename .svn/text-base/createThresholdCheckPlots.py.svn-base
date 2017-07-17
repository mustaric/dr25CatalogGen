# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 13:33:46 2017

Code to create plots to check each threshold of the robovetter.

@author: smullall
"""


import rvIO as io
import numpy as np
import plotRobovetter as prv
import matplotlib.pyplot as plt
import createRVPerfMetrics as cpm

plt.switch_backend('Agg')
#plt.switch_backend('Qt4Agg')
"""
This code is to loop through the different metrics/thresholds and plot
the prv.plotThresholdCheck plot
and write to an output location.
"""

svnId=62273
root="/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r%s/" % (str(svnId))
#Read in OPS data
ops,a,b,c=io.loadRVInOut(svnId)
ops['flags']=ops['flags'].fillna('ABSENT_FLAG')

rawinv,a,b,c=io.loadRVInOut(svnId,type='INV')
inv=cpm.trimINVData(rawinv,'INV')
inv['flags']=inv['flags'].fillna('ABSENT_FLAG')
rawss1,a,b,c=io.loadRVInOut(svnId,type='SS1')
ss1=cpm.trimINVData(rawss1,'SS1')
ss1['flags']=ss1['flags'].fillna('ABSENT_FLAG')


dr=(0,75)
#List of Flags to Consider, metric name and threshold
flagList=[]
flagList.append(('LPP_DV_TOO_HIGH','lpp_dv',2.20))
flagList.append(('LPP_ALT_TOO_HIGH','lpp_alt',2.2))
flagList.append(('RUBBLE','ngood',2))
flagList.append(('DV_ROBO_ODD_EVEN_TEST_FAIL','dv_oesig',1.05))   #RV_OE_DV_THRESH
flagList.append(('SKYE','ngood',3.0))
#flagList.append(('ROCKY','ngood',3.0))
flagList.append(('MARSHALL','ngood',3.0))
flagList.append(('GHOST','GB_htoc',6.8))
flagList.append(('TRANSITS_NOT_CONSISTENT','sestomes',1.0))
flagList.append(('DV_SIG_PRI_OVER_FRED_TOO_LOW','S_pri_dv',1.0))
#flagList.append(('MS_ALT_DMM_FAIL','altmoddmm',1.2))
#flagList.append(('MS_DV_DMM_FAIL','dvmoddmm',1.0))
#flagList.append(('THIS_TCE_IS_A_SEC','depth',1.0))

for v in flagList:
    want=np.array(map(lambda x:x.find(v[0])>=0,ops['flags']),dtype=bool)
    prv.plotThresholdCheck(ops,want,v[0],v[1],v[2])
    name="%s%s-a-thplot.pdf" % (root,v[0])
    plt.savefig(name)
    
    #Plot up the fraction of the FPs     
    prv.plotFracFpPerMetric(ops,inv,ss1,v[0],drange=dr,xBins=[0,100,300,350,390,500])
    name="%s%s-fpgrid-thplot.pdf" % (root,v[0])
    plt.savefig(name)
    
    #Plot up the fraction of the FPs     
    #prv.plotFracFpPerMetric(ops,ss1,v[0],reltype="SS1",drange=dr)
    #name="%s%s-fpss1-thplot.pdf" % (root,v[0])
    #plt.savefig(name)
#%%
    
    
 