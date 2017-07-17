# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 16:43:13 2017

@author: smullall
"""

svnId=62256
root="/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r%s/" % (str(svnId))

#Read in data
ops,a,b,c=io.loadRVInOut(svnId)
ops['flags']=ops['flags'].fillna('ABSENT_FLAG')

rawinv,a,b,c=io.loadRVInOut(svnId,type='INV')
inv=cpm.trimINVData(rawinv,'INV')
inv['flags']=inv['flags'].fillna('ABSENT_FLAG')
rawss1,a,b,c=io.loadRVInOut(svnId,type='SS1')
ss1=cpm.trimINVData(rawss1,'SS1')
ss1['flags']=ss1['flags'].fillna('ABSENT_FLAG')


#----------------------
#---------------------

meswant=ops.mes<5000
want=np.array(map(lambda x:x.find(v[0])>=0,ops['flags']),dtype=bool)


#%%
#--------------------
import plotting as pl
import pandas as p
import numpy as np
import matplotlib.pyplot as plt

#Read in the skye thresholds.

for type in ('ops','inv','ss1'):
    infile='/home/smullall/Kepler/RoboVetter/DR25/skye/skye-%s-susan-3.0sig-1.0day-45days-NoFSP-thresh.txt' % (type)
    #infile='/home/smullall/Kepler/RoboVetter/DR25/skye/skye-%s-nonzero-2.5sig-1.0day-NoFSP-NoDirty-thresh.txt' % (type) 
    #infile='/home/smullall/Kepler/RoboVetter/DR25/skye/skye-%s-susan-3.0sig-1.0day-NoFSP-thresh.txt' % (type)
    data=p.read_csv(infile,names=['skygroup','rate'],skiprows=1,delim_whitespace=True)

    thresholds=np.array(data.rate)
    plt.figure(figsize=(12,10))
    pl.plotKeplerFocalPlane(thresholds,label=True)
    plt.title("Rates for %s by skygroup - susan Skye noFSP" % type)
    plt.savefig("%s.png" % infile)