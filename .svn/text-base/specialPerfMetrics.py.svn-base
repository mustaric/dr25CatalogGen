# -*- coding: utf-8 -*-
"""
Created on Mon May 16 12:44:05 2016

@author: smullall

Create Specialty performance metric plots.
"""

import rvIO as io
import plotRobovetter as prv
import matplotlib.pyplot as plt
import numpy as np



plt.figure(figsize=(10,6.5))
plt.subplots_adjust(hspace=0.05, wspace=0.4)
topdir='/soc/nfs/so-nfs/DR25/'
outplot='/home/smullall/Kepler/RoboVetter/DR25/stats/DR25v1/ConfFPStatsr61459.png'

confirmed='/soc/nfs/so-nfs/DR25/other/keplernames.csv'

rvFile='/home/smullall/Kepler/RoboVetter/DR25/svn/robovet/RoboVetterOut-OPS.txt'
tceFile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
fedFile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'

rvdata=io.createRVResults(rvFile,tceFile)

cdata=io.readConfirmed(confirmed)

feddata=io.readFederation(fedFile)
#Merge using the KOIs. Leave behind only the KOIs and data for all.
#Start by merging cdata and feddata.

clist=cdata.index
isConfirmed=np.zeros(len(feddata))
for i,koi in enumerate(feddata['koi']):
    
    if koi in clist:
        isConfirmed[i]=1


feddata['confirmed']=isConfirmed
alldata=rvdata.merge(feddata,left_index=True,right_index=True,how='left')

isconf=alldata['confirmed']==1

data=alldata[isconf]

title='Confirmed\nDR25 RV PC'

plt.subplot(121)
period=data['period']
mes=data['mes']
passed=data['disp']=='PC'
prv.plotGrid(period,mes,passed,(50,100))
plt.title(title)


#CFPs
fpFile='/soc/nfs/so-nfs/DR25/other/fpwg.csv'
fpdata=io.readConfirmed(fpFile)

fplist=fpdata.index
isFP=np.zeros(len(feddata))
for i,koi in enumerate(feddata['koi']):
    if koi in fplist:
        isFP[i]=1

feddata['fp']=isFP
newdata=rvdata.merge(feddata,left_index=True,right_index=True,how='left')

isfp=newdata['fp']==1

data=newdata[isfp]
climRange=(0,50)
title='CFP in DR25 RV PC'
plt.subplot(122)
period=data['period']
mes=data['mes']
passed=data['disp']=='PC'
prv.plotGrid(period,mes,passed,climRange)
plt.title(title)

#
fid=open(rvFile,'r')
line=fid.readline()
svnId=line.replace('$',' ')

figureTitle="DR25v1\n  %s" %(svnId)
plt.annotate(figureTitle,xy=(0.20,0.9),xycoords='figure fraction',  fontsize=11)

plt.savefig(outplot)