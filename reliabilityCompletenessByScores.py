# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 11:16:43 2016

@author: smullall
"""
import pandas as p
import rvIO as io
import plotRobovetter as prv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import createRVPerfMetrics as crpm
#import sys
#import getopt as getopt

id='r62339'
idnum=62339

# To get dataframe with all of the ops,inj and inversion data.
injfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-INJ-PlanetOn.txt'
opsfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-OPS.txt'
invfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-INV.txt'

tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
injtfile='/soc/nfs/so-nfs/DR25/INJ/DATA/TCEs.txt'
invtfile='/soc/nfs/so-nfs/DR25/INV/DATA/TCEs.txt'

outinj='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsINJ-PlanetOn-'+id+'.csv'
outinv='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsINV-'+id+'.csv'
outops='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsOPS-'+id+'.csv' 
 
opsfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-OPS.txt'
tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'
fedfile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016-DMedit.txt'
rvdata=io.createAllResults(opsfile,tcefile,fedfile,cumfile)
rvmetrics,a,b,c=io.loadRVInOut(idnum)
merged=p.merge(rvdata,rvmetrics,how='inner', right_index=True,left_index=True,suffixes=("","y"))
flags=p.read_csv(outops,comment='#', header=0, delimiter=',',index_col='tce')
OpsFlags=p.merge(merged,flags,how='inner', right_index=True,left_index=True,suffixes=("","fl"))

invdata,a,b,c=io.loadRVInOut(idnum,type="INV")
invtrim=crpm.trimINVData(invdata,'INV')
invdata=invtrim
invmetrics=io.createRVResults(invfile,invtfile)
merged=p.merge(invdata,invmetrics,how='inner', right_index=True,left_index=True,suffixes=("","y"))
flags=p.read_csv(outinv,comment='#', header=0, delimiter=',',index_col='tce')
InvFlags=p.merge(merged,flags,how='inner', right_index=True,left_index=True,suffixes=("","fl"))

injmetrics,a,b,c=io.loadRVInOut(idnum,type="INJ-PlanetOn")
injdata=io.createRVResults(injfile,injtfile)
merged=p.merge(injdata,injmetrics,how='inner', right_index=True,left_index=True,suffixes=("","y"))
flags=p.read_csv(outinj,comment='#', header=0, delimiter=',',index_col='tce')
InjFlags=p.merge(merged,flags,how='inner', right_index=True,left_index=True,suffixes=("","fl"))


#%%
def passes(data,s=[0.0,1.0]):
    """
    Use disposition and scores to determine which ones pass.
    Returns true false array
    """
    
    passed= ((data['disp'] == 'PC') & (data['score']>=s[0])) | \
            ((data['disp']=='FP') & (data['score']>=s[1]))  
    return passed

plt.figure(4)
scores=np.arange(0,.9,.05)
Cvalues=np.zeros([len(scores),9])
Rvalues=np.zeros([len(scores),9])
Ivalues=np.zeros([len(scores),9])
PCvalues=np.zeros([len(scores),9])

fpScoreLimit=0.0   #what to let back in if an FP
xBins = [0, 10,200,500]
yBins = [7, 10, 20, 2000]

for (i,scorelimit) in enumerate(scores):

    #Completeness
    totdata=InjFlags
    data=totdata
    period=data['period']
    mes=data['mes']
    scoresCut=[fpScoreLimit,scorelimit]
    passed=passes(data,s=scoresCut)
    climRange=[0,100]
    plt.figure(1)
    plt.clf()
    C=prv.plotGrid(period,mes,passed,climRange,xBins = xBins,yBins = yBins)
    Cvalues[i,:]=C.reshape(9)
    plt.title('Completeness')
    
    #Ineffectiveness
    totdata=InvFlags
    data=totdata
    period=data['period']
    mes=data['mes']
    passed=(data['disp']=='PC') & (totdata['score']>=scorelimit)
    climRange=[0,100]
    plt.figure(2)
    plt.clf()
    I=prv.plotGrid(period,mes,passed,climRange)
    plt.title('inEffectiveness')
    Ivalues[i,:]=I.reshape(9)
    
    #Reliability
    invdata=InvFlags
    opsdata=OpsFlags
    plt.figure(3)
    plt.clf()
    R=prv.plotOnlyReliability(invdata,opsdata,s=(scorelimit,fpScoreLimit),atitle="Reliability",xBins = xBins,yBins = yBins)
    Rvalues[i,:]=R.reshape(9)
    
    #Num PCs
    Npcs=prv.numPassed(opsdata,xBins=xBins,yBins=yBins,s=(scorelimit,fpScoreLimit))
    PCvalues[i,:]=Npcs.reshape(9)
#%%
    
for box in [0,1,2,3,4,5,6,7,8]:
#for box in [0,1,2,3]:
    
    fig1=plt.figure(4,figsize=(8,5))
    cv=Cvalues[:,box]/100.0
    rv=Rvalues[:,box]/100.0
    iv=1-Ivalues[:,box]/100.0
    plt.clf()
    ax1=fig1.add_subplot(111)
    line1=ax1.plot(cv,rv,'ro-')
    ax1.invert_xaxis()

    plt.ylim([0.4,1])
    
    plt.ylabel('Reliability')
    plt.xlabel('Completeness')
    atitle='Reliability vs. Completeness adjusting Score\n box %s fpscore %3.1f ' % (str(box),fpScoreLimit)
    plt.title(atitle)
    for i in (0,3,6,9,12,15):
        st=scores[i].astype('str')
        plt.annotate(st,xy=(cv[i],rv[i]))
    
    ax2=fig1.add_subplot(111,sharex=ax1, frameon=False)
    ax2.plot(cv,iv,'bx--',ms=3)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    plt.ylabel("Effectiveness",color="blue")
    xmin=np.min(cv)-0.5
    plt.xlim([1,.4])
    
    name="/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/adjScore/compRelScore-feb01-"+str(box)+ '-' +id+ '.png'
    plt.savefig(name)
    
    #A second figure showing SNR, and number of PCs
    plt.figure()
    plt.plot(scores,PCvalues[:,box],'ro',label='Npc')
    plt.xlabel('Score Limit')
    plt.ylabel('Number of PCs')
    plt.plot(scores,PCvalues[:,box]*(rv/cv),'*b',ms=8,label='Np*R/C')
    plt.plot(scores,PCvalues[:,box]*rv, 'sg',ms=7,label='Npc*R',alpha=0.7)
    atitle='Adjusted Number of PCs by Score\n box %s fpscore %3.1f ' % (str(box),fpScoreLimit)
    plt.title(atitle)
    plt.xlim(0.2,1.0)
    plt.legend(fontsize='small')
    plt.annotate(id,xy=(0.1,0.02),xycoords='figure fraction',fontsize=11)
    
    name="/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/adjScore/compRelScore-feb01-pcnum-"+str(box)+ '-' +id+ '.png'
    plt.savefig(name)
    
    #A Third plotting Pcvalues*rv vs cv
    plt.figure()
    plt.plot(cv,PCvalues[:,box]*rv,'sg',ms=7,alpha=0.9)
    atitle='Npc*R vs Completeness for box %s fpscore %3.1f' % (str(box),fpScoreLimit)
    plt.title(atitle)
    plt.xlabel('completeness')
    plt.ylabel('Npc * R')
    plt.annotate(id,xy=(0.1,0.02),xycoords='figure fraction',fontsize=11)
    
    name="/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/adjScore/compRelScore-feb01-snr-"+str(box)+ '-' +id+ '.png'
    plt.savefig(name)