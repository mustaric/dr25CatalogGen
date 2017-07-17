# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:00:21 2016

@author: smullall

Code to play with making cuts on the stellar population.
"""

import rvIO as io
import matplotlib.pyplot as plt

tcefile="/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt"
opsfile="/home/smullall/Kepler/RoboVetter/DR25/svn/robovet/RoboVetterOut-OPS.txt"

opsdata=io.createRVResults(opsfile,tcefile)
want=(opsdata['period']>330) & (opsdata['period']<390);
longopsdata=opsdata[want]

want=opsdata['N']==0
pcrvdata=opsdata[want]

invtcefile="/soc/nfs/so-nfs/DR25/INV/DATA/TCEs.txt"
invfile="/home/smullall/Kepler/RoboVetter/DR25/svn/robovet/RoboVetterOut-INV.txt"

invdata=io.createRVResults(invfile,invtcefile)
want=(invdata['period']>330) & (invdata['period']<390)
longinvdata=invdata[want]


injtcefile="/soc/nfs/so-nfs/DR25/INJ/DATA/TCEs.txt"
injfile="/home/smullall/Kepler/RoboVetter/DR25/svn/robovet/RoboVetterOut-INJ-PlanetOn.txt"

injdata=io.createRVResults(injfile,injtcefile)

stelFile='/soc/nfs/so-nfs/DR25/other/stellarProps-DR25.csv'

steldata=io.readStellarInternal(stelFile)

staropsdata=opsdata.merge(steldata,how='left',left_on='kic',right_index=True,suffixes=("","_stel"))

starinvdata=invdata.merge(steldata,how='left',left_on='kic',right_index=True,suffixes=("","_stel"))

starinjdata=injdata.merge(steldata,how='left',left_on='kic',right_index=True,suffixes=("","_stel"))



#%%
#This takes a very long time
onesteldata=io.addTceToStellar(steldata,longopsdata,'OpsLong')

twosteldata=io.addTceToStellar(onesteldata,opsdata,'Ops')

threesteldata=io.addTceToStellar(twosteldata,invdata,'Inv')

foursteldata=io.addTceToStellar(threesteldata,longinvdata,'InvLong')

fivesteldata=io.addTceToStellar(foursteldata,pcrvdata,'pcLong')

sixsteldata=io.addTceToStellar(fivesteldata,longopsdata,'OpsSection')

stellarCounts=sixsteldata;

#%%

wantinv=stellarCounts['Inv']>0;
wantinvlong=stellarCounts['InvLong']>0;
wantopslong=stellarCounts['OpsLong']>0;
wantlongpcs=stellarCounts['pcLong']>0;

#%%
plt.subplot(311)
plt.hist(stellarCounts['cdppSlope'][wantinvlong],100,histtype='step',color='m',lw=2,normed=True)
plt.hist(stellarCounts['cdppSlope'][~wantinvlong],100,histtype='step',color='b',lw=2,normed=True)
plt.title('long inversion')
plt.subplot(312)
plt.hist(stellarCounts['cdppSlope'][wantops],100,histtype='step',color='m',lw=2,normed=True)
plt.hist(stellarCounts['cdppSlope'][~wantops],100,histtype='step',color='b',lw=2,normed=True)
plt.title('long ops')
plt.subplot(313)
plt.hist(stellarCounts['cdppSlope'][wantlongpcs],50,histtype='step',color='m',lw=2,normed=True)
plt.hist(stellarCounts['cdppSlope'][~wantlongpcs],50,histtype='step',color='b',lw=2,normed=True)
plt.title('N flag = 0 Ops')

#%%

def plotHistStar(want1,want2,want3,stellarCounts,param,bin):
    """
    """
    plt.figure()
    plt.subplot(311)
    n,bins,p=plt.hist(stellarCounts[param][want1],bins=bin,histtype='step',color='m',lw=2,normed=True)
    plt.hist(stellarCounts[param][~want1],bins=bins,histtype='step',color='b',lw=2,normed=True)
    plt.title('long inversion')
    plt.subplot(312)
    n,bins,p=plt.hist(stellarCounts[param][want2],bins=bin,histtype='step',color='m',lw=2,normed=True)
    plt.hist(stellarCounts[param][~want2],bins=bins,histtype='step',color='b',lw=2,normed=True)
    plt.title('long ops')
    plt.subplot(313)
    n,bins,p=plt.hist(stellarCounts[param][want3],bins=bin,histtype='step',color='m',lw=2,normed=True)
    plt.hist(stellarCounts[param][~want3],bins=bins,histtype='step',color='b',lw=2,normed=True)
    plt.title('N flag = 0 Ops')


def getSvnNum(fname):
    """
    Get the first line of file that contains the id string
    """
    fid=open(fname,'r')
    line=fid.readline()
    line=line.replace('$',' ')
    
    return line

#%%

plotHistStar(wantinvlong,wantopslong,wantlongpcs,stellarCounts,'cdppSlope')
plotHistStar(wantinvlong,wantopslong,wantlongpcs,stellarCounts,'kpMag')
bin=np.linspace(0,5000,num=200)
plotHistStar(wantinvlong,wantopslong,wantlongpcs,stellarCounts,'cdpp105',bin)
plotHistStar(wantinvlong,wantopslong,wantlongpcs,stellarCounts,'cdppAve',bin)


#%%

#Run grid for certain cuts.
import plotRobovetter as prv
import numpy as np

svnId=getSvnNum(opsfile)

cdppcut=300;


want=(starinvdata['cdppSlope']<-0.1) & (starinvdata['cdpp105']< cdppcut);
invcut=starinvdata[want]
print "inv: %u" % len(invcut)

want=(starinjdata['cdppSlope']<-0.1) & (starinjdata['cdpp105']< cdppcut);
injcut=starinjdata[want]
print "inj: %u" % len(injcut)


want=(staropsdata['cdppSlope']<-0.1) & (staropsdata['cdpp105']< cdppcut);
opscut=staropsdata[want]
print "ops: %u" % len(opscut)


prv.plotFullPageGrid(opscut,invcut,injcut)
figureTitle="cdppSlope<-0.1 and cdpp10.5 <%s \n %s " % (cdppcut,svnId)
plt.annotate(figureTitle,xy=(0.20,0.9),xycoords='figure fraction', fontsize=11)

outname="/home/smullall/Kepler/RoboVetter/DR25/stats/starAnalysis/gridStellarCut1-r61459.png"
plt.savefig(outname)

#%%
pcs=(opscut['disp']=='PC') ;#& (opscut['score']>0.5)
x=np.log10(opscut.loc[pcs,'period'])
y=np.log10(opscut.loc[pcs,'mes'])
prv.plot2dHist(x,y,[-.35,2.95],[.83,1.7],nxbins=100,nybins=80,xlabel="log(Period)",ylabel="log(MES)")
plt.annotate(figureTitle,xy=(0.20,0.9),xycoords='figure fraction', fontsize=11)

outname="/home/smullall/Kepler/RoboVetter/DR25/stats/starAnalysis/pcStellarCut1-r61459.png"
plt.savefig(outname)

#%%
plt.figure()
climRange=(0,100)
want=(starinvdata['cdppSlope']<-0.1) & (starinvdata['cdpp105']< 500) & (starinvdata['period']>330) & (starinvdata['mes']<=390);
data=starinvdata[want]
title='Inversion, Star Cut'
period=data['period']
mes=data['mes']
passed=data['disp']=='PC'
prv.plotGrid(period,mes,passed,climRange)
plt.title(title)

#%%
plt.figure()
want=(staropsdata['cdppSlope']<-0.1) & (staropsdata['cdpp105']< 500 )& (staropsdata['period']>330) & (staropsdata['mes']<=390);
data=staropsdata[want]
print len(data)

climRange=(0,100)
title='Ops, Star Cut'
period=data['period']
mes=data['mes']
passed=data['disp']=='PC'
prv.plotGrid(period,mes,passed,climRange)
plt.title(title)

#%%


#Ops population cuts for certain periods. uses StellarCounts

#Create 2dhistogram of the ops with TCEs compared to those that do not have TCEs.
wantSect=stellarCounts['OpsSection']>0;
x=(stellarCounts['cdppSlope'][wantSect],stellarCounts['cdppSlope'][~wantSect])
y=(stellarCounts['cdpp105'][wantSect], stellarCounts['cdpp105'][~wantSect])
prv.plot2dMultiNorm(x,y,[-0.8,.8],[0,1000],nxbins=80,nybins=80,xlabel="cdppSlope",ylabel="cdpp 10.5hrs")
plt.annotate("Stars with 330-390 day TCEs (blue and 2dhist) compared to those without (red). \nNormalized Histograms",xy=(0.1,0.95),xycoords='figure fraction',fontsize=12)
outname="/home/smullall/Kepler/RoboVetter/DR25/stats/starAnalysis/HistPeriodSection.png"
plt.savefig(outname)

#%%
wantops=stellarCounts['Ops']>0;
x=(stellarCounts['cdppSlope'][wantops],stellarCounts['cdppSlope'][~wantops])
y=(stellarCounts['cdpp105'][wantops], stellarCounts['cdpp105'][~wantops])
prv.plot2dMultiNorm(x,y,[-0.8,.8],[0,1000],nxbins=80,nybins=80,xlabel="cdppSlope",ylabel="cdpp 10.5hrs")
plt.annotate("Stars with TCEs (blue and 2dhist) compared to those without (red). \nNormalized Histograms",xy=(0.1,0.95),xycoords='figure fraction',fontsize=12)
outname="/home/smullall/Kepler/RoboVetter/DR25/stats/starAnalysis/HistPeriodOpsAll.png"
plt.savefig(outname)


#%%
#Quick look at scores
#Make a cut at scores needs to be greater than some level and create the pc curve.
scorecut=0.9
pcs=(staropsdata['disp']=='PC') & (staropsdata['score']>scorecut)
N=len(pcs[pcs])
x=np.log10(staropsdata.loc[pcs,'period'])
y=np.log10(staropsdata.loc[pcs,'mes'])
prv.plot2dHist(x,y,[-.35,2.95],[.83,1.7],nxbins=100,nybins=100,xlabel="log(Period)",ylabel="log(MES)")
figTitle="PCs with Score > %f   N=%u \n %s" % (scorecut,N,svnId)
plt.annotate(figTitle,xy=(0.20,0.9),xycoords='figure fraction', fontsize=11)

outname="/home/smullall/Kepler/RoboVetter/DR25/stats/starAnalysis/pcScoreCut4-r61459.png"
plt.savefig(outname)

#%%
scorecut=0.2
injpcs=(starinjdata['disp']=='PC') & (starinjdata['score']>scorecut)
prv.plotGrid(starinjdata['period'],starinjdata['mes'],injpcs,drange=(0,100))

#%%
scorecut=0.2
invpcs=(starinvdata['disp']=='PC') & (starinvdata['score']>scorecut)
plt.clf()
prv.plotGrid(starinvdata['period'],starinvdata['mes'],invpcs,drange=(0,100))