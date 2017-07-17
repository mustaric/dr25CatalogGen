# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 09:41:30 2016

@author: smullall

Using the EB list from Villanova, and the TCE lists
Match a KIC,period and within a day of the epoch for things that are likely
things found before.

Print an array of all Inverted TCEs and True if should keep it and False if
should not use it for Effectiveness stats.

Specify Reasons -- EB and morph of match, vs TCE list, vs KOI vs Confirmed
morph=-1 is HB star
"""
import rvIO as io
import pandas as p
import numpy as np
import plotRobovetter as prv
import matplotlib.pyplot as plt
import pdb

def periodEpochMatch((per,ep),periods,epochs,tol,eptol):
    """
    for a particular period (per) and epoch (ep)
    return a true false array if they match the aray of periods and epochs
    tol is an accepted fractional error
    """
    depoch=eptol*per;
    ar=np.zeros(len(periods))==0
    if len(periods)==0:
        ar=[False]
    Pdiff=np.array(np.abs(periods-per))
    Ediff=np.array(np.abs(ep-epochs))   
      
    
    for i,v in enumerate(periods):
        
        if (Pdiff[i]>per*tol) & (Ediff[i]>depoch):
            ar[i]=False
        else:
            ar[i]=True
           
    return ar

def periodMatch(per,periods,tol):
    """
    for a particular period (per)
    return a true false array if they match the aray of periods
    tol is an accepted fractional error
    """
    #pdb.set_trace()
    ar=np.zeros(len(periods))==0   
    if len(periods)==0:
        ar=[False]
    Pdiff=np.array(np.abs(periods-per))
    
    for i,v in enumerate(periods):
        if (Pdiff[i]>per*tol) :
           ar[i]=False
        else:
           ar[i]=True
           
    return ar
  
def markKICs(df,kiclist,colname):
    """
    Mark those tces on a particular star
    """
    
    onStar=np.array(map(lambda x: len(set(kiclist) & set([x]))>0,df['kic']),dtype=bool)
    
    df[colname]=onStar
    
    return df
     
def addCleanColumns(df,plimit,meslimit):
    """
    Add columns saying if it is on a funny star.
    """
    #Create columns of those with certain attributes
    ebkics=ebs['KIC'][ebs['morph']<=0.6]
    df=markKICs(df,ebkics,'onEB')
    
    otherkics=other['kic']
    df=markKICs(df,otherkics,'onOther')    
    
    pulkics=pul['kic']
    df=markKICs(df,pulkics,'onBPulsator')
    
    pul2kics=[]#pul2['kic']
    df=markKICs(df,pul2kics,'onDSct')
    
    pul3kics=pul3['kic'][pul3['amp']>10000]
    df=markKICs(df,pul3kics,'onGamDor')    
    
    #On short period high snr KIC
    cumkics=cumKOI['kepid'][(cumKOI['koi_period']<plimit) & (cumKOI['koi_max_mult_ev']>meslimit)]
    df=markKICs(df,cumkics,'onCum')
    
    #on short period high snr new KIC
    want=(rvData['N']==0) & (rvData['period']<plimit) & (rvData['mes']>meslimit)
    opskics=rvData['kic'][want]
    df=markKICs(df,opskics,'onOps')
    
    df['keep']=~(df['onEB'] | df['onBPulsator'] | df['onDSct'] | df['onCum'] | df['onOps'] | df['onOther'])

    return df    
    
#%

invfile='/soc/nfs/so-nfs/DR25/INV/DATA/TCEs.txt'
invrvfile='/soc/nfs/so-nfs/DR25/INV/RoboVet/RoboVetterOut-INV.txt'
invTCE=io.createRVResults(invrvfile,invfile)
#invTCE=io.readTceInfo(invfile)

ebcatalogfile='/home/smullall/Kepler/RoboVetter/DR25/cleanINV/ebCatalog.csv'
ebs=p.read_csv(ebcatalogfile,header=0,comment='#')
ebs['epoch']=ebs['bjd0']+2400000.0-2454833.0

pulsatorfile='/home/smullall/Kepler/RoboVetter/DR25/cleanINV/pulsators.csv'
pul=p.read_csv(pulsatorfile,header=0,comment='#')

murphyfile='/home/smullall/Kepler/RoboVetter/DR25/cleanINV/murphy_deltaScuti.txt'
pul2=p.read_csv(murphyfile,header=0,comment='#')

gammadorfile='/home/smullall/Kepler/RoboVetter/DR25/cleanINV/gammDor.txt'
pul3=p.read_csv(gammadorfile,header=0,comment='#',delim_whitespace=True)

#rvfile='/soc/nfs/so-nfs/DR25/OPS/RoboVet/RoboVetterOut-OPS.txt'
#tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
#opsTCE=io.readTceInfo(tcefile)
#rvData=io.createRVResults(rvfile,tcefile)
rvData,a,b,c=io.loadRVInOut(62141,'OPS')

cumfile='/soc/nfs/so-nfs/DR25/other/nexsciCumulative.csv'
cumfile='/soc/nfs/so-nfs/DR25/other/cumulative-20170327.csv'
cumKOI=io.readNexsciCumulative(cumfile)

ss1tcefile='/soc/nfs/so-nfs/DR25/SS1/DATA/TCEs.txt'
ss1Data=io.readTceInfo(ss1tcefile)

otherfile='/home/smullall/Kepler/RoboVetter/DR25/cleanINV/Others.csv'
other=p.read_csv(otherfile,header=0,comment='#')
#%

#Create a string of Trues and put in a column in invTCEs.
keep=np.zeros(len(invTCE))==0
#reasons=['' for x in xrange(len(invTCE))]
plimit=900.0
meslimit=7 
invTCE=addCleanColumns(invTCE,plimit,meslimit)
invTCEs=invTCE[invTCE['keep']]
print "--INV--"
print 'OnEB:',  len(invTCE[invTCE['onEB']]), 100  * len(invTCE[invTCE['onEB']])/len(invTCE)
print 'OnPulsator:', len(invTCE[invTCE['onBPulsator']]), 100 * len(invTCE[invTCE['onBPulsator']])/len(invTCE)
print 'onGamDor:', len(invTCE[invTCE['onGamDor']]), 100 * len(invTCE[invTCE['onGamDor']])/len(invTCE)
print 'onCumKoi:', len(invTCE[invTCE['onCum']]), 100 * len(invTCE[invTCE['onCum']])/len(invTCE)
print 'onOpsKoi:', len(invTCE[invTCE['onOps']]), 100 * len(invTCE[invTCE['onOps']])/len(invTCE)
print 'total', len(invTCE), len(invTCE[invTCE['keep']])

#rvData=addCleanColumns(rvData,plimit,meslimit)
#want=rvData['N']==1
#allrvData=rvData
#rvData=rvData[want]
#opsTCE=rvData
#print "--OPS--"
#print 'OnEB:',  len(rvData[rvData['onEB']]), 100  * len(rvData[rvData['onEB']])/len(rvData)
#print 'OnPulsator:', len(rvData[rvData['onBPulsator']]), 100 * len(rvData[rvData['onBPulsator']])/len(rvData)
#print 'onGamDor:', len(rvData[rvData['onGamDor']]), 100 * len(rvData[rvData['onGamDor']])/len(rvData)
#print 'onCumKoi:', len(rvData[rvData['onCum']]), 100 * len(rvData[rvData['onCum']])/len(rvData)
#print 'onOpsKoi:', len(rvData[rvData['onOps']]), 100 * len(rvData[rvData['onOps']])/len(rvData)
#print 'total', len(rvData), len(rvData[rvData['keep']])
#rvData=allrvData
#rvData=rvData[rvData['keep']]

ss1Data=addCleanColumns(ss1Data,plimit,meslimit)
print "--SS1--"
print 'OnEB:',  len(ss1Data[ss1Data['onEB']]), 100  * len(ss1Data[ss1Data['onEB']])/len(ss1Data)
print 'onBPulsator:', len(ss1Data[ss1Data['onBPulsator']]), 100 * len(ss1Data[ss1Data['onBPulsator']])/len(ss1Data)
print 'onGamDor:', len(ss1Data[ss1Data['onGamDor']]), 100 * len(ss1Data[ss1Data['onGamDor']])/len(ss1Data)
print 'onCumKoi:', len(ss1Data[ss1Data['onCum']]), 100 * len(ss1Data[ss1Data['onCum']])/len(ss1Data)
print 'onOpsKoi:', len(ss1Data[ss1Data['onOps']]), 100 * len(ss1Data[ss1Data['onOps']])/len(ss1Data)
print 'total', len(ss1Data), len(ss1Data[ss1Data['keep']])
ss1TCEs=ss1Data[ss1Data['keep']]

invTCE.to_csv('/home/smullall/Kepler/RoboVetter/DR25/cleanINV/invTCEClean-900day-7mes-Mar2017.csv',sep=',',np_rep=0,columns=['keep','onEB','onCum','onOther'])
ss1Data.to_csv('/home/smullall/Kepler/RoboVetter/DR25/cleanINV/ss1TCEClean-900day-7mes-Mar2017.csv',sep=',',np_rep=0,columns=['keep','onEB','onCum','onOther'])
#%%
tag='cleanCompare-20170327.png'
all=p.concat([invTCEs,rvData,ss1TCEs])
allbins=prv.returnBins(all,num=100)
plt.figure(figsize=(8.5,11))
plt.subplots_adjust(hspace=0.5, wspace=0.65)
prv.plotMultiHist(rvData,allbins)
prv.plotMultiHist(invTCEs,allbins)
prv.plotMultiHist(ss1TCEs,allbins)
plt.annotate(xy=(.2,.94),s='Clean-Up Ops,N (blue), INV (green), SS1 (red) 20d KOI cut',xycoords="figure fraction")
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/cleanINV/all-' + tag)
#%
wantrv=rvData['period']>200
wantinv=invTCEs['period']>200
wantss1=ss1TCEs['period']>200
all=p.concat([invTCEs[wantinv],rvData[wantrv],ss1TCEs[wantss1]])
allbins=prv.returnBins(all,num=100)
plt.figure(figsize=(8.5,11))
plt.subplots_adjust(hspace=0.45, wspace=0.65)
prv.plotMultiHist(rvData[wantrv],allbins)
prv.plotMultiHist(invTCEs[wantinv],allbins)
prv.plotMultiHist(ss1TCEs[wantss1],allbins)
plt.annotate(xy=(.2,.94),s='Clean-Up P>200 Ops,N (blue), INV (green), SS1 (red) all KOI cut',xycoords="figure fraction")
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/cleanINV/pgt200-' + tag)

wantrv=rvData['period']<200
wantinv=invTCEs['period']<200
wantss1=ss1TCEs['period']<200
all=p.concat([invTCEs[wantinv],rvData[wantrv],ss1TCEs[wantss1]])
allbins=prv.returnBins(all,num=100)
plt.figure(figsize=(8.5,11))
plt.subplots_adjust(hspace=0.45, wspace=0.65)
prv.plotMultiHist(rvData[wantrv],allbins)
prv.plotMultiHist(invTCEs[wantinv],allbins)
prv.plotMultiHist(ss1TCEs[wantss1],allbins)
plt.annotate(xy=(.2,.94),s='Clean-Up P<200 Ops,N (blue), INV (green), SS1 (red) all KOI cut',xycoords="figure fraction")
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/cleanINV/plt200-' + tag)

wantrv=rvData['mes']<20
wantinv=invTCEs['mes']<20
wantss1=ss1TCEs['mes']<20
all=p.concat([invTCEs[wantinv],rvData[wantrv],ss1TCEs[wantss1]])
allbins=prv.returnBins(all,num=100)
plt.figure(figsize=(8.5,11))
plt.subplots_adjust(hspace=0.45, wspace=0.65)
prv.plotMultiHist(rvData[wantrv],allbins)
prv.plotMultiHist(invTCEs[wantinv],allbins)
prv.plotMultiHist(ss1TCEs[wantss1],allbins)
plt.annotate(xy=(.2,.94),s='Clean-Up mes<20 Ops,N (blue), INV (green), SS1 (red) all KOI cut',xycoords="figure fraction")
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/cleanINV/meslt20-' + tag)

wantrv=(rvData['period']>330) & (rvData['period']<410)
wantinv=(invTCEs['period']>330) & (invTCEs['period']<410)
wantss1=(ss1TCEs['period']>330) & (ss1TCEs['period']<410)
all=p.concat([invTCEs[wantinv],rvData[wantrv],ss1TCEs[wantss1]])
allbins=prv.returnBins(all,num=100)
plt.figure(figsize=(8.5,11))
plt.subplots_adjust(hspace=0.45, wspace=0.65)
prv.plotMultiHist(rvData[wantrv],allbins)
prv.plotMultiHist(invTCEs[wantinv],allbins)
prv.plotMultiHist(ss1TCEs[wantss1],allbins)
plt.annotate(xy=(.2,.94),s='Clean-Up P:330-410 Ops,N (blue), INV (green), SS1 (red) all KOI cut',xycoords="figure fraction")
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/cleanINV/p330-410-' + tag)

#%%
want=(invTCE['mes']>15) & (invTCE['keep']==1) &  (invTCE['N']==0)
invTCEs[want]
len(invTCEs[invTCEs['N']==0])
#%%
#plt.figure(1)
#plt.plot(invTCE['period'][invTCE['keep']],np.log10(invTCE['mes'][invTCE['keep']]),'.')
#plt.plot(invTCE['period'][~invTCE['keep']],np.log10(invTCE['mes'][~invTCE['keep']]),'r.')
#print len(invTCE[invTCE['keep']])
##%%
#plt.figure()
#x=[np.log10(invTCE['period'][invTCE['keep']]),np.log10(rvData['period'][rvData['N']==1]),np.log10(ss1TCEs['period'])]
#y=[np.log10(invTCE['mes'][invTCE['keep']]),np.log10(rvData['mes'][rvData['N']==1]),np.log10(ss1TCEs['mes'])]
#xlim=(np.log10(.5),np.log10(700))
#ylim=(np.log10(7.0),np.log10(500))
#prv.plot2dMulti(x,y,xlim,ylim,nxbins=80,nybins=80,xlabel="log period",ylabel="log MES",makelog=False)
#plt.annotate(xy=(.1,.92),xycoords='figure fraction',s='SS1(blue) compared to DR25 N-Flag (red) after removing EBs and KOIs')
#plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/cleanINV/ss1-KOIEBClean-wOrig.png')
##plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/cleanINV/inv-KOIEBClean.png')
##%%
#plt.figure()
#wantinv=(invTCE['period']>200) & (invTCE['mes']<20)
#wantrv=(rvData['period']>200) & (rvData['mes']<20)
#wantss1=(ss1TCEs['period']>200) & (ss1TCEs['mes']<20)
#
#x=[np.log10(invTCE['period'][invTCE['keep'] & wantinv]),np.log10(rvData['period'][(rvData['N']==1) & wantrv]),np.log10(ss1TCEs['period'][wantss1])]
#y=[np.log10(invTCE['mes'][invTCE['keep'] & wantinv]),np.log10(rvData['mes'][(rvData['N']==1) & wantrv]), np.log10(ss1TCEs['mes'][wantss1])]
#xlim=(np.log10(200),np.log10(700))
#ylim=(np.log10(7.0),np.log10(20))
#prv.plot2dMulti(x,y,xlim,ylim,nxbins=80,nybins=80,xlabel="log period",ylabel="log MES",makelog=False)
#plt.annotate(xy=(.1,.92),xycoords='figure fraction',s='SS1(blue) compared to DR25 N-Flag (red) after removing EBs and KOIs')
#plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/cleanINV/ss1-KOIEBClean-LongPerLowMes.png')
##plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/cleanINV/inv-KOIEBClean-LongPerLowMes.png')
##%%
##What PCs are in this remaining set.
##Start with high MES.
#
#(invdata,a,b,c)=io.loadRVInOut(61801,type="INV")
#
#want=(invTCE['mes']>15) & (invTCE['keep']==1) &  (invdata['N']==0)
#len(want[want])
#print invdata[want]['mes']




#%%
#
#for i,tce in enumerate(invTCE.index):
#    period=invTCE.loc[tce]['period']
#    epoch=invTCE.loc[tce]['epoch']
#    kic=int(tce[:9])   
#    
##    #Match to EBs
##    want=(ebs['KIC'] == kic) & (ebs['morph']<0.7)
##    if len(want[want])>0:
##    #ar=periodMatch(period,ebs[want]['period'],tol)
##    #if any(ar):
##        invTCE.loc[tce,'keep']=False
##        invTCE.loc[tce,'reason']=invTCE.loc[tce,'reason'] +'EB--'
##    #Match to pulsators
##    want=(pul['kic'] == kic)
##    if len(want[want])>0:
##        invTCE.loc[tce,'keep']=False
##        invTCE.loc[tce,'reason']=invTCE.loc[tce,'reason'] +'BPulsator--'    
##    
##    want=(pul2['kic'] == kic)
##    if len(want[want])>0:
##        invTCE.loc[tce,'keep']=False
##        invTCE.loc[tce,'reason']=invTCE.loc[tce,'reason'] +'DScutiType--' 
#    #Match current TCE List
#    want=(rvData['kic'] == kic) & (rvData['N']==0)
#    ar=periodEpochMatch((period,epoch),rvData[want]['period'],rvData[want]['epoch'],tol,eptol)
#    if any(ar):
#        invTCE.loc[tce,'keep']=False
#        invTCE.loc[tce,'reason']=invTCE.loc[tce,'reason'] +'opsKOI--'
#    
#    #Match Cumulative Table
#    want=(cumKOI['kepid'] == kic)   #& (cumKOI['koi_fpflag_nt']==0)
#    ar=periodEpochMatch((period,epoch),cumKOI[want]['koi_period'],cumKOI[want]['koi_time0bk'],tol,eptol)
#    if any(ar):
#        invTCE.loc[tce,'keep']=False
#        invTCE.loc[tce,'reason']=invTCE.loc[tce,'reason'] + 'CummulativeTable--'
#    
#    #Remove all Period <125 (3 transits in a year) day objects from the Cummulative table.
#    want=((cumKOI['kepid'] == kic) & (cumKOI['koi_period']<125)) | ((rvData['kic'] == kic) & (rvData['period']<125))
#    if len(want[want])>0:
#        invTCE.loc[tce,'keep']=False
#        invTCE.loc[tce,'reason']=invTCE.loc[tce,'reason'] +'ShortPeriodKOI--'
