# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 10:05:26 2016

@author: smullall
"""

import rvIO as io
import plotRobovetter as prv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pdb
import createRVPerfMetricsPlots as pmp

#Set Score Cut on next line
def passes(data,s=[0.7,0.7]):
#def passes(data,s=[0.0,1.0]):
    """
    Use disposition and scores to determine which ones pass.
    Returns true false array
    """
    
    passed= ((data['disp'] == 'PC') & (data['score']>=s[0])) | \
            ((data['disp']=='FP') & (data['score']>s[1]))  
    return passed

def fracbox2pc(data):
    
    want=(data.period>200) & (data.period<500) & (data.mes<10)
    
    wantpc=(want) & passes(data)
    
    return np.float(len(wantpc[wantpc]))/np.float(len(want[want]))
    
    
def box2(data,di=''):
    want = (data.period>200) & (data.period<500) & (data.mes<10)
    
    if di == 'PC':
        want= want & passes(data)
    elif di == 'FP':
        want= want & ~passes(data)
    
    return len(want[want])


def box2Rel(ops,inv):
    """
    Return reliaiblity for box 2
    """
    InEff=fracbox2pc(inv)
    Npc=np.float(box2(ops,di='PC'))
    Nfp=np.float(box2(ops,di='FP'))

    R=1.0 - (np.float(InEff)/(1-InEff)) * (Nfp/Npc)    

    return R


def fracboxHzpc(data,hz=(0.25,2.0)):
    
    want=(data.srad>hz[0]) & (data.srad<hz[1]) & (data.rplanet>hz[0]) & (data.rplanet<hz[1]) & (data.logg>=4.0) & (data.tstar>=4000)
    
    wantpc=(want) & passes(data)
    

    return np.float(len(wantpc[wantpc]))/np.float(len(want[want]))
    
    
def boxHz(data,di='',hz=(0.25,2.0)):
    
    want=(data.srad>hz[0]) & (data.srad<hz[1]) & (data.rplanet>hz[0]) & (data.rplanet<hz[1]) & (data.logg>=4.0) & (data.tstar>=4000)
    if di == 'PC':
        want= want & passes(data)
    elif di == 'FP':
        want= want & ~passes(data)
    
    return len(want[want])


def boxHzRel(ops,inv,hz=(0.25,2.0)):
    """
    Return reliaiblity for box 2
    """
    InEff=fracboxHzpc(inv,hz=hz)
    Npc=boxHz(ops,di='PC')
    Nfp=boxHz(ops,di='FP')
  

    R=1 - (InEff/(1-InEff)) * (np.float(Nfp)/np.float(Npc))    
    
    return R
    
def numrange(data,key='period',range=(360,380)):
    
    want=(data[key]>=range[0]) & (data[key]<=range[1]) & passes(data)
    
    return len(want[want])
    
    
#---------
outname='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/compare/compFeb022017-Score.png'
#outname='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/compare/compFeb022017.png'
scores=(0.7,0.7)
#scores=(0.0,1.0)
rvlist=[62338,62339,62340]  #,62094,62095,62113,62107,62116]
colors=['k','#FF0033','#996600']#,'#00CC00','b','c','#FF4500','y','#FFCCFF' ]
col="black, red, brown"#, green, blue, cyan, orange, yellow, lightpink"


C2=np.zeros(len(rvlist))
R2=np.zeros(len(rvlist))
Chz=np.zeros(len(rvlist))
Rhz=np.zeros(len(rvlist))
Nhzpc=np.zeros(len(rvlist))
Nhzinv=np.zeros(len(rvlist))
N370d=np.zeros(len(rvlist))
N340d=np.zeros(len(rvlist))
fConfirm=np.zeros(len(rvlist))
fCFP=np.zeros(len(rvlist))
adjHZpcs=np.zeros(len(rvlist))

confirmed='/soc/nfs/so-nfs/DR25/other/keplernames.csv'
fedFile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
fpFile='/soc/nfs/so-nfs/DR25/other/fpwg.csv'
fpdata=io.readConfirmed(fpFile)
cdata=io.readConfirmed(confirmed)
feddata=io.readFederation(fedFile)
#Run outside of loop because we want all confirmed the same.
#opsdata,a,b,c=io.loadRVInOut(rvlist[0],type="OPS")

#
hz=(0.25,2.5)
for i,id in enumerate(rvlist):
    
    ops,a,b,c=io.loadRVInOut(id,type="OPS")
    inv,a,b,c=io.loadRVInOut(id,type="INV")
    inj,a,b,c=io.loadRVInOut(id,type="INJ-PlanetOn")
    confdata,isconf=io.combineConfirmed(ops,cdata,feddata)
    fpwgdata,isfpwg=io.combineFPWG(ops,fpdata,feddata)
    #pcEye,fpEye=io.loadByEyeVet(id)
   
    C2[i]=fracbox2pc(inj)
    R2[i]=box2Rel(ops,inv)
    
    Chz[i]=fracboxHzpc(inj,hz=hz)
    Rhz[i]=boxHzRel(ops,inv,hz=hz)
    
    Nhzpc[i]=boxHz(ops,di='PC',hz=hz)
    Nhzinv[i]=boxHz(inv,di='PC',hz=hz)
    
    N370d[i]=numrange(ops,range=(360,380))
    N340d[i]=numrange(ops,range=(340,350)) + numrange(ops,range=(390,400))
    
    fConfirm[i]=np.float(len(confdata[passes(confdata)]))/np.float(len(confdata))
    fCFP[i]=np.float(len(fpwgdata[passes(fpwgdata)]))/np.float(len(fpwgdata))
    
adjHZpcs=Nhzpc*Rhz/Chz
    
#%%
ss=np.arange(30+len(colors),30,-1)    
    
plt.figure(figsize=(8.5,9))
plt.subplot(2,2,1)

plt.scatter(100*R2,100*C2,c=colors,marker='o',linewidths=0,alpha=0.8,s=ss,label="box2")
plt.scatter(100*Rhz,100*Chz,c=colors,marker='v',linewidths=0,alpha=0.8,s=ss+2,label="HZ")
plt.xlabel('Reliability percent')
plt.ylabel('Completeness percent')
#plt.xlim((-10,55))
#plt.ylim((50,100))
plt.gca().invert_xaxis()
plt.title("200-500d mes<10 and HZ %s" % str(hz) )
plt.legend(fontsize=8,loc="lower right")


plt.subplot(2,2,2)
plt.scatter(Nhzinv,Nhzpc,c=colors,marker='o',linewidths=0,alpha=0.8,s=ss, label='Actual PCs')
plt.scatter(Nhzinv, adjHZpcs,c=colors, marker='*',linewidths=0,alpha=0.8,s=ss+2,label='Adj. R/C PCs')
plt.xlabel('Number HZ Inverted PCs')
plt.ylabel('Number HZ OPS PCs')
plt.title('Number HZ PCs OPS and INV')
plt.legend(loc=2,fontsize=9,framealpha=0.5)
    
    
plt.subplot(2,2,4)
plt.scatter(N340d,N370d-N340d,c=colors,marker='o',linewidths=0,alpha=0.8,s=ss)
#plt.plot(np.arange(7,16,1),np.arange(7,16,1),'r--')
plt.ylabel('Number PCs 360-380days - Num 340d ')
plt.xlabel('Number PCs 340-350 and 390-400 days')
plt.title('370day spike height')

plt.subplot(2,2,3)
plt.scatter(100*fCFP,100*fConfirm,c=colors,marker='o',linewidths=0,alpha=0.8,s=ss)
plt.ylabel('Percent Confirmed PCs')
plt.xlabel('Percent CFPs PCs')
#plt.xlim(2.2,3.01)
plt.gca().invert_xaxis()
plt.title('Confirmed and CFP PCs')

rev=str(rvlist)
#col=str(colors)
plt.annotate(rev,xy=(0.01,0.01),xycoords='figure fraction')
plt.annotate(col,xy=(0.01,0.03),xycoords='figure fraction')
plt.annotate("Scores "+ str(scores),xy=(0.01,0.95),xycoords='figure fraction')

plt.savefig(outname)