# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 11:16:43 2016

@author: smullall
"""
import rvIO as io
import plotRobovetter as prv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import sys
import getopt as getopt

# To get dataframe with all of the ops,inj and inversion data.
injfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-INJ-PlanetOn.txt'
opsfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-OPS.txt'
invfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-INV.txt'

tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
injtfile='/soc/nfs/so-nfs/DR25/INJ/DATA/TCEs.txt'
invtfile='/soc/nfs/so-nfs/DR25/INV/DATA/TCEs.txt'

outinj='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsINJ-PlanetOn-r61884.csv'
outinv='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsINV-r61884.csv'
outops='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsOPS-r61884.csv' 
 
id='r61884'
opsfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-OPS.txt'
tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'
fedfile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
rvdata=io.createAllResults(opsfile,tcefile,fedfile,cumfile)
rvmetrics,a,b,c=io.loadRVInOut(61884)
merged=p.merge(rvdata,rvmetrics,how='inner', right_index=True,left_index=True,suffixes=("","y"))
flags=p.read_csv(outops,comment='#', header=0, delimiter=',',index_col='tce')
OpsFlags=p.merge(merged,flags,how='inner', right_index=True,left_index=True,suffixes=("","fl"))


id='r61884'
invdata,a,b,c=io.loadRVInOut(61884,type="INV")
invmetrics=io.createRVResults(invfile,invtfile)
merged=p.merge(invdata,invmetrics,how='inner', right_index=True,left_index=True,suffixes=("","y"))
flags=p.read_csv(outinv,comment='#', header=0, delimiter=',',index_col='tce')
InvFlags=p.merge(merged,flags,how='inner', right_index=True,left_index=True,suffixes=("","fl"))

injmetrics,a,b,c=io.loadRVInOut(61884,"INJ-PlanetOn")
injdata=io.createRVResults(injfile,injtfile)
merged=p.merge(injdata,injmetrics,how='inner', right_index=True,left_index=True,suffixes=("","y"))
flags=p.read_csv(outinj,comment='#', header=0, delimiter=',',index_col='tce')
InjFlags=p.merge(merged,flags,how='inner', right_index=True,left_index=True,suffixes=("","fl"))

#%%
#rvmetrics,a,b,c=io.loadRVInOut(61884)
rvmetrics,a,b,c=io.loadRVInOut(61710)
opsdata=rvmetrics
plt.figure(figsize=(8.5,5.5))
pcwant=opsdata['disp'] == 'PC'
pcops=opsdata[pcwant]
cumfile='/soc/nfs/so-nfs/DR25/other/nexsciCumulative.csv'
cumdata=io.readNexsciCumulative(cumfile)
cumpcs=(cumdata['koi_pdisposition'] == 'CANDIDATE') & (~np.isnan(cumdata['koi_max_mult_ev']))    
cumperiod=cumdata.loc[cumpcs,'koi_period']
cummes=cumdata.loc[cumpcs,'koi_max_mult_ev']
x=(np.log10(pcops['period']),np.log10(cumperiod))
y=(np.log10(pcops['mes']),np.log10(cummes))    

axHisty,axHistx = prv.plot2dMulti(x,y,[-.35,2.95],[.83,1.7],clim=[0,13],nxbins=100,nybins=80,xlabel="log(Period)",ylabel="log(MES)")
handles,labels=axHistx.get_legend_handles_labels()
labels[0]='DR25 Candidates'
labels[1]='DR24 Candidates'
axHistx.legend(handles,labels,bbox_to_anchor=(1.6, 1) )
plt.annotate(id,xy=(0.15,0.91),xycoords='figure fraction',fontsize=12) 	
outplot3="/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/rv3replots/histops-r61710.png"
plt.savefig(outplot3)

#%%

plt.figure(4)
scores=np.arange(0,.75,.02)
Cvalues=np.zeros([len(scores),9])
Rvalues=np.zeros([len(scores),9])
Ivalues=np.zeros([len(scores),9])

for (i,scorelimit) in enumerate(scores):

    #Completeness
    totdata=InjFlags
    data=totdata
    period=data['period']
    mes=data['mes']
    passed=(data['disp']=='PC') & (data['score']>=scorelimit)
    climRange=[0,100]
    plt.figure(1)
    plt.clf()
    C=prv.plotGrid(period,mes,passed,climRange)
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
    R=prv.plotOnlyReliability(invdata,opsdata,key='score',limit=scorelimit,atitle="Reliability")
    Rvalues[i,:]=R.reshape(9)
#%%
fig1=plt.figure(4,figsize=(8,5))
cv=Cvalues[:,0]/100.0
rv=Rvalues[:,0]/100.0
iv=1-Ivalues[:,0]/100.0
plt.clf()
ax1=fig1.add_subplot(111)
line1=ax1.plot(cv,rv,'ro-')
ax1.invert_xaxis()

plt.ylim([0.5,1])
plt.ylabel('Reliability')
plt.xlabel('Completeness')
plt.title('Reliability vs. Completeness as you adjust the Score\n 0-10 days, MES<10')
for i in (0,3,5,9,12,16,20):
    st=scores[i].astype('str')
    plt.annotate(st,xy=(cv[i],rv[i]))

ax2=fig1.add_subplot(111,sharex=ax1, frameon=False)
ax2.plot(cv,iv,'bx--',ms=3)
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
plt.ylabel("Effectiveness",color="blue")
plt.xlim([1,0])

name="/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/rv3replots/compRelScore-0-r61884.png"
plt.savefig(name)