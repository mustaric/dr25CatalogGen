# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 10:55:05 2016

@author: smullall
"""


import rvIO as io
import plotRobovetter as prv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
rc('text',usetex=True)
import pandas as p
import numpy as np

tcefile="/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt"
opsfile="/home/smullall/Kepler/RoboVetter/DR25/svn/robovet/RoboVetterOut-OPS.txt"

opsdata=io.createRVResults(opsfile,tcefile)

opsperiod=opsdata['period']
opsmes=opsdata['mes']

dr24tcefile="/soc/nfs/so-nfs/DR24/Q17/DATA/Q1Q17TCEs.txt"

colnames=('tce','kic','pn','period','epoch','mes','depth','duration','rplanet','rstar','tstar','a','rprstar','arstar','snr','teq','secmees','secphase','posmes','posphase','mesmad')
dr24data=p.read_csv(dr24tcefile,skiprows=3,names=colnames,delim_whitespace=True,index_col='tce')

dr24RVfile='/home/smullall/Kepler/RoboVetter/DR24/RoboVetter-Output.txt'
colnames=('tce','disp','N','S','C','E','comment')
dr24rv=p.read_csv(dr24RVfile,skiprows=1,names=colnames,delim_whitespace=True,index_col='tce')

dr24ops=p.merge(dr24data,dr24rv,how='left',left_index=True,right_index=True)

dr24period=dr24data['period']
dr24mes=dr24data['mes']

#%%

cumfile='/soc/nfs/so-nfs/DR25/other/nexsciCumulative.csv'
cumdata=io.readNexsciCumulative(cumfile)

cumpcs=(cumdata['koi_pdisposition'] == 'CANDIDATE') & (~np.isnan(cumdata['koi_max_mult_ev']))

cumperiod=cumdata.loc[cumpcs,'koi_period']
cummes=cumdata.loc[cumpcs,'koi_max_mult_ev']

#%%

#Compare DR24 and DR25 plots

opspcs=opsdata['disp']=='PC'
#opsfgk=(opsdata['logg']>4.0) & (opsdata['tstar']>=4000.0) & (opsdata['tstar']<7000.0)
#opspcs=opspcs & opsfgk
opsperiod=opsdata.loc[opspcs,'period']
opsmes=opsdata.loc[opspcs,'mes']
opsprad=opsdata.loc[opspcs,'rplanet']


#cumpcs=(cumdata['koi_pdisposition'] == 'CANDIDATE') 
cumpcs=(cumdata['koi_pdisposition'] == 'CANDIDATE') & (~np.isnan(cumdata['koi_max_mult_ev']))
#cumfgk=(cumdata['koi_slogg']>4.0) & (cumdata['koi_steff']>=4000.0) & (cumdata['koi_steff']<7000.0)
#cumpcs=cumpcs & cumfgk
cumperiod=cumdata.loc[cumpcs,'koi_period']
cumprad=cumdata.loc[cumpcs,'koi_prad']

x=[np.log10(opsperiod),np.log10(cumperiod)]
y=[np.log10(opsmes),np.log10(cummes)]

ylim=np.log10([7.0,50.0])
xlim=np.log10([.5,600])
yax,xax = prv.plot2dMulti(x,y,xlim,ylim,nxbins=100,nybins=100,xlabel="log(Period)",ylabel="log(MES)",makelog=False)
handles,labels=xax.get_legend_handles_labels()
labels[0]='Preliminary DR25 PCs'
labels[1]='DR24 PCs'
xax.legend(handles,labels,bbox_to_anchor=(1.7, 1) )
plt.sca(xax)
plt.title('Planet Candidates for DR25')

lab=np.array([1,3,10,30,100,300])
locat=np.log10(lab)
plt.sca(xax)
plt.xticks(locat,lab.astype(str))

lab=np.array([7,8,10,15,20,40])
locat=np.log10(lab)
plt.sca(yax)
plt.yticks(locat,lab.astype(str))

yax.axis([0,950, ylim[0],ylim[1]])
xax.axis([xlim[0],xlim[1],0,2050])
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/scheduleExt/planetCandidatesHist-scaled.png')

yax.axis([0,130, ylim[0],ylim[1]])
xax.axis([xlim[0],xlim[1],0,130])
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/scheduleExt/planetCandidatesHist.png')




#%%
#TCE plot
opsperiod=opsdata['period']
opsmes=opsdata['mes']

dr24per=dr24ops['period']
dr24mes=dr24ops['mes']

#cumperiod=cumdata.loc[c'koi_period']
#cumprad=cumdata.loc[cumpcs,'koi_prad']

x=[np.log10(opsperiod),np.log10(dr24per)]
y=[np.log10(opsmes),np.log10(dr24mes)]

ylim=np.log10([7.0,50.0])
xlim=np.log10([.5,600])
yax,xax = prv.plot2dMulti(x,y,xlim,ylim,nxbins=100,nybins=100,xlabel="log(Period)",ylabel="log(MES)",makelog=False)
handles,labels=xax.get_legend_handles_labels()
labels[0]='DR25 TCEs'
labels[1]='DR24 TCEs'
xax.legend(handles,labels,bbox_to_anchor=(1.7, 1) )
yax.axis([0,950, ylim[0],ylim[1]])
xax.axis([xlim[0],xlim[1],0,2050])


plt.sca(xax)
plt.title('Threshold Crossing Events DR25')

lab=np.array([1,3,10,30,100,300])
#lab=np.array([.5,5,10,50,100,500])
locat=np.log10(lab)
plt.sca(xax)
plt.xticks(locat,lab.astype(str))

lab=np.array([7,8,10,15,20,40])
locat=np.log10(lab)
plt.sca(yax)
plt.yticks(locat,lab.astype(str))


plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/scheduleExt/tceHist.png')

#%%
opspcs=(opsdata['disp']=='PC') 
opsperiod=opsdata.loc[opspcs,'period']
opsmes=opsdata.loc[opspcs,'mes']
opsprad=opsdata.loc[opspcs,'rplanet']
opsscore=opsdata.loc[opspcs,'score']


plt.figure()
plt.plot(np.log10(opsperiod),np.log10(opsprad),'.')
plt.xlim(-0.6,2.9)
plt.ylim(-.5,1.8)
lab=np.array([.4,1,4,10,40])
locat=np.log10(lab)
plt.yticks(locat,lab.astype(str))
plt.ylabel(r"Planet Radius (Earth radii)")

lab=np.array([0.3,1,3,10,30,100,300])
locat=np.log10(lab)
plt.xticks(locat,lab.astype(str))
plt.xlabel('Period (days)')

plt.title('Planet Candidates from DR25 r61710')

plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/scheduleExt/pc-radii-period.png')

#%%
#Get population of new KOIs sorted by MES for periods between 350 and 390
fedfile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
alldr25Data=io.createAllResults(opsfile,tcefile,fedfile,cumfile)

opsdata=alldr25Data
want=(opsdata['period']>=350) & (opsdata['period']<390) & (opsdata['disp']=='PC')

len(want[want])

these=opsdata[want].sort(columns='mes',ascending=False)
outfile='/home/smullall/Kepler/RoboVetter/DR25/stats/yearperiod/pcs.txt'

#these.to_csv(outfile,sep=',',index_label='tce')

want2=(opsdata['period']>=400) & (opsdata['period']<500) & (opsdata['disp']=='PC')
those=opsdata[want2]
plt.hist(those['duration'],bins=20)