# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 14:17:51 2016

@author: smullall

code to create some plots for the AAS.
"""

import rvIO as io
import plotRobovetter as prv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as p
import numpy as np

tcefile="/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt"
opsfile="/home/smullall/Kepler/RoboVetter/DR25/svn/robovet/RoboVetterOut-OPS.txt"

opsdata=io.createRVResults(opsfile,tcefile)

opsperiod=opsdata['period']
opsmes=opsdata['mes']

dr24tcefile="/soc/nfs/so-nfs/DR24/Q17/DATA/Q1Q17TCEs.txt"

#dr24data=io.readTceInfo(dr24tcefile)

colnames=('tce','kic','pn','period','epoch','mes','depth','duration','rplanet','rstar','tstar','a','rprstar','arstar','snr','teq','secmees','secphase','posmes','posphase','mesmad')
dr24data=p.read_csv(dr24tcefile,skiprows=3,names=colnames,delim_whitespace=True,index_col='tce')

dr24period=dr24data['period']
dr24mes=dr24data['mes']

cumfile='/soc/nfs/so-nfs/DR25/other/nexsciCumulative.csv'
cumdata=io.readNexsciCumulative(cumfile)

cumpcs=(cumdata['koi_pdisposition'] == 'CANDIDATE') & (~np.isnan(cumdata['koi_max_mult_ev']))

cumperiod=cumdata.loc[cumpcs,'koi_period']
cummes=cumdata.loc[cumpcs,'koi_max_mult_ev']


x=(np.log10(opsperiod),np.log10(dr24period),np.log10(cumperiod))
y=(np.log10(opsmes),np.log10(dr24mes),np.log10(cummes))

xlim=[-.34,3.0]
ylim=[.84,1.6]
yax,xax = prv.plot2dMulti(x,y,xlim,ylim,nxbins=120,nybins=140,xlabel="log(Period)",ylabel="log(MES)")

handles,labels=xax.get_legend_handles_labels()
labels[0]='DR25 TCEs'
labels[1]='DR24 TCEs'
labels[2]='DR24 Candidates'
xax.legend(handles,labels,bbox_to_anchor=(1.6, 1) )

plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/aasplots/2dhistCompare.png')

#%%
#Now want the inversion data
invFile="/home/smullall/Kepler/RoboVetter/DR25/svn/robovet/RoboVetterOut-INV.txt"
invTceFile="/soc/nfs/so-nfs/DR25/INV/DATA/TCEs.txt"

invdata=io.createRVResults(invFile,invTceFile)

invwant=(invdata['E']==0) & (invdata['C']==1) | (invdata['S']==1) | (invdata['N']==1)
invperiod=invdata.loc[invwant,'period']
invmes=invdata.loc[invwant,'mes']

opswant=(opsdata['E']==0) & (opsdata['C']==0) & (opsdata['S']==0) & (opsdata['N']==1)
opsperiod=opsdata.loc[opswant,'period']
opsmes=opsdata.loc[opswant,'mes']

x=(np.log10(invperiod),np.log10(opsperiod))
y=(np.log10(invmes), np.log10(opsmes))

yax,xax = prv.plot2dMulti(x,y,xlim,ylim,nxbins=120,nybins=140,xlabel="log(Period)",ylabel="log(MES)")
handles,labels=xax.get_legend_handles_labels()
labels[0]='DR25 INV TCEs'
labels[1]='DR25 not-transit-like TCEs'
xax.legend(handles,labels,bbox_to_anchor=(.8, 1) )

plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/aasplots/2dhistInvCompare.png')

#%%

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

x=(np.log10(opsperiod),np.log10(cumperiod))
y=(opsprad,cumprad)

ylim=[.5,10]

yax,xax = prv.plot2dMulti(x,y,xlim,ylim,nxbins=60,nybins=60,xlabel="log(Period)",ylabel="Planet Radius (R$_{Earth}$)")
handles,labels=xax.get_legend_handles_labels()
labels[0]='Preliminary DR25 PCs'
labels[1]='DR24 PCs'
xax.legend(handles,labels,bbox_to_anchor=(1.7, 1) )

plt.sca(xax)
plt.title('All Candidates')
#plt.title('FGK Dwarf Candidates')
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/aasplots/planetCandidatesHist.png')

#%%

opspcs=opsdata['disp']=='PC'
opsfgk=(opsdata['logg']>4.0) & (opsdata['tstar']>=4000.0) & (opsdata['tstar']<7000.0)
opspcs=opspcs & opsfgk
opsperiod=opsdata.loc[opspcs,'period']
opsmes=opsdata.loc[opspcs,'mes']
opsprad=opsdata.loc[opspcs,'rplanet']


#cumpcs=(cumdata['koi_pdisposition'] == 'CANDIDATE') 
cumpcs=(cumdata['koi_pdisposition'] == 'CANDIDATE') & (~np.isnan(cumdata['koi_max_mult_ev']) )
cumfgk=(cumdata['koi_slogg']>4.0) & (cumdata['koi_steff']>=4000.0) & (cumdata['koi_steff']<7000.0)
cumpcs=cumpcs & cumfgk
cumperiod=cumdata.loc[cumpcs,'koi_period']
cumprad=cumdata.loc[cumpcs,'koi_prad']
cummes=cumdata.loc[cumpcs,'koi_max_mult_ev']

x=(np.log10(opsperiod),np.log10(cumperiod))
y=(np.log10(opsmes),np.log10(cummes))

#ylim=[.5,10]
ylim=[.84,1.9]
yax,xax = prv.plot2dMulti(x,y,xlim,ylim,nxbins=60,nybins=60,xlabel="log(Period)",ylabel="log(MES)")
handles,labels=xax.get_legend_handles_labels()
labels[0]='Preliminary DR25 PCs'
labels[1]='DR24 PCs'
xax.legend(handles,labels,bbox_to_anchor=(1.7, 1) )

plt.sca(xax)
plt.title('All Candidates')
plt.title('FGK Dwarf Candidates')
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/aasplots/planetCandidatesHistFGKmes.png')