# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 13:12:24 2017

@author: smullall

Collect both DV and MCMC fits and try to plot the differences in the fits.


"""




import numpy as np
import pandas as p
import matplotlib.pyplot as plt
import testHammer as th
import rvIO as io
import plotRobovetter as prv
import createRVPerfMetricsPlots as pmp

id=62353
ops,inj,inv=io.getAllData(id,tcedata="DATA-SUP")

hammerfile='/home/smullall/Kepler/Nexsci/deliveryWorkspace/DR25KOI/testHammer-complete.txt'
dfham=th.readHammerDF2(hammerfile)


#%%

#Merge the two dataframes on the KOI index.

kois=ops[ops.isdr25koi]
both=p.merge(kois,dfham,how='left',left_on='dr25koiName',right_index=True,suffixes=('','_ham'))
both.set_index('dr25koiName',inplace=True)
#%%

pcs=both.disp=='PC'

plt.figure()
plt.plot(both.period[pcs],both.rplanet[pcs],'r.',ms=6)
plt.plot(both.koi_period[pcs],both.koi_prad_ham[pcs],'b.',ms=5)
plt.xlabel('period')
plt.ylabel('planet radius')

#%%

xname1='koi_insol_ham'
yname1='koi_prad_ham'
xname2='srad'
yname2='rplanet'

pcs=(both.disp=='PC') & (both.koi_insol != 0)

plt.figure()
plt.plot(both.srad[pcs],both.rplanet[pcs],'rs',ms=5,label='DV Fits Sup')
plt.plot(both.koi_insol_ham[pcs],both.koi_prad_ham[pcs],'bo',ms=5,label='MCMC (201702221)')
plt.plot(both.koi_insol[pcs],both.koi_prad[pcs],'g^',ms=5,label='Cumulative Fits')
plt.xlabel('Insolation Flux')
plt.ylabel('Planet Radius (earth)')
plt.xlim((0,3))
plt.ylim((0,3))
plt.legend(loc='lower left')


for v in both.index[pcs]:
    x=(both.loc[v][[xname1]].astype(float),both.loc[v][[xname2]].astype(float))
    y=(both.loc[v][[yname1]].astype(float),both.loc[v][[yname2]].astype(float))
    
    plt.plot(x,y,'-k')

plt.savefig('/home/smullall/Kepler/Nexsci/deliveryWorkspace/DR25KOI/plots/mcmc-dvFitsSup-InsolPrad2.png')

#%%

xname=('period','koi_period_ham')
yname=('depth','koi_depth_ham')
labels=('DV Fits Sup','MCMC (201702221)')
want=(both.disp=='PC') & (both.koi_fittype!='none')

prv.plotChanges(both,xname,yname,want,labels)
plt.xlim((0.5,50))
plt.ylim((0,1000))
#%%
plt.figure()
diff=both.depth-both.koi_depth_ham
plt.plot(both.koi_period_ham,diff,'.')

plt.figure()
plt.plot(both.koi_depth_ham,both.depth,'.')

#%%
#Determine whihch have a rprstar that differs from koi_ror_ham by more than 1 sigma

error=(both.koi_ror_err1 + -1 * both.koi_ror_err2)/2
diff=np.abs(both.rprstar - both.koi_ror_ham)

want=(diff > 5.0*error)
pcs=both.disp=='PC'

print len(want[want & pcs])
both[want & pcs][['kic','koi_period_ham','rprstar','koi_ror_ham','koi_ror_err1','mes']].to_csv('/home/smullall/Kepler/Nexsci/deliveryWorkspace/DR25KOI/rprstar-comparison-mcmcVdvsup.csv')

error=(both.koi_dor_err1 + -1 * both.koi_dor_err2)/2
diff=np.abs(both.arstar - both.koi_dor)

want=(diff > 5.0*error)
pcs=(both.disp=='PC') & (both.koi_dor >0)

print len(want[want & pcs])
both[want & pcs][['kic','koi_period_ham','arstar','koi_dor','koi_dor_err1','mes']].to_csv('/home/smullall/Kepler/Nexsci/deliveryWorkspace/DR25KOI/arstar-comparison-mcmcVdvsup.csv')

#%%
plt.figure()
pcs=(both.disp=='PC') & (both.koi_dor>0)
n,b,p=plt.hist(both.arstar[pcs],bins=300,histtype='step')
plt.hist(both.koi_dor[pcs],bins=b,histtype='step')

plt.figure()
diffs=(both.arstar[pcs]-both.koi_dor[pcs])
plt.plot(diffs,both.arstar[pcs],'.')

#%%
#Get just the fields that we need to send through some common plots for ops.

params=[u'score', u'disp', u'N', u'S', u'C', u'E', u'flags', 'koi_period_ham','koi_time0bk_ham',\
        'koi_prad_ham','koi_insol_ham','koi_max_mult_ev_ham','koi_srad_ham','koi_steff_ham','koi_slogg_ham','koi_smass','koi_sma','isdr25koi']
hamops=both[params]

#Convert the names to those I used before.
oldnames=['koi_period_ham','koi_time0bk_ham',\
        'koi_prad_ham','koi_insol_ham','koi_max_mult_ev_ham','koi_srad_ham','koi_steff_ham','koi_slogg_ham','koi_sma']
newnames=['period','epoch','rplanet','srad','mes','rstar','tstar','logg','a']

change={}
for i,v in enumerate(oldnames):
    change[v]=newnames[i]
        

newhamops=hamops.rename(columns=change,inplace=False)
#Note this only contains KOIs.

#%%
#Now plot period radius diagram for this catalog.

outplot6="/home/smullall/Kepler/Nexsci/deliveryWorkspace/DR25KOI/plots/pradiusMCMC.png" 
pmp.plotPeriodRadius(newhamops,outplot6, "MCMC", label="OPS",colorkey="score",colorlim=[0,1])  

outplot6b="/home/smullall/Kepler/Nexsci/deliveryWorkspace/DR25KOI/plots/pradius-DVsup.png" 
pmp.plotPeriodRadius(ops,outplot6b, "DVsup", label="OPS",colorkey="score",colorlim=[0,1]) 

outplot7="/home/smullall/Kepler/Nexsci/deliveryWorkspace/DR25KOI/plots/pradiusMCMC-tstar.png" 
pmp.plotPeriodRadius(newhamops,outplot7, "MCMC", label="OPS",colorkey="tstar",colorlim=[4000,6000]) 