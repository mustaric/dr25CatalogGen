# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 11:59:17 2017

Evaluate Kelsey's Fits.
compare to DV.
Look for errant values.

@author: smullall
"""


import numpy as np
import pandas as p
import rvIO as io
import matplotlib.pyplot as plt

id=62353
mcmcFile='/home/smullall/Kepler/RoboVetter/DR25/mcmc/bftable_022117.dat'

ops,inj,inv=io.getAllData(id)

mcmcdata=p.read_csv(mcmcFile,delimiter=',')

koiNames =map( lambda x: "K%08.2f" % (x), mcmcdata.KOI)
mcmcdata["koiName"]=koiNames
mcmcdata.set_index("koiName",inplace=True)

#
#Merge to ops for comparison. Names being the same is a problem.

both=p.merge(ops,mcmcdata,how='left', left_on='dr25koiName',right_index=True,suffixes=('','_mcmc'))


#Merge to the cumulative table...keeping the whole cumulative table around.
svnfile='/home/smullall/Kepler/cumulativeTable/dr25/dr25-status.csv'
svndata=io.readcumulativeTable2(svnfile)

allkois=p.merge(svndata,mcmcdata,how='left',left_index=True,right_index=True)


#%%
#Compare the DV fits with the Kelsey Fits.

isdr25koi=both.isdr25koi

#How many KOIs have no period?

good=np.isreal(both.Period_mcmc)
negerr=both.Periodsig<0
zeroerr=both.Periodsig==0

print "Not Real Period"
print both[~good & isdr25koi].index

print "Negative Error"
print both[negerr & isdr25koi]

print "Zero Error"
print both[zeroerr & isdr25koi]

#%%
#do we have fits for all KOIs.  Len should be 9584

invisible=allkois.dr25disp=='Invisible'
donotdel=allkois.dr25disp=='DoNotDeliver'

hasperiod=allkois.Period>0

good=hasperiod & ~invisible & ~donotdel

len(good[good])

#%%
#compare periods and epochs.

cat=both[isdr25koi & ~zeroerr]

plt.figure(figsize=(12,8))
plt.subplot(221)
plt.xlabel('MCMC Period error')
plt.ylabel('DV Period - MCMC Period')
plt.plot(cat.Periodsig,cat.period-cat.Period_mcmc,'r.')
plt.ylim([-.1,.1])
plt.xlim([0,.05])

plt.subplot(222)
plt.xlabel('log MCMC Tdepth')
plt.ylabel('log(DV depth - MCMC Tdepth)')
plt.plot(np.log10(cat.Tdepth),(cat.dv_depth-cat.Tdepth),'b.')

plt.subplot(223)
plt.xlabel('log MCMC Period')
plt.ylabel('dv_impact - MCMC b')
plt.plot(np.log10(cat.Period_mcmc),cat.Impact - cat.b,'g.')

plt.subplot(224)
plt.xlabel('log MCMC Period')
plt.ylabel('dv_duration -MCMC Tdur')
plt.plot(np.log10(cat.Period_mcmc),cat.duration - cat.Tdur,'m.')
