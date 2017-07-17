# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 13:54:34 2016

@author: smullall
"""

import rvIO as io
import plotRobovetter as prv
import numpy as np
import matplotlib.pylab as plt
import multiCuts as mc
#%%
rvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/RoboVetterOut-INJ-PlanetOn.txt'
tcefile='/soc/nfs/so-nfs/DR25/INJ/DATA/TCEs.txt'
injDataOnT=io.createRVResults(rvfile,tcefile)

allrvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/RoboVetterOut-INJ-All.txt'
injDataAll=io.createRVResults(allrvfile,tcefile)

pnStart=3

injDataOnT['MultiCut'] = mc.cutMultis(injDataOnT,injDataAll,pnStart)
injCut=(injDataOnT['N']==0) & injDataOnT['MultiCut']

#MultiPC cut 0==pass
#MultiCut=np.zeros((len(injDataOnT['period']),1))
#
#for i,tce in enumerate(injDataOnT.index):
#    #print i,tce,injDataOnT.loc[tce]['kic']
#    if injDataOnT.loc[tce]['pn'] >= pnStart:
#        #Find those with the same kic number
#        kic=injDataOnT.loc[tce]['kic']
#        pn=injDataOnT.loc[tce]['pn']
#        samekic=injDataAll['kic']==kic
#        lowerpn=injDataAll['pn']<pn
#        Ns=injDataAll['N'][lowerpn & samekic]
#        if len(Ns)>0:
#            if ~np.any(Ns==0):
#                MultiCut[i]=1
#            
#injDataOnT['MultiCut']=MultiCut

#%%
bins=[1,2,3,4,5,6,7,8,9,10]
plt.figure()
plt.subplot(211)
fullhist=plt.hist(injDataOnT['pn'],bins)
ax=plt.gca()
ax.set_yscale('log')



cut=injDataOnT['MultiCut']
cuthist=plt.hist(injDataOnT['pn'][cut],bins,color='red')
plt.ylabel('Number of Injected TCEs')
plt.xlabel('Planet Number')
plt.title('Injected TCEs by Planet Number and proposed Cuts (red)')

plt.subplot(212)
plt.plot(bins[0:9],cuthist[0]/fullhist[0],'o--')
plt.plot(bins[0:9],cuthist[0]/injDataOnT['period'].count(),'rs-')
ax=plt.gca()
ax.set_yscale('log')
plt.ylim((0,1))

fileout='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/multi-cuts-INJ.png'
plt.savefig(fileout)

#%%
#Do this for OPS
rvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/RoboVetterOut-OPS.txt'
tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
opsData=io.createRVResults(rvfile,tcefile)
pnStart=3
opsData['MultiCut']=mc.cutMultis(opsData,opsData,pnStart)

#%%
bins=[1,2,3,4,5,6,7,8,9,10]
plt.figure(figsize=(11,7))
plt.subplot(311)
fullhist=plt.hist(opsData['pn'],bins)
ax=plt.gca()
ax.set_yscale('log')

cut=opsData['MultiCut']
cuthist=plt.hist(opsData['pn'][cut],bins,color='red')
plt.ylabel('Number of Injected TCEs')
plt.xlabel('Planet Number')
plt.title('Inv TCEs (blue) by Planet Number and proposed Cuts (red)')


plt.subplot(312)
plt.plot(bins[0:9],cuthist[0]/fullhist[0],'o--',label="fraction of bin")
plt.plot(bins[0:9],cuthist[0]/opsData['period'].count(),'rs-',label="fraction of total")
ax=plt.gca()
ax.legend(loc="center right")
#ax.set_yscale('log')
plt.ylim((0,1))


koisCut=(opsData['N']==0) & (opsData['MultiCut'])
len(koisCut[koisCut])
opskois=(opsData['N']==1)

plt.subplot(313)
plt.hist(opsData['pn'][opskois],bins,color='blue')
plt.hist(opsData['pn'][koisCut],bins,color='red')
ax=plt.gca()
ax.set_yscale('log')
plt.title('invKOIs and proposed cut invKOIs (red)')
plt.ylabel('number of KOIs')

fileout='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/multi-cuts-OPS.png'
plt.savefig(fileout)

#%%
#Do this for INV
rvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/RoboVetterOut-INV.txt'
tcefile='/soc/nfs/so-nfs/DR25/INV/DATA/TCEs.txt'
invData=io.createRVResults(rvfile,tcefile)
pnStart=3
invData['MultiCut']=mc.cutMultis(invData,invData,pnStart)

#%%
bins=[1,2,3,4,5,6,7,8,9,10]
plt.figure(figsize=(11,7))
plt.subplot(311)
fullhist=plt.hist(invData['pn'],bins)
ax=plt.gca()
ax.set_yscale('log')

cut=invData['MultiCut']
cuthist=plt.hist(invData['pn'][cut],bins,color='red')
plt.ylabel('Number of Injected TCEs')
plt.xlabel('Planet Number')
plt.title('Ops TCEs by Planet Number and proposed Cuts (red)')


plt.subplot(312)
plt.plot(bins[0:9],cuthist[0]/fullhist[0],'o--',label="fraction of bin")
plt.plot(bins[0:9],cuthist[0]/invData['period'].count(),'rs-',label="fraction of total")
ax=plt.gca()
ax.legend(loc="center right")
#ax.set_yscale('log')
plt.ylim((0,1))


invkoisCut=(invData['N']==0) & (invData['MultiCut'])
len(invkoisCut[invkoisCut])
invkois=(invData['N']==1)

plt.subplot(313)
plt.hist(invData['pn'][invkois],bins,color='blue')
plt.hist(invData['pn'][invkoisCut],bins,color='red')
ax=plt.gca()
ax.set_yscale('log')
plt.title('InvKOIs and proposed multi-cut InvKOIs')
plt.ylabel('number of invKOIs')

fileout='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/multi-cuts-INV.png'
plt.savefig(fileout)

#%%

plt.figure(figsize=(10,10))
plt.plot(invData['period'][invkoisCut],invData['mes'][invkoisCut],'bs')
plt.plot(opsData['period'][koisCut],opsData['mes'][koisCut],'m.')
plt.plot(injDataOnT['period'][injCut],injData['mes'][injCut],'go')
ax=plt.gca()
#ax.set_yscale('log')
ax.set_xscale('log')
plt.ylabel('mes')
plt.ylim((6,20))
plt.xlim((0,700))
plt.xlabel('Period')
plt.ylabel('MES')
plt.title('KOIs cut for Inv(blue square) and Ops(magenta dot) using Multi cut')

fileout='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/multi-cutKois-ops-inv.png'
plt.savefig(fileout)


#%%
#Examine the Reliability if we only consider Inverted things that were found
#For pn=1 or 2.
import statsRobovetter as srv

rvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/RoboVetterOut-INV.txt'
tcefile='/soc/nfs/so-nfs/DR25/INV/DATA/TCEs.txt'
invData=io.createRVResults(rvfile,tcefile)


rvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/RoboVetterOut-OPS.txt'
tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
opsData=io.createRVResults(rvfile,tcefile)

pnlimit=1

opswant=srv.tcesWithLowPn(opsData,pnlimit)
invwant=srv.tcesWithLowPn(invData,pnlimit)

prv.plotReliability(invData[invwant],opsData[opswant])


rvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/RoboVetterOut-INJ-PlanetOn.txt'
tcefile='/soc/nfs/so-nfs/DR25/INJ/DATA/TCEs.txt'
injDataOnT=io.createRVResults(rvfile,tcefile)

allrvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/RoboVetterOut-INJ-All.txt'
injDataAll=io.createRVResults(allrvfile,tcefile)


injwant=srv.tcesWithLowPn(injDataAll,pnlimit)


