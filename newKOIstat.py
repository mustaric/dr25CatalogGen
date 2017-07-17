# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 15:46:13 2016

@author: smullall
"""

#Script to look at the new KOI population
#


import rvIO as io
import plotRobovetter as prv
import numpy as np
import matplotlib.pylab as plt


rvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61794/RoboVetterOut-OPS.txt'
tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
fedfile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'

opsData=io.createAllResults(rvfile,tcefile,fedfile,cumfile)

opsData['iskoi']=opsData['koi'].notnull()

#Specify not a KOI and does not have the N flag set.
#want=~opsData['iskoi'] & (opsData['N']==0)
#x=np.log10(opsData['period'][want])
#y=np.log10(opsData['mes'][want])

#prv.plot2dHist(x,y,[-.35,2.95],[.83,1.7],nxbins=70,nybins=70,xlabel='log(Period)',ylabel='log(MES)')

fgk=(opsData['logg']>4.0) & (opsData['tstar']>=4000.0) & (opsData['tstar']<7000.0);

wantnew= (~opsData['iskoi']) & (opsData['N']==0)
x=np.log10(opsData['period'][wantnew & fgk])
y=np.log10(opsData['mes'][wantnew & fgk])
plt.figure(figsize=(9,5))
plt.subplot(211)
plt.plot(x,y,'ro')

wantpc= (~opsData['iskoi']) & (opsData['disp']=='PC')
x=np.log10(opsData['period'][wantpc & fgk])
y=np.log10(opsData['mes'][wantpc & fgk])
plt.plot(x,y,'k.')
plt.ylim((.8,1.6))
titstr='New KOIs (%u) and PCs(black) (%i) on FGK dwarfs' % (len(wantnew[wantnew & fgk]),len(wantpc[wantpc & fgk]))
plt.title(titstr)
plt.ylabel('log(MES)')

#%

plt.subplot(212)
wantnew= (~opsData['iskoi']) & (opsData['N']==0)
x=np.log10(opsData['period'][wantnew & ~fgk])
y=np.log10(opsData['mes'][wantnew & ~fgk])
plt.plot(x,y,'ro')

wantpc= (~opsData['iskoi']) & (opsData['disp']=='PC')
x=np.log10(opsData['period'][wantpc & ~fgk])
y=np.log10(opsData['mes'][wantpc & ~fgk])
plt.plot(x,y,'k.')
plt.ylim((.8,1.6))
titstr='New KOIs (%u) and PCs(black) (%i) on rest of stars' % (len(wantnew[wantnew & ~fgk]),len(wantpc[wantpc & ~fgk]))
plt.title(titstr)
plt.xlabel('log(Period)')
plt.ylabel('log(MES)')

#%%
outfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/newKOIs.png'
plt.savefig(outfile)

#%%

fgk=(opsData['logg']>4.0) & (opsData['tstar']>=4000.0) & (opsData['tstar']<7000.0);

wantnew= (~opsData['iskoi']) & (opsData['N']==0)
x=np.log10(opsData['period'][wantnew & fgk])
y=opsData['score'][wantnew & fgk]
plt.figure(figsize=(9,5))
plt.subplot(211)
plt.plot(x,y,'ro')

wantpc= (~opsData['iskoi']) & (opsData['disp']=='PC')
x=np.log10(opsData['period'][wantpc & fgk])
y=opsData['score'][wantpc & fgk]
plt.plot(x,y,'k.')
plt.ylim((0,1))
titstr='New KOIs (%u) and PCs(black) (%i) on FGK dwarfs' % (len(wantnew[wantnew & fgk]),len(wantpc[wantpc & fgk]))
plt.title(titstr)
plt.ylabel('(Score)')

#%

plt.subplot(212)
wantnew= (~opsData['iskoi']) & (opsData['N']==0)
x=np.log10(opsData['period'][wantnew & ~fgk])
y=opsData['score'][wantnew & ~fgk]
plt.plot(x,y,'ro')

wantpc= (~opsData['iskoi']) & (opsData['disp']=='PC')
x=np.log10(opsData['period'][wantpc & ~fgk])
y=opsData['score'][wantpc & ~fgk]
plt.plot(x,y,'k.')
plt.ylim((0,1))
titstr='New KOIs (%u) and PCs(black) (%i) on rest of stars' % (len(wantnew[wantnew & ~fgk]),len(wantpc[wantpc & ~fgk]))
plt.title(titstr)
plt.xlabel('log(Period)')
plt.ylabel('(Score)')

#%%
#Creat histogram of new KOIs and old KOIs that were found but no longer deemed worthy.

rvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/RoboVetterOut-OPS.txt'
tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
fedfile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'

opsData=io.createAllResults(rvfile,tcefile,fedfile,cumfile)

opsData['iskoi']=opsData['koi'].notnull()


newkoiwant=(~opsData['iskoi']) & (opsData['N']==0)
newpcswant=(~opsData['iskoi']) & (opsData['disp']=='PC')
koisnomore=(opsData['iskoi']) & (opsData['N']==1)
pcsnomore=(opsData['iskoi']) & (opsData['dr24disp']=='PC') & (opsData['disp']=='FP')
plt.figure(figsize=(9,8))
plt.plot((opsData['period'][newkoiwant]),(opsData['mes'][newkoiwant]),
         'ro',label='new KOIs')

plt.plot((opsData['period'][koisnomore]),(opsData['mes'][koisnomore]),
         'xb',label='KOIs that now have N==1')
ax=plt.gca()
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('log Period')
plt.ylabel('log MES')
plt.xlim((.5,700))
plt.ylim((5,100000))
outfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/newAndRemovedKOIs.png'
plt.savefig(outfile)

newpcswant=(~opsData['iskoi']) & (opsData['disp']=='PC')
pcsnomore=(opsData['iskoi']) & (opsData['dr24disp']=='PC') & (opsData['disp']=='FP')
plt.figure(figsize=(9,8))
plt.plot((opsData['period'][newpcswant]),(opsData['mes'][newpcswant]),
         'ro',label='PCs on new KOIs')

plt.plot((opsData['period'][pcsnomore]),(opsData['mes'][pcsnomore]),
         'xb',label='dr24-PCs that now are FPs')
ax=plt.gca()
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('log Period')
plt.ylabel('log MES')
plt.xlim((.5,700))
plt.ylim((5,100000))
outfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/newPCsAndRemovedPCs.png'
plt.savefig(outfile)