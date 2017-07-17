# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 13:13:37 2016

@author: smullall
"""

#Script to start looking at stats coming from DV for cuts.

import dvio
import dr25
import rvIO as rvio
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#%%
def plot2pops(ops,inj,xkey,ykey):
    """
    """
    plt.plot(inj[xkey],(inj[ykey]),'bs')
    plt.plot(ops[xkey],(ops[ykey]),'ro')
    ax=plt.gca()
    ax.set_yscale('log')
    plt.xlabel(xkey)
    plt.ylabel(ykey)    
    
#%%

opsdvData=dvio.loadDvAsDataFrame(dr25.opsDvPath+'/dvOutputMatrix.mat')
opsdvData.index=dr25.tceIdToStr(opsdvData.index)

injdvData1=dvio.loadDvAsDataFrame(dr25.injDvPath1+'/dvOutputMatrix.mat')
injdvData2=dvio.loadDvAsDataFrame(dr25.injDvPath2+'/dvOutputMatrix.mat')
injdvData3=dvio.loadDvAsDataFrame(dr25.injDvPath3+'/dvOutputMatrix.mat')
injdvData4=dvio.loadDvAsDataFrame(dr25.injDvPath4+'/dvOutputMatrix.mat')

frames=[injdvData1,injdvData2,injdvData3,injdvData4]
injdvData=pd.concat(frames)
injdvData.index=dr25.tceIdToStr(injdvData.index)


rvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/RoboVetterOut-OPS.txt'
tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
fedfile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'

opsData=rvio.createAllResults(rvfile,tcefile,fedfile,cumfile)


rvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61710/RoboVetterOut-INJ-PlanetOn.txt'
tcefile='/soc/nfs/so-nfs/DR25/INJ/DATA/TCEs.txt'
injDataOnT=rvio.createRVResults(rvfile,tcefile)

#%%
#Now merge together the robovetter and the interesting columns.
#

bsfa='bootstrap_falseAlarmRate'
bsthresh='planetCandidate.bootstrapThresholdForDesiredPfa'
chi1='planetCandidate.chiSquare1'
chi1dof='planetCandidate.chiSquareDof1'
chi2='planetCandidate.chiSquare2'
chi2dof='planetCandidate.chiSquareDof2'
bsmesstd='planetCandidate.bootstrapMesStd'
chigof='planetCandidate.chiSquareGof'
chigofdof='planetCandidate.chiSquareGofDof'
tcount='planetCandidate.observedTransitCount'
robust='planetCandidate.weakSecondaryStruct.robustStatistic'

#OPS
ops=pd.merge(opsData,opsdvData,left_index=True,right_index=True,how='left')
opslong=(ops['period']>200) & (ops['mes']<8) & (ops['N']==0);

ops['redchi1']=ops[chi1]/ops[chi1dof]
ops['redchi2']=ops[chi2]/ops[chi2dof]
ops['redgof']=ops[chigof]/ops[chigofdof]

#%%
#INJ
inj=pd.merge(injDataOnT,injdvData,left_index=True,right_index=True,how='left')
injkoi=(inj['period']>200) & (inj['mes']<8) & (inj['N']==0);
inj['redchi1']=inj[chi1]/inj[chi1dof]
inj['redchi2']=inj[chi2]/inj[chi2dof]
inj['redgof']=inj[chigof]/inj[chigofdof]



plt.figure()
plt.subplot(221)
plot2pops(ops[opslong],inj[injkoi],'period','mes')
plt.subplot(224)
plot2pops(ops[opslong],inj[injkoi],'period','redchi1')
plt.subplot(222)
plot2pops(ops[opslong],inj[injkoi],'period','redchi2')
plt.subplot(223)
plot2pops(ops[opslong],inj[injkoi],'period','redgof')


plt.figure()
plt.subplot(221)
plot2pops(ops[opslong],inj[injkoi],'mes',bsfa)
plt.subplot(224)
plot2pops(ops[opslong],inj[injkoi],'mes',chi1)
plt.subplot(222)
plot2pops(ops[opslong],inj[injkoi],'mes',chi2)
plt.subplot(223)
plot2pops(ops[opslong],inj[injkoi],'mes',robust)

#%%
plt.figure()
plt.subplot(221)
plot2pops(ops[opslong],inj[injkoi],bsmesstd,bsfa)
plt.subplot(224)
plot2pops(ops[opslong],inj[injkoi],'redchi2',bsfa)
plt.subplot(222)
plot2pops(ops[opslong],inj[injkoi],'redchi1',bsfa)
plt.subplot(223)
plot2pops(ops[opslong],inj[injkoi],bsthresh,bsfa)

#%%
plt.figure()
plt.subplot(221)
plot2pops(ops[opslong],inj[injkoi],bsmesstd,'redchi1')
plt.subplot(224)
plot2pops(ops[opslong],inj[injkoi],'redchi2','redchi1')
plt.subplot(222)
plot2pops(ops[opslong],inj[injkoi],'redgof','redchi1')
plt.subplot(223)
plot2pops(ops[opslong],inj[injkoi],'score','redchi1')