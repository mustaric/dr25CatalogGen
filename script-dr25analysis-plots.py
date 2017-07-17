# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 15:51:25 2017

@author: smullall
"""

import pandas as p
import numpy as np
import plotRobovetter as prv
import createRVPerfMetrics as crpm
import matplotlib.pyplot as plt
import createRVPerfMetricsPlots as pmp
import rvIO as io

id=62353
ops,inj,inv=io.getAllData(id,tcedata="DATA-SUP")
#Remove those that are banned. -- these were not turned into KOIs and should not be considered for Reliability measures.
w=(~ops.isdr25koi) & (ops.N==0)
ops.drop(ops[w].index)

#%%
#Create list of confirmed exoplanets that failed the robovetter.
#Create a histogram by period and MES

isconf=ops.koi_disposition=='CONFIRMED'
isfp=ops.disp=='FP'


myparamlist=['period','mes','disp','ngood','koi_disposition','dr25koiName','flag','flags']
print ops[isconf & isfp][myparamlist]
ops[isconf & isfp][myparamlist].to_csv('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/confirmed-failures.csv')

plt.figure()
plt.plot(ops[isconf & isfp].mes,ops[isconf & isfp].ngood,'ob')

#%%
#
opssub=(ops.mes<12) & (ops.mes>7)
invsub=(inv.mes<12) & (inv.mes>7)
prv.plot1DReliability(ops[opssub],inv[invsub],'ngood',(2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5),s=[0.0,1.0])
plt.title('MES: 7--12')
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/ntransitReliability2.png')

#%%
ntrans=ops['ngood']
mes=ops['mes']
passed=crpm.passes(ops,(0,1)) & (ops.isdr25koi)
climRange=(0,1)
xBins=(3,4,5,10,100)
yBins=(7,10,20,100)
prv.plotGrid(ntrans,mes,passed,climRange,xBins=xBins,yBins=yBins)

#%%

new=(ops.disp=='PC') & (~ops.isdr24koi) & (ops.isdr25koi)
len(new[new])
outname='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/periodRadiusNew.png'
pmp.plotPeriodRadius(ops[new],outname, 'New PCs\n Score Cut (0.0,1.0)',colorkey="score",s=[0,1.0])


#%%
#Reliability of the catalog by showing points of the inv/ss1 over the candidates.
#
want=(ops.tstar>=4000) & (ops.logg>4.0) & (ops.tstar<=7000)
iwant=(inv.tstar>=4000) & (inv.logg>4.0) & (inv.tstar<=7000)
pcs=ops.disp=='PC'
ipcs=inv.disp=='PC'
w=want & pcs
iw=iwant & ipcs

plt.figure()

plus, =plt.plot(ops[w].period,ops[w].rplanet,'k+',ms=10,lw=2,label="")
circle, =plt.plot(ops[w].period,ops[w].rplanet,'ko',fillstyle="none",mec='k',markeredgewidth=1.3,ms=10)

dot, =plt.plot(inv[iw].period,inv[iw].rplanet,'r.',ms=9)
plt.xlim((150,525))
plt.ylim((0.5,6))
plt.xlabel('Period (days)',fontsize=15)
plt.ylabel('Planet Radius (Earth radii)',fontsize=15)

pa=(384.846,1.21)
pb=(267.281,1.41)
gcircle,= plt.plot(pb[0],pb[1],'ko',fillstyle="none",mec='g',markeredgewidth=2.0,ms=10)
plt.plot(pa[0],pa[1],'ko',fillstyle="none",mec='g',markeredgewidth=2.0,ms=10)

plt.legend([(circle,plus), dot, gcircle], ["DR25 Candidates", "False Alarms", "This Proposal"],loc="lower right",fontsize=14,framealpha=0.5)



plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/candidatesFalseAlarms-rockyHZ.png')

#%%
#Calculate a model to represent the 
import plotRobovetter as prv
rpops=ops.drop(ops[np.isnan(ops.rplanet)].index)
rpinv=inv.drop(inv[np.isnan(inv.rplanet)].index)
want=(rpops.tstar>=4500) & (rpops.logg>4.0) & (rpops.tstar<=7000)
iwant=(rpinv.tstar>=4500) & (rpinv.logg>4.0) & (rpinv.tstar<=7000)
pwant=rpops.period>200;
ipwant=rpinv.period>200
plt.figure()
prv.plot1DReliability(rpops[want & pwant],rpinv[iwant & ipwant],'rplanet',(0.5,2.5,3,4,5,6,8),s=[0.0,1.0])
pwant=rpops.rplanet<2
ipwant=rpinv.rplanet<2
plt.figure()
prv.plot1DReliability(rpops[want & pwant],rpinv[iwant & ipwant],'period',(0,100,200,300,400,500),s=[0.0,1.0])

#%%
from scipy import interpolate
xbins=(100,200,300,400,500)
ybins=(0.5,3,5,7)
plt.figure(figsize=(8.5,4))

R=prv.plotOnlyReliability(rpinv,rpops,xmetric='period',ymetric='rplanet',xBins=xbins,yBins=ybins,drange=(1,99),key='',limit=0.0,s=[0.0,1.0])
x=np.arange(50,550,2)
y=np.arange(0.1,7,.1)
xx,yy = np.meshgrid(x,y)

xv=(150,250,350,450)
yv=(1.75,3.5,5.5)
f=interpolate.interp2d(xv,yv,R,kind='linear')

Rnew=f(x,y)
plt.figure()
CS=plt.contour(xx,yy,Rnew,linewidths=3)
plt.clabel(CS,fmt='%2.0f%%',inline=1,fonsize=11)
line,=plt.plot([],[],'b-',lw=3)

want=(ops.tstar>=4000) & (ops.logg>4.0) & (ops.tstar<=7000)
iwant=(inv.tstar>=4000) & (inv.logg>4.0) & (inv.tstar<=7000)
pcs=ops.disp=='PC'
ipcs=inv.disp=='PC'
w=want & pcs
iw=iwant & ipcs

plus, =plt.plot(ops[w].period,ops[w].rplanet,'k+',ms=10,lw=2,label="")
circle, =plt.plot(ops[w].period,ops[w].rplanet,'ko',fillstyle="none",mec='k',markeredgewidth=1.3,ms=10)

dot, =plt.plot(inv[iw].period,inv[iw].rplanet,'r.',ms=8)

plt.xlabel('Period (days)',fontsize=15)
plt.ylabel('Planet Radius (Earth radii)',fontsize=15)

pa=(384.846,1.21)
pb=(267.281,1.41)
gcircle,= plt.plot(pb[0],pb[1],'ko',fillstyle="none",mec='g',markeredgewidth=2.0,ms=12)
plt.plot(pa[0],pa[1],'ko',fillstyle="none",mec='g',markeredgewidth=2.0,ms=10)
plt.annotate('Kepler-452b',xy=pa,xytext=(368,1.0),fontsize=12,color='g')
plt.annotate('Kepler-62f',xy=pb,xytext=(258,1.5),fontsize=12,color='g')

plt.legend([(circle,plus), gcircle,dot, line], ["DR25 Candidates","This Proposal", "False Alarms", "Reliability"],loc="lower left",fontsize=14,framealpha=0.5,numpoints=1)

plt.xlim((150,500))
plt.ylim((0.2,4.1))


plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/candidatesFalseAlarms-rockyHZ-Rel2.png')
#%%
#
#Creating a table of the best candidates
