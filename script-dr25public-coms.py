# -*- coding: utf-8 -*-
"""
Created on Thu May 25 10:10:40 2017

Plots for the Coms disucssion.
Also a playground before it ends up in the public release script for the paper.

@author: smullall
"""

import publicIO as io
import pandas as p
import numpy as np
import plotRobovetter as prv
import createRVPerfMetrics as crpm
import matplotlib.pyplot as plt
import createRVPerfMetricsPlots as pmp
import publicPaperPlots as ppp

#dr25=io.DR25_IO(ddir="/Users/sthompson/kepler/DR25/public0512/KeplerPublicProducts/")
#For my linux box.
dr25=io.DR25_IO(ddir="/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/publicData/")

koi=dr25.loadKOITable()
ops=dr25.loadOps()
inv=dr25.loadInv()
scr=dr25.loadScr()
fas=dr25.loadBothFA()
inj=dr25.loadInj()

ops['mes']=ops.MES
fas['mes']=fas.MES
inj['mes']=inj.MES
#
outdir='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/Coms/'

#%%
#Plot a histogram, coarse, of the new and the old planet candidates
cumkoi=dr25.loadCumKOITable()
real=~np.isnan(cumkoi.Rp)
newones=cumkoi[cumkoi.New]
prevones=cumkoi[~cumkoi.New]
longP=cumkoi.period>100
pcs=(cumkoi.disp=='CANDIDATE') & longP

bins=np.linspace(0.5,15,num=25)
bins=[0,1.3,2,6,20]

data=[cumkoi[~cumkoi.New & pcs & real].Rp,cumkoi[cumkoi.New & pcs & real].Rp]
h1,b1=np.histogram(cumkoi[~cumkoi.New & pcs & real].Rp,bins=bins)
h2,b2=np.histogram(cumkoi[cumkoi.New & pcs & real].Rp,bins=bins)


plt.figure()
#data=[cumkoi[cumkoi.New & pcs],cumkoi[~cumkoi.New & pcs]]
#plt.hist(data,bins=bins,histtype='bar',stacked=True,color=['blue','yellow'])

ind=np.arange(len(h1))
p1=plt.bar(ind,h1,0.5,color='blue',label='Prev. Known')
p2=plt.bar(ind,h2,0.5,color='yellow',bottom=h1,label='New Candidates')
plt.xticks(ind+.25,('Earth','Super\nEarth','Neptune','Jupiter'))

plt.legend(loc='best')
plt.title('Orbital Periods longer than 100 days')

plt.savefig('%s/coms-bigBinsNewCandidatesLongP.png' % outdir)

#%%
#For DR25 catalog what we found
#Then correct for completeness.
real=~np.isnan(koi.Rp)
pcs=koi.disp=='CANDIDATE'
longP=(koi.period>100.0) & (koi.period<400)

bins=[0,1.3,2,6,20]
h,b=np.histogram(koi[real & pcs & longP].Rp,bins=bins)
#NOTE-- THIS CORRECTION IS WRONG...JUST AN EXAMPLE
hcor=h*[100,8,1.5,1.0]
ind=np.arange(len(bins)-1)

plt.figure()
plt.subplot(211)
p1=plt.bar(ind,h,0.5,color='blue')
plt.xticks(ind+.25,('Earth','Super Earth','Neptune','Jupiter'),fontsize=11)
plt.yticks([])
plt.ylim([0,np.max(hcor)+25])
plt.annotate('Very Hard\n  to Find',xy=(0.05,40))
plt.annotate('Hard\nto Find',xy=(1.1,60))
plt.annotate('Kinda Easy \nto Find',xy=(2.0,180))
plt.annotate('Easy\n to Find', xy=(3.05,80))
plt.title('DR25 catalog for Period 100--400 days \n Number of Planets Found')
plt.subplot(212)
p1=plt.bar(ind,hcor,0.5)
plt.xticks(ind+.25,('Earth','Super Earth','Neptune','Jupiter'))
plt.yticks([])
plt.ylim([0,np.max(hcor)+25])

plt.title('Actual Population of Transiting Planets')

plt.savefig('%s/coms-bigBinsLongPeriodCorr.png' % outdir)


#%%
#Plot the cumulative table candidates for HZ.
#
#This just shows that we have a lot of crap left in the HZ.
plt.figure(figsize=(8.5,11))

cumkoi['kepoi_name']=cumkoi.index
for v in cumkoi[np.isnan(cumkoi.koi_score)].index:
    cumkoi.set_value(v,'koi_score',1.0)


cumhzframe=ppp.HZplot(cumkoi,scut=0.5,addConf=False,mfactor=10,ebars=False)



#%%
#Plot the cumulative table new and old as period vs. radius.
#Just DR25 candidates with new ones in yellow
koipcs=(koi.disp=='CANDIDATE') & (koi.koi_score>0.5) #& (koi.koi_duration<15)
cumpc=(cumkoi.disp=='CANDIDATE') & (cumkoi.koi_score>0.5) #& (cumkoi.koi_duration<15)
logrp=np.log10(koi[koipcs].Rp)
logP=np.log10(koi[koipcs].period)

newlogrp=np.log10(cumkoi[cumpc& cumkoi.New].Rp)
newlogP=np.log10(cumkoi[cumpc & cumkoi.New].period)

plt.figure()
plt.plot(logP,logrp,'bo',label='Previous Candidates',mew=.1)
plt.plot(newlogP,newlogrp,'ro',label='New Candidates',mew=.1)
plt.xlim([-.5,2.8])
plt.ylim([-0.55,1.5])
plt.xlabel('Period (days)')
plt.ylabel('Planet Radius')
lab=np.array([0.4,1,3,10,100,300])
locat=np.log10(lab)
plt.xticks(locat,lab.astype(str))
lab=np.array([0.4,1.0,2.0,4.0,10.0])
locat=np.log10(lab)
plt.yticks(locat,lab.astype(str))
plt.legend(loc='lower right',numpoints=1)

plt.savefig('%s/coms-dr25periodradius-scatter.png' % outdir)

#%%
#Plot the cumulative table new and old as period vs. radius.
#Just DR25 candidates with new ones in red
#Overplot shaded regions of low reliability and low completeness.
#
koipcs=(koi.disp=='CANDIDATE') & (koi.koi_score>0.6) #& (koi.koi_duration<15)
cumpc=(cumkoi.disp=='CANDIDATE') & (cumkoi.koi_score>0.6) #& (cumkoi.koi_duration<15)
logrp=np.log10(koi[koipcs].Rp)
logP=np.log10(koi[koipcs].period)

newlogrp=np.log10(cumkoi[cumpc& cumkoi.New].Rp)
newlogP=np.log10(cumkoi[cumpc & cumkoi.New].period)

plt.figure()
plt.plot(logP,logrp,'bo',label='Previous Candidates',mew=.1)
plt.plot(newlogP,newlogrp,'ro',label='New Candidates',mew=.1)
plt.xlim([-.5,2.8])
plt.ylim([-0.55,1.5])
plt.xlabel('Period (days)')
plt.ylabel('Planet Radius')
lab=np.array([0.4,1,3,10,100,300])
locat=np.log10(lab)
plt.xticks(locat,lab.astype(str))
lab=np.array([0.4,1.0,2.0,4.0,10.0])
locat=np.log10(lab)
plt.yticks(locat,lab.astype(str))
plt.legend(loc='lower right',numpoints=1)

plt.add_patch

plt.savefig('%s/coms-dr25periodradius-shaded.png' % outdir)


#%%
#Calculating up some Stats
#This is a broad HZ for the catalog completeness and effectiveness
#Used to Verfy D.2 Tcert documentation Numbers
metric1='Sp'
range1=[0.2,1.8]
metric2='Rp'
range2=[0.4,1.8]

#metric2='mes'
#range2=[7,10]

opshz=ppp.inBox(ops,metric1,range1,metric2,range2)
invhz=ppp.inBox(inv,metric1,range1,metric2,range2)
scrhz=ppp.inBox(scr,metric1,range1,metric2,range2)
fashz=ppp.inBox(fas,metric1,range1,metric2,range2)
injhz=ppp.inBox(inj,metric1,range1,metric2,range2)

scl=0.5
templim=5000
glim=4.0
C=ppp.fracPass(inj[injhz & (inj.Ts>templim) & (inj.logg>glim)],s=[scl,1.0])
Einv=1-ppp.fracPass(inv[invhz & (inv.Ts>templim) & (inv.logg>glim)],s=[scl,1.0])
Escr=1-ppp.fracPass(scr[scrhz& (scr.Ts>templim) & (scr.logg>glim)],s=[scl,1.0])
Efas=1-ppp.fracPass(fas[fashz& (fas.Ts>templim) & (fas.logg>glim)],s=[scl,1.0])
R=pmp.estimateReliability(fas[fashz & (fas.Ts>templim) & (fas.logg>glim)],ops[opshz& (ops.Ts>templim) & (ops.logg>glim)],s=[scl,1.0])

numops=len(ops[opshz& (ops.Ts>templim)  & (ops.logg>glim) & (ops.score>scl) & (ops.disp=='PC')])
Tot=numops*R[0]/C
print numops,Tot
print ops[opshz& (ops.Ts>templim) & (ops.logg>glim)&(ops.score>scl) & (ops.disp=='PC')][['mes','period','ntrans','Rp','Sp','Ts','score']]
print len(fas[fashz& (fas.Ts>templim) & (fas.logg>glim)&(fas.score>scl) & (fas.disp=='PC')][['mes','period','ntrans','Rp','Sp','Ts','score']])
print np.mean(inj[injhz& (inj.Ts>templim) & (inj.logg>glim)][['Expected_MES']])
print np.mean(ops[opshz& (ops.Ts>templim) & (ops.logg>glim)][['period']])

print "\def \SRcompleteness {%5.2f}" % (C*100)
print "\def \SRinveffect {%5.2f}" % (Einv*100)
print "\def \SRscreffect {%5.2f}" % (Escr*100)
print "\def \SRfaseffect {%5.2f}" % (Efas*100)
print "\def \SRreliability {%5.2f}" % (R[0]*100)
print C,Einv,Escr, Efas,R

#%%
#Habitable Zone
cumkoi=dr25.loadCumKOITable()
koi=dr25.loadKOITable()
#Plot the HZ candidates
plt.figure(figsize=(11,11))
hzframe=ppp.HZplotSimple(koi, scut=0.5,pradcut=1.8, mfactor=23, addConf=False, ylimits=(3200,6400), xlimits=(2.5,0.0),ebars=False)
plt.savefig('%s/coms-hzTstarInsolsimple.png' % outdir)


#%%
import hzCalc
plt.figure()
small=(koi.koi_prad<1.5) & (koi.koi_score>0.0) & (koi.koi_pdisposition=='CANDIDATE')
shperiod=(koi.period<=200) & (koi.koi_steff<7000)
plt.plot(np.log10(koi[small].koi_insol), koi[small].koi_steff,'ob', label='Not Counted')
plt.plot(np.log10(koi[small & shperiod].koi_insol), koi[small & shperiod].koi_steff,'go',label='Counted',mew=0.1)
plt.legend(loc='best',numpoints=1)
plt.xlim(4.2,-1.2)
plt.ylim(2500,7200)
plt.xlabel('Period (days)',fontsize=14)
plt.xlabel('Insolation Flux',fontsize=14)
plt.ylabel('Stellar Temperature (K)',fontsize=14)

#Draw a HZ across the plot.
teff=np.linspace(2000,7000,50)
hzBorders=hzCalc.hzBorder(teff)
Shot=hzBorders[:,0]
Scool=hzBorders[:,3]
plt.plot(np.log10(Shot),teff,'-b',zorder=2)
plt.plot(np.log10(Scool),teff,'-r',zorder=2)
    
plt.title('')
#plt.annotate('\n2.5\nplanets\n',xy=(30,2750),alpha=0.6,\
#     color='k',fontsize=14,horizontalalignment='center',bbox=dict(facecolor='red', alpha=0.8))
#plt.annotate('\n\n\n       1/4      \n  planets\n   per star   \n\n',xy=(50,4450),alpha=0.6,\
#     color='k',fontsize=14,horizontalalignment='center',bbox=dict(facecolor='yellow', alpha=0.8))
#plt.annotate('\n\n\n\n      1/3 -- 1.9     \n  planets per star      \n\n\n\n',xy=(190,4250),alpha=0.6,\
#     color='k',fontsize=14,horizontalalignment='center',bbox=dict(facecolor='orange', alpha=0.85))

plt.savefig('%s/coms-occurRateStudies.png' % outdir)

#%%
import hzCalc
plt.figure()
small=(koi.koi_prad<4) & (koi.koi_score>0.5) & (koi.koi_pdisposition=='CANDIDATE')

plt.plot(koi[small].koi_period, koi[small].koi_steff,'ob')

#plt.legend(loc='best',numpoints=1)
plt.xlim(0,500)
plt.ylim(2500,7200)
plt.xlabel('Period (days)',fontsize=14)

plt.ylabel('Stellar Temperature (K)',fontsize=14)

#Draw a HZ across the plot.
teff=np.linspace(2000,7000,50)
hzBorders=hzCalc.hzBorder(teff)
Shot=hzBorders[:,0]
Scool=hzBorders[:,3]
plt.plot(np.log10(Shot),teff,'-b',zorder=2)
plt.plot(np.log10(Scool),teff,'-r',zorder=2)
    
plt.title('')
plt.annotate('\n2.5\nplanets\n',xy=(30,2750),alpha=0.6,\
     color='k',fontsize=14,horizontalalignment='center',bbox=dict(facecolor='red', alpha=0.8))
plt.annotate('\n\n\n\n       0.25      \n  planets\n   per star   \n\n',xy=(50,4450),alpha=0.6,\
     color='k',fontsize=14,horizontalalignment='center',bbox=dict(facecolor='yellow', alpha=0.8))
plt.annotate('\n\n\n\n\n      0.3 -- 1.9     \n  planets per star      \n\n\n\n',xy=(190,4250),alpha=0.6,\
     color='k',fontsize=14,horizontalalignment='center',bbox=dict(facecolor='orange', alpha=0.85))

plt.savefig('%s/coms-occurRateStudies2.eps' % outdir)