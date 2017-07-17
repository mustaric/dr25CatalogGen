# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 16:57:00 2017
Play with public data

@author: sthompson
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
#%%
#outdir='/Users/sthompson/kepler/DR25/public0512/KeplerPublicProducts/plots'
outdir='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/Paper/'

#%%
#GET SUMMARY CATALOG NUMBERS
cumkoi=dr25.loadCumKOITable()
koi=dr25.loadKOITable()

C=ppp.fracPass(inj,s=[0,1.0])
Einv=1-ppp.fracPass(inv,s=[0,1.0])
Escr=1-ppp.fracPass(scr,s=[0,1.0])
Efas=1-ppp.fracPass(fas,s=[0,1.0])
R=pmp.estimateReliability(fas,ops)

print "\def \\nkois {%u}" % len(koi)
print "\def \\ncand {%u}" % len(koi[koi.disp=='CANDIDATE'])
print "\def \\newkois {%u}"  % len(cumkoi[cumkoi.New])
print "\def \\newcand {%u}" % len(cumkoi[cumkoi.New & (cumkoi.disp=='CANDIDATE')])
print "\def \\completeness {%5.2f}" % (C*100)
print "\def \\reliability {%5.2f}" % (R[0]*100)


#%%
#Get Reliability for different boxes.
#
metric1='period'
range1=[200,500]
metric2='mes'
range2=[7,10]
wops=ppp.inBox(ops,metric1,range1,metric2,range2)
wfas=ppp.inBox(fas,metric1,range1,metric2,range2)
winv=ppp.inBox(inv,metric1,range1,metric2,range2)
wscr=ppp.inBox(scr,metric1,range1,metric2,range2)

Rfas=pmp.estimateReliability(fas[wfas],ops[wops])
Rinv=pmp.estimateReliability(inv[winv],ops[wops])
Rscr=pmp.estimateReliability(scr[wscr],ops[wops])
print(Rfas)
print(Rinv)
print(Rscr)

#%%
#Randomly select half of the false alarm sets 100 times.
npts=100
rels=np.zeros(npts)
for i in np.arange(0,npts,1):
    wrand=(np.round(np.random.rand(len(fas)))==1)
    newfas=fas[wrand]
    wfas=ppp.inBox(newfas,metric1,range1,metric2,range2)
    R=pmp.estimateReliability(newfas[wfas],ops[wops])
    rels[i]=R[0]

print np.mean(rels)
print np.std(rels)    
    
#%%
#TCE HISTOGRAMS
#Plot Histogram of the obsTCE list and compare to DR24
#
outname=outdir + "fig-obstcePeriods.png"

dr24tcefile="/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/publicData/other/Q1Q17TCEs.txt"
colnames=('tce','kic','pn','period','epoch','mes','depth','duration','rplanet','rstar','tstar','a','rprstar','arstar','snr','teq','secmees','secphase','posmes','posphase','mesmad')
dr24ops=p.read_csv(dr24tcefile,skiprows=5,names=colnames,delim_whitespace=True,index_col='tce')

b=np.linspace(-0.301,2.83,180); #Bins

plt.figure(figsize=(8.5,4))

plt.hist(np.log10(ops.period),histtype='step',bins=b,label='DR25 obsTCE',color='black',lw=2.5,zorder=2)
plt.hist(np.log10(dr24ops.period),histtype='step',bins=b,label='DR24 obsTCE',color='green',lw=1.5,zorder=1)
plt.xlim([-.4,2.9])
plt.legend(loc='upper left',fontsize=13)
plt.ylabel('Counts',fontsize=14)
plt.savefig(outname)

#%%
#Calculating up some Stats
#This is a broad HZ for the catalog completeness and effectiveness
#Used to Verfy D.2 Tcert documentation Numbers
metric1='period'
range1=[200,400.0]
metric2='Rp'
range2=[0.5,2.0]

opshz=ppp.inBox(ops,metric1,range1,metric2,range2)
invhz=ppp.inBox(inv,metric1,range1,metric2,range2)
scrhz=ppp.inBox(scr,metric1,range1,metric2,range2)
fashz=ppp.inBox(fas,metric1,range1,metric2,range2)
injhz=ppp.inBox(inj,metric1,range1,metric2,range2)

scl=0.5
templim=4000
glim=4.0
C=ppp.fracPass(inj[injhz & (inj.Ts>templim) & (inj.logg>glim)],s=[scl,1.0])
Einv=1-ppp.fracPass(inv[invhz & (inv.Ts>templim) & (inv.logg>glim)],s=[scl,1.0])
Escr=1-ppp.fracPass(scr[scrhz& (scr.Ts>templim) & (scr.logg>glim)],s=[scl,1.0])
Efas=1-ppp.fracPass(fas[fashz& (fas.Ts>templim) & (fas.logg>glim)],s=[scl,1.0])
R=pmp.estimateReliability(fas[fashz& (fas.Ts>templim) & (fas.logg>glim)],ops[opshz& (ops.Ts>templim) & (ops.logg>glim)],s=[scl,1.0])

numopspc=len(ops[opshz& (ops.Ts>templim)  & (ops.logg>glim) & (ops.score>scl) & (ops.disp=='PC')])
numopsfp=len(ops[opshz& (ops.Ts>templim)  & (ops.logg>glim) & ((ops.score<scl) | (ops.disp=='FP'))])
numfaspc=len(fas[fashz& (fas.Ts>templim) & (fas.logg>glim) & (ops.score>scl) & (ops.disp=='PC')]) 
numfas=len(fas[fashz& (fas.Ts>templim) & (fas.logg>glim)])
sigR=pmp.reliabiltiyErrorBinomial(R[0],Efas,numopsfp,numopspc,numfas,numfaspc)

Tot=numopspc*R[0]/C
print numopspc,Tot
#print ops[opshz& (ops.Ts>templim) & (ops.logg>glim)&(ops.score>scl) & (ops.disp=='PC')][['mes','period','ntrans','Rp','Sp','Ts']]
#print (fas[fashz& (fas.Ts>templim) & (fas.logg>glim)&(fas.score>scl) & (fas.disp=='PC')][['mes','period','ntrans','Rp','Sp','Ts']])


print "\def \SRcompleteness {%5.2f}" % (C*100)
print "\def \SRinveffect {%5.2f}" % (Einv*100)
print "\def \SRscreffect {%5.2f}" % (Escr*100)
print "\def \SRfaseffect {%5.2f}" % (Efas*100)
print "\def \SRreliability {%5.2f +- %5.2f}" % (R[0]*100, sigR[0]*100)
print C,Einv,Escr, Efas,R


#%%
outname=outdir + "fig-simulTcePeriods.png"
b=np.linspace(-0.301,2.83,180); #Bins
plt.figure(figsize=(8.5,10))

plt.subplot(311)
plt.hist(np.log10(inv.period),histtype='step',bins=b,label='invTCE',color='red',lw=2.4,zorder=2)
plt.hist(np.log10(ops.period),histtype='step',bins=b,label='DR25 obsTCE',color='black',lw=1.8,zorder=1)
plt.legend(loc='upper left',fontsize=12)
plt.xlim([-.4,2.9])

plt.subplot(312)
plt.hist(np.log10(scr.period),histtype='step',bins=b,label='scrTCE',color='green',lw=2.4,zorder=2)
plt.hist(np.log10(ops.period),histtype='step',bins=b,label='DR25 obsTCE',color='black',lw=1.8,zorder=1)
plt.legend(loc='upper left',fontsize=12)
plt.xlim([-.4,2.9])
plt.ylabel('Counts',fontsize=14)

plt.subplot(313)

plt.hist(np.log10(inj.period),histtype='step',bins=b,label='injTCE',color='black',lw=2.4)
plt.xlabel('log$_{10}$(Period (days))')
plt.legend(loc='upper left',fontsize=12)
plt.xlim([-.4,2.9])


plt.savefig(outname)
#hfas,bedge=np.histogram(np.log10(fas.period),bins=b)
#hinv,bedge=np.histogram(np.log10(inv.period),bins=b)
#hscr,bedge=np.histogram(np.log10(scr.period),bins=b)
#hops,bedge=np.histogram(np.log10(ops.period),bins=b)
#plt.step(bedge[:-1],hops.astype(float)-hscr,color='red',lw=2.0,label='(invTCE+scrTCE)/opsTCE')
#plt.step(bedge[:-1],hops.astype(float)-hinv,color='green',lw=2.0,label='(invTCE+scrTCE)/opsTCE')
#plt.legend(loc='upper left', fontsize=12)

#%%
#Summary of the Catalog Plot
#Plot of the Rp vs Period with the marginalized distributions on the side. Score also shown.
plt.figure(figsize=(8.5,7))
#plt.clf()
want=~np.isnan(koi.Rp)
y=np.log10(koi[want].Rp)
x=np.log10(koi[want].period)
pcs=(koi[want].disp=='CANDIDATE')
scorevalues=koi[want].koi_score
xlim=(-0.4,2.9)
ylim=(-0.5,1.66)
nxbins=115
nybins=130

axHisty,axHistx,axT=ppp.plot2dHistPaper(x,y,pcs,xlim,ylim,scorevalues,\
            scorelimit=0.7,nxbins=nxbins,nybins=nybins,xlabel="Period (days)",\
            ylabel="Rp (R$_{\oplus}$)",\
            makelog=False,midtype='scatter',colormap='plasma',msize=14)

plt.savefig('%s/fig-radiusPeriodScore-hist.png' % outdir)
#%%
#Period Radius with a star type cut.

plt.figure(figsize=(8.5,7))
#plt.clf()
want1=~np.isnan(koi.Rp)
want2=ppp.inBox(koi,'koi_steff',[5300,6000],'koi_slogg',[4.0,6.0])
want=want1&want2
y=np.log10(koi[want].Rp)
x=np.log10(koi[want].period)
pcs=(koi[want].disp=='CANDIDATE')
scorevalues=koi[want].koi_score
xlim=(-0.4,2.9)
ylim=(-0.5,1.66)
nxbins=115
nybins=130

axHisty,axHistx,axT=ppp.plot2dHistPaper(x,y,pcs,xlim,ylim,scorevalues,\
            scorelimit=0.7,nxbins=nxbins,nybins=nybins,xlabel="Period (days)",\
            ylabel="Rp (R$_{\oplus}$)",\
            makelog=False,midtype='scatter',colormap='plasma',msize=14)
plt.title('G Dwarf Stars')
plt.savefig('%s/fig-radiusPeriodScore-hist-Gdwarf.png' % outdir)
#%%
#Summary of the Catalog Plot
#Plot of the Rp vs Period with the marginalized distributions on the side. Score also shown.
#Create the plot for the false alarms. A bit of a change from above
plt.figure(figsize=(8.5,7))

want=~np.isnan(fas.Rp)
y=fas[want].Rp
x=fas[want].period
pcs=(fas[want].disp=='PC') 
scorevalues=fas[want].score
#xlim=(1.1,2.9)
xlim=(70,650)
#ylim=(-0.5,1.66)
ylim=(.5,8)
nxbins=30
nybins=30

axHisty,axHistx,axT=ppp.plot2dHistPaper(x,y,pcs,xlim,ylim,scorevalues,\
            scorelimit=0.5,nxbins=nxbins,nybins=nybins,xlabel="Period (days)",\
            ylabel="Rp (R$_{\oplus}$)",\
            makelog=False,midtype='scatter',colormap='plasma',msize=14,logticks=False)
#
plt.title('False Alarms')
plt.savefig('%s/fig-radiusPeriodScore-hist-invscr.png' % outdir)


#%%
#Linear long period of the ops dataset
#Summary of the Catalog Plot
#Plot of the Rp vs Period with the marginalized distributions on the side. Score also shown.
plt.figure(figsize=(8.5,7))
#plt.clf()
want=~np.isnan(koi.Rp)
y=koi[want].Rp
x=koi[want].period
pcs=(koi[want].disp=='CANDIDATE')
scorevalues=koi[want].koi_score
#xlim=(1.1,2.9)
xlim=(70,650)
#ylim=(-0.5,1.66)
ylim=(.5,8)
nxbins=30
nybins=30

axHisty,axHistx,axT=ppp.plot2dHistPaper(x,y,pcs,xlim,ylim,scorevalues,\
            scorelimit=0.5,nxbins=nxbins,nybins=nybins,xlabel="Period (days)",\
            ylabel="Rp (R$_{\oplus}$)",\
            makelog=False,midtype='scatter',colormap='plasma',msize=14,logticks=False)

plt.savefig('%s/fig-radiusPeriodScore-hist-linear.png' % outdir)

#%%
#Plot the Habitable zone for scores > 0.5
cumkoi=dr25.loadCumKOITable()
koi=dr25.loadKOITable()
#Plot the HZ candidates
plt.figure(figsize=(8.5,11))
hzframe=ppp.HZplot2(koi, scut=0.5,pradcut=1.8, mfactor=18, addConf=False, ylimits=(3000,6600), xlimits=(2.5,0.0),mycolormap="GnBu")
plt.savefig('%s/fig-hzTstarInsol.png' % outdir)

#%%
#Plot the number of candidates as a function of score.
metric1='period'
metric2='mes'
range1=[200,500]
range2=[7,10]
scores=np.arange(0.05,0.95,.1)
ppp.plotNpcScore(ops,fas,inj,metric1,metric2,range1,range2,\
               scores=scores,annotate="Period: 200--500\nMES:7--10")
               
plt.savefig('%s/fig-varyScoreNcandidatesBox2.png' % outdir)
#%%
metric1='Rp'
metric2='Sp'
range1=[0.5,2.0]
range2=[0.5,2.0]
scores=np.arange(0.05,0.95,.1)
ppp.plotNpcScore(ops,fas,inj,metric1,metric2,range1,range2,\
               scores=scores,annotate="Rp:0.5--2.0\nSp:0.5--2.0")
               
plt.savefig('%s/fig-varyScoreNcandidatesHZ.png' % outdir)

#%%
metric1='ngoodt'
metric2='MES'
range1=[0,6]
range2=[7,10]
scores=np.arange(0.05,1.0,.1)
ppp.plotNpcScore(ops,fas,inj,metric1,metric2,range1,range2,\
               scores=scores,annotate="Ntransits:0--6\nMES:7--10")
               
plt.savefig('%s/fig-varyScoreNcandidatesNtransits.png' % outdir)
#%%
#Calculate the Reliability and Complenetess for certain boxes.
#This is the calculation for period 200-500 and MES 7.1--10
metric1='period'
metric2='mes'
range1=(200,500)
range2=(7.0,10)
opsbox=ppp.inBox(ops,metric1,range1,metric2,range2)
fasbox=ppp.inBox(fas,metric1,range1,metric2,range2)
R,E=pmp.estimateReliability(fas[fasbox],ops[opsbox],s=[0.5,1.0])
print R,E

#%%
##%%
#Now some finer detail. Show period and mes dependence of completeness, reliability 
#
plt.figure(figsize=(8.5,11))
plt.subplot(421)
#bins=np.linspace(-0.3,2.7,num=23)
bins=np.linspace(0,550,num=12)
opsgroup=[ops.mes<=10,ops.mes>10]
bothgroup=[fas.mes<=10,fas.mes>10]
prv.plot1DReliabilityGroups(ops,fas,'period',bins,opsgroup,bothgroup,s=[0.0,1.0],xlabel='Period (d)',labels=['mes$\leq$10','mes>10'])
plt.ylim(0,1.1)
plt.xlim(0,500)
plt.legend(loc='best',fontsize=12)
#lab=np.array([1,3,10,30,100,300])
#locat=np.log10(lab)
#plt.xticks(locat,lab.astype(str))

ntransname='ntrans'
plt.subplot(423)
bins=np.linspace(7.0,12,num=11)
opsgroup=[ops.period<=100,ops.period>100]
bothgroup=[fas.period<=100,fas.period>100]
prv.plot1DReliabilityGroups(ops,fas,'mes',bins,opsgroup,bothgroup,s=[0.0,1.0],xlabel='MES',labels=['P$\leq$100 d','P>100 d'])
plt.ylim(0,1.1)
plt.legend(loc='best',fontsize=12)

plt.subplot(425)
bins=np.linspace(3,16,num=14)
bins=np.linspace(3,17,num=8)
opsgroup=[ops.mes<=10,ops.mes>10]
bothgroup=[fas.mes<=10,fas.mes>10]
prv.plot1DReliabilityGroups(ops,fas,ntransname,bins,opsgroup,bothgroup,s=[0.0,1.0],xlabel='N Transits',labels=['mes$\leq$10','mes$>$10'])
plt.ylim(0,1.1)
plt.xlim(3,15)
plt.legend(loc='best',fontsize=12)

plt.subplot(427)
bins=np.linspace(0,40,num=17)
opsgroup=[ops.mes<=10,ops.mes>10]
bothgroup=[fas.mes<=10,fas.mes>10]
#opsgroup=[ops.period<=100,ops.period>100]
#bothgroup=[fas.period<=100,fas.period>100]
prv.plot1DReliabilityGroups(ops,fas,'duration',bins,opsgroup,bothgroup,s=[0.0,1.0],xlabel='Duration (hrs)',labels=['mes$\leq$10','mes$>$10'])
plt.ylim(0,1.1)
plt.xlim(0,25)
plt.legend(loc='best',fontsize=12)

plt.subplot(422)
#bins=np.linspace(-0.3,2.7,num=46)
bins=np.linspace(0,500,num=80)
injgroup=[inj.Expected_MES<=10,inj.Expected_MES>10]
prv.plot1DCompletenessGroups(inj,'period',bins,injgroup,s=[0,1.0],xlabel='Period (d)',labels=['Exp_mes$\leq$10','Exp_mes>10'])
plt.ylim(0,1.1)
plt.xlim(0,500)
plt.legend(loc='best',fontsize=12)
#lab=np.array([1,3,10,30,100,300])
#locat=np.log10(lab)
#plt.xticks(locat,lab.astype(str))

plt.subplot(424)
bins=np.linspace(7.0,12,num=11)
injgroup=[inj.period<=100,inj.period>100]
prv.plot1DCompletenessGroups(inj,'Expected_MES',bins,injgroup,s=[0.0,1.0],xlabel='EXP MES',labels=['P$\leq$100 d','P>100 d'])
plt.ylim(0,1.1)
plt.legend(loc='best',fontsize=12)

plt.subplot(426)
bins=np.linspace(3,17,num=15)
bins=np.linspace(3,17,num=8)
injgroup=[inj.Expected_MES<=10,inj.Expected_MES>10]
prv.plot1DCompletenessGroups(inj,ntransname,bins,injgroup,s=[0.0,1.0],xlabel='N Transits',labels=['Exp_MES$\leq$10','Exp_MES$>$10'])
plt.ylim(0,1.1)
plt.xlim(3,15)
plt.legend(loc='best',fontsize=12)

plt.subplot(428)
#bins=np.linspace(-0.3,2.7,num=46)
bins=np.linspace(0,40,num=25)
#injgroup=[inj.Expected_MES<=10,inj.Expected_MES>10]
#injgroup=[inj.period<=100,inj.period>100]
prv.plot1DCompletenessGroups(inj,'duration',bins,injgroup,s=[0,1.0],xlabel='Duration (hrs)',labels=['Exp_mes$\leq$10','Exp_mes>10'])
plt.ylim(0,1.1)
plt.xlim(0,25)
plt.legend(loc='best',fontsize=12)




plt.tight_layout()

plt.savefig(outdir + 'fig-compRel1D-PerMes.png')


#%%
#Completeness and REliability by Duration
plt.figure()
plt.subplot(122)
#bins=np.linspace(-0.3,2.7,num=46)
bins=np.linspace(0,40,num=17)
#injgroup=[inj.Expected_MES<=10,inj.Expected_MES>10]
#injgroup=[inj.period<=100,inj.period>100]
prv.plot1DCompletenessGroups(inj,'duration',bins,injgroup,s=[0,1.0],xlabel='Duration (hrs)',labels=['Exp_mes$\leq$10','Exp_mes>10'])
plt.ylim(0,1.1)
plt.xlim(0,25)
plt.legend(loc='best',fontsize=12)

plt.subplot(121)
bins=np.linspace(0,40,num=17)
opsgroup=[ops.mes<=10,ops.mes>10]
bothgroup=[fas.mes<=10,fas.mes>10]
#opsgroup=[ops.period<=100,ops.period>100]
#bothgroup=[fas.period<=100,fas.period>100]
prv.plot1DReliabilityGroups(ops,fas,'duration',bins,opsgroup,bothgroup,s=[0.0,1.0],xlabel='Duration (hrs)',labels=['mes$\leq$10','mes$>$10'])
plt.ylim(0,1.1)
plt.xlim(0,25)
plt.legend(loc='best',fontsize=12)


#%%
#Radius Insolation Flux Score

pcs=(koi.koi_pdisposition=='CANDIDATE') & (koi.koi_steff>3500) & (koi.koi_steff<6000) & (koi.koi_slogg>4.0)

prad=koi[pcs].koi_prad
praderr1=koi[pcs].koi_prad_err1
praderr2=koi[pcs].koi_prad_err2

period=koi[pcs].koi_period
insol=koi[pcs].koi_insol
insoler1=koi[pcs].koi_insol_err1
teff=koi[pcs].koi_steff
score=koi[pcs].koi_score
mes=koi[pcs].koi_max_mult_ev
ebarwant=(insol<3) & (prad<3)
mederbar1=np.median(praderr1[ebarwant])
mederbar2=np.median(praderr2[ebarwant])
medinsoler1=np.median(insoler1[ebarwant])

#norm=colors.Normalize(vmin=3500,vmax=9500,clip=True)
#mapper = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap('jet'))

plt.figure(figsize=(8.5,6))
#plt.scatter	(np.log10(period),np.log10(prad),s=score*30+8,c=teff,edgecolors='k',linewidth=0.2,cmap="rainbow_r",marker='o')
plt.errorbar(0,0,yerr=[.1487],ecolor='gray',fmt='none')
plt.scatter	(np.log10(insol),np.log10(prad),s=score*30+8,c=teff,edgecolors='k',linewidth=0.2,cmap="rainbow_r",marker='o')

#plt.plot([1],[1],'y*',ms=14)
plt.clim([4000,7500])
cbar=plt.colorbar()
cbar.set_label('Stellar Temperature (K)',fontsize=13)
plt.ylim([-.46,1.1])
#plt.xlim([-.5,3.0])
plt.xlim([-1.2,1.8]) #For insolation flux
lab=np.array([0.4,1,2,3,10,20,100])
locat=np.log10(lab)
plt.xticks(locat,lab.astype(str))
lab=np.array([0.4,1.0,2.0,4.0,10.0])
locat=np.log10(lab)
plt.yticks(locat,lab.astype(str))
#plt.xlabel('Period (days)',fontsize=14)
plt.xlabel('Insolation Flux ($\oplus$)',fontsize=14)
plt.ylabel('Planet Radius (R$_{\oplus}$)',fontsize=14)

#Need a legend for the score.
import matplotlib.lines as mlines
bigdot=mlines.Line2D([],[],color='red',marker='.',markersize=12,label='Score=1.0',lw=0)
smalldot=mlines.Line2D([],[],color='red',marker='.',markersize=5,label='Score=0.1',lw=0)
plt.legend(handles=[bigdot,smalldot],numpoints=1,framealpha=0.8,loc='upper left')

plt.savefig('%s/fig-CatalogRadiusInsolScore.png' % outdir)
#%%
#Plot the linear version of insolation flux and radius with error bars
plt.figure(figsize=(8.5,6))

plt.errorbar(1,1,yerr=[mederbar1],xerr=[medinsoler1],fmt='none',ecolor='grey')
plt.scatter	(insol,prad,s=score*30+8,c=teff,edgecolors='k',linewidth=0.2,cmap="rainbow_r",marker='o')
plt.clim([4000,7500])
cbar=plt.colorbar()
cbar.set_label('Stellar Temperature (K)',fontsize=13)
plt.ylim([0,3.5])
plt.xlim([3.0,0])
plt.xlabel('Insolation Flux ($\oplus$)',fontsize=14)
plt.ylabel('Planet Radius (R$_{\oplus}$)',fontsize=14)
bigdot=mlines.Line2D([],[],color='red',marker='.',markersize=12,label='Score=1.0',lw=0)
smalldot=mlines.Line2D([],[],color='red',marker='.',markersize=5,label='Score=0.1',lw=0)
plt.legend(handles=[bigdot,smalldot],numpoints=1,framealpha=0.8,loc='lower right')
plt.savefig('%s/fig-CatalogRadiusInsolScoreLinear.png' % outdir)
#%%
#Plot the period vs score for low MES objects in the catalog 
#For both the PCs and the FAs
#
opsw=(ops.MES<=10) & (ops.disp=='PC')
fasw=(fas.MES<=10) & (fas.disp=='PC')

plt.figure(figsize=(8.5,8.5))
plt.plot(np.log10(ops[opsw].period),ops[opsw].score,'.',label='obsPC',ms=13)
plt.plot(np.log10(fas[fasw].period),fas[fasw].score,'r*',label='falPC')
#plt.plot(ops[opsw].NTran,ops[opsw].score,'.',label='obsPC',ms=13)
#plt.plot(fas[fasw].NTran,fas[fasw].score,'r*',label='falPC')

#%%
#Period Radius Reliabiltiy Plots.
#For FGK Dwarf Stars.
opsFgkWant=ppp.fgkwant(ops)
injFgkWant=ppp.fgkwant(inj)
fasFgkWant=ppp.fgkwant(fas)
climrange=(0,100)
outname=outdir + 'fig-FgkReliabilityPR.png'

xbins=np.arange(0,600,100)
xbins=[0,5,100,200,300,400,500]
ybins=np.arange(0,13,2)
plt.figure(1,figsize=(5,5))
plt.clf()
Rfgk=ppp.plotReliability(fas[fasFgkWant],ops[opsFgkWant],xmetric='period',ymetric='Rp',xBins=xbins,yBins=ybins,drange=(1,99),\
                        atitle="FKG Dwarf Reliability", s=[0,1.0],xlabel='Period (d)',ylabel='R$_P$(R$_{\oplus}$)')
plt.savefig(outname)

plt.figure(2,figsize=(5,5))
plt.clf()
Rall=ppp.plotReliability(fas,ops,xmetric='period',ymetric='Rp',xBins=xbins,yBins=ybins,drange=(1,99),\
                        atitle="All Star Reliability", s=[0,1.0],xlabel='Period (d)',ylabel='R$_P$(R$_{\oplus}$)')
plt.savefig(outdir + 'fig-AllReliabilityPR.png')
#pmp.completenessAndReliability(ops,fas,inj,opsFgkWant,bothFgkWant,injFgkWant,climrange,outname,figTitle,s=[0,1.0])

#%
#Planet Radius Completeness Plots
#For FGK Dwarf Stars
injFgkWant=ppp.fgkwant(inj)
#xbins=np.arange(0,600,100)
#xbins=[0,10,50,100,300,400,500]
#ybins=np.arange(0,13,2)
outname=outdir + 'fig-FgkCompletePR.png'
plt.figure(3,figsize=(5,5))
plt.clf()
Cfgk=ppp.plotCompleteness(inj[injFgkWant],xmetric='period',ymetric='Rp',xBins=xbins,yBins=ybins,drange=(1,99),\
                        atitle="FGK Dwarf Completeness",s=[0.0,1.0],xlabel='Period (d)',ylabel='R$_p$(R$_{\oplus}$)')
plt.savefig(outname)
outname=outdir + 'fig-AllCompletePR.png'
plt.figure(4,figsize=(5,5))
plt.clf()
Call=ppp.plotCompleteness(inj,xmetric='period',ymetric='Rp',xBins=xbins,yBins=ybins,drange=(1,99),\
                        atitle="All Stars Completeness",s=[0.0,1.0],xlabel='Period (d)',ylabel='R$_p$(R$_{\oplus}$)')

plt.savefig(outname)

#%%
#%Plot Effectiveness.
fasFgkWant=ppp.fgkwant(fas)
#xbins=np.arange(0,600,100)
#xbins=[0,10,50,100,300,400,500]
#ybins=np.arange(0,13,2)
outname=outdir + 'fig-FgkEffect-PR.png'
plt.figure(3,figsize=(5,5))
plt.clf()
Efgk=ppp.plotCompleteness(fas[fasFgkWant],xmetric='period',ymetric='Rp',xBins=xbins,yBins=ybins,drange=(1,99),\
                        atitle="FGK Dwarf Effectiveness",s=[0.0,1.0],xlabel='Period (d)',ylabel='R$_p$(R$_{\oplus}$)')
plt.savefig(outname)
outname=outdir + 'fig-AllEffect-PR.png'
plt.figure(4,figsize=(5,5))
plt.clf()
Eall=ppp.plotCompleteness(fas,xmetric='period',ymetric='Rp',xBins=xbins,yBins=ybins,drange=(1,99),\
                        atitle="All Stars Effectiveness",s=[0.0,1.0],xlabel='Period (d)',ylabel='R$_p$(R$_{\oplus}$)')

plt.savefig(outname)


#%%
#Period MES Reliabiltiy and Completeness Plots.
#For FGK Dwarf Stars.
opsFgkWant=ppp.fgkwant(ops)
injFgkWant=ppp.fgkwant(inj)
fasFgkWant=ppp.fgkwant(fas)
climrange=(0,100)
outname=outdir + 'fig-FgkReliabilityPmes.png'

xbins=[0,10,200,500]
ybins=[0,10,20,2000]

plt.figure(1,figsize=(6.4,5))
plt.clf()
Rfgk=ppp.plotReliability(fas[fasFgkWant],ops[opsFgkWant],xmetric='period',ymetric='mes',xBins=xbins,yBins=ybins,drange=(0,98),\
                        atitle="FKG Dwarf Reliability", s=[0,1.0],xlabel='Period (d)',ylabel='MES')
plt.savefig(outname)

plt.figure(2,figsize=(6.4,5))
plt.clf()
Rall=ppp.plotReliability(fas,ops,xmetric='period',ymetric='mes',xBins=xbins,yBins=ybins,drange=(0,98),\
                        atitle="All Star Reliability", s=[0,1.0],xlabel='Period (d)',ylabel='MES')
plt.savefig(outdir + 'fig-AllReliabilityPmes.png')
#pmp.completenessAndReliability(ops,fas,inj,opsFgkWant,bothFgkWant,injFgkWant,climrange,outname,figTitle,s=[0,1.0])

#%
#Planet Radius Completeness Plots
#For FGK Dwarf Stars
injFgkWant=ppp.fgkwant(inj)
#xbins=np.arange(0,600,100)
#xbins=[0,10,50,100,300,400,500]
#ybins=np.arange(0,13,2)
outname=outdir + 'fig-FgkCompletePmes.png'
plt.figure(3,figsize=(6.4,5))
plt.clf()
Cfgk=ppp.plotCompleteness(inj[injFgkWant],xmetric='period',ymetric='mes',xBins=xbins,yBins=ybins,drange=(0,98),\
                        atitle="FGK Dwarf Completeness",s=[0.0,1.0],xlabel='Period (d)',ylabel='MES')
plt.savefig(outname)
outname=outdir + 'fig-AllCompletePmes.png'
plt.figure(4,figsize=(6.4,5))
plt.clf()
Call=ppp.plotCompleteness(inj,xmetric='period',ymetric='mes',xBins=xbins,yBins=ybins,drange=(0,98),\
                        atitle="All Stars Completeness",s=[0.0,1.0],xlabel='Period (d)',ylabel='MES')

plt.savefig(outname)

#%
#%Effectiveness
fasFgkWant=ppp.fgkwant(fas)
#xbins=np.arange(0,600,100)
#xbins=[0,10,50,100,300,400,500]
#ybins=np.arange(0,13,2)
outname=outdir + 'fig-FgkEffect-Pmes.png'
plt.figure(3,figsize=(6.4,5))
plt.clf()
Cfgk=ppp.plotCompleteness(fas[fasFgkWant],xmetric='period',ymetric='mes',xBins=xbins,yBins=ybins,drange=(0,98),\
                        atitle="FGK Dwarf Ineffectiveness",s=[0.0,1.0],xlabel='Period (d)',ylabel='MES')
plt.savefig(outname)
outname=outdir + 'fig-AllEffect-Pmes.png'
plt.figure(4,figsize=(6.2,5))
plt.clf()
Call=ppp.plotCompleteness(fas,xmetric='period',ymetric='mes',xBins=xbins,yBins=ybins,drange=(0,98),\
                        atitle="All Stars Ineffectiveness",s=[0.0,1.0],xlabel='Period (d)',ylabel='MES')

plt.savefig(outname)


#%%
#Duration MES Completeness Plots
#For FGK Dwarf Stars
injFgkWant=ppp.fgkwant(inj)
xbins=np.arange(0,40,5)
#xbins=[0,10,50,100,300,400,500]
ybins=np.arange(7,20,2)
#outname=outdir + 'fig-FgkCompletePR.png'
plt.figure()
Cfgk=ppp.plotCompleteness(inj[injFgkWant],xmetric='duration',ymetric='mes',xBins=xbins,yBins=ybins,drange=(1,99),\
                        atitle="FGK Dwarf Completeness",s=[0.0,1.0],xlabel='Duration (hrs)',ylabel='MES')
#plt.savefig(outname)
#outname=outdir + 'fig-AllCompletePR.png'
plt.figure()
Call=ppp.plotCompleteness(inj,xmetric='duration',ymetric='mes',xBins=xbins,yBins=ybins,drange=(1,99),\
                        atitle="All Stars Completeness",s=[0.0,1.0],xlabel='Duration (hrs)',ylabel='MES')

#plt.savefig(outname)

#%%
#Create the Banned List Latex Table.
banned=dr25.loadBannedTCEs()
bpandas=p.DataFrame(banned,columns=['TCE_ID'])

outfile=outdir + 'banned.tex'
bpandas.to_latex(outfile,columns=['TCE_ID'],index=False,header=False)

#%%
#Create the list of INV TCE-Ids that we used to do analysis.
#
outtab=outdir + 'tab-inv.tex'
stub=inv[['period','mes','disp']][:]
stub.to_latex(outtab,header=True,index=True)

#%%
#Create the list of SCR TCE-Ids that we used to do analysis.
#
outtab=outdir + 'tab-scr.tex'
stub=scr[['period','mes','disp']][:]
stub.to_latex(outtab,header=True,index=True)
