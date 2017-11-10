# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 13:50:23 2017

Paper Plots

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
type='DATA-SUP'
ops,inj,inv=io.getAllData(id,tcedata=type,rtype='INV')
ops,inj,scr=io.getAllData(id,tcedata=type,rtype='SS1')
ops,inj,both=io.getAllData(id,tcedata=type)
#Remove those that are banned. -- these were not turned into KOIs and should not be considered for Reliability measures.
#w=(~ops.isdr25koi) & (ops.N==0)
#ops.drop(ops[w].index)
ops['logperiod']=np.log10(ops.period)
both['logperiod']=np.log10(both.period)
inj['logperiod']=np.log10(inj.period)
outdir="/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/Paper/"
plt.switch_backend('Qt4Agg')
#%%
#Plot histogram of the populations
#This needs DATA for the type above
outname=outdir + "fig-tcePeriods.png"
b=np.linspace(-0.301,2.83,250); #Bins

dr24ops=io.loadDR24Catalog()

plt.figure(figsize=(8.5,5.4))

plt.subplot(311)
print (len(ops))
plt.hist(np.log10(ops.period),histtype='step',bins=b,label='DR25 OBS-TCE',color='black',lw=2)
plt.hist(np.log10(dr24ops.period),histtype='step',bins=b,label='DR24 OBS-TCE',color='green',lw=1.3)
plt.xlim([-.4,2.9])
plt.legend(loc='upper left',fontsize=11)

plt.subplot(312)
plt.hist(np.log10(inv.period),histtype='step',bins=b,label='INV-TCE',color='blue',lw=2)
plt.hist(np.log10(scr.period),histtype='step',bins=b,label='SCR-TCE',color='red',lw=2)
plt.legend(loc='upper left',fontsize=11)
plt.xlim([-.4,2.9])
plt.ylabel('Counts',fontsize=14)
plt.subplot(313)
plt.hist(np.log10(inj.period),histtype='step',bins=b,label='INJ-TCE',color='black',lw=2)
plt.xlabel('log$_{10}$(Period (days))')
plt.legend(loc='upper left',fontsize=11)
plt.xlim([-.4,2.9])
plt.savefig(outname)

#%%
#I want to plot the fraction of certain types of fails for large Period bins
#Just the fraction of false positives for OPS.
#For this I need the binary file that says why things fail.
ss1file='/soc/nfs/so-nfs/DR25/minorFlags/IndividualFlagsSS1-r62353.csv'
opsfile='/soc/nfs/so-nfs/DR25/minorFlags/IndividualFlagsOPS-r62353.csv'
invfile='/soc/nfs/so-nfs/DR25/minorFlags/IndividualFlagsINV-r62353.csv'

ss1rf=p.read_csv(ss1file,comment='#',index_col='tce')
invrf=p.read_csv(invfile,comment='#',index_col='tce')
opsrf=p.read_csv(opsfile,comment='#',index_col='tce')

#Probably a good idea to group up these metrics into one figure.
#Then it is a nice concise plot for the paper.
#metric='g_LPP_HIGH'
#metric='g_MARSHALL'
#metric='g_RUBBLE'
metric='g_ITRANS_ALL'
#metric='g_SKYE'
#metric='g_SIG_PRI_OVER_FRED'

metrics=['g_LPP_HIGH','g_ITRANS_ALL','g_SKYE','g_SIG_PRI_OVER_FRED']
names=['LPP','INDIVID.\nTRANS','SKYE','MS$_1$']

plt.figure(figsize=[8.5,5.5])

for i,m in enumerate(metrics):
    metric=metrics[i]
    name=names[i]    
    
    fps=ops.N==1
    opslist=ops[fps].index
    opsfps=ops[fps]
    metricfail=opsrf.loc[opslist][metric]==1
    
    opsfpsN,opsfpsB=np.histogram(np.log10(opsfps.period),bins=40)
    opsmetN,opsfpsB=np.histogram(np.log10(opsfps[metricfail].period),bins=opsfpsB)
    
    ss1list=scr.index
    metricfail=ss1rf.loc[ss1list][metric]==1
    ss1fpsN,ss1fpsB=np.histogram(np.log10(scr.period),bins=opsfpsB)
    ss1metN,ss1fpsB=np.histogram(np.log10(scr[metricfail].period),bins=ss1fpsB)
    
    invlist=inv.index
    metricfail=invrf.loc[invlist][metric]==1
    invfpsN,invfpsB=np.histogram(np.log10(inv.period),bins=opsfpsB)
    invmetN,invfpsB=np.histogram(np.log10(inv[metricfail].period),bins=invfpsB)
    
    plt.subplot(2,2,i+1)
    plt.plot(opsfpsB[:-1],opsmetN.astype(float)/opsfpsN,'-k',label='OPS NTL',lw=2)
    plt.plot(ss1fpsB[:-1],ss1metN.astype(float)/ss1fpsN,'--b',label='SCR',lw=1.8)
    plt.plot(ss1fpsB[:-1],invmetN.astype(float)/invfpsN,':r',label='INV',lw=1.8)
    plt.legend(loc='best',fontsize=10)
    plt.xlabel('log$_{10}$(Period (d))',fontsize=12)
    plt.ylabel('Fraction of FPs',fontsize=12)
    plt.ylim(0,1)
    plt.annotate(name,(1.25,0.4),xycoords="data",fontsize=11)
    plt.tight_layout()

plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/Paper/fig-fractionFailsByMetric.png')
#%%
#Same as above but for MES

plt.figure(figsize=[8.5,5.5])

for i,m in enumerate(metrics):
    metric=metrics[i]
    name=names[i]    
    
    fps=ops.N==1
    opslist=ops[fps].index
    opsfps=ops[fps]
    metricfail=opsrf.loc[opslist][metric]==1
    
    mesbins=np.linspace(7,40,num=40)
    opsfpsN,opsfpsB=np.histogram(opsfps.mes,bins=mesbins)
    opsmetN,opsfpsB=np.histogram(opsfps[metricfail].mes,bins=opsfpsB)
    
    ss1list=scr.index
    metricfail=ss1rf.loc[ss1list][metric]==1
    ss1fpsN,ss1fpsB=np.histogram(scr.mes,bins=opsfpsB)
    ss1metN,ss1fpsB=np.histogram(scr[metricfail].mes,bins=ss1fpsB)
    
    invlist=inv.index
    metricfail=invrf.loc[invlist][metric]==1
    invfpsN,invfpsB=np.histogram(inv.mes,bins=opsfpsB)
    invmetN,invfpsB=np.histogram(inv[metricfail].mes,bins=invfpsB)
    
    plt.subplot(2,2,i+1)
    plt.plot(opsfpsB[:-1],opsmetN.astype(float)/opsfpsN,'-k',label='OPS NTL',lw=2)
    plt.plot(ss1fpsB[:-1],ss1metN.astype(float)/ss1fpsN,'--b',label='SCR',lw=1.5)
    plt.plot(ss1fpsB[:-1],invmetN.astype(float)/invfpsN,':r',label='INV',lw=1.5)
    plt.legend(loc='best',fontsize=10)
    plt.xlabel('MES',fontsize=12)
    plt.ylabel('Fraction of FPs',fontsize=12)
    plt.ylim(0,1)
    plt.xlim(7,30)
    plt.annotate(name,(20,.01),xycoords="data",fontsize=11)
    plt.tight_layout()

plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/Paper/fig-fractionFailsByMetricMes.png')
#%%
#Create the grids for Completeness and Reliability.
def fgkwant(data):
    want=(data['logg']>=4.0) & (data['tstar']>=4000.0) & (data['tstar']<7000.0);
    return want

def getHZpercentS(data,Srange=(0.75,1.25),pRadius=(0.75,1.25),s=(0.0,1.0)):
    """
    Return the percentage of PCs in the HZ with small radius.
    Use the input ranges.
    """
    
    hzwant=(data['srad'] >=Srange[0]) & (data['srad']<=Srange[1]) & \
           (data['rplanet'] >=pRadius[0]) & (data['rplanet']<=pRadius[1])

    hzdata=data.loc[hzwant]
    pcwant=(hzdata['disp'] == 'PC') & (hzdata['score']>s[0])
    fpwant=(hzdata['disp']== 'FP') & (hzdata['score']>s[1])
    passed=pcwant | fpwant    
    
    total=float(len(hzdata))
    pcs=float(len(hzdata[passed]))    
    percent=100 * (pcs/total) 
    
    return (percent,total,pcs,hzwant)
#%%
#Continue working
opsFgkWant=fgkwant(ops)
injFgkWant=fgkwant(inj)
bothFgkWant=fgkwant(both)
climrange=(0,100)
outname='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/FinalAnalysis/Paper/fig-completeReliabilityCard.png'
figTitle=""
pmp.completenessAndReliability(ops,both,inj,opsFgkWant,bothFgkWant,injFgkWant,climrange,outname,figTitle,s=[0.0,1.0])
#%%
#Same again, but use a score cut of 0.6
opsFgkWant=fgkwant(ops)
injFgkWant=fgkwant(inj)
bothFgkWant=fgkwant(both)
climrange=(0,100)
outname=outdir + 'fig-completeReliabilityCard-score.png'
figTitle=""
pmp.completenessAndReliability(ops,both,inj,opsFgkWant,bothFgkWant,injFgkWant,climrange,outname,figTitle,s=[0.6,0.6])
#%%
#Create the score card for reference
figureTitle=""
outname=outdir + 'fig-scorecard.png'
pRadius=[0.5,1.7]
Srange=[0.5,1.7]
scores=[0,1.0]
hzpercent=np.ones((3,2),dtype=float)
hztotal=np.ones((3,2),dtype=int)
hzpcs=np.ones((3,2),dtype=int)
(hzpercent[0,0],hztotal[0,0],hzpcs[0,0],hzopswantS)=getHZpercentS(ops,\
                                        Srange=Srange,pRadius=pRadius,s=scores)  
(hzpercent[0,1],hztotal[0,1],hzpcs[0,1],hzopsfgkwantS)=getHZpercentS(ops[opsFgkWant],\
                                        Srange=Srange,pRadius=pRadius,s=scores)
(hzpercent[1,0],hztotal[1,0],hzpcs[1,0],hzinvwantS)=getHZpercentS(both,\
                                        Srange=Srange,pRadius=pRadius,s=scores)
(hzpercent[1,1],hztotal[1,1],hzpcs[1,1],hzinvfgkwantS)=getHZpercentS(both[bothFgkWant],\
                                        Srange=Srange,pRadius=pRadius,s=scores)
(hzpercent[2,0],hztotal[2,0],hzpcs[2,0],hzinjwantS)=getHZpercentS(inj,\
                                        Srange=Srange,pRadius=pRadius, s=scores)
(hzpercent[2,1],hztotal[2,1],hzpcs[2,1],hzinjwantS)=getHZpercentS(inj[injFgkWant],\
                                        Srange=Srange,pRadius=pRadius,s=scores)  

pmp.createScoreCard(ops,ops[opsFgkWant],both,both[bothFgkWant],inj,inj[injFgkWant],Srange,pRadius,\
                        hzpercent,hztotal,hzpcs,climrange,figureTitle,outname)
#%%
#Now some finer detail. Show period and mes dependence of completeness, reliability 
#
plt.figure(figsize=(8.5,8))
plt.subplot(321)
#bins=np.linspace(-0.3,2.7,num=23)
bins=np.linspace(0,550,num=12)
opsgroup=[ops.mes<=10,ops.mes>10]
bothgroup=[both.mes<=10,both.mes>10]
prv.plot1DReliabilityGroups(ops,both,'period',bins,opsgroup,bothgroup,s=[0.0,1.0],xlabel='Period (d)',labels=['mes$\leq$10','mes>10'])
plt.ylim(0,1.1)
plt.xlim(0,500)
plt.legend(loc='best',fontsize=12)
#lab=np.array([1,3,10,30,100,300])
#locat=np.log10(lab)
#plt.xticks(locat,lab.astype(str))

plt.subplot(323)
bins=np.linspace(7.0,12,num=11)
opsgroup=[ops.period<=100,ops.period>100]
bothgroup=[both.period<=100,both.period>100]
prv.plot1DReliabilityGroups(ops,both,'mes',bins,opsgroup,bothgroup,s=[0.0,1.0],xlabel='MES',labels=['P$\leq$100 d','P>100 d'])
plt.ylim(0,1.1)
plt.legend(loc='best',fontsize=12)

plt.subplot(325)
bins=np.linspace(3,11,num=9)
opsgroup=[ops.mes<=10,ops.mes>10]
bothgroup=[both.mes<=10,both.mes>10]
prv.plot1DReliabilityGroups(ops,both,'ngood',bins,opsgroup,bothgroup,s=[0.0,1.0],xlabel='N Transits',labels=['mes$\leq$10','mes$>$10'])
plt.ylim(0,1.1)
plt.xlim(3,10)
plt.legend(loc='best',fontsize=12)

plt.subplot(322)
#bins=np.linspace(-0.3,2.7,num=46)
bins=np.linspace(0,500,num=80)
injgroup=[inj.mes<=10,inj.mes>10]
prv.plot1DCompletenessGroups(inj,'period',bins,injgroup,s=[0,1.0],xlabel='Period (d)',labels=['mes$\leq$10','mes>10'])
plt.ylim(0,1.1)
plt.xlim(0,500)
plt.legend(loc='best',fontsize=12)
#lab=np.array([1,3,10,30,100,300])
#locat=np.log10(lab)
#plt.xticks(locat,lab.astype(str))

plt.subplot(324)
bins=np.linspace(7.0,12,num=11)
injgroup=[inj.period<=100,inj.period>100]
prv.plot1DCompletenessGroups(inj,'mes',bins,injgroup,s=[0.0,1.0],xlabel='MES',labels=['P$\leq$100 d','P>100 d'])
plt.ylim(0,1.1)
plt.legend(loc='best',fontsize=12)

plt.subplot(326)
bins=np.linspace(3,11,num=9)
injgroup=[inj.mes<=10,inj.mes>10]
prv.plot1DCompletenessGroups(inj,'ngood',bins,injgroup,s=[0.0,1.0],xlabel='N Transits',labels=['mes$\leq$10','mes$>$10'])
plt.ylim(0,1.1)
plt.xlim(3,10)
plt.legend(loc='best',fontsize=12)


plt.tight_layout()

plt.savefig(outdir + 'fig-compRel1D-PerMes.png')
#plt.subplot(323)
#bins=np.linspace(2.5,10.5,num=9)
#opsgroup=[ops.mes<=100]
#bothgroup=[both.mes<=100]
#prv.plot1DReliabilityGroups(ops,both,'ngood',bins,opsgroup,bothgroup,s=[0.0,1.0],xlabel=' Number of Good Transits',labels=['mes$\leq$100'])
#plt.ylim(0,1.1)
#Next Add in the 1d Completeness plots.  These can be finer grids.

#%%
#Adjust score and show how completeness and reliability change.
metric1='period'
metric2='mes'
range1=(200,500)
range2=(7,10)
prv.plotRelCompScore(ops,both,inj,metric1,metric2,range1,range2,scores=np.arange(0.2,1,.1),\
                    Rlim=(.2,1.04),Clim=(0.8,.2))
#plt.annotate("Period:100-500 days\nMES:7--10",(0.4,0.988),xycoords="data",fontsize=13)
#plt.title("Completeness and Reliablity for different Score Thresholds")

plt.savefig('%s/fig-CRadjustScore-DR25.png' % outdir)

#%%
#Plot the Pradius vs Period, colored by stellar temp and sized by score (1.0 is big, 0.1 is small)
#mcmcFile='/home/smullall/Kepler/RoboVetter/DR25/mcmc/bftable_022117.dat'
mcmcFile='/home/smullall/Kepler/RoboVetter/DR25/mcmc/bftable_040717.dat'

mcmcdata=p.read_csv(mcmcFile,delimiter=',')

koiNames =map( lambda x: "K%08.2f" % (x), mcmcdata.KOI)
mcmcdata["koiName"]=koiNames
mcmcdata.set_index("koiName",inplace=True)
pcs=ops.disp=='PC'
mcmcops=p.merge(ops,mcmcdata,how='left', left_on='dr25koiName',right_index=True,suffixes=('','_mcmc'))
#%%
#Continue this plot
import matplotlib.colors as colors
import matplotlib.cm as cm

pcs=ops.disp=='PC'
prad=mcmcops[pcs].Rp_mcmc
period=mcmcops[pcs].Period_mcmc
insol=mcmcops[pcs].Srad
teff=ops[pcs].tstar
score=ops[pcs].score
mes=ops[pcs].mes

#norm=colors.Normalize(vmin=3500,vmax=9500,clip=True)
#mapper = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap('jet'))

plt.figure(figsize=(8.5,6))
#plt.scatter	(np.log10(period),np.log10(prad),s=score*30+8,c=teff,edgecolors='k',linewidth=0.2,cmap="rainbow_r",marker='o')
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
#Histograms of planet radii.
#Still uses the MCMC fits (mcmcops)
plimit=100
plimit2=300
radbins=np.logspace(np.log10(.5),np.log10(12),num=20,base=10.0)
want=mcmcops.disp=='PC'
pwant=(mcmcops.Period_mcmc>plimit) & (mcmcops.Period_mcmc<plimit2)
injwant=inj.disp=='PC'
injpwant=inj.period>plimit & (inj.period<plimit2)
bothwant=both.disp=='PC'
bothpwant=both.period>plimit & (both.period<plimit2)

npcinj,bb=np.histogram(inj[injwant & injpwant].rplanet,bins=radbins)
ninj,bb=np.histogram(inj[injpwant].rplanet,bins=radbins)
C=npcinj.astype(float)/ninj.astype(float)
R,E=prv.get1DReliability(ops[pwant],both[bothpwant],'rplanet',bins=radbins,s=[0,1])

plt.figure(8)
plt.clf()
numops,bb,patch=plt.hist(mcmcops.Rp[want & pwant],bins=radbins,histtype='step',lw=3,label='Candidate Counts')
plt.annotate('Period>%4.1f <300 d' % plimit,[.6,.7],xycoords='figure fraction',fontsize=17)
adjNum=numops*R/C

plt.step(bb[:-1],adjNum,'r-',where='pre',lw=3,label='Adjusted for C&R')
plt.xlim([0,12])
plt.xlabel('Planet Radius ($\oplus$)',fontsize=13)
plt.ylabel('Counts',fontsize=13)
plt.legend(loc='best')
plt.savefig('%s/pradius-histogram2.png' % outdir)

#%%
#Plot the supplemental DV (D.2 product) vs MCMC (KOI Table) planet radii 

plt.figure(8)
plt.clf()
pcs=(mcmcops.disp=='PC') & (mcmcops.mes>10)
pcs2=(mcmcops.disp=='PC') & (mcmcops.mes<=10)
plt.plot(mcmcops[pcs].Rp_mcmc,mcmcops[pcs].rplanet,'.',ms=5,label='DR25 PCs')
#plt.plot(mcmcops[pcs2].Rp_mcmc,mcmcops[pcs2].rplanet,'.r',ms=5,label='DR25 PCs MES$\leq$10')
plt.xlim((0,8))
plt.ylim((0,8))
plt.plot([0,20],[0,20],'k-',lw=1)
plt.xlabel('MCMC Planet Radius ($\oplus$)',fontsize=15)
plt.ylabel('DV SUPP Planet Radius ($\oplus$)',fontsize=15)
plt.legend(fontsize=14,numpoints=1,loc='best')

plt.savefig('%s/fig-comparePradius-mcmcSup.png' % outdir)

#%%
#Plot histograms for different MES regions.
plt.figure(8)
plt.clf()
n=False
pcs=(mcmcops.disp=='PC')
want1=mcmcops.mes<10 
want2=(mcmcops.mes>=10) & (mcmcops.mes<20)
want3=mcmcops.mes>20
mcmcops['diffrp']=(mcmcops.Rp_mcmc-mcmcops.rplanet)/mcmcops.Rp_mcmc
want=~np.isnan(mcmcops.diffrp)
b=np.linspace(-0.9,0.9,num=70)
plt.hist(mcmcops[pcs & want1 & want].diffrp,bins=b,histtype='step',label='MES<10',normed=n,lw=2)
plt.hist(mcmcops[pcs & want2 & want].diffrp,bins=b,histtype='step',label='10<MES<20',normed=n,lw=2)
plt.hist(mcmcops[pcs & want3 & want].diffrp,bins=b,histtype='step',label='MES>20',normed=n,lw=2)
plt.vlines(0,0,600,linestyles='dotted')
plt.legend(loc='best',fontsize=14)
plt.xlabel('Rp$_{MCMC}$/Rp$_{DV Supp}$ - 1',fontsize=15)
plt.ylabel('Counts',fontsize=15)
plt.ylim(0,420)
plt.xlim([-.7,.7])

plt.savefig('%s/fig-comparePradius-histogram.png' % outdir)

print np.median(mcmcops[pcs & want1 & want].diffrp)
print np.median(mcmcops[pcs & want2 & want].diffrp)
print np.median(mcmcops[pcs & want3 & want].diffrp)
print np.median(mcmcops[pcs & want].diffrp)
print '-----'
print np.std(mcmcops[pcs & want1 & want].diffrp)
print np.std(mcmcops[pcs & want2 & want].diffrp)
print np.std(mcmcops[pcs & want3 & want].diffrp)
print np.std(mcmcops[pcs & want].diffrp)
print '-----'
print np.mean(mcmcops[pcs & want1 & want].diffrp)
print np.mean(mcmcops[pcs & want2 & want].diffrp)
print np.mean(mcmcops[pcs & want3 & want].diffrp)
#%%
# Plot of the Rp vs Period with the marginalized distributions on the side. Score also shown.
plt.figure(figsize=(8.5,5.1))
#plt.clf()
want=~np.isnan(mcmcops.Rp_mcmc)
y=np.log10(mcmcops[want].Rp_mcmc)
x=np.log10(mcmcops[want].Period_mcmc)
pcs=mcmcops[want].disp=='PC'
scorevalues=mcmcops[want].score
xlim=(-0.4,2.9)
ylim=(-0.5,1.66)
nxbins=120
nybins=100

axHisty,axHistx,axT=prv.plot2dHistPaper(x,y,pcs,xlim,ylim,scorevalues,\
            scorelimit=0.8,nxbins=nxbins,nybins=nybins,xlabel="Period (days)",ylabel="Rp (R$_{\oplus}$)",\
            makelog=False,midtype='scatter',colormap='plasma_r',msize=19)

plt.savefig('%s/fig-radiusPeriodScore-hist.png' % outdir)

#%%
#HZ plot
#Remove things with large error bars or low SNR plots
#pcs=(mcmcops.disp == 'PC')
import hzCalc
scut=0.6
pcs=(mcmcops.score > scut)
highsnr=mcmcops['SNR']>0
smallerrors=(-1*mcmcops['Rp-sig']/mcmcops['Rp_mcmc'] < 2.0 ) & (-1*mcmcops['Srad-sig']/mcmcops['Srad']<2.0)
want=mcmcops['Rp_mcmc']+mcmcops['Rp-sig'] < 2.0
want2=np.array(map( lambda x,y:hzCalc.inHZ(x,y), mcmcops['tstar'], mcmcops['Srad']+mcmcops['Srad-sig'] ))
#want2=mcmcops.Srad + mcmcops['Srad-sig'] < 2.2
#want3=mcmcops.Srad + mcmcops['Srad+sig'] > 0.2  #This doesn't actually cut anything.

rearth=1
rsupearth=2
mfactor=15
pars=['KOI','kepler_name','disp','score','period','Rp','Rp_mcmc','Rp-sig','Rp+sig','tstar','SNR','Srad','Srad-sig','Srad+sig','mes','ngood','logg']
hzframe=mcmcops[pcs & want & want2 & want3 & highsnr & smallerrors][pars]

#Add in those three Confirmed ones we are missing.
missing=["009002278-04","008845205-01"]#,"010604335-02","005640085-02"]  #[62f,283c,159c]
for tce in missing:
    mydf=mcmcops[mcmcops.index==tce][pars]
    hzframe=hzframe.append(mydf)
    print len(hzframe)

hzframe['RpBig']=hzframe.Rp_mcmc+hzframe['Rp+sig']
hzframe['RpSmall']=hzframe.Rp_mcmc+hzframe['Rp-sig']

print len(hzframe)

plt.figure(figsize=(8.5,11))

for i,hz in enumerate(hzframe.index):
    Srad=hzframe.loc[hz].Srad
    Smin=Srad+hzframe.loc[hz]['Srad-sig']
    Smax=Srad+hzframe.loc[hz]['Srad+sig']
    rp=hzframe.loc[hz].Rp_mcmc
    rpBig=hzframe.loc[hz].RpBig
    rpSmall=hzframe.loc[hz].RpSmall
    tstar=hzframe.loc[hz].tstar
    sc=hzframe.loc[hz].score
    kepname="%s" % (str(hzframe.loc[hz].KOI))
#plt.scatter(hzframe.Srad,hzframe.tstar,s=(hzframe.Rp*mfactor)**2,marker='o',c=hzframe.score,cmap="plasma",linewidth=0.8)
    plt.hlines(tstar,Smin,Smax,colors='steelblue',label='.',lw=1.2,zorder=1)
    plt.scatter(Srad,tstar,s=(rpBig*mfactor)**2,marker='o',c=sc,cmap="brg_r",linewidth=0.1,vmin=0,vmax=1,zorder=2)
    plt.scatter(Srad,tstar,s=(rpSmall*mfactor)**2,marker='o',c="white",linewidth=0.1,zorder=3)
    #plt.scatter(Srad,tstar,s=(rp*mfactor)**2,marker='o',c=sc,cmap="brg_r",linewidth=.4,vmin=0,vmax=1,zorder=2)
    if hzframe.loc[hz].KOI>7621.01:
        plt.annotate(s=kepname,xy=(Srad-.1,tstar+8),xycoords='data',color="grey",fontsize=10)
    #elif p.notnull(hzframe.loc[hz].kepler_name):
    #    plt.annotate(s=hzframe.loc[hz].kepler_name,xy=(Srad-.1,tstar+8),xycoords='data',color="grey",fontsize=10)
    
#plt.scatter(hzframe.Srad,hzframe.tstar,s=(hzframe.RpBig*mfactor)**2,marker='o',c=hzframe.score,cmap="coolwarm_r",linewidth=0.1,vmin=0,vmax=1)
#plt.scatter(hzframe.Srad,hzframe.tstar,s=(hzframe.RpSmall*mfactor)**2,marker='o',c='white',linewidth=0.1)

small=plt.scatter(1,5800,s=(rearth*mfactor)**2,marker='o',c='grey',linewidth=0.2,label='R$_\oplus$')
big=plt.scatter(1.0,7300,s=(rsupearth*mfactor)**2,marker='o',c='grey',linewidth=0.2,label='2R$_\oplus$')
half=plt.scatter(1.0,6500,s=(1.5*rearth*mfactor)**2,marker='o',c='grey',linewidth=0.2,label='1.5R$_\oplus$')
plt.xlabel('Insolation Flux')
plt.ylabel('Effective Temperature (K)')

plt.legend(handles=[small,half],scatterpoints=1,framealpha=0.8,loc='upper left',fontsize=14)

#Draw a HZ across the plot.
teff=np.linspace(2000,7000,50)
hzBorders=hzCalc.hzBorder(teff)
Shot=hzBorders[:,0]
Scool=hzBorders[:,3]
#hzhot=7000*sx-7500
#hzcool=17500*sx-500
plt.plot(Shot,teff,'--b')
plt.plot(Scool,teff,'--r')
plt.title('Score cut of %f' % scut)

plt.xlim((2.5,0))
plt.ylim((3200,6400))
plt.savefig('%s/fig-hzTeffInsol.png' %outdir)

#plt.xlim((2.5,0))
#plt.ylim((2500,7000))
#plt.savefig('%s/fig-hzTeffInsol2.png' %outdir)
#%%
plt.figure()
hzframe['tstar'].hist(bins=16)
plt.xlabel('Tstar (K)')
plt.savefig('%s/hzTeffGapHist.png' % outdir)
#%%
#Plot radius vs effective temperature of these that are nominally in the habitable zone
#

plt.figure()
#plt.plot(hzframe.Rp_mcmc,hzframe.tstar,'or')
plt.xlabel('Rp_mcmc')
plt.ylabel('Tstar (K)')
plt.scatter(hzframe.Rp_mcmc,hzframe.tstar,s=30,c=hzframe.score,cmap="brg_r",marker='o',linewidth=0.2,vmin=0,vmax=1)
plt.errorbar(hzframe.Rp_mcmc,hzframe.tstar,fmt='none',xerr=hzframe['Rp-sig'],linewidth=0.2,marker='none',linewidth=1)
plt.xlim((0.2,2.4))
plt.ylim((2500,7000))
plt.title('Near or in HZ')
#plt.errorbar(hzframe.Rp_mcmc,hzframe.tstar,xerr=hzframe.Rp-sig)
plt.savefig('%s/fig-hzTeffRadius.png' % outdir)