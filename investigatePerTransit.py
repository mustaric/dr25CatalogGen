# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 13:48:37 2016

@author: smullall

Code to read in individual transit statistics to investigate the SES or MES region that Marshall functions on.

"""
import numpy as np
import pandas as p
import matplotlib.pyplot as plt
import rvIO as io
import pdb



#%%
#This takes a while

def addMarshallValues(injTrOn,marshData):

    #injTrOn['marshall']=np.zeros(len(injTrOn))
    pwant=injTrOn['period']>=150
    injTrOnL=injTrOn[pwant]
    injTrOnL['marshall']=np.zeros(len(injTrOnL))
    injTrOnL['sigbic']=np.zeros(len(injTrOnL))
    injTrOnL['model']=np.zeros(len(injTrOnL))
    
    for idx in injTrOnL.index.values:
        print idx
        tce=injTrOnL.loc[idx,'tce']
        #period=injTrOnL.loc[idx,'period']    
        want=(marshData['tce']==tce)
        #print tce,period
        #pdb.set_trace()
        try:
            minidx=np.argmin(np.abs(marshData[want]['time']-injTrOnL.loc[idx,'time']))
            injTrOnL.loc[idx,'marshall']=np.float(marshData.loc[minidx,'marshall'])  
            injTrOnL.loc[idx,'sigbic']=marshData.loc[minidx,'sigbic'] 
            injTrOnL.loc[idx,'model']=marshData.loc[minidx,'model'] 
        except ValueError:
            injTrOnL.loc[idx,'model']='None Known' 
    
        #injTrOnL.loc[idx,'marshall']        
        

        #for jdx in marshData[want].index.values:
            #if (np.abs(marshData.loc[jdx,'time']-injTrOnL.loc[idx,'time']) < 0.25*period):
                #injTrOnL.loc[idx,'marshall']=np.float(marshData.loc[jdx,'marshall'])
                #break
            
    return injTrOnL
    
#%%
    
injTransitFile='/soc/nfs/so-nfs/DR25/INJ/DATA/TransitTimes.txt'

injTrData=p.read_csv(injTransitFile,sep='\s+',header=None,names=('tce','time','ses','corr','norm'))

injMarshallFile='/home/smullall/Kepler/RoboVetter/DR25/stats/marshallAnalysis/INJ-2016-10-24.txt'

marshData=p.read_csv(injMarshallFile,sep='\s+',header=None,names=('tce','time','marshall','sigbic','model'),comment='#')

injData,a,b,c=io.loadRVInOut(61884,type="INJ-PlanetOn")
tces=injData.index.values

injTrOn=p.merge(injTrData,injData,how='inner',left_on='tce',right_index=True,left_index=False)

injTrOnL=addMarshallValues(injTrOn,marshData)
plt.figure();plt.plot(injTrOnL['ses'],injTrOnL['marshall'],'r.',ms=1)

#%%
            
#now we need to look at the ss data to see where it lies.
invTransitFile = '/soc/nfs/so-nfs/DR25/INV/DATA/TransitTimes.txt'
invMarshalFile= '/soc/nfs/so-nfs/DR25/INV/DATA/Marshall-INV-2016-09-23.txt'
invMarshalFile='/soc/nfs/so-nfs/DR25/INv/DATA/Marshall/inv-2016-09-08-MoreInfo.txt'
invTrData=p.read_csv(invTransitFile,sep='\s+',header=None,names=('tce','time','ses','corr','norm'))
invmarshData=p.read_csv(invMarshallFile,sep='\s+',header=None,names=('tce','time','marshall','sigbic','model'), comment='#')

invtce='/soc/nfs/so-nfs/DR25/INV/DATA/TCEs.txt'
invData=io.readTceInfo(invtce)
tces=invData.index.values
invTrOn=p.merge(invTrData,invData,how='inner',left_on='tce',right_index=True,left_index=False)

invTrOnMar=addMarshallValues(invTrOn,invmarshData)

plt.plot(invTrOnMar['ses'],invTrOnMar['marshall'],'b.',ms=1)


#%
#OPS
opsTransitFile='/soc/nfs/so-nfs/DR25/OPS/DATA/TransitTimes.txt'
#opsMarshalFile='/soc/nfs/so-nfs/DR25/OPS/DATA/Marshall-OPS-2016-09-23.txt'
opsMarshalFile='/soc/nfs/so-nfs/DR25/OPS/Marshall/ops-2016-09-07-moreInfo.txt'
opsTrData=p.read_csv(opsTransitFile,sep='\s+',header=None,names=('tce','time','ses','corr','norm'))
opsmarshData=p.read_csv(opsMarshalFile,sep='\s+',header=None,names=('tce','time','marshall','sigbic','model'),comment='#')
opstce='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
opsData=io.readTceInfo(opstce)
tces=opsData.index.values
opsTrOn=p.merge(opsTrData,opsData,how='inner',left_on='tce',right_index=True,left_index=False)
opsTrOnMar=addMarshallValues(opsTrOn,opsmarshData)

plt.figure()
plt.plot(opsTrOnMar['ses'],opsTrOnMar['marshall'],'g.',ms=1)
#%%


invmar=invTrOnMar['marshall']
invses=invTrOnMar['ses']
injmar=injTrOnL['marshall']
injses=injTrOnL['ses']


steps=np.arange(-2,30,6)
steps=(-2,0,1,2,3,4,6,20)
pl=len(steps)
bins=np.arange(-150,150,3)
plt.figure(figsize=(8.5,11))
plt.subplots_adjust(hspace=0.75, wspace=0.75)
for i,v in enumerate(steps):
    
    plt.subplot(pl,1,i+1)
    if i==pl-1:
        end=50
    else:
        end=steps[i+1]
    start=v
    wantinv = (invses>start) & (invses<end)
    wantinj = (injses>start) & (injses<end)
    plt.hist(injmar[wantinj],histtype='step',color='red',bins=bins)
    plt.hist(invmar[wantinv],histtype='step',color='blue',bins=bins)
    range=plt.gca().get_ylim()
    plt.vlines(10,0,range[1],color='g')
    plt.title('bin=%4.1f -- %4.1f ses' % (start,end))
    if i==0:
        plt.legend(('inj','inv'))


#%%
steps=np.arange(1,40,.25)
threshold=12.0
ses=invses
mar=invmar

mypass=np.zeros(len(steps))
for i,v in enumerate(steps):
    wantN=(ses<v) & (ses!=0.0)& (mar != 2.0)
    wantPass=wantN & (mar<threshold) & (mar != 2.0)
    mypass[i]=np.float(len(ses[wantPass]))/np.float(len(ses[wantN])) * 100.0

invpass=mypass

ses=injses
mar=injmar
mypass=np.zeros(len(steps))
for i,v in enumerate(steps):
    wantN=(ses<v) & (ses!=0.0)& (mar != 2.0)
    wantPass=wantN & (mar<threshold) & (mar != 2.0)
    mypass[i]=np.float(len(ses[wantPass]))/np.float(len(ses[wantN])) * 100.0

injpass=mypass

plt.figure()
plt.plot(steps,injpass,'r')
plt.plot(steps,invpass,'b')
plt.xlabel('SES')
plt.ylabel('Fraction that Pass Threshold')
plt.title('Fraction below given SES that pass Marshall Score Threshold=%3.1f' % (threshold))
plt.legend(('inj','inv'))

#%%

#Look at high signal to noise Ops
#
ebfile='/home/smullall/Kepler/RoboVetter/DR25/cleanINV/ebCatalog.csv'
ebs=p.read_csv(ebfile,header=0,comment='#')

opsEBs=opsTrOnMar.merge(ebs,how='inner',left_on='kic',right_on='KIC',suffixes=('','_eb'))

#%%
#Can we look at the median values or the standard deviation (but it is only 3 measurements)? 
#Seems like if they are all high then it is very likely an EB
#
data=opsEBs
uniqTCEs=unique(opsEBs['tce'])

num=np.zeros(len(uniqTCEs))
numbad=np.zeros(len(uniqTCEs))
fracbad=np.zeros(len(uniqTCEs))
med=np.zeros(len(uniqTCEs))
avg=np.zeros(len(uniqTCEs))
std=np.zeros(len(uniqTCEs))
mes=np.zeros(len(uniqTCEs))
ses=np.zeros(len(uniqTCEs))

for i,tce in enumerate(uniqTCEs):
    want=data['tce']==tce
    goodwant=(data['ses']!=0) & want
    num[i]=len(data[goodwant])
    badwant=(data[goodwant]['marshall']>10) | (data[goodwant]['marshall']<-100)
    numbad[i]=len(data[goodwant][badwant])
    fracbad[i]=np.float(numbad[i])/np.float(num[i])
    med[i]=np.median(data[goodwant]['marshall'])
    avg[i]=np.mean(data[goodwant]['marshall'])
    std[i]=np.std(data[goodwant]['marshall'])
    mes[i]=np.median(data[goodwant]['mes'])
    ses[i]=np.median(data[goodwant]['ses'])
    
ebMarshallComb=p.DataFrame({'tce':uniqTCEs,'mes':mes,'ses':ses,'num':num,'numbad':numbad,'fracbad':fracbad,'med':med,'std':std,'avg':avg})




#%

data=invTrOnMar
uniqTCEs=unique(data['tce'])
num=np.zeros(len(uniqTCEs))
numbad=np.zeros(len(uniqTCEs))
fracbad=np.zeros(len(uniqTCEs))
med=np.zeros(len(uniqTCEs))
avg=np.zeros(len(uniqTCEs))
std=np.zeros(len(uniqTCEs))
mes=np.zeros(len(uniqTCEs))
ses=np.zeros(len(uniqTCEs))

for i,tce in enumerate(uniqTCEs):
    want=data['tce']==tce
    goodwant=(data['ses']!=0) & want
    num[i]=len(data[goodwant])
    badwant=(data[goodwant]['marshall']>10) | (data[goodwant]['marshall']<-100)
    numbad[i]=len(data[goodwant][badwant])
    fracbad[i]=np.float(numbad[i])/np.float(num[i])
    med[i]=np.median(data[goodwant]['marshall'])
    avg[i]=np.mean(data[goodwant]['marshall'])
    std[i]=np.std(data[goodwant]['marshall'])
    mes[i]=np.median(data[goodwant]['mes'])
    ses[i]=np.median(data[goodwant]['ses'])

invMarshallComb=p.DataFrame({'tce':uniqTCEs,'mes':mes,'ses':ses,'num':num,'numbad':numbad,'fracbad':fracbad,'med':med,'std':std,'avg':avg})

#%%
plt.figure()
plt.plot(np.log10(invMarshallComb['mes']),invMarshallComb['med'],'r.',ms=2)
plt.plot(np.log10(ebMarshallComb['mes']),ebMarshallComb['med'],'b.',ms=2)

plt.xlabel('mes')
plt.ylabel('median marshall')

#%%
plt.figure()
plt.plot(np.log10(invMarshallComb['ses']),invMarshallComb['fracbad'],'r.',ms=1)
plt.plot(np.log10(ebMarshallComb['ses']),ebMarshallComb['fracbad'],'b.',ms=2)
plt.xlabel('ses')
plt.ylabel('fracbad')

#%%
plt.figure()
v,b,patch=plt.hist(invMarshallComb['std'],bins=100002,histtype='step',color='red',normed=True)
v2,b2,patch=plt.hist(ebMarshallComb['std'],bins=b,histtype='step',color='blue',normed=True)
#%%
plt.figure()
plt.plot(np.log10(invMarshallComb['std']),invMarshallComb['numbad'],'r.',ms=2)
plt.plot(np.log10(ebMarshallComb['std']),ebMarshallComb['numbad'],'b.',ms=2)

plt.xlabel('std')
plt.ylabel('median marshall')

#%%
plt.figure()
b=np.arange(-40,150,2)
lowbic=opsTrOnMar['sigbic']<700000;
lowses=(opsTrOnMar['ses']<3) & (opsTrOnMar['ses']>0)
want=(opsTrOnMar['model']=='offsetOnly')
plt.hist(opsTrOnMar[want & lowses & lowbic]['marshall'],bins=b,color='red',histtype='step',label='offsetOnly')
plt.xlabel('marshall score')
want=(opsTrOnMar['model']=='spsd')
plt.hist(opsTrOnMar[want & lowses & lowbic]['marshall'],bins=b,color='blue',histtype='step',label='spsd')
want=(opsTrOnMar['model']=='stepUp')
plt.hist(opsTrOnMar[want & lowses & lowbic]['marshall'],bins=b,color='green',histtype='step',label='stepUp')
want=(opsTrOnMar['model']=='stepDown')
plt.hist(opsTrOnMar[want & lowses & lowbic]['marshall'],bins=b,color='magenta',histtype='step',label='stepDown')
want=(opsTrOnMar['model']=='sigmoidBox')
plt.hist(opsTrOnMar[want & lowses & lowbic]['marshall'],bins=b,color='black',histtype='step',label='Transit')
plt.title('OPS prefered model. Marshall Score for SES<3')
plt.legend()
plt.savefig("/home/smullall/Kepler/RoboVetter/DR25/stats/marshallAnalysis/preferredModelOpsLowLowSes.png")

#%%
#As a function of SES number that pass.
plt.figure()
b=np.arange(-4,30,.1)
lowbic=opsTrOnMar['sigbic']<1000000;
lowses=(opsTrOnMar['ses']<3) & (opsTrOnMar['ses']>0)
failmarsh=opsTrOnMar['marshall']>=10;
want=(opsTrOnMar['model']=='offsetOnly')
plt.hist(opsTrOnMar[want & lowbic & failmarsh]['ses'],bins=b,color='red',histtype='step',label='offsetOnly')
plt.xlabel('ses')
want=(opsTrOnMar['model']=='spsd')
plt.hist(opsTrOnMar[want &  lowbic & failmarsh]['ses'],bins=b,color='blue',histtype='step',label='spsd')
want=(opsTrOnMar['model']=='stepUp')
plt.hist(opsTrOnMar[want &  lowbic & failmarsh]['ses'],bins=b,color='green',histtype='step',label='stepUp')
want=(opsTrOnMar['model']=='stepDown')
plt.hist(opsTrOnMar[want &  lowbic & failmarsh]['ses'],bins=b,color='magenta',histtype='step',label='stepDown')
want=(opsTrOnMar['model']=='sigmoidBox')
plt.hist(opsTrOnMar[want &  lowbic & failmarsh]['ses'],bins=b,color='black',histtype='step',label='Transit')
plt.title('OPS prefered model. Number of transits that fail Marshall by method.')
plt.legend()


#%%
#SES distribution
from scipy import stats
data=injTrOnL;
ngood=data.ngood;
ngood[ngood==0]==1;

nses=data.ses-(data.new_mes/(ngood)**(0.5))

plt.figure()
n,bins,patches=plt.hist(nses,bins=np.arange(-10,10,.1),normed=True,histtype='step',label='Injections')


m,s=stats.norm.fit(nses)
print m,s
pdf_g=stats.norm.pdf(bins,m,s)
plt.plot(bins,pdf_g,label="Norm")

mu,stdev = stats.norm.fit(nses[(nses>-3) & (nses<3)])
print mu,stdev
pdf_gamma =stats.norm.pdf(bins,mu,stdev)
plt.plot(bins-mu,pdf_gamma,label="Norm Fit")
plt.legend()
plt.show()
plt.title('SES-<ses> for Injections')
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/marshallAnalysis/sesDistribution.png')

plt.figure()
plt.plot(nses,data.marshall,'.')
want=(data['sigbic']<1e6) & (data['model']=='offsetOnly')
plt.figure();plt.plot(nses[want],data.marshall[want],'.')
plt.xlabel('ses-<ses>')
plt.ylabel('marshall offset only scores')
plt.title('Injections')
plt.vlines(-2*stdev,0,40,colors='g')
plt.vlines(-1*stdev,0,40,colors='g')
plt.vlines(-3*stdev,0,40,colors='g')
plt.hlines(10,-6,4,colors='g')
plt.xlim((-8,8))
plt.ylim((0,40))
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/marshallAnalysis/sesMarshallScore.png')