# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 15:40:49 2016

@author: smullall
"""

#Script to read in the output of the DR25 ops robovetter file.
#And output lists of TCEs for vetting.

import rvIO as rvIO
import pandas as p

rvfile='/soc/nfs/so-nfs/DR25/OPS/RoboVet/RoboVetterOut-OPS.txt'

rvdata=rvIO.readRVDispFile(rvfile)

tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/Q1Q17DR25TCEs.txt'
tcedata=rvIO.readTceInfo(tcefile)
result=p.concat([rvdata,tcedata],axis=1,join='outer')


fedfile='/home/smullall/Kepler/RoboVetter/DR25/OPS/koimatch_DR25_03032016.txt'
cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'

data=rvIO.createAllResults(rvfile,tcefile,fedfile,cumfile)


#
import numpy as np
import csv

outfile='/home/smullall/Kepler/RoboVetter/DR25/Sample/sampleMesPer.txt'
fid=open(outfile,'w')
writer=csv.writer(fid)

tces=data.index;
fullmergdata=p.DataFrame()

Nset=data['N']==1
SCset=(data['S']==1) | (data['C']==1) & (~Nset)
PCset=data['disp']=='PC'

pcuts=(0,10,300,700)
mescuts=(7,10,20,1000)

counter=(0,1,2)

neach=10  #Number of each type

Ncount=np.zeros((3,3))
SCcount=np.zeros((3,3))
PCcount=np.zeros((3,3))

for i in counter:
    for j in counter:
        
        want=(data['period']>pcuts[i]) & (data['period']<pcuts[i+1]) & \
             (data['mes']>mescuts[j]) & (data['mes']<mescuts[j+1])
             
        
        Ncount[i,j]=len(want[want & Nset])
        SCcount[i,j]=len(want[want & SCset])
        PCcount[i,j]=len(want[want & PCset])
        
        NTces=tces[want & Nset]
        SCTces=tces[want & SCset]
        PCTces=tces[want & PCset]        
        
        whN=np.random.choice(np.where(NTces)[0],neach,replace=False)
        whSC=np.random.choice(np.where(SCTces)[0],neach,replace=False)        
        whPC=np.random.choice(np.where(PCTces)[0],neach,replace=False)
        
        dataN=data.loc[NTces[whN]]
        dataSC=data.loc[SCTces[whSC]]
        dataPC=data.loc[PCTces[whPC]]
        mergdata=p.concat([dataN,dataSC,dataPC])
        
        fullmergdata=fullmergdata.append(mergdata)
        

print Ncount
print SCcount
print PCcount

finalData=fullmergdata.drop(['kic','pn','rprstar','arstar','snr','duration','a','pratio','fed','dr24disp','dr24flag','dr24comment'],axis=1)

outfile='/home/smullall/Kepler/RoboVetter/DR25/Sample/sample1.csv'
finalData.to_csv(outfile)


#%%
#List all those in the most interesting realm
outfile='/home/smullall/Kepler/RoboVetter/DR25/Sample/sample-hz-2.csv'

want=(data['srad']<2) & (data['srad']>0.5) & (data['rplanet']>0.5) & (data['rplanet']<2) & (data['tstar']>4000)


tcelist=tces[want & (data['N']==0)]
#| ((data['N']==1) & (data['C']==1))

hzdata=data.loc[tcelist]

finalHz=hzdata.drop(['kic','pn','rprstar','arstar','snr','duration','a','pratio','fed','dr24disp','dr24flag','dr24comment'],axis=1)

finalHz.to_csv(outfile)
 
len(tcelist)

