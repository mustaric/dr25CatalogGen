# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 14:40:33 2016

@author: smullall
"""

import rvIO as io
import numpy as np
import pandas as p

#Get list of new KOIs
#Or other intersting parameter spaces

#Get list of new KOIs, preferably not on Dirty Multis
tfile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
cumfile='/soc/nfs/so-nfs/DR25/other/nexsciCumulative.csv'
fedFile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
opsfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61884/RoboVetterOut-OPS.txt'
opsDataCum=io.createAllResults(opsfile,tfile,fedFile,cumfile)

dr25cumfile='/home/smullall/Kepler/cumulativeTable/dr25/dr25-status.csv'
cumData=io.readcumulativeTable(dr25cumfile)

onkoi=np.array( map(lambda x: len( set(cumData['kicCum']) & set([x]) ) > 0,opsDataCum['kic']), dtype=bool)

opsDataCum['onkoi']=onkoi

want=(opsDataCum['N']==0) & (opsDataCum['onkoi']==0)

opsDataCum[want]
 opsDataCum[want].sort('mes').to_csv("/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61884/newNewKois.csv")
 
 
 
#%%
# To get dataframe with all of the ops,inj and inversion data.
 id='r61884'
injfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-INJ-PlanetOn.txt'
opsfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-OPS.txt'
invfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-INV.txt'

tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
injtfile='/soc/nfs/so-nfs/DR25/INJ/DATA/TCEs.txt'
invtfile='/soc/nfs/so-nfs/DR25/INV/DATA/TCEs.txt'

outinj='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsINJ-PlanetOn-r61884.csv'
outinv='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsINV-r61884.csv'
outops='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsOPS-r61884.csv' 
 
id='r61884'
opsfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-OPS.txt'
tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'
fedfile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
rvdata=io.createAllResults(opsfile,tcefile,fedfile,cumfile)
rvmetrics,a,b,c=io.loadRVInOut(61884)
merged=p.merge(rvdata,rvmetrics,how='inner', right_index=True,left_index=True,suffixes=("","y"))
flags=p.read_csv(outops,comment='#', header=0, delimiter=',',index_col='tce')
OpsFlags=p.merge(merged,flags,how='inner', right_index=True,left_index=True,suffixes=("","fl"))


id='r61884'
invdata,a,b,c=io.loadRVInOut(61884,type="INV")
invmetrics=io.createRVResults(invfile,invtfile)
merged=p.merge(invdata,invmetrics,how='inner', right_index=True,left_index=True,suffixes=("","y"))
flags=p.read_csv(outinv,comment='#', header=0, delimiter=',',index_col='tce')
InvFlags=p.merge(merged,flags,how='inner', right_index=True,left_index=True,suffixes=("","fl"))

injmetrics,a,b,c=io.loadRVInOut(61884,"INJ-PlanetOn")
injdata=io.createRVResults(injfile,injtfile)
merged=p.merge(injdata,injmetrics,how='inner', right_index=True,left_index=True,suffixes=("","y"))
flags=p.read_csv(outinj,comment='#', header=0, delimiter=',',index_col='tce')
InjFlags=p.merge(merged,flags,how='inner', right_index=True,left_index=True,suffixes=("","fl"))

#%%
#Create some new columns for the sigpri/Fred
OF=OpsFlags
OF['sigpriFred_alt']=OF['σ_pri_alt']/OF['F_red_alt'] - OF['F_red_alt']
OF['sigpriFred_dv']=OpsFlags['σ_pri_dv']/OpsFlags['F_red_dv'] - OF['F_red_dv']
OF['sig-priminuster_alt']=OF['σ_pri_alt']-OF['σ_ter_alt']-OF['σ_fa2_alt']
OF['sig-priminuster_dv']=OF['σ_pri_dv']-OF['σ_ter_dv']-OF['σ_fa2_dv']
OF['sig-secminuspos_alt']=OF['σ_sec_alt']-OF['σ_pos_alt']-OF['σ_fa2_alt']
OF['sig-secminuspos_dv']=OF['σ_sec_dv']-OF['σ_pos_dv']-OF['σ_fa2_dv']
OpsFlags=OF

OF=InvFlags
OF['sigpriFred_alt']=OF['σ_pri_alt']/OF['F_red_alt'] - OF['F_red_alt']
OF['sigpriFred_dv']=OpsFlags['σ_pri_dv']/OpsFlags['F_red_dv'] - OF['F_red_dv']
OF['sig-priminuster_alt']=OF['σ_pri_alt']-OF['σ_ter_alt']-OF['σ_fa2_alt']
OF['sig-priminuster_dv']=OF['σ_pri_dv']-OF['σ_ter_dv']-OF['σ_fa2_dv']
OF['sig-secminuspos_alt']=OF['σ_sec_alt']-OF['σ_pos_alt']-OF['σ_fa2_alt']
OF['sig-secminuspos_dv']=OF['σ_sec_dv']-OF['σ_pos_dv']-OF['σ_fa2_dv']
InvFlags=OF

OF=InjFlags
OF['sigpriFred_alt']=OF['σ_pri_alt']/OF['F_red_alt'] - OF['F_red_alt']
OF['sigpriFred_dv']=OpsFlags['σ_pri_dv']/OpsFlags['F_red_dv'] - OF['F_red_dv']
OF['sig-priminuster_alt']=OF['σ_pri_alt']-OF['σ_ter_alt']-OF['σ_fa2_alt']
OF['sig-priminuster_dv']=OF['σ_pri_dv']-OF['σ_ter_dv']-OF['σ_fa2_dv']
OF['sig-secminuspos_alt']=OF['σ_sec_alt']-OF['σ_pos_alt']-OF['σ_fa2_alt']
OF['sig-secminuspos_dv']=OF['σ_sec_dv']-OF['σ_pos_dv']-OF['σ_fa2_dv']
InjFlags=OF

#%%
from glue import qglue

app=qglue(ops=OpsFlags,inv=InvFlags,inj=InjFlags)

#%%


