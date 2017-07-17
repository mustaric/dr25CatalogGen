# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 11:24:06 2016

@author: smullall
"""

import rvIO as io
import pandas as p
import matplotlib.pyplot as plt


file1='/home/smullall/Kepler/cumulativeTable/q1q12/pars-xml.txt';
file2='/home/smullall/Kepler/cumulativeTable/q1q16/pars-xml.txt';
file3='/home/smullall/Kepler/cumulativeTable/q1q17/pars-xml.txt';

rvfile='/soc/nfs/so-nfs/DR25/OPS/RoboVet/RoboVetterOut-OPS.txt'
tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/Q1Q17DR25TCEs.txt'
fedfile='/home/smullall/Kepler/RoboVetter/DR25/OPS/koimatch_DR25_03032016.txt'
cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'

datadr25=io.createAllResults(rvfile,tcefile,fedfile,cumfile)
iskoi=datadr25['koi'].isnull()

tce25=datadr25[~iskoi]
koi25=tce25.set_index('koi')


data12=p.read_csv(file1,sep='|',header=0,comment='#',index_col='koiName')
data12.rename(columns={'koi_max_mult_ev':'mes12'},inplace=True)
data16=p.read_csv(file2,sep='|',header=0,comment='#',index_col='koiName')
data16.rename(columns={'koi_max_mult_ev':'mes16'},inplace=True)
data17=p.read_csv(file3,sep='|',header=0,comment='#',index_col='koiName')
data17.rename(columns={'koi_max_mult_ev':'mes17'},inplace=True)


merge1=koi25.merge(data17,how='inner',right_index=True,left_index=True)
merge2=merge1.merge(data16,how='inner',right_index=True,left_index=True)
merge3=merge2.merge(data12,how='inner',right_index=True,left_index=True)

mesdf=merge3[['mes','mes17','mes16','mes12']]

messtd=mesdf.std(axis=1,skipna=True)
mesavg=mesdf.mean(axis=1,skipna=True)
lowmes=mesavg<30;
relvar=messtd[lowmes]/mesavg[lowmes]

plt.plot(mesavg,messtd,'b.')
 

 