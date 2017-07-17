# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 15:52:30 2016

@author: smullall
"""
import matplotlib.pyplot as plt
import createRVPerfMetrics as pm 
import rvIO as io

svnid=62163
ss1data,a,b,c=io.loadRVInOut(svnid,type='SS1')
cleanSS1=pm.trimINVData(ss1data,'SS1')
invdata,a,b,c=io.loadRVInOut(svnid,type='INV')
cleanInv=pm.trimINVData(invdata,'INV')


#%%

plt.figure(1)
plt.clf()
want=cleanSS1.N==0
plt.hist(cleanSS1.score[want],histtype='step',label='SS1 N=0')

want=cleanInv.N==0
plt.hist(cleanInv.score[want],histtype='step',color='r',label='INV N=0')
plt.legend

plt.figure(2)
plt.clf()
want=cleanSS1.disp=='PC'
plt.hist(cleanSS1.score[want],histtype='step',color='b',label='SS1 PC')

want=cleanInv.disp=='PC'
plt.hist(cleanInv.score[want],histtype='step',color='r',label='INV PC')

plt.legend()
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r'+ str(62163) + '/FalseAlarmAnalysis/pcHist.png')

#%%

