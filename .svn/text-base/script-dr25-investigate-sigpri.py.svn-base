# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 16:23:31 2016

@author: smullall
"""

import rvIO as io
import matplotlib.pylab as plt

rev=61710
#Read in a particular revion of the OPS data
opsRVAll,dfin,dfout,tcedf=io.loadRVInOut(rev)

#Read in a particular revision of the INJ data
injRVAll,dfin,dfout,tcedf=io.loadRVInOut(rev,type="INJ-PlanetOn")

#Read in a particular revision of the INV data
invRVAll,dfin,dfout,tcedf=io.loadRVInOut(rev,type="INV")

#%%
opspcs=opsRVAll['disp']=='PC'
injpcs=injRVAll['disp']=='PC'
invpcs=invRVAll['disp']=='PC'

#opsdvspri=(opsRVAll['σ_pri_dv']<=0) & (opsRVAll['σ_pri_alt']<=0)
#injdvspri=(injRVAll['σ_pri_dv']<=0) & (injRVAll['σ_pri_alt']<=0)
#invdvspri=(invRVAll['σ_pri_dv']<=0) & (invRVAll['σ_pri_alt']<=0)

opsdvspri=(opsRVAll['σ_pri_alt']<=0)
injdvspri=(injRVAll['σ_pri_alt']<=0)
invdvspri=(invRVAll['σ_pri_alt']<=0)

plt.figure()
plt.plot(opsRVAll[opspcs & opsdvspri]['period'],opsRVAll[opspcs & opsdvspri]['mes'],'.')

plt.figure()
plt.plot(invRVAll[invpcs & invdvspri]['period'],invRVAll[invpcs & invdvspri]['mes'],'.')
print len(invRVAll[invpcs & invdvspri]['period'])
print len(injRVAll[injpcs & injdvspri]['period'])
print len(opsRVAll[opspcs & opsdvspri]['period'])