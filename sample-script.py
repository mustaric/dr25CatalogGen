# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 11:42:48 2016

@author: smullall
"""

import rvIO as io
import numpy as np
import pandas as p


def passes(data,s=[0.5,0.4]):
    """
    Use disposition and scores to determine which ones pass.
    Returns true false array
    """
    
    passed= ((data['disp'] == 'PC') & (data['score']>=s[0])) | \
            ((data['disp']=='FP') & (data['score']>s[1])) 



    return passed



opsfile='

num=62076

rvdata=io.createAllResults(opsfile,tcefile,fedfile,cumfile)
rvmetrics,a,b,c=io.loadRVInOut(num)
opsdata=p.merge(rvdata,rvmetrics,how='inner', right_index=True,left_index=True,suffixes=("","y"))
#opsdata,a,b,c=io.loadRVInOut(num)
invdata,a,b,c=io.loadRVInOut(num,type='INV')
injdta,a,b,c=io.loadRVInOut(num,type='INJ-PlanetOn')


scores=(0.4,0.3)
hz=(0,2)

want=passes(opsdata,s=scores) & (opsdata.rplanet<2.5) & (opsdata.srad<2.5)

hzones=opsdata[want][['disp','flags','score','tstar','rplanet','mes','period','koi']]

hzones.to_csv('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r62076/hzops.csv')

print hzones

want=passes(invdata,s=scores) & (invdata.rplanet<2.5) & (invdata.srad<2.5)

hzones=invdata[want][['disp','flags','score','tstar','rplanet','mes','period']]

hzones.to_csv('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r62076/hzinv.csv')

print hzones