# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 14:17:43 2016

@author: smullall
"""

#Code to start playing wtih pulling many samples and 
#thereby estimate the error on the reliability

import numpy as np
import rvIO as io
import statsRobovetter as sr
import pandas as p

svnid=62163

opsdata,a,b,c=io.loadRVInOut(svnid,'OPS')
invdata,a,b,c=io.loadRVInOut(svnid,'INV')
ss1data,a,b,c=io.loadRVInOut(svnid,'SS1')

trimInv=sr.trimData(invdata,'INV')
trimSs1=sr.trimData(ss1data,'SS1')

fpData=p.concat([trimInv,trimSs1],ignore_index=True,copy=True)



#%%
#Draw from the population N times and calculate the reliability

N=10  ;

