# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 11:29:39 2016

@author: smullall
"""
import numpy as np
import pandas as p
import matplotlib.pylab as plt

#Script to see how Inversion has changed after rerunning.



#Determine if it is a bad transit
def isGoodTransit(bic,score):
    """
    Take a one line data frame
    Return pass/fail for Marshall
    
    """
    
    threshold=10.0
    if np.float(bic)>1.5e6:
        return True
    elif score<threshold:
        return True
    else:
        return False

def giveName(tce,midpt):
    
    s="%12s-%04u" % (tce,np.int(midpt))
    return s
    
#Read in new Inversion File.
invfile='/soc/nfs/so-nfs/fergal/marshall2.0/inv-2016-12-11p.txt'
inv2=p.read_csv(invfile,delim_whitespace=True,names=['tce','midpt','score','bic','model'],comment='#',index_col=False)

#Read old inversion File.
invfile='/soc/nfs/so-nfs/DR25/INV/DATA/Marshall-INV-2016-09-23-MoreInfo.txt'
inv1=p.read_csv(invfile,delim_whitespace=True,names=['tce','midpt','score', 'bic','model'],comment='#',index_col=False)


marshPass1=np.array(map(lambda x,y: isGoodTransit(x,y),inv1.bic,inv1.score),dtype=bool)
marshPass2=np.array(map(lambda x,y: isGoodTransit(x,y),inv2.bic,inv2.score),dtype=bool)

inv1['pass']=marshPass1
inv2['pass']=marshPass2

name1=map(lambda x,y: giveName(x,y),inv1.tce,inv1.midpt)
name2=map(lambda x,y: giveName(x,y),inv2.tce,inv2.midpt)

inv1['name']=name1
inv2['name']=name2

invtr=p.merge(inv1,inv2,how='inner',left_on='name',right_on='name',suffixes=('_1','_2'))


#Read in new SS1 File.
invfile='/soc/nfs/so-nfs/fergal/marshall2.0/ss1-2016-12-09p.txt'
inv2=p.read_csv(invfile,delim_whitespace=True,names=['tce','midpt','score','bic','model'],comment='#',index_col=False)

#Read old SS1 File.
invfile='/soc/nfs/so-nfs/DR25/SS1/DATA/Marshall-SS1-2016-11-21-MoreInfo.txt'
inv1=p.read_csv(invfile,delim_whitespace=True,names=['tce','midpt','score', 'bic','model'],comment='#',index_col=False)


marshPass1=np.array(map(lambda x,y: isGoodTransit(x,y),inv1.bic,inv1.score),dtype=bool)
marshPass2=np.array(map(lambda x,y: isGoodTransit(x,y),inv2.bic,inv2.score),dtype=bool)

inv1['pass']=marshPass1
inv2['pass']=marshPass2

name1=map(lambda x,y: giveName(x,y),inv1.tce,inv1.midpt)
name2=map(lambda x,y: giveName(x,y),inv2.tce,inv2.midpt)

inv1['name']=name1
inv2['name']=name2

ss1tr=p.merge(inv1,inv2,how='inner',left_on='name',right_on='name',suffixes=('_1','_2'))