# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 12:41:50 2016

@author: smullall
"""
import numpy as np

def cutMultis(sampleDF,allTcesDF,pnStart):
    """

    Make a cut on multi-planets. But don't remove those
    whose earlier    
    
    Parameters
    -----
    sample = tce data frame of just those TCEs you want the array returned for
    allTces == data frame of all TCEs found in that run
    pnStart == the planet number at which you should start making the cut at
    
    """

    print sampleDF.index
    
    MultiCut=np.zeros(sampleDF['period'].count())
    
    for i,tce in enumerate(sampleDF.index):
        #print i,tce,sampleDF.loc[tce]['kic']
        if sampleDF.loc[tce]['pn'] >= pnStart:
            #Find those with the same kic number
            kic=sampleDF.loc[tce]['kic']

            samekic=allTcesDF['kic']==kic
            lowerpn=allTcesDF['pn']<pnStart
            Ns=allTcesDF['N'][lowerpn & samekic]
            if len(Ns)>0:
                if ~np.any(Ns==0):
                    MultiCut[i]=1
        
    tfMultiCut = MultiCut==1
        
    return tfMultiCut