# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 08:55:34 2016

@author: smullall
"""

import numpy as np
#import pandas as p
#import rvIO as io
import matplotlib.pyplot as plt
import plotRobovetter as prv
import pandas as p

def passes(data,s=[0.0,1.0]):
    """
    Use disposition and scores to determine which ones pass.
    Returns true false array
    """
    
    passed= ((data['disp'] == 'PC') & (data['score']>=s[0])) | \
            ((data['disp']=='FP') & (data['score']>=s[1]))  
            
    return passed

def getCounts(data,s=[0.0,1.0]):
    """
    Return the number of planet candidates.
    Return PCs with a score cut.
    Return number of KOIs.
    Return the number of EBs (S and C flags)
    
    """
    pcs=passes(data,s=[0.0,1.0])
    spcs=passes(data,s=s)
    koi=data.N==0
    ebs=((data.S==1) | (data.C==1)) & koi
    
    npcs=len(pcs[pcs])
    nspcs=len(spcs[spcs])
    nkoi=len(koi[koi])
    nebs=len(ebs[ebs])

    print "Number PCs: %u" % (npcs)
    print "Number score PCs: %u" % (nspcs)
    print "Number KOIs: %u" % (nkoi)
    print "Number EBs (KOI with S|C flag): %u" % (nebs)    
    
    return npcs,nspcs,nkoi,nebs


    

def tcePlots(result):
    """
    Create a couple plots based on information read from the
    robovetter results file and the tce information file
    """
    
    period=result['period'].as_matrix()
    mes=result['mes'].as_matrix()

    makeKOI=result['N'] == 0
    print len(makeKOI)
    
    fig1=plt.figure()
    prv.plotGrid(period,mes,makeKOI.as_matrix())
    plt.title('Fraction of TCEs that we make into KOIs')
    
    
    period=result['period'][makeKOI]
    mes=result['mes'][makeKOI]
    isPC=result['disp'][makeKOI] == 'PC'

    fig2=plt.figure()
    prv.plotGrid(period,mes,isPC.as_matrix())
    plt.title('Fraction of KOIs that are PCs')    


    period=result['period'][makeKOI]
    mes=result['mes'][makeKOI]
    isPC=result['score'][makeKOI] > 0.95

    fig3=plt.figure()
    prv.plotGrid(period,mes,isPC.as_matrix())
    plt.title('Fraction of KOIs that are >0.95 PCs')

    return fig1,fig2,fig3
    
    
def plotFlagFails(result,flagName):
    """
    Create plots based on certain flags
    """
    
    hasflag=hasMinorFlag(result,flagName)
    
    period=result['period']
    mes=result['mes']
    
    plt.figure()
    prv.plotGrid(period,mes,hasflag.as_matrix())
    plt.title('Fraction of TCEs that have %s flag set' % flagName)
 

def plotPcPc(result):
    """
    Plot based on previous PCs
    """
    
    prevPC=result['dr24disp']=='PC'
    print len(prevPC[prevPC])
    
    period=result['period'][prevPC]
    mes=result['mes'][prevPC]
    isPC=result['disp'][prevPC]=='PC'
    
    plt.figure()
    prv.plotGrid(period,mes,isPC) 
    plt.title('Fraction of dr24-PCs that are made into dr25-PCs')       
    
    
def plotFpFp(result):
    """
    Plot based on previous PCs
    """
    
    prevFP=result['dr24disp']=='FP'
    print len(prevFP[prevFP])
    
    period=result['period'][prevFP]
    mes=result['mes'][prevFP]
    isFP=result['disp'][prevFP]=='FP'
    
    plt.figure()
    prv.plotGrid(period,mes,isFP) 
    plt.title('Fraction of FPs(that federated) that are made into dr25-FPs') 
    
def hasMinorFlag(data,flagName):
    """
    Return true false data frame if the 
    tce has that flagName
    this function is mostly here to remind you how to do it.
    Sets to false if no information
    """
    hasFlag=data['flags'].str.contains(flagName, na=False)
    
    return hasFlag

def tcesWithLowPn(data,pnLimit):
    """
    Return a true false array indicating if the TCE 
    Is on a star with fewer than (or equal to) pnLimit tces
    requried tce as index, pn for planet number and kic columns
    """
    
    #zero indicates that you do not want to count this TCE.
    want=np.zeros(len(data['pn']))
    
    if len(data.index.unique()) != len(data.index):
        print "TCE indicies are not unique."
        return want
    
    for i,tce in enumerate(data.index):
        
        if (data.loc[tce]['pn']) <= pnLimit:
            kic=data.loc[tce]['kic']
            samekic=data['kic']==kic
            maxpn=np.max(data[samekic]['pn'])
            if maxpn <= pnLimit:
                want[i]=1
                
    tfwant=want==1
    
    return tfwant

def estimateReliability(invdata,opsdata,s=(0.0,1.0)):
    """
    Given an inversion run cut down to the population you want to consider.
    Return a reliability number.
    1-Reliability = Nfp/Npc * (1-E)/E
    E is the effectiveness as measured by Inversion.
    """
    
    invwant=((invdata['disp'] == 'PC') & (invdata['score']>s[0])) |\
            ((invdata['disp']=='FP') & (invdata['score']>s[1]))
            
    opsPC=((opsdata['disp'] == 'PC') & (opsdata['score']>s[0])) |\
            ((opsdata['disp']=='FP') & (opsdata['score']>s[1]))   
    
    if len(invwant)>0:
        E=len(invwant[~invwant])/float(len(invwant))
        if (E !=0) & (len(opsPC[opsPC])!=0):
            R= 1-len(opsPC[~opsPC])/len(opsPC[opsPC]) * (1.0-E)/E
        else:   
            R=-1
    else:
        E=-1
        R=-1
        
    return R,E

def trimData(invdata,rtype):
    """
    Trim the inverted/SeasonScrambled data of the known astrophysical events
    Note the hard coded files for the list of things to keep.
    """
    
    if rtype == "INV":
        keepfile="/soc/nfs/so-nfs/DR25/INV/DATA/invTCEClean-20161114.csv"
    elif rtype == "SS1":
        #This is temporary file and should be changed.
        keepfile="/home/smullall/Kepler/RoboVetter/DR25/cleanINV/ss1TCEClean-300day-7mes-Nov2016.csv"
    #keepfile="/soc/nfs/so-nfs/DR25/INV/DATA/invTCEClean-300day-09282016.csv"
    keepinv=p.read_csv(keepfile,header=0,index_col='tce',comment='#')
    newdata=p.merge(invdata,keepinv,how="left",left_index=True,right_index=True)
    returnData=newdata[newdata['keep']]    
    
    return returnData

def arrayReliability(fps,pcs,eff):
    """
    Given an array of number of fps, number of pcs, and effectiveness
    for the same bins,
    Return the array of reliailbity.
    """
    
    U=fps.astype(float)/pcs.astype(float) * ((1.0-eff)/eff)
    
    R=1.0-U
    
    return R
    
    