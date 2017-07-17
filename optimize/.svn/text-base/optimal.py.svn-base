# -*- coding: utf-8 -*-
"""

Compute optimal values for RBV thresholds

Experimental code.

Created on Mon Apr 18 14:22:52 2016

@author: fmullall
"""

__version__ = "$Id$"
__URL__ = "$URL$"

import matplotlib.pyplot as mp
import scipy.optimize as spOpt
import pandas as pd
import numpy as np

import logging

import dr25


def printScores(injData, invData, thresholds=None):
    """Print scores by adding in tests one at a time"""

    colList = """Marshall LPP_DV LPP_ALT SES/MES
                modshift1 modshift2 modshift3
                modshift1alt modshift2alt modshift3alt
                """.split()
    sgnList = np.ones(len(colList))

    if thresholds is None:
        thresholds = [10, -2.5, -2.6, .9, +1, +1, +1, +1, +1, +1]

    if True:
        for nTest in range(len(colList)):

            cl = colList[:nTest+1]
            sl = sgnList[:nTest+1]
            tl = thresholds[:nTest+1]

            def lambdaFunc(thresholds):
                return costFunc(injData, invData, cl, sl, tl)

            comp = 1 - computeFailFraction(injData, cl, sl, tl)
            eff =  computeFailFraction(invData, cl, sl, tl)
            args = (nTest, lambdaFunc(tl), comp, eff)
            print "%i %.4f %.4f %.4f" %(args)

        #Fit at this threshold
        bestFit = fit(lambdaFunc, cl, sl, tl, injData, invData)
        comp = 1 - computeFailFraction(injData, cl, sl, tl)
        eff =  computeFailFraction(invData, cl, sl, tl)
        args = (bestFit.fun, comp, eff)
        print "B %.4f %.4f %.4f" %(args)

        return bestFit


def optimiseThresholds(inj, inv, wgt=1,nTrial=200):

    colList = """Marshall LPP_DV LPP_ALT SES/MES
                modshift1 modshift2 modshift3
                modshift1alt modshift2alt modshift3alt
                """.split()
    #colList = """Marshall SES/MES modshift1 modshift1alt""".split()
    sgnList = np.ones(len(colList))
    
   # tRange = [ [5, 15], [0.2,1.7], [.1, 2], [.1, 2] ]

    tRange= [ [5, 10], [.0015, .005], [.0015, .005], [.6,1.2],
              [.1, 2] , [0, 2],  [0, 2],
              [.1, 2] , [0, 2],  [0, 2]
            ]

    mp.clf()
    mp.axis([0,1,0,1])

    ax = mp.gca()
    for r in [.05, .1, .15]:
        circle = mp.Circle((0, 1), r, color='k', alpha=.2)
        ax.add_artist(circle)

    def lambdaFunc(thresholds):
        return costFunc(inj, inv, colList, sgnList, thresholds, weight=wgt)

    #Create some dummy values for initial condidions
    threshold = np.zeros(len(colList))
    bestFit = fit(lambdaFunc, colList, sgnList, threshold, inj, inv)

    trials=np.zeros([nTrial,len(colList)])
    
    for i in range(nTrial):
        for j in range(len(colList)):
            lwr, upr = tRange[j]
            threshold[j] = lwr + (upr-lwr)*np.random.rand()

        res = fit(lambdaFunc, colList, sgnList, threshold, inj, inv)
        optThres = res.x
        comp = 1 - computeFailFraction(inj, colList, sgnList, optThres)
        eff =  computeFailFraction(inv, colList, sgnList, optThres)
        trials[i,:]=res.x
        if res.fun < bestFit.fun:
            bestFit= res
        mp.plot(1-eff, comp, 'ko')

        if i % 10 == 0:
            #Update the plot
            mp.pause(.01)
            print i, bestFit.fun, bestFit.x
    
    bestOptThres=bestFit.x
    bestComp=1 - computeFailFraction(inj,colList,sgnList,bestOptThres)
    bestEff=computeFailFraction(inv, colList, sgnList, bestOptThres)
    mp.plot(1-bestEff, bestComp, 'rs')

    return bestFit,trials,colList


def fit(func, colList, sgnList, thresholds, injData, invData):

    logger = logging.getLogger("Fit")
    logger.setLevel(logging.INFO)

    res = spOpt.minimize(func, thresholds, method='Nelder-Mead')
    return res

def load():
    injFile= "/home/smullall/Kepler/RoboVetter/DR25/svn/robovet/RoboVetter-Input-INJ-PlanetOn.txt"
    invFile="/home/smullall/Kepler/RoboVetter/DR25/svn/robovet/RoboVetter-Input-INV.txt"
    #injFile = "/net/spmdb-fs/so-nfs/DR25/INJ/RoboVet/RoboVetter-Input-INJ-PlanetOn.txt"
    #invFile = "/net/spmdb-fs/so-nfs/DR25/INV/RoboVet/RoboVetter-Input-INV-vDR24.txt"
        
    injStatus = '/net/spmdb-fs/so-nfs/DR25/INJ/RoboVet/RoboVetterOut-INJ-PlanetOn.txt'
    invStatus = '/home/smullall/Kepler/RoboVetter/DR25/svn/robovet/RoboVetterOut-INV.txt'

    inj = readFile(injFile)
    injFed = dr25.loadFederation(dr25.injFedFile)
    inj = dr25.federate(inj, injFed)
    #inj = trimFalsePositives(inj, injStatus)
    inj = editColumns(inj)


    inv = readFile(invFile)
    invFed = dr25.loadFederation(dr25.invFedFile)
    inv = dr25.antiFederate(inv, invFed)
    #inv = trimFalsePositives(inv, invStatus)
    inv = editColumns(inv)

    return inj, inv


def trimFalsePositives(data, statusFile):
    """Remove TCEs that fail because of S, C or E flags"""

    status = np.loadtxt(statusFile, dtype=str, usecols=(0,1,2,3,4,5,6))
    tceId = dr25.tceStrToId(status[:,0])

    status = pd.DataFrame(status, index=tceId)
    colList = "TCE Score Disposition N S C E".split()
    status.columns = colList

    idx = (status['S'] == "0") & (status['C'] == "0") & (status['E'] == "0")

    return pd.merge(data, status[idx], left_index=True, right_index=True)





def editColumns(data):

    #data['LPP_DV'] = np.log10(data['LPP_DV'])
    #data['LPP_ALT'] = np.log10(data['LPP_ALT'])

    #Add MES/SES column, but only for P>90 days
    data['SES/MES'] = data['SES_max']/ data['MES']
    idx  = data['Period'] < 90
    data.loc[idx, 'MES/SES'] = -1  #Will always pass

    idx = data['Period'] < 150
    data.loc[idx, 'Marshall'] = -100  #Will always pass


    #Modshift cols
    sigPri = data['σ_pri_dv']
    sigTer = data['σ_ter_dv']
    sigPos = data['σ_pos_dv']
    sigFa =  data['σ_fa1_dv']
    sigFa2 = data["σ_fa2_dv"]
    fRed = data['F_red_dv']


    eps = 1e-99   # A very small number
    #DR24-RoboVetter.cpp:326, :333 and :340
    data['modshift1'] = fRed * sigFa / (sigPri + eps) * (sigPri > 0)
    data['modshift2'] = sigFa2 / (sigPri - sigTer + eps) * (sigPri > 0) * (sigTer > 0)
    data['modshift3'] = sigFa2 / (sigPri - sigPos) * (sigPri > 0) * (sigPos > 0)


    #Modshift alternative detrending cols
    sigPri = data['σ_pri_alt']
    sigTer = data['σ_ter_alt']
    sigPos = data['σ_pos_alt']
    sigFa =  data['σ_fa1_alt']
    sigFa2 = data["σ_fa2_alt"]
    fRed = data['F_red_alt']

    #DR24-RoboVetter.cpp:326, :333 and :340
    data['modshift1alt'] = fRed * sigFa / (sigPri + eps)    * (sigPri > 0)
    data['modshift2alt'] = sigFa2 / (sigPri - sigTer + eps) * (sigPri > 0) * (sigTer > 0)
    data['modshift3alt'] = sigFa2 / (sigPri - sigPos + eps) * (sigPri > 0) * (sigPos > 0)



    return data


def readFile(filename):

    data = np.loadtxt(filename, dtype=str)
    fp = open(filename)
    hdr = fp.readline()
    hdr = fp.readline()
    data[:,0] = dr25.tceStrToId(data[:,0])
    fp.close()


    hdr = hdr[1:].split()
    df = pd.DataFrame(data.astype(float), columns=hdr)
    df.index = df.iloc[:,0]
    return df





def costFunc(injData, invData, colList, sgnList, thresholdList, weight=1):
    """
    Compute the cost of choosing a given set of thresholds.
    Cost is defined as the distance from perfect completeness and perfect
    reliability according to some distance metric

    Inputs:
    ----------
    injData, invData:
        Pandas Dataframe, one TCE per row, of injection TCEs and inversion TCEs
    colList
        List of column names used to make pass fail decision. Must be the
        same for injData and invData
    sgnList
        List of signs. If sgn[i] == +1 a TCE fails if its value is greater
        than threshold. If sng[i] == -1, the TCE fails if the value is less than
        threshold.
    thresholdList
        List of values establishing the pass fail threshold.

    Optional Inputs:
    -------------------
    weight
        (float) How much reliability should be weighted over completness.
        1 means they are equally weighted. .5 implies we'd rather have
        two false alarms in the final catalogue than lose a single real
        transit.

    Returns:
    -----------
    sqrt( c**2 + weight*r**2)
    where c is the fraction of injected TCEs passed, and r is the
    fraction of inversion TCEs failed.
    """

#    comp = computeFailFraction(injData, colList, sgnList, thresholdList)
#    eff = 1 - computeFailFraction(invData, colList, sgnList, thresholdList)
    comp = computeFailFraction(injData, colList, sgnList, thresholdList)
    eff = 1 - computeFailFraction(invData, colList, sgnList, thresholdList)

#    print "Cp", comp, eff, np.hypot(comp, eff)
    w = np.sqrt(weight)
    return np.hypot(comp, w*eff)


def computeFailFraction(data, colList, sgnList, thresholdList):
    """Compute fraction of TCEs that fail.

    Inputs:
    ----------
    data:
        Pandas Dataframe, one TCE per row.
    colList
        List of column names used to make pass fail decision
    sgnList
        List of signs. If sgn[i] == +1 a TCE fails if its value is greater
        than threshold. If sng[i] == -1, the TCE fails if the value is less than
        threshold.
    thresholdList
        List of values establishing the pass fail threshold.

    Returns:
    ----------
    Fraction of TCEs that fail, a value between [0,1]
    """

    assert len(colList) == len(sgnList)
    assert len(colList) == len(thresholdList)
    assert( np.all( np.fabs(sgnList) == 1) )

    nTce = len(data)
    nFeatures = len(colList)
    fpIdx = np.zeros( nTce, dtype=bool)

    for i in range(nFeatures):
        y = data[colList[i]]
        fpIdx |= sgnList[i] * (y - thresholdList[i]) > 0

    return np.sum(fpIdx) / float(nTce)



