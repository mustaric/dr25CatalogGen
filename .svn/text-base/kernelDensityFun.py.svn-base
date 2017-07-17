# -*- coding: utf-8 -*-
"""
Created on Sat Dec 24 14:07:23 2016

@author: sthompson
"""

import rvIO as io
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

from astroML.density_estimation import KNeighborsDensity
from astroML.plotting import hist

from matplotlib.colors import LogNorm


def kNeighbors(rawdata,k, nbins, method='bayesian', mins=None, maxs=None,sctype='linear'):
    """
    Do the k nearest neighbors version of a kernel density estimation on the data.
    First I normalize the N dimensions so that they are on similar scales.
    rawdata is (ndatapts, Ndimensions)
    mins and maxs is array the length of nDim of rawdata
    They are used to rescale the data.
    nbins has same number of dimensions and gives the number of bins in each dimension.
    """
    Ndim=len(rawdata[0,:])
    if mins ==None:
        mins=np.zeros(Ndim)
        for i in np.arange(0,Ndim):
            mins[i]=np.min(rawdata[:,i])

    if maxs==None:    
        maxs=np.ones(Ndim)
        for i in np.arange(0,Ndim):
            maxs[i]=np.max(rawdata[:,i])
    
    data=rescale(rawdata,mins,maxs,sctype=sctype)
    #Create list of arrays to give mesh grid to return.
    bins=list()    
    for i in np.arange(0,Ndim):
        onebin=np.linspace(mins[i],maxs[i],nbins[i])
        bins.append(onebin)
    
    thebins=np.vstack(map(np.ravel,np.meshgrid(*bins))).T
    print thebins.shape
    
    nbrs = KNeighborsDensity(method, n_neighbors=k).fit(data)
    dens=nbrs.eval(thebins)/(data.size)
    
    dens_n=dens.reshape(*nbins)
    
    #dens_unscale=unrescale(dens_n,mins,maxs,sctype='linear')    
    
    
    return dens_n,bins

def rescale(data,mins,maxs,sctype='linear'):
    """
    Rescale data from 0-1 either linear or log
    """
    newdata=data    
    Ndim=len(data[0,:])    

    if sctype=='linear':    
        for i in np.arange(0,Ndim):
            wid=maxs[i]-mins[i]
            newdata[:,i]= (data[:,i]-mins[i])/wid
    else:
        print "warning: sctype %s not supported by rescale." % sctype
    
    return newdata
        
def unrescale(data,mins,maxs,sctype='linear'):
    """
    Do the inverse of rescale.
    """

    newdata=data    
    Ndim=len(data[0,:])    

    if sctype=='linear':    
        for i in np.arange(0,Ndim):
            wid=maxs[i]-mins[i]
            newdata[:,i]=data[:,i]*wid +mins[i]
    else:
        print "warning: sctype %s not supported by unrescale." % sctype
    
    return newdata

def plot_2dkde(dens,mins,maxs):
    
    plt.imshow(dens,origin='lower',norm=LogNorm(), \
     extent=[mins[0],maxs[0],mins[1],maxs[1]],aspect='auto',cmap='rainbow')