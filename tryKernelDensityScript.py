# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 08:33:05 2016

@author: sthompson
"""
import rvIO as io
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

from astroML.density_estimation import KNeighborsDensity
from astroML.plotting import hist

# Scikit-learn 0.14 added sklearn.neighbors.KernelDensity, which is a very
# fast kernel density estimator based on a KD Tree.  We'll use this if
# available (and raise a warning if it isn't).
try:
    from sklearn.neighbors import KernelDensity
    use_sklearn_KDE = True
except:
    import warnings
    warnings.warn("KDE will be removed in astroML version 0.3.  Please "
                  "upgrade to scikit-learn 0.14+ and use "
                  "sklearn.neighbors.KernelDensity.", DeprecationWarning)
    from astroML.density_estimation import KDE
    use_sklearn_KDE = False
    
#----------------------------------------------------------------------
# This function adjusts matplotlib settings for a uniform feel in the textbook.
# Note that with usetex=True, fonts are rendered with LaTeX.  This may
# result in an error if LaTeX is not installed on your system.  In that case,
# you can set usetex to False.
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=False)

#Get the data from the OPS TCEs
#
tcefile='/Users/sthompson/kepler/DR25/Robovetter/Versions/OPS/TCEs.txt'
opsdata=io.readTceInfo(tcefile)

periods=np.log10(opsdata.period)
#periods=opsdata.period
mes=opsdata.mes

nbins=1000
bins=np.linspace(np.log10(0.5),np.log10(800),1000)
#bins=np.linspace(0.5,700,nbins)
kde = KernelDensity(0.01, kernel='gaussian')
kde.fit(np.array(periods).reshape(-1,1))
dens_kde = np.exp(kde.score_samples(bins[:,None])) 

# Compute density with Bayesian nearest neighbors
k=100
N=len(periods)
nbrs = KNeighborsDensity('bayesian', n_neighbors=k).fit(periods[:,None])
dens_nbrs = nbrs.eval(bins[:,None]) / N

plt.figure()
plt.hist(periods,bins=nbins,histtype='step',color='green', normed=True)
plt.plot(bins,dens_kde,'b-',label='KDE')
plt.plot(bins,dens_nbrs,'r-',label='KNN')
plt.legend()
#%%
Nx=200
Ny=100
xmin=0.5
xmax=700
ymin=7
ymax=40
pbins=np.linspace(np.log10(xmin),np.log10(xmax),Nx)
mbins=np.linspace(np.log10(ymin),np.log10(ymax),Ny)

#pbins=np.linspace(0.5,800,Nx)
#mbins=np.linspace(7,1000,Ny)

data=np.array([periods,np.log10(mes)]).transpose()
bin2=np.vstack(map(np.ravel,np.meshgrid(pbins,mbins))).T
kde = KernelDensity(0.02, kernel='gaussian')
log_dens=kde.fit(data).score_samples(bin2)
dens_kde = data.shape[0]*np.exp(log_dens).reshape(Ny,Nx)
#%%
k=10
nbrs = KNeighborsDensity('bayesian', n_neighbors=k).fit(data)
dens2_nbrs=nbrs.eval(bin2)/(data.size)
#%
from matplotlib.colors import LogNorm
plt.figure(1)
plt.clf()
plt.subplot(211)
plt.imshow(dens_kde, origin='lower',norm=LogNorm(), \
           extent=np.log10([xmin,xmax,ymin,ymax]), aspect='auto',cmap='rainbow')
plt.subplot(212)
plt.scatter(data[:,0],data[:,1],s=1,lw=0,c=u'r')
plt.ylim(np.log10([ymin,ymax]))    
plt.xlim(np.log10([xmin,xmax])) 

plt.figure(2)
plt.clf()
plt.subplot(211)
plt.imshow(dens2_nbrs.reshape(Nx,Ny),origin='lower',norm=LogNorm(), \
     extent=np.log10([xmin,xmax,ymin,ymax]),aspect='auto',cmap='rainbow')
plt.subplot(212)
plt.scatter(data[:,0],data[:,1],s=1,lw=0,c=u'k')   
plt.ylim(np.log10([ymin,ymax]))    
plt.xlim(np.log10([xmin,xmax]))  