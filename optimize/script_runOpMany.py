# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 16:31:41 2016

@author: smullall
"""

import optimal as op
import matplotlib.pyplot as plt
import corner
import numpy as np

mesbin=[7,10,20]
mesend=[10,20,1000]
perbin=[.5,10,200,300,400]
perend=[10,200,500,400,500]

inj,inv=op.load()

out="/home/smullall/Kepler/RoboVetter/DR25/optimize/bestFits-Aug4.txt"
fid=open(out,"w")

x=np.zeros([9,10])
v=0

for i,m in enumerate(mesbin):
    for j,p in enumerate(perbin):
    

        wantinj=(inj['Period']>perbin[j]) & (inj['Period']<perend[j]) & (inj['MES']>mesbin[i]) & (inj['MES']<=mesend[i])
        wantinv=(inv['Period']>perbin[j]) & (inv['Period']<perend[j]) & (inv['MES']>mesbin[i]) & (inv['MES']<=mesend[i])

        name='/home/smullall/Kepler/RoboVetter/DR25/optimize/picks-%i-%i' %(int(m),int(p))    
        bestFit,trials,cols=op.optimiseThresholds(inj[wantinj],inv[wantinv],wgt=1.0,nTrial=700)
        plt.title(name)
        plt.savefig(name)      

        print m,p
        print bestFit.x
        fid.write(str(bestFit.x))

        x[v,:]=bestFit.x

        
        v=v+1
        
        try:
            corner.corner(bestFit,labels=cols,show_titles=True)
        
            name='/home/smullall/Kepler/RoboVetter/DR25/optimize/corner-%i-%i.png' % (int(m),int(p))
            plt.savefig(name)
        except:
            print "Corner Failed"

        
#for i in np.linspace(0,10):
#    plt.imshow(np.resize(x[:,i],(len(perbin),len(mesbin))),interpolation="nearest",cmap=plt.cm.YlGnBu_r, origin="bottom") 
#    ax = plt.gca()
#    xLabels = map(lambda x: "%i" %(x), perbin)
#    ax.xaxis.set_ticks(np.arange(len(perbin)))
#    ax.xaxis.set_ticklabels(xLabels)
#    plt.xlabel("Period (days)")
#
#    yLabels = map(lambda x: "%i" %(x), mesbin)
#    ax.yaxis.set_ticks(np.arange(len(mesbin)))
#    ax.yaxis.set_ticklabels(yLabels)
#    plt.ylabel("MES")
        
        
        
