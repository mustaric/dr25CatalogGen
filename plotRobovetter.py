# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 08:55:34 2016

@author: smullall
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter, MaxNLocator   #You need this below
import statsRobovetter as rvs


def plotGrid(period, mes, pcIdx, drange=(0,100),xBins = [0, 10, 200, 500],yBins = [7, 10, 20, 2000],clabel="Percent Passed"):
    """
    Originally Fergal's code to plot 2d hist of period and mes.
    Use plt.clim(lwr, upr) to change the colour scheme
    pcIdx is a true false array where True means it 'passed'
    """
    
    assert(len(period) == len(mes))
    assert(len(period) == len(pcIdx))

    allHist = compute2dHist(period, mes, xBins, yBins)
    pcHist = compute2dHist(period[pcIdx], mes[pcIdx], xBins, yBins)

    #Add a little to denominate to prevent Nans
    percent = 100*pcHist/(allHist + 1e-7)

    #plt.clf()
    extent = [0, len(xBins)-1, 0, len(yBins)-1]
    plt.imshow(percent, interpolation="nearest", cmap=plt.cm.YlGnBu_r, \
        origin="bottom", extent=extent)
    cb = plt.colorbar()
    cb.set_label(clabel)

    writeLabels(percent, allHist, drange)

    ax = plt.gca()
    xLabels = map(lambda x: "%i" %(x), xBins)
    ax.xaxis.set_ticks(np.arange(len(xBins)))
    ax.xaxis.set_ticklabels(xLabels)
    plt.xlabel("Period (days)")

    yLabels = map(lambda x: "%i" %(x), yBins)
    ax.yaxis.set_ticks(np.arange(len(yBins)))
    ax.yaxis.set_ticklabels(yLabels)
    plt.ylabel("MES") 
    
    plt.clim(drange)
    
    return percent
    
    
    
def writeLabels(percent, number, drange):
    nR, nC = percent.shape
    assert(percent.shape == number.shape)

    midVal = .5*(np.max(drange) + np.min(drange))
    for r in range(nR):
        for c in range(nC):
            text = "%.1f%%\n%i" %(percent[r,c], number[r,c])

            cType = 'k'
            if percent[r, c] < midVal:
                cType='w'
            plt.text(c+.5, r+.5, text, va="center", ha="center", color=cType)

def writeLabels2(percent, number,number2, drange):
    nR, nC = percent.shape
    assert(percent.shape == number.shape)

    midVal = .5*(np.max(drange) + np.min(drange))
    for r in range(nR):
        for c in range(nC):
            text = "%.1f%%\n%i/%i" %(percent[r,c], number[r,c],number2[r,c])

            cType = 'k'
            if percent[r, c] < midVal:
                cType='w'
            plt.text(c+.5, r+.5, text, va="center", ha="center", color=cType,fontsize=10)


def compute2dHist(x, y, xBins, yBins):
    h, xedege, yedge = np.histogram2d(x,y, [xBins, yBins])
    h = h.transpose()
    return h 



def passes(data,s=[0.0,1.0]):
    """
    Use disposition and scores to determine which ones pass.
    Returns true false array
    """
    
    passed= ((data['disp'] == 'PC') & (data['score']>=s[0])) | \
            ((data['disp']=='FP') & (data['score']>=s[1]))  
    return passed

def plotReliability(invdata,opsdata,xBins=[0, 10, 200, 500],yBins=[7, 10, 20, 2000],drange=(1,99),key='mes',limit=0.0,s=[0.0,1.0]):
    """
    Plot grid of reliability numbers
    Similar to plotgrid above, but needs that to run twice to get the numbers 
    we need the dataframes for inversion and ops.
    """
    
    plt.figure(figsize=(10,6.5))
    plt.subplots_adjust(hspace=0.55, wspace=0.3)
    opsPeriod=opsdata['period']
    opsMes=opsdata['mes']
    opsDisp= passes(opsdata,s)  & (opsdata[key]>limit)

    allHistops = compute2dHist(opsPeriod, opsMes, xBins, yBins)
    pcHistops = compute2dHist(opsPeriod[opsDisp], opsMes[opsDisp], xBins, yBins)

    invPeriod=invdata['period']
    invMes=invdata['mes']
    invDisp=passes(invdata,s) & (invdata[key]>limit)

    allHistinv = compute2dHist(invPeriod, invMes, xBins, yBins)
    pcHistinv = compute2dHist(invPeriod[invDisp], invMes[invDisp], xBins, yBins)

    E = 1 - pcHistinv/allHistinv
    
    R = 100.0 *(1.0 - ((allHistops-pcHistops)/pcHistops) * (1.0-E)/E)
    print E
    #import pdb
    #pdb.set_trace()

    plt.subplot(121)
    extent = [0, len(xBins)-1, 0, len(yBins)-1]
    plt.imshow(R, interpolation="nearest", cmap=plt.cm.YlGnBu_r, \
        origin="bottom", extent=extent)
    cb = plt.colorbar()
    cb.set_label("Catalog Reliability")
    
    writeLabels2(R, pcHistops,allHistops-pcHistops, drange)
    
    ax = plt.gca()
    xLabels = map(lambda x: "%i" %(x), xBins)
    ax.xaxis.set_ticks(np.arange(len(xBins)))
    ax.xaxis.set_ticklabels(xLabels)
    plt.xlabel("Period (days)")

    yLabels = map(lambda x: "%i" %(x), yBins)
    ax.yaxis.set_ticks(np.arange(len(yBins)))
    ax.yaxis.set_ticklabels(yLabels)
    plt.ylabel("MES") 
    plt.title("Reliability\n#PCs/#FPs")
    plt.clim(drange)
    
    plt.subplot(122)
    extent = [0, len(xBins)-1, 0, len(yBins)-1]
    plt.imshow(E*100, interpolation="nearest", cmap=plt.cm.YlGnBu_r, \
        origin="bottom", extent=extent)
    cb = plt.colorbar()
    cb.set_label("Catalog Effectiveness using Inversion")
    
    writeLabels(E*100, allHistinv, drange)
    
    ax = plt.gca()
    xLabels = map(lambda x: "%i" %(x), xBins)
    ax.xaxis.set_ticks(np.arange(len(xBins)))
    ax.xaxis.set_ticklabels(xLabels)
    plt.xlabel("Period (days)")

    yLabels = map(lambda x: "%i" %(x), yBins)
    ax.yaxis.set_ticks(np.arange(len(yBins)))
    ax.yaxis.set_ticklabels(yLabels)
    plt.ylabel("MES") 
    plt.title('Effectiveness\n#InvTces')
    
    plt.clim(drange)

def plotOnlyReliability(invdata,opsdata,xmetric='period',ymetric='mes',xBins=[0, 10, 200, 500],yBins=[7, 10, 20, 2000],drange=(1,99),\
                        atitle="Reliability\n#PCs/#FPs",key='score',limit=0.0,s=[0.0,1.0],xlabel='metric',ylabel='metric'):
    """
    Plot grid of reliability numbers
    Similar to plotgrid above, but needs that to run twice to get the numbers 
    we need the dataframes for inversion and ops.
    """
    
    opsPeriod=opsdata[xmetric]
    opsMes=opsdata[ymetric]
    opsDisp=passes(opsdata,s) #  & (opsdata[key]>limit) -- Commented out, why here?

    allHistops = compute2dHist(opsPeriod, opsMes, xBins, yBins)
    pcHistops = compute2dHist(opsPeriod[opsDisp], opsMes[opsDisp], xBins, yBins)

    invPeriod=invdata[xmetric]
    invMes=invdata[ymetric]
    invDisp=passes(invdata,s) # & (invdata[key]>limit) -- commented out, why here?

    allHistinv = compute2dHist(invPeriod, invMes, xBins, yBins)
    pcHistinv = compute2dHist(invPeriod[invDisp], invMes[invDisp], xBins, yBins)

    E = 1 - pcHistinv/allHistinv
    
    R = 100.0 *(1.0 - ((allHistops-pcHistops)/pcHistops) * (1.0-E)/E)
    #print E
    #import pdb
    #pdb.set_trace()

    extent = [0, len(xBins)-1, 0, len(yBins)-1]
    plt.imshow(R, interpolation="nearest", cmap=plt.cm.YlGnBu_r, \
        origin="bottom", extent=extent)
    cb = plt.colorbar()
    cb.set_label("Catalog Reliability")
    
    writeLabels2(R, pcHistops,allHistops-pcHistops, drange)
    
    ax = plt.gca()
    xLabels = map(lambda x: "%i" %(x), xBins)
    ax.xaxis.set_ticks(np.arange(len(xBins)))
    ax.xaxis.set_ticklabels(xLabels)
    if xlabel=='metric':
        plt.xlabel(xmetric)
    else:
        plt.xlabel(xlabel)
        
    yLabels = map(lambda x: "%i" %(x), yBins)
    ax.yaxis.set_ticks(np.arange(len(yBins)))
    ax.yaxis.set_ticklabels(yLabels)
    if ylabel=='metric':
        plt.ylabel(ymetric) 
    else:
        plt.ylabel(ylabel)
    plt.title(atitle,fontsize=11)
    plt.clim(drange)
    
    return R

def numPassed(opsdata,xBins=[0, 10, 200, 500],yBins=[7, 10, 20, 2000],s=[0.0,1.0]):
    """
    Return 2d array of the number of OPS passed in each bin.
    """
    opsPeriod=opsdata['period']
    opsMes=opsdata['mes']
    opsDisp=passes(opsdata,s) #  & (opsdata[key]>limit) -- Commented out, why here?

    pcHistops = compute2dHist(opsPeriod[opsDisp], opsMes[opsDisp], xBins, yBins)
    
    return pcHistops
    

import matplotlib.gridspec as gridspec
import pdb
def plot2dHist(x,y,xlim,ylim,nxbins=50,nybins=50,xlabel="",ylabel="",makelog=False,clim=None,showPoints=False):
    """
    Create a complicated 2d histogram plot with projected histograms onthe side
    """
    
    # Define the locations for the axes
    #left, width = 0.12, 0.55
    #bottom, height = 0.12, 0.55
    #bottom_h = left_h = left+width+0.02
    
    #rect_temperature = [left, bottom, width, height] # dimensions of temp plot
    #rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
    #rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram
    
    gs=gridspec.GridSpec(3,3)    
    plt.figure(1, figsize=(9,9))
    
    axTemperature = plt.subplot(gs[1:3,0:2])  #plt.axes(rect_temperature) # temperature plot
    axHistx = plt.subplot(gs[0,0:2],sharex=axTemperature) #plt.axes(rect_histx) # x histogram ,
    axHisty = plt.subplot(gs[1:3,2], sharey=axTemperature) #plt.axes(rect_histy) # y histogram
    if makelog:
        axHistx.set_xscale=('log')
        axHisty.set_xscale=('log')
        axTemperature.set_xscale('log')
        axTemperature.set_yscale('log')
    #pdb.set_trace()
    
    xbins = np.linspace(start = xlim[0], stop = xlim[1], num = nxbins)
    ybins = np.linspace(start = ylim[0], stop = ylim[1], num = nybins)

    #Plot the histograms
    axHistx.hist(x, bins=xbins, color = 'blue',histtype='step',lw=2)
    axHisty.hist(y, bins=ybins, orientation='horizontal', color = 'blue',histtype='step',lw=2)    
    
    #Set up the histogram bins   
    nbins=(nxbins+nybins)/4
    print nxbins,nybins,nbins
    x2bins = np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/nbins)
    y2bins = np.arange(ylim[0], ylim[1], (ylim[1]-ylim[0])/nbins)
    
    #xcenter = (xbins[0:-1]+xbins[1:])/2.0
    #ycenter = (ybins[0:-1]+ybins[1:])/2.0
    #aspectratio = 1.0*(xlim[1] - xlim[0])/(1.0*(ylim[1] - ylim[0])
    
    H, xedges,yedges = np.histogram2d(y,x,bins=(y2bins,x2bins))
    if clim==None:
        clim=[0,np.max(H)]
    
    # Plot the temperature data
    axTemperature.imshow(H, extent=[xlim[0],xlim[1],ylim[0],ylim[1]],
          interpolation='nearest', origin='lower',aspect="auto", cmap=plt.cm.YlGnBu_r, clim=clim)
          
    #Plot points if requested
    if showPoints:
        print "showing Points"
        axTemperature.plot(x,y,'y.',ms=1)

    #Set up the plot limits
    axTemperature.set_xlim(xlim)
    axTemperature.set_ylim(ylim)  
    
    #Plot the axes labels
    axTemperature.set_xlabel(xlabel,fontsize=18)
    axTemperature.set_ylabel(ylabel,fontsize=18)
    #plt.setp(axTemperature.get_xticklabels(), fontsize=14)
     
    #Set up the histogram limits
    #axHistx.set_xlim( xlim[0], xlim[1] )
    #axHisty.set_ylim( ylim[0], ylim[1] )
    
    #plt.setp(axTemperature.get_yticklabels(),fontsize=12,visible=True)
        
    #Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))
    
    return axHisty,axHistx,ybins,xbins
    
def plot2dHistNorm(x,y,xlim,ylim,nxbins=50,nybins=50,xlabel="",ylabel="",clim=None,showPoints=False):
    """
prv.compute2dHist(x,y,50,50)
    Create a complicated 2d histogram plot with projected histograms onthe side
    """
    
    # Define the locations for the axes
    #left, width = 0.12, 0.55
    #bottom, height = 0.12, 0.55
    #bottom_h = left_h = left+width+0.02
    
    #rect_temperature = [left, bottom, width, height] # dimensions of temp plot
    #rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
    #rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram
    
    gs=gridspec.GridSpec(3,3)    
    plt.figure(figsize=(9,9))
    
    axTemperature = plt.subplot(gs[1:3,0:2])  #plt.axes(rect_temperature) # temperature plot
    axHistx = plt.subplot(gs[0,0:2],sharex=axTemperature) #plt.axes(rect_histx) # x histogram ,
    axHisty = plt.subplot(gs[1:3,2], sharey=axTemperature) #plt.axes(rect_histy) # y histogram
    #pdb.set_trace()
    
    xbins = np.linspace(start = xlim[0], stop = xlim[1], num = nxbins)
    ybins = np.linspace(start = ylim[0], stop = ylim[1], num = nybins)

    #Plot the histograms
    axHistx.hist(x, bins=xbins, color = 'blue',histtype='step',normed=True,lw=2)
    axHisty.hist(y, bins=ybins, orientation='horizontal', color = 'blue',histtype='step',normed=True,lw=2)    
    
    #Set up the histogram bins   
    nbins=(nxbins+nybins)/4
    x2bins = np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/nbins)
    y2bins = np.arange(ylim[0], ylim[1], (ylim[1]-ylim[0])/nbins)
    
    #xcenter = (xbins[0:-1]+xbins[1:])/2.0
    #ycenter = (ybins[0:-1]+ybins[1:])/2.0
    #aspectratio = 1.0*(xlim[1] - xlim[0])/(1.0*(ylim[1] - ylim[0])
    
    H, xedges,yedges = np.histogram2d(y,x,bins=(y2bins,x2bins))
    
    if clim==None:
        clim=[0,np.max(H)]
    
    # Plot the temperature data
    axTemperature.imshow(np.log10(H), extent=[xlim[0],xlim[1],ylim[0],ylim[1]],
          interpolation='nearest', origin='lower',aspect="auto", cmap=plt.cm.YlGnBu_r, clim=clim)

    #Plot points if requested
    if showPoints:
        print "showing Points"
        axTemperature.plot(x,y,'y.',ms=1)
        
    #Set up the plot limits
    axTemperature.set_xlim(xlim)
    axTemperature.set_ylim(ylim)  
    
    #Plot the axes labels
    axTemperature.set_xlabel(xlabel,fontsize=18)
    axTemperature.set_ylabel(ylabel,fontsize=18)
    #plt.setp(axTemperature.get_xticklabels(), fontsize=14)
     
    #Set up the histogram limits
    #axHistx.set_xlim( xlim[0], xlim[1] )
    #axHisty.set_ylim( ylim[0], ylim[1] )
    
    #plt.setp(axTemperature.get_yticklabels(),fontsize=12,visible=True)
        
    #Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))
    
    return axHisty,axHistx,ybins,xbins
    
def plot2dMulti(x,y,xlim,ylim,nxbins=80,nybins=80,\
        xlabel="log(Period)",ylabel="log(MES)",makelog=False,clim=None,showPoints=False):
    """
    Create a 2dhistogram with other histograms overlaid on the axes.
    The zeroth, element in the x and y list will be shown in the 2dpart.   
    """
    color=["red","green","magenta","black","orange"]
    
    (axHisty,axHistx,ybins,xbins)=plot2dHist(x[0],y[0],xlim,ylim,nxbins=nxbins,\
        nybins=nybins,xlabel=xlabel,ylabel=ylabel,makelog=makelog,clim=clim,showPoints=showPoints)
    
    for i in np.arange(1,len(x)):
        
        print i
        
        axHisty.hist(y[i],bins=ybins,orientation='horizontal', color = color[i-1],histtype='step',lw=1)    
    
        axHistx.hist(x[i],bins=xbins,color=color[i-1],histtype='step',lw=1)  
    
    return axHisty,axHistx
    


def plot2dMultiNorm(x,y,xlim,ylim,nxbins=80,nybins=80,xlabel="log(Period)",ylabel="log(MES)"):
    """
    Create a 2dhistogram with other histograms overlaid on the axes.
    The zeroth, element in the x and y list will be shown in the 2dpart.   
    """
    color=["red","green","magenta","black","orange"]
    
    (axHisty,axHistx,ybins,xbins)=plot2dHistNorm(x[0],y[0],xlim,ylim,nxbins=nxbins,nybins=nybins,xlabel=xlabel,ylabel=ylabel)
    
    for i in np.arange(1,len(x)):
        
        print i
        
        axHisty.hist(y[i],bins=ybins,orientation='horizontal', color = color[i-1],histtype='step',normed=True, lw=1)    
    
        axHistx.hist(x[i],bins=xbins,color=color[i-1],histtype='step',normed=True, lw=1)  
    
    return axHisty,axHistx
    
def plot2dDiff(x,y,xlim,ylim,nxbins=80,nybins=80,xlabel="log(Period)",ylabel="log(MES)"):
    """
    Take histogram of both, but plot difference of histograms in the middle.
    """
    
def plot2dHistPaper(x,y,pcs,xlim,ylim,scorevalues,scorelimit=0.75,nxbins=50,nybins=50,xlabel="",ylabel="",\
                    makelog=False,clim=None,midtype='scatter',msize=10,colormap='plasma'):
    """
    Create a complicated 2d histogram plot with projected histograms onthe side
    Give the option for either a 2d histogram in the middle or a scatter - colored by score.
    msize is the marker size to use
    Show marginalized distribution for with and without a score cut.
    """
    gs=gridspec.GridSpec(3,3)    
    plt.figure(1, figsize=(9,9))
    
    
    axTemperature = plt.subplot(gs[1:3,0:2])  #plt.axes(rect_temperature) # temperature plot
    axHistx = plt.subplot(gs[0,0:2],sharex=axTemperature) #plt.axes(rect_histx) # x histogram ,
    axHisty = plt.subplot(gs[1:3,2], sharey=axTemperature) #plt.axes(rect_histy) # y histogram
    if makelog:
        axHistx.set_xscale=('log')
        axHisty.set_xscale=('log')
        axTemperature.set_xscale('log')
        axTemperature.set_yscale('log')
   
    spcs=scorevalues>scorelimit    
    
    
    xbins = np.linspace(start = xlim[0], stop = xlim[1], num = nxbins)
    ybins = np.linspace(start = ylim[0], stop = ylim[1], num = nybins)
    

    #Plot the histograms
    axHistx.hist(x[pcs], bins=xbins, color = 'blue',histtype='step',lw=2,label='PCs')
    axHistx.hist(x[spcs], bins=xbins, color = 'red',histtype='step',lw=1.5,label='Score > %4.2f'%scorelimit)
    axHistx.legend(loc=(1.1,0),fontsize=12)
    
    axHisty.hist(y[pcs], bins=ybins, orientation='horizontal', color = 'blue',histtype='step',lw=2)
    axHisty.hist(y[spcs], bins=ybins, orientation='horizontal', color = 'red',histtype='step',lw=1.5)    
    
    #Set up the histogram bins   
    nbins=(nxbins+nybins)/4
    print nxbins,nybins,nbins
    x2bins = np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/nbins)
    y2bins = np.arange(ylim[0], ylim[1], (ylim[1]-ylim[0])/nbins)
    
    #xcenter = (xbins[0:-1]+xbins[1:])/2.0
    #ycenter = (ybins[0:-1]+ybins[1:])/2.0
    #aspectratio = 1.0*(xlim[1] - xlim[0])/(1.0*(ylim[1] - ylim[0])
    
    if midtype=='hist':
        H, xedges,yedges = np.histogram2d(y,x,bins=(y2bins,x2bins))
        if clim==None:
            clim=[0,np.max(H)]
        
        # Plot the temperature data
        
        axTemperature.imshow(H, extent=[xlim[0],xlim[1],ylim[0],ylim[1]],\
                                 interpolation='nearest', origin='lower',
                                 aspect="auto", cmap=colormap, clim=clim)

    elif midtype=='scatter':
         cax=axTemperature.scatter(x[pcs],y[pcs],s=msize,c=scorevalues[pcs],edgecolors='k',\
                                 linewidth=0.2,cmap=colormap,marker='o',label='DR25 PCs')
         cbar=plt.colorbar(cax)
         plt.clim=clim
         cbar.set_label('Score',fontsize=14)
         #plt.legend(loc='best')
    
    
    lab=np.array([0.4,1.0,2.0,4.0,10.0,40])
    locat=np.log10(lab)
    plt.yticks(locat,['{:.0f}'.format(l) for l in lab])
    lab=np.array([0.4,1,4,10,40,100,400])
    locat=np.log10(lab)
    
    #plt.xticks(locat,lab.astype(str)) 
    axTemperature.set_xticks(locat)
    axTemperature.set_xticklabels(['{:.0f}'.format(l) for l in lab])
    #pdb.set_trace()

    #Set up the plot limits
    axTemperature.set_xlim(xlim)
    axTemperature.set_ylim(ylim)  
    
    #Plot the axes labels
    axTemperature.set_xlabel(xlabel,fontsize=16)
    axTemperature.set_ylabel(ylabel,fontsize=16)
    #plt.setp(axTemperature.get_xticklabels(), fontsize=14)
           
    #Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))
    

    
    return axHisty,axHistx,axTemperature



def plotConfirmedPCGrid(rvdata,cdata,feddata):
    """
    Plot fraction of confirmed planets that we continue to make into PCs
    Takes a robovetter dataframe with KOIs, mes and period
    Need the confirmed csv file, KeplerNames = cdata
    Needs the federation between the TCEs and the KOI = feddata
    """
    #Merge using the KOIs. Leave behind only the KOIs and data for all.
    #Start by merging cdata and feddata.
    
    clist=cdata.index
    isConfirmed=np.zeros(len(feddata))
    for i,koi in enumerate(feddata['koi']):
        
        if koi in clist:
            isConfirmed[i]=1

    feddata['confirmed']=isConfirmed
    alldata=rvdata.merge(feddata,left_index=True,right_index=True,how='left')
    
    isconf=alldata['confirmed']==1
    
    data=alldata[isconf]
    
    title='Confirmed\nDR25 RV PC'

    period=data['period']
    mes=data['mes']
    passed=data['disp']=='PC'
    
    plotGrid(period,mes,passed,(50,100))
    plt.title(title)
    
def plotFpwgPCGrid(rvdata,fpdata,feddata):
    """
    Plot fraction of FPWG planets that are passed by robovetter.
    """
    
    fplist=fpdata.index
    isFP=np.zeros(len(feddata))
    for i,koi in enumerate(feddata['koi']):
        if koi in fplist:
            isFP[i]=1

    feddata['fp']=isFP
    newdata=rvdata.merge(feddata,left_index=True,right_index=True,how='left')
    
    isfp=newdata['fp']==1
    
    data=newdata[isfp]

    title='CFP\n DR25 RV PC'
    period=data['period']
    mes=data['mes']
    passed=data['disp']=='PC'
    plotGrid(period,mes,passed,(0,50))
    plt.title(title)


def plotFullPageGrid(opsdata,invdata,injdata,climrange=(0,100)):
    """
    Create a full page diagnostic plot of the simple grid.
    """
    
    fig=plt.figure(figsize=(8.5,11))
    gs=gridspec.GridSpec(7,2)
    plt.subplots_adjust(hspace=0.75, wspace=0.75)
    
        
    #Create Plot
    ax=plt.subplot(gs[0:2,0])
    plotMPGrid(opsdata,'All OPS PC Rate',climrange)
    
    #Create for Cut Dwarf Plot
    fgkRvWant=(opsdata['logg']>4.0) & (opsdata['tstar']>=4000.0) & (opsdata['tstar']<7000.0);
    fgkData=opsdata[fgkRvWant]
    ax=plt.subplot(gs[0:2,1])
    plotMPGrid(fgkData,'FGK Dwarf OPS PC Rate',climrange)
    
    #(hzpercent[0,1],hztotal[0,1],hzpcs[0,1],hzopsfgkwant)=getHZpercent(fgkData)
    
    #INV
    #(hzpercent[1,0],hztotal[1,0],hzpcs[1,0],hzinvwant)=getHZpercent(invdata,Srange=Srange,pRadius=pRadius)
    
    #Create Plot
    ax=plt.subplot(gs[2:4,0])
    plotMPGrid(invdata,'All INV PC Rate',(climrange[0],climrange[1]/2))
    
    #Create for FGK Dwarf Plot
    fgkRvWant=(invdata['logg']>4.0) & (invdata['tstar']>=4000.0) & (invdata['tstar']<7000.0);
    fgkData=invdata[fgkRvWant]
    #(hzpercent[1,1],hztotal[1,1],hzpcs[1,1],hzinvfgkwant)=getHZpercent(fgkData)
    
    ax=plt.subplot(gs[2:4,1])
    plotMPGrid(fgkData,'FGK Dwarf INV PC Rate',(climrange[0],climrange[1]/2))
    ax.plot([-2,1.1],[1.15,1.15],transform=ax.transAxes, clip_on=False,color='k')
     
    #INJ
    #(hzpercent[2,0],hztotal[2,0],hzpcs[2,0],hzinjwant)=getHZpercent(injdata,Srange=Srange,pRadius=pRadius)
    
    #Create Plot
    plt.subplot(gs[4:6,0])
    plotMPGrid(injdata,'All INJ PC Rate',(climrange[1]/2,climrange[1]))
    
    #Create for FGK Dwarf Plot
    fgkRvWant=(injdata['logg']>4.0) & (injdata['tstar']>=4000.0) & (injdata['tstar']<7000.0);
    fgkData=injdata[fgkRvWant]
    #(hzpercent[2,1],hztotal[2,1],hzpcs[2,1],hzinjwant)=getHZpercent(fgkData)
    
    ax=plt.subplot(gs[4:6,1])
    plotMPGrid(fgkData,'FGK Dwarf INJ PC Rate',(climrange[1]/2,climrange[1]))
    ax.plot([-2,1.1],[1.15,1.15],transform=ax.transAxes, clip_on=False,color='k')   
    
    return fig,ax
    
      
    
def plotMPGrid(data,title,climRange=(0,100)):
    """
    Plot the prv.PlotGrid of MES vs Period with given title and range.
    """
    period=data['period']
    mes=data['mes']
    passed=data['disp']=='PC'
    plotGrid(period,mes,passed,climRange)
    plt.title(title)
  
def changeLogXLabel(lab):
    """
    Plotting commands to create new period with linear numbers
    lab=np.array([.3,1,3,10])
    """
    #Period
    locat=np.log10(lab)
    plt.xticks(locat,lab.astype(str))
    
def changeLogYLabel(lab):
    """
    Plotting commands to create new period with linear numbers
    lab=np.array([.3,1,3,10])
    """
    #MES
    locat=np.log10(lab)
    plt.yticks(locat,lab.astype(str))

    
    
def plotNewKoisPcs(opsData):
    """Plot the new KOIs and the new PCs as a function of period and mes
    rvfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/r61794/RoboVetterOut-OPS.txt'
    tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
    fedfile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
    cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'

    opsData=io.createAllResults(rvfile,tcefile,fedfile,cumfile)
    """
    xlabel=np.array([0.3,1,3,10,30,100,372])
    ylabel=np.array([7,10,15,20,30,40])
    
    opsData['iskoi']=opsData['koi'].notnull()
    wantnewkoi=(~opsData['iskoi']) & (opsData['N']==0)
    wantnewpcs=((~opsData['iskoi']) | (opsData['dr24disp']=='FP')) & (opsData['disp']=='PC') 
    missedKOIs=(opsData['iskoi']) & (opsData['N']==1)
    missedPCs=((opsData['iskoi']) & (opsData['dr24disp']=='PC')) & (opsData['disp']=='FP')

    numNewKoi=len(wantnewkoi[wantnewkoi])
    numNewPCs=len(wantnewpcs[wantnewpcs])
    numMissedKOIs=len(missedKOIs[missedKOIs])   
    numMissedPCs=len(missedPCs[missedPCs])
    
    plt.subplots_adjust(hspace=0.40, wspace=0.25)
    plt.subplot(221)
    plt.plot(np.log10(opsData['period'][wantnewkoi]),np.log10(opsData['mes'][wantnewkoi]),'r.')
    plt.title("New KOIs: %u" % (numNewKoi),fontsize=11)
    plt.xlabel('Period (d)')
    plt.ylabel('MES')   
    plt.ylim((.8,1.6))
    changeLogXLabel(xlabel)
    changeLogYLabel(ylabel)
    
    plt.subplot(222)
    plt.plot(np.log10(opsData['period'][wantnewpcs]),np.log10(opsData['mes'][wantnewpcs]),'b.')
    plt.title("New PCs (on new KOI or FP in DR24): %u" % (numNewPCs),fontsize=11)
    plt.xlabel('Period (d)')
    plt.ylabel('MES') 
    plt.ylim((.8,1.6))
    changeLogXLabel(xlabel)
    changeLogYLabel(ylabel)
    
    plt.subplot(224)
    plt.plot(0,0,'r*',ms=14)
    plt.plot(np.log10(opsData['srad'][wantnewpcs]),np.log10(opsData['rplanet'][wantnewpcs]),'b.')
    plt.title("New PCs (on new KOI or FP in DR24)",fontsize=11)
    plt.xlabel('log Insolation Flux')
    plt.ylabel('log Planet Radius')
    
    
    plt.subplot(223)
    plt.plot(np.log10(opsData['period'][missedKOIs]),np.log10(opsData['mes'][missedKOIs]),'g.',ms=6)
    plt.plot(np.log10(opsData['period'][missedPCs]),np.log10(opsData['mes'][missedPCs]),'k.',ms=4)    
    
    plt.ylim((.5,5.1))
    plt.title('DR25 FPs that were PC in DR24 (black,%u)\n Missed KOIs (green, %u)' % (numMissedPCs,numMissedKOIs),fontsize=11)
    plt.xlabel('Period (d)')
    plt.ylabel('MES') 
    changeLogXLabel(xlabel)
    changeLogYLabel(np.array([10,30,100,300,1000,10000]))
    
    
def returnBins(df,columns=["period","mes","duration","depth","pn"],logged=[True,True,False,True,False],num=100):
    """
    Bins for plotMultiHist
    """    
    allbins=list(list(columns))
    L=len(columns)

    for i,c in enumerate(columns):
        
        plt.subplot(L,1,i+1)
        if (columns != 'pn'):
            if logged[i]:
                xmin=np.log10(np.percentile(df[c],2.0))
                xmax=np.log10(np.percentile(df[c],98.0))
                bins=np.linspace(xmin,xmax,num=num)
                allbins[i]=bins
            else:
                xmin=np.percentile(df[c],1.0)
                xmax=np.percentile(df[c],99.0)
                bins=np.linspace(xmin,xmax,num=num)
                allbins[i]=bins
        else:
            allbins[i]=[1,2,3,4,5,6,7,8,9,10]
        
        allbins[i]=(bins)
    
    return allbins
    
    
def plotMultiHist(df,allbins,columns=["period","mes","duration","pn"],logged=[True,True,False,False]):
    """
    for a list of columns, plot a multipanel plot of the distributions
    """
    L=len(columns)
    for i,c in enumerate(columns):
        print c
        
        plt.subplot(L,1,i+1)
        if logged[i]:
            plt.hist(np.log10(df[c]),bins=allbins[i],histtype='step')
        else:
            plt.hist(df[c],bins=allbins[i],histtype='step')

        plt.xlabel(c)
    
    return allbins
    
def plotThresholdCheck(df,mFail, title, metric,thresh,params=['period','mes','score'],logScale=True):
    """
    metricFail is the array of T/F where T indicates it has failed that metric
    title should contain the name of the comment field that is considered a metricFail
    metrics contains the name of the metric and thresh is the threshold for it.
    params should contain three strings with the names of the parameters to plot against each other
    The first two will also be plotted against each other.
    """

    #Need a string of wants that fail for the specified reason

    p2=df[params[2]]
    if logScale:
        p0=np.log10(df[params[0]])
        p1=np.log10(df[params[1]])
        m=np.log10(df[metric])
        keep=np.array(np.isfinite(m))
        #print type(keep)        
        p0=p0[keep]
        p1=p1[keep]
        p2=p2[keep]
        m=m[keep]
        metricFail=mFail[keep]
        th=np.log10(thresh)
        #print "mFail: %u, keep: %u, p2: %u" % (len(mFail),len(keep[keep]),len(p2))
        df=df[keep]

    else:
        p0=df[params[0]]
        p1=df[params[1]]
        m=df[metric] 
        th=thresh
        metricFail=mFail
        
    #print "metricFail-len one %u" % len(metricFail)
    pcs=df.disp=='PC'
    plt.figure(figsize=(8.5,11))
    plt.subplot(411)
    x=p0
    y=p1

    #plt.plot(x[~pcs & ~metricFail],y[~pcs & ~metricFail],'r.',label='FP mPass')
    #plt.plot(x[pcs & ~metricFail],y[pcs & ~metricFail],'r*',label='PC mPass')
    plt.plot(x[~pcs & metricFail],y[~pcs & metricFail],'k.',ms=3,label='FP mFail')
    plt.plot(x[pcs & metricFail],y[pcs & metricFail],'rv',ms=4,label='PC mFail')
    
    plt.xlabel(params[0])
    plt.ylabel(params[1])
    plt.legend(fontsize='x-small')
    plt.title(title)
    plt.ylim(np.percentile(y,.1),np.percentile(y,99.9))
    plt.xlim(np.percentile(x,.01),np.percentile(x,99.9))
    
    plt.subplot(412)
    x=p0
    y=m
    plt.plot(x[~pcs & ~metricFail],y[~pcs & ~metricFail],'r.',ms=3,label='FP mPass')
    plt.plot(x[pcs & ~metricFail],y[pcs & ~metricFail],'m.',ms=3,label='PC mPass')
    plt.plot(x[~pcs & metricFail],y[~pcs & metricFail],'k.',ms=3,label='FP mFail')
    plt.plot(x[pcs & metricFail],y[pcs & metricFail],'kv',ms=4,label='PC mFail')
    plt.hlines(th,np.min(x)-.1,np.max(x)+.1,colors='blue',linestyles='dashed')

    plt.xlabel(params[0])
    plt.ylabel(metric)
    
    plt.ylim(np.percentile(y,.1),np.percentile(y,99.9))
    
    
    plt.subplot(413)
    x=p1
    y=m
    plt.plot(x[~pcs & ~metricFail],y[~pcs & ~metricFail],'r.',ms=3,label='FP mPass')
    plt.plot(x[pcs & ~metricFail],y[pcs & ~metricFail],'m.',ms=3,label='PC mPass')
    plt.plot(x[~pcs & metricFail],y[~pcs & metricFail],'k.',ms=3,label='FP mFail')
    plt.plot(x[pcs & metricFail],y[pcs & metricFail],'kv',ms=4,label='PC mFail')
    plt.hlines(th,np.min(x),np.max(x),colors='blue',linestyles='dashed')
    
    plt.xlabel(params[1])
    plt.ylabel(metric)
    plt.legend(fontsize='x-small')
    plt.ylim(np.percentile(y,.1),np.percentile(y,99.9))
    
    plt.subplot(414)
    x=p2
    y=m
    
    plt.hist(x[metricFail],bins=200,color='black',log=True,histtype='step',label='mFail')  
    plt.hist(x[~metricFail],bins=200,color='red',log=True,histtype='step',label='mPass')
    plt.xlim(-.02,1.02)
    plt.xlabel(params[2])
    plt.legend()
    
    #plt.plot(x[~pcs & metricFail],y[~pcs & metricFail],'k.',label='FP mFail')
    #plt.plot(x[pcs & metricFail],y[pcs & metricFail],'kv',label='PC mFail')
    #plt.plot(x[~pcs & ~metricFail],y[~pcs & ~metricFail],'r.',label='FP mPass')
    #plt.plot(x[pcs & ~metricFail],y[pcs & ~metricFail],'rv',label='PC mPass')
    
    #plt.xlabel(params[2])
    #plt.ylabel(metric)

def get1DReliability(ops,inv,metric,bins,s=[0,1]):
    """
    Return a 1D reliability array given a ops and inversion set.
    """
    bins=np.array(bins)
    opspcs=passes(ops,s=s)
    opsfps=~opspcs
    invpcs=passes(inv,s=s)
    invfps=~invpcs
    
    nopspcs,bb,patch=plt.hist(ops[metric][opspcs],bins=bins)
    nopsfps,bb,patch=plt.hist(ops[metric][opsfps],bins=bins)
    ninvpcs,bb,patch=plt.hist(inv[metric][invpcs],bins=bins)
    ninvfps,bb,patch=plt.hist(inv[metric][invfps],bins=bins)
    
    eff=ninvfps.astype(float)/(ninvfps+ninvpcs).astype(float)
    
    rel=rvs.arrayReliability(nopsfps.astype(float),nopspcs.astype(float),eff)
    
    return rel,eff

def plot1DReliability(ops,inv,metric,bins,s=[0.0,1.0],xlabel='metric'):
    """
    Given an ops and inv data frame
    plot the reliability as a function of that metric for the bins specified.
    s is a score cut for consideration of the 
    """
    bins=np.array(bins)
    opspcs=passes(ops,s=s)
    opsfps=~opspcs
    invpcs=passes(inv,s=s)
    invfps=~invpcs
    
    nopspcs,bb,patch=plt.hist(ops[metric][opspcs],bins=bins)
    nopsfps,bb,patch=plt.hist(ops[metric][opsfps],bins=bins)
    ninvpcs,bb,patch=plt.hist(inv[metric][invpcs],bins=bins)
    ninvfps,bb,patch=plt.hist(inv[metric][invfps],bins=bins)
    
    eff=ninvfps.astype(float)/(ninvfps+ninvpcs).astype(float)
    plt.clf()

    midbins=(bins[:-1]+bins[1:])/2.0    
    
    rel=rvs.arrayReliability(nopsfps.astype(float),nopspcs.astype(float),eff)
    print rel
    #plt.figure()
    plt.plot(midbins,rel,'-ko',lw=2.5,ms=5)
    if xlabel=='metric':
        plt.xlabel(metric)
    else:
        plt.xlabel(xlabel)
    plt.ylabel('Reliability')
    #plt.show()

def plot1DReliabilityGroups(ops,inv,metric,bins,opsgroup,invgroup,s=[0.0,1.0],xlabel='metric',labels=['a','b','c','d']):
    """
    Given an ops and inv data frame
    plot the reliability as a function of that metric for the bins specified.
    s is a score cut for consideration of the 
    """
    bins=np.array(bins)
    opspcs=passes(ops,s=s)
    opsfps=~opspcs
    invpcs=passes(inv,s=s)
    invfps=~invpcs
    
    for i in np.arange(0,len(opsgroup),1):
   
        nopspcs,bb=np.histogram(ops[metric][opspcs & opsgroup[i]],bins=bins)
        nopsfps,bb=np.histogram(ops[metric][opsfps & opsgroup[i]],bins=bins)
        ninvpcs,bb=np.histogram(inv[metric][invpcs & invgroup[i]],bins=bins)
        ninvfps,bb=np.histogram(inv[metric][invfps & invgroup[i]],bins=bins)
        
        eff=ninvfps.astype(float)/(ninvfps+ninvpcs).astype(float)
    
        #midbins=(bins[:-1]+bins[1:])/2.0    
        
        rel=rvs.arrayReliability(nopsfps.astype(float),nopspcs.astype(float),eff)

        plt.step(bb[:-1],rel,'-',lw=2.1-i*0.1,ms=5,label=labels[i],where='post')
    

    if xlabel=='metric':
        plt.xlabel(metric)
    else:
        plt.xlabel(xlabel)
    plt.ylabel('Reliability')
    plt.legend(loc='best')

def plot1DCompletenessGroups(inj,metric,bins,injgroup,s=[0.0,1.0],xlabel='metric',labels=['a','b','c','d']):
    """
    Plot the completeness against the metric binned by bins. 
    injgroups contains n arrays of t/F of ways you want to group against any other metric.
    labels should give the text to go into a legend for the plot
    """
    bins=np.array(bins)
    injpcs=passes(inj,s=s)
    #injfps=~injpcs
    
    for i in np.arange(0,len(injgroup),1):
        
        ninj,bb=np.histogram(inj[metric][injgroup[i]],bins=bins)
        npc,bb=np.histogram(inj[metric][injgroup[i] & injpcs], bins=bins)
        C=npc.astype(float)/ninj.astype(float)        
        
        plt.step(bb[:-1],C,lw=2.1-i*0.1,label=labels[i],where='post')
    
    if xlabel=='metric':
        plt.xlabel(metric)
    else:
        plt.xlabel(xlabel)
    
    plt.ylabel('Completeness')
    plt.legend(loc='best')
    
def plotFracFpPerMetric(ops,inv,ss1, metricflag,reltype=("INV","SS1"),\
        xBins = [0, 10, 200, 500],yBins = [7, 10, 20, 2000],\
        drange=(0,100),s=[0.0,1.0]):
    """
    using plotGrid
    plot the percent of FPs in OPS and INV data frames that fail that metric.
    metricflag is the name to match in the FP flags.
    """
    
    plt.figure(figsize=(14,4))
    
    plt.subplot(132)
    plt.title('OPS %s' % metricflag)
    
    fps=~passes(ops,s=s)
    per=np.array(ops.loc[fps]['period'])
    mes=np.array(ops.loc[fps]['mes'])
    pcIdx=np.array(map(lambda x:x.find(metricflag)>=0,ops.loc[fps]['flags']),dtype=bool)
    plotGrid(per, mes, pcIdx, drange=drange,xBins = xBins,yBins = yBins,clabel="Percent FPs Failed by Metric")
    
    
    plt.subplot(131)
    plt.title('%s %s' % (reltype[0],metricflag))
    
    fps=~passes(inv,s=s)
    peri=np.array(inv.loc[fps]['period'])
    mesi=np.array(inv.loc[fps]['mes'])
    pcIdxi=np.array(map(lambda x:x.find(metricflag)>=0,inv.loc[fps]['flags']),dtype=bool)
    plotGrid(peri, mesi, pcIdxi, drange=drange,xBins = xBins,yBins = yBins,clabel="Percent FPs Failed by Metric")
    
    
    plt.subplot(133)
    plt.title('%s %s' % (reltype[1],metricflag))
    
    fps=~passes(ss1,s=s)
    peri=np.array(ss1.loc[fps]['period'])
    mesi=np.array(ss1.loc[fps]['mes'])
    pcIdxi=np.array(map(lambda x:x.find(metricflag)>=0,ss1.loc[fps]['flags']),dtype=bool)
    plotGrid(peri, mesi, pcIdxi, drange=drange,xBins = xBins,yBins = yBins,clabel="Percent FPs Failed by Metric")


def inBox(data,metric1,range1,metric2,range2):
    """
    get true false array of those things in the box defined by metric1 and metic2
    and range1 and range2
    """
    want=(data[metric1]>=range1[0]) & (data[metric1]<range1[1]) & \
        (data[metric2]>=range2[0]) & (data[metric2]<range2[1])
    
    return want


import createRVPerfMetricsPlots as pmp
def plotRelCompScore(ops,inv,inj,metric1,metric2,range1,range2,scores=np.arange(0,1,.05),\
    fpscore=1.0,Rlim=(.2,1.02),Clim=(1.0,.2)):
    """
    ops data frame
    inv dataframe
    inj dataframe with metric1, metric2, disposition and score.
    a false alarm population data frame
    metric1 with range1 as a min, max to cut on
    metric2 with range2 as a min, max to cut on
    scores is the scores to calculate it across.
    fpscore is the score to let fps back in.
    """
    
    opsbox=inBox(ops,metric1,range1,metric2,range2)
    invbox=inBox(inv,metric1,range1,metric2,range2)
    injbox=inBox(inj,metric1,range1,metric2,range2)
    
    C=np.zeros([len(scores)])
    R=np.zeros([len(scores)])
    E=np.zeros([len(scores)])

    for (i,scorev) in enumerate(scores):
        
        #opspc=passes(ops,s=(scorev,fpscore))
        #injpc=passes(inj,s=(scorev,fpscore))  - use this if you want a fixed fpscore
        injpc=passes(inj,s=(scorev,scorev))
        
        C[i]=np.float(len(inj[injbox & injpc])) / np.float(len(inj[injbox]))
        
        #R[i],E[i]=pmp.estimateReliability(inv[invbox],ops[opsbox],s=(scorev,fpscore)) - use this if you want a fixed fpscore
        R[i],E[i]=pmp.estimateReliability(inv[invbox],ops[opsbox],s=(scorev,scorev))
    
    fig1=plt.figure(figsize=(5,5.4),dpi=150)
    ax1=fig1.add_subplot(111)
    line1=ax1.plot(C,R,'ro-',label='Reliability')
    
    plt.ylim(Rlim)
    ax1.invert_xaxis()
    
    plt.ylabel('Reliability',fontsize=16)
    plt.yticks(fontsize=14)
    plt.xlabel('Completeness',fontsize=16)
    #atitle='Reliability vs. Completeness adjusting Score\n box %s fpscore %3.1f ' % (str(box),fpScoreLimit)
    #plt.title(atitle)
    for i in np.arange(0,len(scores),2):
        st=scores[i].astype('str')
        plt.annotate(st,xy=(C[i]-.005,R[i]-.032),fontsize=16)
    
    ax2=fig1.add_subplot(111,sharex=ax1, frameon=False)
    line2=ax2.plot(C,E,'bs--',ms=4, lw=2,label='Effectiveness')
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    plt.ylabel("Effectiveness",fontsize=16)
    plt.yticks(fontsize=12)
    plt.xlim(Clim)
    plt.ylim([.982,1.0009])
    import matplotlib.lines as mlines
    red_line = mlines.Line2D([],[],color='red', label='Reliability',lw=2)
    blue_line = mlines.Line2D([],[],color='blue',label='Effectiveness',lw=2)
    #bline=plt.plot([-1],[-1],'bx-')
    plt.legend(handles=[red_line,blue_line],loc="lower right",fontsize=14)
    #plt.legend([line1, line2], ["Reliability","Effectiveness"],loc="best",fontsize=12,framealpha=0.75,numpoints=1)
    
    print C,R
    
    
def plotHistResults(ops,nb=150):
    """
    Plot used for koi documentation of the TCEs and the N==1 and PCs.
    """
    logperiods=np.log10(ops.period)    
    tlike=(ops.isdr25koi) #(ops.N==0) & (ops.E==0)
    pclike=tlike  & (ops.disp=='PC')        
    
    plt.figure(figsize=(10,7))
    n,bins,p=plt.hist(logperiods,bins=nb,color='firebrick',label='TCEs',linewidth=.5) 
    plt.hist(logperiods[tlike],bins=bins,color='goldenrod',label='KOIs',linewidth=.5)
    plt.hist(logperiods[pclike],bins=bins,color='powderblue', label='Planet Candidates',linewidth=.4)
    plt.title('Kepler Q1--Q17 DR25')
    plt.legend(loc='upper left')
    
    changeLogXLabel(np.array([0.5,1.0,3,10,30,100,500]))
    plt.xlabel('Period (days)')
    plt.ylabel('Number')


def plotChanges(both,xname,yname,want,labels):
    #xnames, ynames and labels are arrays
    #plots x[0] vs y[0] and x[1] vs y[1] etc.
    #Then connects a line between the first set and the second set.

    shapes=('o','s','^','v','*')

    plt.figure()
    
    for i,v in enumerate(labels):
        plt.plot(both[want][xname[i]],both[want][yname[i]],shapes[i],ms=5,label=labels[i],lw=0)
    
    plt.xlabel(xname[0])
    plt.ylabel(yname[0])
    plt.legend(loc='lower left')
    
    
    for i,v in enumerate(both.index[want]):
        x=(both.loc[v][[xname[0]]].astype(float),both.loc[v][[xname[1]]].astype(float))
        y=(both.loc[v][[yname[0]]].astype(float),both.loc[v][[yname[1]]].astype(float))
        
        plt.plot(x,y,'-k')