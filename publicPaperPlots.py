# -*- coding: utf-8 -*-
"""
Created on Fri May 12 11:03:22 2017

Paper Plots from Public Data

@author: sthompson
@Date: May 12, 2017
"""

import matplotlib.pyplot as plt
import numpy as np
import hzCalc as hzCalc
import createRVPerfMetricsPlots as pmp
import plotRobovetter as prv
import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter, MaxNLocator   #You need this below
import pdb
    #pdb.set_trace()

def passes(data,s=[0.0,1.0]):
    """
    Use disposition and scores to determine which ones pass.
    Returns true false array
    """
    
    passed= ((data['disp'] == 'PC') & (data['score']>=s[0])) | \
            ((data['disp']=='FP') & (data['score']>=s[1]))  
    return passed

def inBox(data,metric1,range1,metric2,range2):
    """
    get true false array of those things in the box defined by metric1 and metic2
    and range1 and range2
    """
    want=(data[metric1]>=range1[0]) & (data[metric1]<range1[1]) & \
        (data[metric2]>=range2[0]) & (data[metric2]<range2[1])
    
    return want

def fracPass(inj,s=[0.0,1.0]):
    """Return the fractional completeness given the injected data frame
    """
    passed=passes(inj,s)
    C=len(inj[passed])/np.float(len(inj))
    
    return C
    
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
            plt.text(c+.5, r+.5, text, va="center", ha="center", color=cType,fontsize=11)

def writeLabels(percent, drange):
    nR, nC = percent.shape

    midVal = .5*(np.max(drange) + np.min(drange))
    for r in range(nR):
        for c in range(nC):
            if np.isnan(percent[r,c]):
                text = ""
            else:
                text = "%.0f%%" %(percent[r,c])

            cType = 'k'
            if percent[r, c] < midVal:
                cType='w'
            plt.text(c+.5, r+.5, text, va="center", ha="center", color=cType,fontsize=11,fontweight='bold')

def fgkwant(data):
    want=(data['logg']>=4.0) & (data['Ts']>=4000.0) & (data['Ts']<7000.0);
    return want

def compute2dHist(x, y, xBins, yBins):
    h, xedege, yedge = np.histogram2d(x,y, [xBins, yBins])
    h = h.transpose()
    return h 

def HZplot(koi, scut=0.5,pradcut=1.8, mfactor=15, addConf=True, \
            ylimits=(3200,6400), xlimits=(2.5,0.0), ebars=True):

    """
    Input koi dataframe
    Input score cut
    input scale factor for radii.
    addConf = True means add back in the confirmed planets we lots.
    
    Plots the habitable zone objects considering error bars of 1 sigma.
    as Insolation flux vs. Stellar temperature 
    It is colored by Score.
    
    """
    pcs=(koi.koi_score >= scut) & (koi.koi_pdisposition=='CANDIDATE')
    smallerrors=(-1*koi.koi_prad_err2/koi.koi_prad < 2.0 ) & (-1*koi.koi_insol_err2/koi.koi_insol<2.0)
    want=koi.koi_prad+koi.koi_prad_err2 <= pradcut
    wantlow=np.array(map( lambda x,y:hzCalc.inHZ(x,y), koi.koi_steff, koi.koi_insol+koi.koi_insol_err2))
    wantup=np.array(map( lambda x,y:hzCalc.inHZ(x,y), koi.koi_steff, koi.koi_insol+koi.koi_insol_err1))
    wantin=np.array(map( lambda x,y:hzCalc.inHZ(x,y), koi.koi_steff, koi.koi_insol))
    
    rearth=1
    rsupearth=2
    mfactor=15
    pars=['kepoi_name','kepler_name','koi_disposition','koi_score','koi_period',\
          'koi_prad','koi_prad_err2','koi_prad_err1','koi_steff','koi_steff_err1',\
          'koi_steff_err2','koi_insol','koi_insol_err1','koi_insol_err2','koi_max_mult_ev','koi_slogg']
    hzframe=koi[pcs & want & (wantlow | wantup | wantin) &smallerrors][pars]
    
    #Add in those three Confirmed ones we are missing.
    if addConf:
        missing=["009002278-04","008845205-01", "8142787-01"]#,"010604335-02","005640085-02"]  #[62f,283c,159c]
        for tce in missing:
            mydf=koi[koi.index==tce][pars]
            hzframe=hzframe.append(mydf)
    
    hzframe['RpBig']=hzframe.koi_prad+hzframe.koi_prad_err1
    hzframe['RpSmall']=hzframe.koi_prad+hzframe.koi_prad_err2
    
    print len(hzframe)

    plt.figure(figsize=(8.5,11))
    
    for i,hz in enumerate(hzframe.index):
        Srad=hzframe.loc[hz]['koi_insol']
        Smin=Srad+hzframe.loc[hz]['koi_insol_err2']
        Smax=Srad+hzframe.loc[hz]['koi_insol_err1']
        rp=hzframe.loc[hz]['koi_prad']
        rpBig=hzframe.loc[hz]['RpBig']
        rpSmall=hzframe.loc[hz]['RpSmall']
        tstar=hzframe.loc[hz]['koi_steff']
        tmin=tstar+hzframe.loc[hz]['koi_steff_err2']
        tmax=tstar+hzframe.loc[hz]['koi_steff_err1']
        sc=hzframe.loc[hz]['koi_score']
        kepname="%s" % (str(hzframe.loc[hz].kepoi_name))

        if ebars:
            plt.hlines(tstar,Smin,Smax,colors='steelblue',label='.',lw=1.3,zorder=1)
            plt.vlines(Srad,tmin,tmax,colors='steelblue',label='.',lw=1.3,zorder=1)
        plt.scatter(Srad,tstar,s=(rpBig*mfactor)**2,marker='o',c=sc,cmap="PiYG",linewidth=0.2,vmin=-.3,vmax=1,zorder=2)
        plt.scatter(Srad,tstar,s=(rpSmall*mfactor)**2,marker='o',c="white",linewidth=0.2,zorder=3,alpha=0.95)
        #plt.scatter(Srad,tstar,s=(rp*mfactor)**2,marker='o',c=sc,cmap="brg_r",linewidth=.4,vmin=0,vmax=1,zorder=2)
        if hzframe.loc[hz].kepoi_name>7621.01:
            plt.annotate(s=kepname,xy=(Srad-.1,tstar+8),xycoords='data',color="grey",fontsize=10)
        #elif p.notnull(hzframe.loc[hz].kepler_name):
        #    plt.annotate(s=hzframe.loc[hz].kepler_name,xy=(Srad-.1,tstar+8),xycoords='data',color="grey",fontsize=10)
        
    #plt.scatter(hzframe.Srad,hzframe.tstar,s=(hzframe.RpBig*mfactor)**2,marker='o',c=hzframe.score,cmap="coolwarm_r",linewidth=0.1,vmin=0,vmax=1)
    #plt.scatter(hzframe.Srad,hzframe.tstar,s=(hzframe.RpSmall*mfactor)**2,marker='o',c='white',linewidth=0.1)
    
    small=plt.scatter(1,5777,s=(rearth*mfactor)**2,marker='o',edgecolor='blue',linewidth=1.5,label='R$_\oplus$',facecolors='none')
    big=plt.scatter(1.0,8300,s=(rsupearth*mfactor)**2,marker='o',c='grey',linewidth=1.4,label='2R$_\oplus$')
    half=plt.scatter(1.0,8500,s=(1.5*rearth*mfactor)**2,marker='o',edgecolor='blue',linewidth=1.5,label='1.5R$_\oplus$',facecolors='none')
    plt.xlabel('Insolation Flux ($S_{\oplus}$)', fontsize=12)
    plt.ylabel('Stellar Effective Temperature (K)', fontsize=12)
    
    plt.legend(handles=[small,half],scatterpoints=1,framealpha=0.8,loc='upper right',fontsize=14)
    
    #Draw a HZ across the plot.
    teff=np.linspace(2000,7000,50)
    hzBorders=hzCalc.hzBorder(teff)
    Shot=hzBorders[:,0]
    Scool=hzBorders[:,3]
    #hzhot=7000*sx-7500
    #hzcool=17500*sx-500
    plt.plot(Shot,teff,'--b')
    plt.plot(Scool,teff,'--r')
    #plt.title('Score cut of %f' % scut)
    
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    
    return hzframe
    #plt.xlim((2.5,0))
    #plt.ylim((2500,7000))
    #plt.savefig('%s/fig-hzTeffInsol2.png' %outdir)

#%%
#Plot the Score vs. the Npc*R*C
def plotNpcScore(ops,inv,inj,metric1,metric2,range1,range2,xlim=(0.05,1.0),scores=np.arange(0,1,.05),annotate=""):
    """
    For a particular box in parameter space (metric1,metric2,range1,range2), calculate the 
    Number of Real Pcs (Npc*R/C)
    Then plot that as a function of score cut for that parameter space.
    fpscore is a way of letting FPs back in
    """

    opsbox=inBox(ops,metric1,range1,metric2,range2)
    invbox=inBox(inv,metric1,range1,metric2,range2)
    injbox=inBox(inj,metric1,range1,metric2,range2)
    
    C=np.zeros([len(scores)])
    R=np.zeros([len(scores)])
    E=np.zeros([len(scores)])
    Nopspc=np.zeros([len(scores)])
    
    for (i,scorev) in enumerate(scores):
        opsPCs=passes(ops,s=[scorev,1.0])
        Nopspc[i]=len(opsPCs[opsPCs & opsbox])   
        #injpc=passes(inj,s=(scorev,fpscore))  - use this if you want a fixed fpscore
        injpc=passes(inj,s=(scorev,1.0)) 
        C[i]=np.float(len(inj[injbox & injpc])) / np.float(len(inj[injbox]))
        
        #R[i],E[i]=pmp.estimateReliability(inv[invbox],ops[opsbox],s=(scorev,fpscore)) - use this if you want a fixed fpscore
        R[i],E[i]=pmp.estimateReliability(inv[invbox],ops[opsbox],s=(scorev,1.0))
    
    Nreal=R*Nopspc/C
    
    plt.figure(figsize=(5,5),dpi=140)
    
    plt.plot(scores,Nreal,'bo-',label='Corrected Candidates',linewidth=2)
    plt.errorbar(scores,Nreal,yerr=np.sqrt(Nopspc)*R/C,color='b')
    plt.plot(scores,Nopspc,'ro--',label='Observed Candidates',linewidth=1.5)
    plt.xlabel('Disposition Score Threshold',fontsize=14)
    plt.ylabel('Number of Candidates ',fontsize=14)
    plt.legend(loc='upper center',fontsize=13)
    plt.annotate(annotate,xy=(0.15,0.14),xycoords='figure fraction',fontsize=14)
    plt.xlim(xlim)
#%%    
##-------------------------------------------------------
def plot2dHistPaper(x,y,pcs,xlim,ylim,scorevalues,scorelimit=0.7,nxbins=50,nybins=50,xlabel="",ylabel="",\
                    makelog=False,clim=None,midtype='scatter',msize=14,colormap='plasma_r',logticks=True):
    """
    Create a complicated 2d histogram plot with projected histograms onthe side
    Give the option for either a 2d histogram in the middle or a scatter - colored by score.
    msize is the marker size to use
    Show marginalized distribution for with and without a score cut.
    """
    gs=gridspec.GridSpec(3,3)    
    #plt.figure(1, figsize=(9,9))
    
    
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
    axHistx.hist(x[spcs], bins=xbins, color = 'red',histtype='step',lw=2.1,label='Score > %4.2f'%scorelimit)
    axHistx.hist(x[pcs], bins=xbins, color = 'blue',histtype='step',lw=1.3,label='PCs')
    axHistx.legend(loc=(1.1,0),fontsize=12)

    
    axHisty.hist(y[pcs], bins=ybins, orientation='horizontal', color = 'blue',histtype='step',lw=1.3,zorder=2)
    axHisty.hist(y[spcs], bins=ybins, orientation='horizontal', color = 'red',histtype='step',lw=2.1,zorder=1)    
    plt.xlabel('counts')
    plt.yticks(fontsize=9)
    
    #Set up the histogram bins   
    nbins=(nxbins+nybins)/4
    print nxbins,nybins,nbins
    x2bins = np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/nbins)
    y2bins = np.arange(ylim[0], ylim[1], (ylim[1]-ylim[0])/nbins)
    
    #xcenter = (xbins[0:-1]+xbins[1:])/2.0
    #ycenter = (ybins[0:-1]+ybins[1:])/2.0
    #aspectratio = 1.0*(xlim[1] - xlim[0])/(1.0*(ylim[1] - ylim[0])
    #pdb.set_trace()
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
                                 linewidth=0.03,cmap=colormap,marker='o',label='DR25 PCs',vmin=0,vmax=1.06)
         cbar=plt.colorbar(cax)
         #cbar.set_clim(0,1.0)
         cbar.set_label('Score',fontsize=14)
         cbar.set_ticks([0.0,0.2,0.4,0.6,0.8,1.0])
         #plt.legend(loc='best')
         
    if logticks:
        lab=np.array([0.4,1.0,2.0,3.0,10.0,30])
        locat=np.log10(lab)
        plt.yticks(locat,['{:.1f}'.format(l) for l in lab])
        lab=np.array([0.4,1,3,10,30,100,300])
        locat=np.log10(lab)
        
        #plt.xticks(locat,lab.astype(str)) 
        axTemperature.set_xticks(locat)
        axTemperature.set_xticklabels(['{:.1f}'.format(l) for l in lab])
        #pdb.set_trace()

    #Set up the plot limits
    axTemperature.set_xlim(xlim)
    axTemperature.set_ylim(ylim)  
    
    axTemperature.set_xlabel(xlabel,fontsize=16)
    axTemperature.set_ylabel(ylabel,fontsize=16)
    #plt.setp(axTemperature.get_xticklabels(), fontsize=14)
      
    #Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))
    #plt.show()
    
    return axHisty,axHistx,axTemperature

def plotReliability(invdata,opsdata,xmetric='period',ymetric='mes',xBins=[0, 10, 200, 500],yBins=[7, 10, 20, 2000],drange=(1,99),\
                        atitle="Reliability\n#PCs/#FPs",s=[0.0,1.0],xlabel='period',ylabel='mes'):
    """
    Plot grid of reliability numbers
    Similar to plotgrid above, but needs that to run twice to get the numbers 
    we need the dataframes for inversion and ops.
    """
    
    opsPeriod=opsdata[xmetric]
    opsMes=opsdata[ymetric]
    opsDisp=passes(opsdata,s) 

    allHistops = compute2dHist(opsPeriod, opsMes, xBins, yBins)
    pcHistops = compute2dHist(opsPeriod[opsDisp], opsMes[opsDisp], xBins, yBins)

    invPeriod=invdata[xmetric]
    invMes=invdata[ymetric]
    invDisp=passes(invdata,s) 

    allHistinv = compute2dHist(invPeriod, invMes, xBins, yBins)
    pcHistinv = compute2dHist(invPeriod[invDisp], invMes[invDisp], xBins, yBins)
    
    opssmall=(pcHistops<=3) | (allHistinv<20);
    
    E = 1 - pcHistinv/allHistinv
    
    R = 100.0 *(1.0 - ((allHistops-pcHistops)/pcHistops) * (1.0-E)/E)
    R[opssmall]=np.nan
        
    #print E
    #import pdb
    #pdb.set_trace()

    extent = [0, len(xBins)-1, 0, len(yBins)-1]
    plt.imshow(R, interpolation="nearest", cmap=plt.cm.YlGnBu_r, \
        origin="bottom", extent=extent)
    cb = plt.colorbar()
    cb.set_label("Percent Reliable")
    
    writeLabels(R, drange)
    
    ax = plt.gca()
    xLabels = map(lambda x: "%i" %(x), xBins)
    ax.xaxis.set_ticks(np.arange(len(xBins)))
    ax.xaxis.set_ticklabels(xLabels,fontsize=14)
    if xlabel=='metric':
        plt.xlabel(xmetric,fontsize=14)
    else:
        plt.xlabel(xlabel,fontsize=14)
        
    yLabels = map(lambda x: "%i" %(x), yBins)
    ax.yaxis.set_ticks(np.arange(len(yBins)))
    ax.yaxis.set_ticklabels(yLabels,fontsize=14)
    if ylabel=='metric':
        plt.ylabel(ymetric,fontsize=14) 
    else:
        plt.ylabel(ylabel,fontsize=14)
    plt.title(atitle,fontsize=14)
    plt.clim(drange)
    
    return R
def plotCompleteness(injdata,xmetric='period',ymetric='mes',xBins=[0, 10, 200, 500],yBins=[7, 10, 20, 2000],drange=(1,99),\
                        atitle="Reliability\n#PCs",s=[0.0,1.0],xlabel='metric',ylabel='metric'):
    """
    Plot grid of reliability numbers
    Similar to plotgrid above, but needs that to run twice to get the numbers 
    we need the dataframes for inversion and ops.
    """
    
    injPeriod=injdata[xmetric]
    injMes=injdata[ymetric]
    injDisp=passes(injdata,s) 

    allHistinj = compute2dHist(injPeriod, injMes, xBins, yBins)
    pcHistinj = compute2dHist(injPeriod[injDisp], injMes[injDisp], xBins, yBins)

    injsmall=allHistinj<=10;
    
    C = 100.0*pcHistinj/allHistinj

    C[injsmall]=np.nan
        
    #print E
    #import pdb
    #pdb.set_trace()

    extent = [0, len(xBins)-1, 0, len(yBins)-1]
    plt.imshow(C, interpolation="nearest", cmap=plt.cm.YlGnBu_r, \
         origin="bottom",extent=extent)
    cb = plt.colorbar()
    cb.set_label("Percent Candidates")
    
    writeLabels(C, drange)
    
    ax = plt.gca()
    xLabels = map(lambda x: "%i" %(x), xBins)
    ax.xaxis.set_ticks(np.arange(len(xBins)))
    ax.xaxis.set_ticklabels(xLabels,fontsize=14)
    if xlabel=='metric':
        plt.xlabel(xmetric,fontsize=14)
    else:
        plt.xlabel(xlabel,fontsize=14)
        
    yLabels = map(lambda x: "%i" %(x), yBins)
    ax.yaxis.set_ticks(np.arange(len(yBins)))
    ax.yaxis.set_ticklabels(yLabels,fontsize=14)
    if ylabel=='metric':
        plt.ylabel(ymetric,fontsize=14) 
    else:
        plt.ylabel(ylabel,fontsize=14)
    plt.title(atitle,fontsize=14)
    plt.clim(drange)
    
    return C
    
    #---
def HZplotSimple(koi, scut=0.5,pradcut=1.8, mfactor=15, addConf=True, \
            ylimits=(3200,6400), xlimits=(2.5,0.0), ebars=False):

    """
    Input koi dataframe
    Input score cut
    input scale factor for radii.
    addConf = True means add back in the confirmed planets we lots.
    
    Plots the habitable zone objects considering error bars of 1 sigma.
    as Insolation flux vs. Stellar temperature 
    It is colored by Score.
    
    """
    pcs=(koi.koi_score >= scut) & (koi.koi_pdisposition=='CANDIDATE')
    smallerrors=(-1*koi.koi_prad_err2/koi.koi_prad < 2.0 ) & (-1*koi.koi_insol_err2/koi.koi_insol<2.0)
    want=koi.koi_prad+koi.koi_prad_err2 <= pradcut
    wantlow=np.array(map( lambda x,y:hzCalc.inHZ(x,y), koi.koi_steff, koi.koi_insol+koi.koi_insol_err2))
    wantup=np.array(map( lambda x,y:hzCalc.inHZ(x,y), koi.koi_steff, koi.koi_insol+koi.koi_insol_err1))
    wantin=np.array(map( lambda x,y:hzCalc.inHZ(x,y), koi.koi_steff, koi.koi_insol))
    
    rearth=1
    rsupearth=2
    mfactor=15
    pars=['kepoi_name','kepler_name','koi_disposition','koi_score','koi_period',\
          'koi_prad','koi_prad_err2','koi_prad_err1','koi_steff','koi_steff_err1',\
          'koi_steff_err2','koi_insol','koi_insol_err1','koi_insol_err2','koi_max_mult_ev','koi_slogg']
    hzframe=koi[pcs & want & (wantlow | wantup | wantin) &smallerrors][pars]
    
    #Add in those three Confirmed ones we are missing.
    if addConf:
        missing=["009002278-04","008845205-01", "8142787-01"]#,"010604335-02","005640085-02"]  #[62f,283c,159c]
        for tce in missing:
            mydf=koi[koi.index==tce][pars]
            hzframe=hzframe.append(mydf)
    
    hzframe['RpBig']=hzframe.koi_prad+hzframe.koi_prad_err1
    hzframe['RpSmall']=hzframe.koi_prad+hzframe.koi_prad_err2
    
    print len(hzframe)

    plt.figure(figsize=(8.5,11))
    
    for i,hz in enumerate(hzframe.index):
        Srad=hzframe.loc[hz]['koi_insol']
        Smin=Srad+hzframe.loc[hz]['koi_insol_err2']
        Smax=Srad+hzframe.loc[hz]['koi_insol_err1']
        rp=hzframe.loc[hz]['koi_prad']
        rpBig=hzframe.loc[hz]['RpBig']
        rpSmall=hzframe.loc[hz]['RpSmall']
        tstar=hzframe.loc[hz]['koi_steff']
        tmin=tstar+hzframe.loc[hz]['koi_steff_err2']
        tmax=tstar+hzframe.loc[hz]['koi_steff_err1']
        sc=hzframe.loc[hz]['koi_score']
        kepname="%s" % (str(hzframe.loc[hz].kepoi_name))

        if ebars:
            plt.hlines(tstar,Smin,Smax,colors='steelblue',label='.',lw=1.3,zorder=1)
            plt.vlines(Srad,tmin,tmax,colors='steelblue',label='.',lw=1.3,zorder=1)
        
        plt.scatter(Srad,tstar,s=(rp*mfactor)**2,marker='o',c='cornflowerblue',linewidth=0.1,vmin=-.3,vmax=1,zorder=2)
        #plt.scatter(Srad,tstar,s=(rpSmall*mfactor)**2,marker='o',c="white",linewidth=0.2,zorder=3,alpha=0.95)
        #plt.scatter(Srad,tstar,s=(rp*mfactor)**2,marker='o',c=sc,cmap="brg_r",linewidth=.4,vmin=0,vmax=1,zorder=2)
        if np.float(hzframe.loc[hz].kepoi_name[-8:])>7621.01:
            plt.scatter(Srad,tstar,s=(rp*mfactor)**2,marker='o',c='darkseagreen',linewidth=0.1,vmin=-.3,vmax=1,zorder=2)
            plt.annotate(s=kepname,xy=(Srad-.1,tstar+8),xycoords='data',color="k",fontsize=10)
        #elif p.notnull(hzframe.loc[hz].kepler_name):
        #    plt.annotate(s=hzframe.loc[hz].kepler_name,xy=(Srad-.1,tstar+8),xycoords='data',color="grey",fontsize=10)
        
    #plt.scatter(hzframe.Srad,hzframe.tstar,s=(hzframe.RpBig*mfactor)**2,marker='o',c=hzframe.score,cmap="coolwarm_r",linewidth=0.1,vmin=0,vmax=1)
    #plt.scatter(hzframe.Srad,hzframe.tstar,s=(hzframe.RpSmall*mfactor)**2,marker='o',c='white',linewidth=0.1)
    
    small=plt.scatter(1,5777,s=(rearth*mfactor)**2,marker='o',edgecolor='blue',linewidth=1.5,label='R$_\oplus$',facecolors='none')
    big=plt.scatter(1.0,8300,s=(rsupearth*mfactor)**2,marker='o',c='grey',linewidth=1.4,label='2R$_\oplus$')
    half=plt.scatter(1.0,8500,s=(1.5*rearth*mfactor)**2,marker='o',edgecolor='blue',linewidth=1.5,label='1.5R$_\oplus$',facecolors='none')
    plt.xlabel('Insolation Flux ($S_{\oplus}$)', fontsize=12)
    plt.ylabel('Stellar Effective Temperature (K)', fontsize=12)
    
    plt.legend(handles=[small,half],scatterpoints=1,framealpha=0.8,loc='upper right',fontsize=14)
    
    #Draw a HZ across the plot.
    teff=np.linspace(2000,7000,50)
    hzBorders=hzCalc.hzBorder(teff)
    Shot=hzBorders[:,0]
    Scool=hzBorders[:,3]
    #hzhot=7000*sx-7500
    #hzcool=17500*sx-500
    plt.plot(Shot,teff,'--b')
    plt.plot(Scool,teff,'--r')
    #plt.title('Score cut of %f' % scut)
    
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    
    return hzframe
    #plt.xlim((2.5,0))
    #plt.ylim((2500,7000))
    #plt.savefig('%s/fig-hzTeffInsol2.png' %outdir)


def HZplot2(koi, scut=0.5,pradcut=1.8, mfactor=15, addConf=False, \
            ylimits=(3200,6400), xlimits=(2.5,0.0), ebars=True,mycolormap="YlGn"):

    """
    Input koi dataframe
    Input score cut
    input scale factor for radii.
    addConf = True means add back in the confirmed planets we lots.
    
    Plots the habitable zone objects considering error bars of 1 sigma.
    as Insolation flux vs. Stellar temperature 
    It is colored by Score.
    
    The circles are simple annulus with no information on radius error.
    Need to indicate which are New
    Need a color bar.
    
    """
    pcs=(koi.koi_score >= scut) & (koi.koi_pdisposition=='CANDIDATE')
    smallerrors=(-1*koi.koi_prad_err2/koi.koi_prad < 2.0 ) & (-1*koi.koi_insol_err2/koi.koi_insol<2.0)
    want=koi.koi_prad+koi.koi_prad_err2 <= pradcut
    wantlow=np.array(map( lambda x,y:hzCalc.inHZ(x,y), koi.koi_steff, koi.koi_insol+koi.koi_insol_err2))
    wantup=np.array(map( lambda x,y:hzCalc.inHZ(x,y), koi.koi_steff, koi.koi_insol+koi.koi_insol_err1))
    wantin=np.array(map( lambda x,y:hzCalc.inHZ(x,y), koi.koi_steff, koi.koi_insol))
    
    rearth=1
    rsupearth=2
    mfactor=12
    pars=['kepoi_name','kepler_name','koi_disposition','koi_score','koi_period',\
          'koi_prad','koi_prad_err2','koi_prad_err1','koi_steff','koi_steff_err1',\
          'koi_steff_err2','koi_insol','koi_insol_err1','koi_insol_err2','koi_max_mult_ev','koi_slogg']
    hzframe=koi[pcs & want & (wantlow | wantup | wantin) &smallerrors][pars]
    
    #Add in the one KOI that gets added after the stellar table erratum
    hzframe=hzframe.append(koi[koi.index=='012885212-02'])
    
    #Add in those three Confirmed ones we are missing.
    if addConf:
        missing=["009002278-04","008845205-01", "008142787-01","010604335-02","005640085-02"]  #[62f,283c,159c]etc
        for tce in missing:
            mydf=koi[koi.index==tce][pars]
            hzframe=hzframe.append(mydf)
    
   # hzframe['RpBig']=hzframe.koi_prad+hzframe.koi_prad_err1
   # hzframe['RpSmall']=hzframe.koi_prad+hzframe.koi_prad_err2
    
    #print len(hzframe)

    
    Srad=hzframe['koi_insol']
    Smin=Srad+hzframe['koi_insol_err2']
    Smax=Srad+hzframe['koi_insol_err1']
    rp=hzframe['koi_prad']
    tstar=hzframe['koi_steff']
    tmin=tstar+hzframe['koi_steff_err2']
    tmax=tstar+hzframe['koi_steff_err1']
    sc=hzframe['koi_score']
    kepnum=np.zeros(len(Srad))
    for i,v in enumerate(hzframe.kepoi_name):
        kepnum[i]=np.float(v[-8:])  #convert koi to a float.

    if ebars:
        plt.hlines(tstar,Smin,Smax,colors='steelblue',label='.',lw=1.3,zorder=1)
        plt.vlines(Srad,tmin,tmax,colors='steelblue',label='.',lw=1.3,zorder=1)

    #Plot a circle with gradient giving the score.    
    plt.scatter(Srad,tstar,s=(rp*mfactor)**2,marker='o',c=sc,cmap=mycolormap,linewidth=0.2,vmin=scut-.15,vmax=1,zorder=2)
    cb=plt.colorbar()
    cb.set_label("Disposition Score",fontsize=13)
    cb.set_ticks([1.0,0.9,0.8,0.7,0.6,0.5,0.4])
    new = (kepnum>7621.01) | (kepnum == 238.03)
    #Put a ring around the confirmed ones.
    plt.scatter(Srad[new],tstar[new],s=(rp[new]*mfactor)**2, marker='o', facecolors='none', edgecolors='r',linewidth=1.5,zorder=2)
    #plt.scatter(Srad[new],tstar[new],s=(mfactor**2)/8, marker='*', c='m',zorder=2)
    small=plt.scatter(1,9777,s=(rearth*mfactor)**2,marker='o',linewidth=0.2,label='R$_\oplus$',facecolors='b')
    big=plt.scatter(1.0,8300,s=(rsupearth*mfactor)**2,marker='o',c='b',linewidth=0.2,label='2R$_\oplus$',facecolors='b')
    ring=plt.scatter(1.0,9000,s=(rearth*mfactor)**2,marker='o',linewidth=2,facecolor='none',label='New in DR25',edgecolor='r')
    #half=plt.scatter(1.0,8500,s=(1.5*rearth*mfactor)**2,marker='o',edgecolor='blue',linewidth=1.5,label='1.5R$_\oplus$',facecolors='none')
    plt.xlabel('Insolation Flux ($S_{\oplus}$)', fontsize=12)
    plt.ylabel('Stellar Effective Temperature (K)', fontsize=12)
    
    plt.legend(handles=[small,big,ring],scatterpoints=1,framealpha=0.8,loc='upper right',fontsize=13)
    
    #Draw a HZ across the plot.
    teff=np.linspace(2000,7000,50)
    hzBorders=hzCalc.hzBorder(teff)
    Shot=hzBorders[:,0]
    Scool=hzBorders[:,3]

    plt.plot(Shot,teff,'--r')
    plt.plot(Scool,teff,'--b')

        
    
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    
    return hzframe
