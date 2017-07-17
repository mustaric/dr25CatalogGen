#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 15:42:58 2016

@author: smullall

This code is intended to create a one page report of the current performance of
the robovetter.

"""

import rvIO as io
import plotRobovetter as prv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import sys
import getopt as getopt
import pdb
import createRVPerfMetricsPlots as pmp

#plt.switch_backend('Agg') 

def plotMPGrid(data,title,climRange=(0,100)):
    """
    Plot the prv.PlotGrid of MES vs Period with given title and range.
    """
    period=data['period']
    mes=data['mes']
    passed=data['disp']=='PC'
    prv.plotGrid(period,mes,passed,climRange)
    plt.title(title)
    

    
def plotMPGrid2(data,title,climRange=(0,100),xBins=[0,10,200,500],yBins=[7,10,20,2000]):
    """
    Plot the prv.PlotGrid of MES vs Period with given title and range.
    """
    period=data['period']
    mes=data['mes']
    passed=data['disp']=='PC'
    prv.plotGrid(period,mes,passed,climRange,xBins=xBins,yBins=yBins)
    plt.title(title)
    
def plotMPGrid3(data,title,climRange=(0,100),xBins=[0,10,200,500],yBins=[7,10,20,2000]):
    """
    Plot the prv.PlotGrid of MES vs Period with given title and range.
    """
    period=data['period']
    mes=data['mes']
    passed=data['N']==0
    prv.plotGrid(period,mes,passed,climRange,xBins=xBins,yBins=yBins)
    plt.title(title)
    
def getData(topdir,rvfile,tcefile):
    """
    Return a pandas dataframe of the robovetter results.
    topdir is through the OPS or INJ, but does not include DATA
    """    
    rvfile='%s/%s' % (topdir,rvfile)
    tcefile='%s/%s' % (topdir,tcefile)

    rvdata=io.createRVResults(rvfile,tcefile)

    return rvdata

def getHZpercent(data,Srange=(0.75,1.25),pRadius=(0.75,1.25)):
    """
    Return the percentage of PCs in the HZ with small radius.
    Use the input ranges.
    """
    
    hzwant=(data['srad'] >=Srange[0]) & (data['srad']<=Srange[1]) & \
           (data['rplanet'] >=pRadius[0]) & (data['rplanet']<=pRadius[1])

    hzdata=data.loc[hzwant]    
    pcwant=hzdata['disp'] == 'PC'
    
    total=float(len(hzdata))
    pcs=float(len(hzdata[pcwant]))    
    percent=100 * (pcs/total) 
    
    return (percent,total,pcs,hzwant)
    
    
def getHZpercentS(data,Srange=(0.75,1.25),pRadius=(0.75,1.25),s=(0.0,1.0)):
    """
    Return the percentage of PCs in the HZ with small radius.
    Use the input ranges.
    """
    
    hzwant=(data['srad'] >=Srange[0]) & (data['srad']<=Srange[1]) & \
           (data['rplanet'] >=pRadius[0]) & (data['rplanet']<=pRadius[1])

    hzdata=data.loc[hzwant]
    pcwant=(hzdata['disp'] == 'PC') & (hzdata['score']>s[0])
    fpwant=(hzdata['disp']== 'FP') & (hzdata['score']>s[1])
    passed=pcwant | fpwant    
    
    total=float(len(hzdata))
    pcs=float(len(hzdata[passed]))    
    percent=100 * (pcs/total) 
    
    return (percent,total,pcs,hzwant)

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

def getSvnNum(fname):
    """
    Get the first line of file that contains the id string
    """
    fid=open(fname,'r')
    line=fid.readline()
    line=line.replace('$',' ')
    
    return line

import pandas as p
def deleteKOIs(invdata,feddata):
    """
    Given the federation frame that is listed by TCE number
    and, the inverted data set, remove those rows that 
    are associated with a KOI.
    """
    newdata=invdata.merge(feddata,right_index=True,left_index=True,how='left')
    
    nokoi=p.isnull(newdata['koi'])
    
    returndata=newdata[nokoi]
    
    return returndata
    
def trimINVData(invdata,rtype):
    """
    Trim the inverted data of the known astrophysical events
    """
    
    if rtype == "INV":
        keepfile="/soc/nfs/so-nfs/DR25/INV/DATA/invTCEClean-20161114.csv"
    elif rtype == "SS1":
        #This is temporary file and should be changed.
        keepfile="/soc/nfs/so-nfs/DR25/SS1/DATA/ss1TCEClean-900day-7mes-Dec2016.csv"
    #keepfile="/soc/nfs/so-nfs/DR25/INV/DATA/invTCEClean-300day-09282016.csv"
    keepinv=p.read_csv(keepfile,header=0,index_col='tce',comment='#')
    newdata=p.merge(invdata,keepinv,how="left",left_index=True,right_index=True)
    returnData=newdata[newdata['keep']]    
    
    return returnData

def usage():
    """
    Description of code
    """
    print "create robovetter Preformance Metrics"
    print " Need to specify top directory, injection, inversion, ops and outputfile"
    print "ht:j:v:o:f:d:i:"
    print "Default dir is /soc/nfs/so-nfs/DR25 for -t and -i \n\n"
    

def main():
    """
    Trying to make this plot generator command line
    """
 
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ht:j:v:o:f:d:i:r:e:", ["help", "inj=","inv=","output=","ops=","title=","dir=","tcedir=","rtype="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print err
        usage()
        sys.exit()
        
    outplot="PerformanceMetrics.png"
    topdir='/soc/nfs/so-nfs/DR25/'
    tcedir='/soc/nfs/so-nfs/DR25/'
    figTitle='DR25 RV Analysis'
    tcefile='DATA/TCEs.txt'
    opsfile='RoboVetterOut-OPS.txt'
    invfile='RoboVetterOut-INV.txt'
    injfile='RoboVetterOut-INJ-PlanetOn.txt'
    ebsfile='RoboVetterOut-INJ-EB.txt'
    rtype='INV'
    revisionnum=62300  #Currently unused
    invfedfile='/soc/nfs/so-nfs/DR25/INV/koimatch_DR25INV_03172016.txt'
    
    for o, a in opts:
        if o in ("-f","--output"):
            outplot = a
            print "Output File is: %s\n" % outplot 
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-o", "--ops"):
            opsfile= a
            print "OPS file is %s%s\n" % (topdir,opsfile)
        elif o in ("-j", "--inj"):
            injfile= a
            print "INJ file is %s%s\n" % (topdir,injfile)
        elif o in ("-v", "--inv"):
            invfile= a
            print "INV file is %s%s\n" % (topdir,invfile)
        elif o in ("-d", "--dir"):
            topdir=a
        elif o in ("-t", "--title"):
            figTitle=a
        elif o in ("-i","--tcedir"):
            tcedir=a
        elif o in ("-r","--revision"):
            revisionnum=a
        elif o in ("-e","--rtype"):   #Reliability Type SS1 or INV
            rtype=a
        else:
            assert False, "Unhandled option"
            sys.exit()   
    
    invfile='RoboVetterOut-%s.txt' % (rtype)  #Get invfile from either INV or SS1    
    
    opsfile="%s%s" % (topdir,opsfile)
    invfile="%s%s" % (topdir,invfile)
    injfile="%s%s" % (topdir,injfile)
    ebsfile="%s%s" % (topdir,ebsfile)
    
    
    Srange=(0.25,1.7)
    pRadius=(0.25,1.7)
#    Srange=(0.25,2.0)
#    pRadius=(0.25,2.0)

    #outplot='/home/smullall/Kepler/RoboVetter/DR25/stats/dr25-dr24rv-performace.pdf'

    #opsfile='OPS/RoboVet/RoboVetterOut-OPS-vDR24.txt'
    #invfile='INV/RoboVet/RoboVetterOut-INV-vDR24.txt'
    #injfile='INJ/RoboVet/RoboVetterOut-INJ-PlanetOn.txt'
    #tcefile='DATA/TCEs.txt'

    
    hzpercent=np.ones((3,2),dtype=float)
    hztotal=np.ones((3,2),dtype=int)
    hzpcs=np.ones((3,2),dtype=int)
    
    svnId=getSvnNum(opsfile)    
    
    #OPS
    tfile="%s/OPS/%s" % (tcedir,tcefile)
    rvdata=io.createRVResults(opsfile,tfile)
    opsdata=rvdata
    opsdata2,a,b,c=io.loadRVInOut(revisionnum,type='OPS')

    (hzpercent[0,0],hztotal[0,0],hzpcs[0,0],hzopswant)=getHZpercent(opsdata,Srange=Srange,pRadius=pRadius)
    
    #Create for FGK Dwarf Plot
    fgkRvWant=(opsdata['logg']>=4.0) & (opsdata['tstar']>=4000.0) & (opsdata['tstar']<7000.0);
    opsfgkData=opsdata[fgkRvWant]
    opsFgkWant=fgkRvWant    
    (hzpercent[0,1],hztotal[0,1],hzpcs[0,1],hzopsfgkwant)=getHZpercent(opsfgkData,Srange=Srange,pRadius=pRadius)
    
    #INV
    if rtype == 'ALL':
        tfile="%s/%s/%s" % (tcedir,"INV",tcefile)
        invfile='RoboVetterOut-%s.txt' % ("INV")
        invdata=io.createRVResults(invfile,tfile)
        invrvdata=trimINVData(invdata,"INV")
        tfile="%s/%s/%s" % (tcedir,"SS1",tcefile)
        invfile='RoboVetterOut-%s.txt' % ("SS1")
        invdata=io.createRVResults(invfile,tfile)
        ss1rvdata=trimINVData(invdata,"SS1")
        newrvdata=p.concat((invrvdata,ss1rvdata),ignore_index=True)
        
        invdata2,a,b,c=io.loadRVInOut(revisionnum,type='INV')
        invrvdata2=trimINVData(invdata2,"INV")
        ss1data2,a,b,c=io.loadRVInOut(revisionnum,type='SS1')
        ss1rvdata2=trimINVData(ss1data2,"SS1")
        newrvdata2=p.concat((invrvdata2,ss1rvdata2),ignore_index=True)
        
    else:
        tfile="%s/%s/%s" % (tcedir,rtype,tcefile)
        invdata=io.createRVResults(invfile,tfile)
        newrvdata=trimINVData(invdata,rtype)
        
        invdata2,a,b,c=io.loadRVInOut(revisionnum,type=rtype)
        newrvdata2=trimINVData(invdata2,rtype)
 
    invdata=newrvdata
    invdata2=newrvdata2   #This contains the RV input and output.. Columns are a superset of invdata
    print "INV data Len: %f" % len(invdata)
    
    (hzpercent[1,0],hztotal[1,0],hzpcs[1,0],hzinvwant)=getHZpercent(invdata,Srange=Srange,pRadius=pRadius)
        
    #Create for FGK Dwarf Plot
    fgkRvWant=(invdata['logg']>=4.0) & (invdata['tstar']>=4000.0) & (invdata['tstar']<7000.0);
    invfgkData=invdata[fgkRvWant]
    invFgkWant=fgkRvWant
    (hzpercent[1,1],hztotal[1,1],hzpcs[1,1],hzinvfgkwant)=getHZpercent(invfgkData,Srange=Srange,pRadius=pRadius)
    
    #INJ
    tfile="%s/INJ/%s" % (tcedir,tcefile)
    injdata=io.createRVResults(injfile,tfile)
    (hzpercent[2,0],hztotal[2,0],hzpcs[2,0],hzinjwant)=getHZpercent(injdata,Srange=Srange,pRadius=pRadius)    
    
    #INJ FGK Dwarf Plot
    injFgkRvWant=(injdata['logg']>=4.0) & (injdata['tstar']>=4000.0) & (injdata['tstar']<7000.0);
    injfgkData=injdata[injFgkRvWant]
    injFgkWant=injFgkRvWant
    (hzpercent[2,1],hztotal[2,1],hzpcs[2,1],hzinjwant)=getHZpercent(injfgkData,Srange=Srange,pRadius=pRadius)    
     
    
    #Create Plot
    climrange=(0,100)
    figureTitle="%s\n  %s relType=%s" %(figTitle,svnId,rtype)
    outname="performance-%s" % (outplot)        
    pmp.createScoreCardScore(opsdata,opsfgkData,invdata,invfgkData,injdata,injfgkData,Srange,pRadius,\
                        hzpercent,hztotal,hzpcs,hzinvwant,hzopswant,climrange,figureTitle,outname)
    
    #Create Plot if do a score cut
    scores=(0.7,0.5)  #PC score cut, FP scorecut

    climrange=(0,100)
    figureTitleScore="ScoreCut ( %f %f) \n %s\n  %s relType=%s" %(scores[0],scores[1],figTitle,svnId,rtype)
    outname="performanceScore-%s" % (outplot) 

    (hzpercent[0,0],hztotal[0,0],hzpcs[0,0],hzopswantS)=getHZpercentS(opsdata,\
                                        Srange=Srange,pRadius=pRadius,s=scores)  
    (hzpercent[0,1],hztotal[0,1],hzpcs[0,1],hzopsfgkwantS)=getHZpercentS(opsfgkData,\
                                        Srange=Srange,pRadius=pRadius,s=scores)
    (hzpercent[1,0],hztotal[1,0],hzpcs[1,0],hzinvwantS)=getHZpercentS(invdata,\
                                        Srange=Srange,pRadius=pRadius,s=scores)
    (hzpercent[1,1],hztotal[1,1],hzpcs[1,1],hzinvfgkwantS)=getHZpercentS(invfgkData,\
                                        Srange=Srange,pRadius=pRadius,s=scores)
    (hzpercent[2,0],hztotal[2,0],hzpcs[2,0],hzinjwantS)=getHZpercentS(injdata,\
                                        Srange=Srange,pRadius=pRadius, s=scores)
    (hzpercent[2,1],hztotal[2,1],hzpcs[2,1],hzinjwantS)=getHZpercentS(injfgkData,\
                                        Srange=Srange,pRadius=pRadius,s=scores)  
                    
    pmp.createScoreCardScore(opsdata,opsfgkData,invdata,invfgkData,injdata,injfgkData,Srange,pRadius,\
                        hzpercent,hztotal,hzpcs,hzinvwantS,hzopswantS,climrange,figureTitleScore,outname,scores) 
       
    
    #pdb.set_trace()
    #Reliability Plot
    #Get information about the HZ zone for inv and ops for HZ.
    R,E = estimateReliability(invdata[hzinvwant],opsdata[hzopswant])    
    hzRelStr="HZ (%4.2f, %4.2f): E= %5.2f%%  R=%5.2f%%" % (Srange[0],Srange[1],E*100,R*100)
    
    prv.plotReliability(invdata,opsdata,xBins=[0, 10, 200, 500],yBins=[7, 10, 20, 2000],drange=(1,99))
    plt.annotate(figureTitle,xy=(0.20,0.9),xycoords='figure fraction',  fontsize=11)
    plt.annotate(hzRelStr,xy=(0.11,0.06),xycoords='figure fraction', fontsize=11)
    outplot2="rel-%s" % (outplot)
    plt.savefig(outplot2)    
    
    #Reliability plot with score cut
    R,E = estimateReliability(invdata[hzinvwantS],opsdata[hzopswantS],s=scores)    
    hzRelStr="HZ (%4.2f, %4.2f): E= %5.2f%%  R=%5.2f%%" % (Srange[0],Srange[1],E*100,R*100)
    
    prv.plotReliability(invdata,opsdata,xBins=[0, 10, 200, 500],yBins=[7, 10, 20, 2000],drange=(1,99),s=scores)
    plt.annotate(figureTitleScore,xy=(0.20,0.9),xycoords='figure fraction',  fontsize=11)
    plt.annotate(hzRelStr,xy=(0.11,0.06),xycoords='figure fraction', fontsize=11)
    outplot2="relScore-%s" % (outplot)
    plt.savefig(outplot2)   

    #One Page
    pmp.createOnePageStats(opsdata,injdata,invdata,scores=[0.0,1.0],figureTitle=figureTitle,hz=Srange)
    outplot2b="Aonepage-%s" % (outplot)
    plt.savefig(outplot2b)
    
    #One Page
    pmp.createOnePageStats(opsdata,injdata,invdata,scores=scores,figureTitle=figureTitle,hz=Srange)
    outplot2b="Asonepage-%s" % (outplot)
    plt.savefig(outplot2b)
    
    #Marginalized Plots
    pmp.marginalizedRunPlots(opsdata2,invdata2,scores=scores,figureTitle=figureTitle,hz=Srange)
    outplot2b="Arunmarginalized-%s" % (outplot)
    plt.savefig(outplot2b)
    
    pmp.marginalizedRunPlots(opsdata2,invdata2,scores=scores,figureTitle=figureTitle,hz=Srange,Effect=True)
    outplot2b="ArunmarginalizedEff-%s" % (outplot)
    plt.savefig(outplot2b)

    #TCE Histogram
    pmp.plotTceDist(opsdata,invdata,bins=200,metrics=('period','mes'),labels=('OPS',rtype))
    plt.annotate(svnId,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10)
    outplot3e="tceHist-%s" % (outplot)
    plt.savefig(outplot3e)

    #Histogram Figure
    plt.figure(figsize=(8.5,5.5))
    pcwant=opsdata['disp'] == 'PC'
    pcops=opsdata[pcwant]
    cumfile='/soc/nfs/so-nfs/DR25/other/nexsciCumulative.csv'
    cumdata=io.readNexsciCumulative(cumfile)
    cumpcs=(cumdata['koi_pdisposition'] == 'CANDIDATE') & (~np.isnan(cumdata['koi_max_mult_ev']))    
    cumperiod=cumdata.loc[cumpcs,'koi_period']
    cummes=cumdata.loc[cumpcs,'koi_max_mult_ev']
    x=(np.log10(pcops['period']),np.log10(cumperiod))
    y=(np.log10(pcops['mes']),np.log10(cummes))    
    
    axHisty,axHistx = prv.plot2dMulti(x,y,[-.35,2.95],[.83,1.7],nxbins=100,nybins=80,\
                            xlabel="log(Period)",ylabel="log(MES)",showPoints=False)
    handles,labels=axHistx.get_legend_handles_labels()
    labels[0]='DR25 Candidates'
    labels[1]='DR24 Candidates'
    axHistx.legend(handles,labels,bbox_to_anchor=(1.6, 1) )
    plt.annotate(svnId,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10) 	
    outplot3="histops-%s" % (outplot)
    plt.savefig(outplot3)
    
    #Histogram Figure  for good scores
    plt.figure(figsize=(8.5,5.5))
    pcwant=(opsdata['disp'] == 'PC') & (opsdata['score']>scores[0])
    fpwant=(opsdata['disp']== 'FP') & (opsdata['score']>scores[1])
    pcops=opsdata[pcwant | fpwant]
    cumfile='/soc/nfs/so-nfs/DR25/other/nexsciCumulative.csv'
    cumdata=io.readNexsciCumulative(cumfile)
    cumpcs=(cumdata['koi_pdisposition'] == 'CANDIDATE') & (~np.isnan(cumdata['koi_max_mult_ev']))    
    cumperiod=cumdata.loc[cumpcs,'koi_period']
    cummes=cumdata.loc[cumpcs,'koi_max_mult_ev']
    x=(np.log10(pcops['period']),np.log10(cumperiod))
    y=(np.log10(pcops['mes']),np.log10(cummes))    
    
    axHisty,axHistx = prv.plot2dMulti(x,y,[-.35,2.95],[.83,1.7],nxbins=100,nybins=80,
                            xlabel="log(Period)",ylabel="log(MES)",showPoints=False)
    handles,labels=axHistx.get_legend_handles_labels()
    labels[0]='DR25 Candidates \n(w score >%3.1f)\nFPs with score>%3.1f' % (scores)
    labels[1]='DR24 Candidates'
    axHistx.legend(handles,labels,bbox_to_anchor=(1.6, 1) )
    plt.annotate(svnId,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10) 	
    outplot3b="histopsScore-%s" % (outplot)
    plt.savefig(outplot3b)    

    #Histogram Figure  for number of transits >=4
    plt.figure(figsize=(8.5,5.5))
    pcwant=(opsdata2['disp'] == 'PC') & (opsdata2['ngood']>=4)
    pcops=opsdata2[pcwant]
    cumfile='/soc/nfs/so-nfs/DR25/other/nexsciCumulative.csv'
    cumdata=io.readNexsciCumulative(cumfile)
    cumpcs=(cumdata['koi_pdisposition'] == 'CANDIDATE') & (~np.isnan(cumdata['koi_max_mult_ev']))    
    cumperiod=cumdata.loc[cumpcs,'koi_period']
    cummes=cumdata.loc[cumpcs,'koi_max_mult_ev']
    x=(np.log10(pcops['period']),np.log10(cumperiod))
    y=(np.log10(pcops['mes']),np.log10(cummes))    
    
    axHisty,axHistx = prv.plot2dMulti(x,y,[-.35,2.95],[.83,1.7],nxbins=100,nybins=80,
                            xlabel="log(Period)",ylabel="log(MES)",showPoints=False)
    handles,labels=axHistx.get_legend_handles_labels()
    labels[0]='DR25 Candidates \n(w ntransits>=4)'
    labels[1]='DR24 Candidates'
    axHistx.legend(handles,labels,bbox_to_anchor=(1.6, 1) )
    plt.annotate(svnId,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10) 	
    outplot4b="histopsNgood-%s" % (outplot)
    plt.savefig(outplot4b) 
    
    #Histogram Figure of Srad Prad for full population
    plt.figure(figsize=(8.5,5.5))
    plt.subplot(211)
    these=(opsdata['rplanet'] == 0) | (opsdata['rplanet'] == np.NaN)
    opsdata.rplanet[these]= 1e-10
    these=opsdata['srad'] == 0
    opsdata.srad[these] = 1e-10
    pcwant=opsdata['disp'] == 'PC'
    pcops=opsdata[pcwant]
   
    y=(np.log10(pcops['rplanet']),np.log10(opsdata['rplanet']))
    x=(np.log10(pcops['srad']),np.log10(opsdata['srad']))    

    axHisty,axHistx = prv.plot2dMulti(x,y,[-1.5,5.0],[-.6,2.5],nxbins=80,nybins=80,xlabel="log(insolation Flux)",ylabel="log(Planet Radius)")
    handles,labels=axHistx.get_legend_handles_labels()
    labels[0]='DR25 Candidates'
    labels[1]='DR25 TCEs'
    axHistx.legend(handles,labels,bbox_to_anchor=(1.6, 1) )    
    
    plt.annotate(svnId,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10) 	
    outplot3b="histopssRS-%s" % (outplot)
    plt.savefig(outplot3b)
    
    
    #Compare to Confirmed and FPWG
    #Hardwired these files
    confirmed='/soc/nfs/so-nfs/DR25/other/keplernames.csv'
    fedFile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
    fpFile='/soc/nfs/so-nfs/DR25/other/fpwg.csv'
    
    fpdata=io.readConfirmed(fpFile)
    cdata=io.readConfirmed(confirmed)
    feddata=io.readFederation(fedFile)
    
    plt.figure(figsize=(8.5,4))
    plt.subplots_adjust(hspace=0.2, wspace=0.2)    
    plt.subplot(121)
    prv.plotFpwgPCGrid(opsdata,fpdata,feddata)
    plt.subplot(122)
    prv.plotConfirmedPCGrid(opsdata,cdata,feddata)
    plt.annotate(svnId,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10) 
    outplot4="confp-%s" % (outplot)
    plt.savefig(outplot4)
    
    #Also do the By Eye Vetting
    pcEye,fpEye=io.loadByEyeVetRvin(rvfile=opsfile)
    pcperiod=pcEye['period']
    pcmes=pcEye['mes']    
    pcpassed=pcEye['disp']=='PC'
    fpperiod=fpEye['period']
    fpmes=fpEye['mes']
    fppassed=fpEye['disp']=='PC'  
    
    plt.figure(figsize=(8.5,4))
    plt.subplots_adjust(hspace=0.2, wspace=0.2)    
    plt.subplot(121)
    prv.plotGrid(pcperiod,pcmes,pcpassed,(0,100))
    plt.title('By Eye PCs Vetted by RV')
    
    plt.subplot(122)
    prv.plotGrid(fpperiod,fpmes,fppassed,(0,100))
    plt.title('By Eye FPs Vetted by RV')
    
    plt.annotate(svnId,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10) 
    outplot4b="byEye-%s" % (outplot)
    plt.savefig(outplot4b)    
    
    
    
    #Report Card -- Reliability and Completeness for FGK and all
#    plt.figure(figsize=(8.5,8))
#    plt.subplots_adjust(hspace=0.25, wspace=0.45)    
#    plt.subplot(221)
#    #Plot reliability
#    prv.plotOnlyReliability(invdata,opsdata,atitle="All Reliability\n#PCs/#FPs")
#    
#    plt.subplot(222)
#    #Plot completeness
#    plotMPGrid(injdata,'Completeness\nAll INJ PC Rate',(climrange[1]/2,climrange[1]))
#    
#    plt.subplot(223)
#    #FGK reliability
#    prv.plotOnlyReliability(invdata[invFgkWant],opsdata[opsFgkWant],atitle="FGK Reliability\n#PCs/#FPs")
#
#    plt.subplot(224)    
#    #FGK Completeness
#    plotMPGrid(injdata[injFgkWant],'FGK Completeness\nINJ PC Rate',(climrange[1]/2,climrange[1]))
#    
#    plt.annotate(svnId,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10) 
    outplot5="relComp-%s" % (outplot)
    pmp.completenessAndReliability(opsdata,invdata,injdata,opsFgkWant,invFgkWant,injFgkWant,climrange,outplot5,figureTitle)
    outplot5="relCompScore-%s" % (outplot)
    pmp.completenessAndReliability(opsdata,invdata,injdata,opsFgkWant,invFgkWant,injFgkWant,climrange,outplot5,figureTitleScore,s=scores)
    
    #Fine Grid Reliability and Completeness
    outplot5b="relCompFineGrid-%s" % (outplot)
    pmp.fineGridRelCompleteness(opsdata,invdata,injdata,opsFgkWant,invFgkWant,injFgkWant,climrange,outplot5b,figureTitle)
    
    outplot5c="relCompFineGridScore-%s" % (outplot)
    pmp.fineGridRelCompleteness(opsdata,invdata,injdata,opsFgkWant,invFgkWant,injFgkWant,climrange,outplot5c,figureTitleScore,s=scores)
    
    
    xb=[0,50,150,350,410,700]
    #Fine Grid Reliability and Completeness
    outplot5d="relCompFineGrid2-%s" % (outplot)
    pmp.fineGridRelCompleteness(opsdata,invdata,injdata,opsFgkWant,invFgkWant,injFgkWant,climrange,outplot5d,figureTitle,xbins=xb)
    
    outplot5e="relCompFineGrid2Score-%s" % (outplot)
    pmp.fineGridRelCompleteness(opsdata,invdata,injdata,opsFgkWant,invFgkWant,injFgkWant,climrange,outplot5e,figureTitleScore,s=scores,xbins=xb)
    
    
    
    #Fine Grid Reliability and Completeness removing complicated federations from injection
    #injfedfile='/soc/nfs/so-nfs/DR25/INJ/injmatch_DR25_03182016.txt'
    #injfed=io.readFederationInj(injfedfile)
    
    #wantnot1=(injfed['pratio']<0.995) | (injfed['pratio']>1.005)
    #injcompfed=injfed[wantnot1].index
    #shortinjdata=injdata.drop(injcompfed,inplace=False,errors='ignore')    
    
    #outplot5d="relCompFineGridCPure-%s" % (outplot)
    #pmp.fineGridRelCompleteness(opsdata,invdata,shortinjdata,opsFgkWant,invFgkWant,injFgkWant,climrange,outplot5d,figureTitle +'/nNo Pure Inj Federation')
  
    
    #Grid Completeness for making KOIs (N flag not set)
    plt.figure(figsize=(12,6))
    plt.subplots_adjust(hspace=0.25, wspace=0.45)    
    
    xbins=[0,10,200,500]
    ybins=[7,10,20,2000]
    
    plt.subplot(121)
    #Plot completeness
    plotMPGrid3(injdata,'Completeness\nAll INJ Transit-like (N=0) Rate',(climrange[1]/2,climrange[1]),xBins=xbins,yBins=ybins)
    
    plt.subplot(122)    
    #FGK Completeness
    plotMPGrid3(injdata[injFgkWant],'FGK Completeness\nINJ Transit-like (N=0) Rate',(climrange[1]/2,climrange[1]),xBins=xbins,yBins=ybins)
    
    plt.annotate(svnId,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10) 
    outplot5c="reCompNflag-%s" % (outplot)
    plt.savefig(outplot5c)
    
    
    
    
    
    #Plot pradius against Period for all PCs.
    
    outplot6="pradiusA-%s" % (outplot)
    pmp.plotPeriodRadius(opsdata,outplot6, figureTitle, label="OPS",colorkey="mes",colorlim=[7,20])    
    outplot6="pradiusAScore-%s" % (outplot)
    pmp.plotPeriodRadius(opsdata,outplot6, figureTitleScore, label="OPS",colorkey="mes",colorlim=[7,20],s=scores)
    outplot6="pradiusB-%s" % (outplot)
    pmp.plotPeriodRadius(opsdata,outplot6, figureTitle, label="OPS",colorkey="score",colorlim=[0,1])
    outplot6="pradiusBScore-%s" % (outplot)
    pmp.plotPeriodRadius(opsdata,outplot6, figureTitleScore, label="OPS",colorkey="score",colorlim=[0,1],s=scores)
    outplot6="pradiusC-%s" % (outplot)
    pmp.plotPeriodRadius(invdata,outplot6, figureTitle, label=rtype,colorkey="mes",colorlim=[7,20])
    outplot6="pradiusCScore-%s" % (outplot)
    pmp.plotPeriodRadius(invdata,outplot6, figureTitleScore, label=rtype,colorkey="mes",colorlim=[7,20],s=scores)
    
    tfile="%s/OPS/%s" % (tcedir,tcefile)
    fedfile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016.txt'
    cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'
    opsDataCum=io.createAllResults(opsfile,tfile,fedfile,cumfile)    
    want=((opsDataCum['iskoi']==0) | (opsDataCum['dr24disp']=='FP'))
    data=opsDataCum[want]
    outplot6="pradiusD-%s" % (outplot)
    pmp.plotPeriodRadius(data,outplot6, figureTitle, label="New Candidates ",colorkey="score",colorlim=[0,1])

    outplot7="insolA-%s" % (outplot)    
    pmp.plotInsolationRadius(opsdata,outplot7, figureTitle,label="OPS",colorkey="tstar",colorlim=[3000,7000],s=[0.0,1.0])
    outplot8="insolAScore-%s" % (outplot)
    pmp.plotInsolationRadius(opsdata,outplot8, figureTitleScore,label="OPS",colorkey="tstar",colorlim=[3000,7000],s=scores)
    outplot9="insolAnew-%s" % (outplot)    
    pmp.plotInsolationRadius(data,outplot9, figureTitle,label="OPS New",colorkey="tstar",colorlim=[3000,7000],s=[0.0,1.0])
    
    pmp.plotPCs(opsdata,invdata,'period','mes',xlim=(10,600),ylim=(7,20),figTitle=figureTitle)
    outplot10="comparePCs-%s" % (outplot)
    plt.savefig(outplot10)
    
    pmp.plotPCs(opsdata,invdata,'period','mes',score=scores,xlim=(10,600),ylim=(7,20),figTitle=figureTitleScore)
    outplot10="comparePCscore-%s" % (outplot)
    plt.savefig(outplot10)

    #INJ - EBs -- Plot in teh grid the efficiency creating these as EBs.
    tfile="%s/INJ/%s" % (tcedir,tcefile)
    ebsdata=io.createRVResults(ebsfile,tfile)   
    outname="ebEffic-%s" % (outplot)
    pmp.ebEfficiency(ebsdata,outname,figureTitle)
    
    
    range1=(200,500)
    range2=(7,10)    
    prv.plotRelCompScore(opsdata,invdata,injdata,'period','mes',range1,range2,scores=np.arange(0,1,.05),\
    fpscore=scores[1],Rlim=(.2,1.01),Clim=(1.0,.2))
    plt.title('Adjust Score for Box 2, fpscore=%.2f' % scores[1])
    plt.annotate(svnId,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10)
    outname="B-adjScore-%s" % (outplot)
    plt.savefig(outname)


if __name__ == "__main__":
    main()
