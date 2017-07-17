# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 13:16:39 2016

@author: smullall

Cleaning up the creatRVPerfMetrics code by putting the plotting
in a separate file in functions so that they are easier to use.

"""
import matplotlib.pyplot as plt
import pandas as p
import numpy as np
import matplotlib.gridspec as gridspec
import plotRobovetter as prv


def passes(data,s=[0.0,1.0]):
    """
    Use disposition and scores to determine which ones pass.
    Returns true false array
    """
    
    passed= ((data['disp'] == 'PC') & (data['score']>=s[0])) | \
            ((data['disp']=='FP') & (data['score']>=s[1]))  
    return passed

def plotMPGrid(data,title,climRange=(0,100)):
    """
    Plot the prv.PlotGrid of MES vs Period with given title and range.
    """
    period=data['period']
    mes=data['mes']
    passed=data['disp']=='PC'
    prv.plotGrid(period,mes,passed,climRange)
    plt.title(title)
    

    
def plotMPGrid2(data,title,climRange=(0,100),xBins=[0,10,200,500],yBins=[7,10,20,2000],s=[0.0,1.0]):
    """
    Plot the prv.PlotGrid of MES vs Period with given title and range.
    """
    period=data['period']
    mes=data['mes']
    passed=passes(data,s)
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
    plt.title(title,fontsize=11)
    
    
def plotMPGridS(data,title,climRange=(0,100),s=(0.0,1.0)):
    """
    Plot the prv.PlotGrid of MES vs Period with given title and range.
    Consider Score in the process
    """
    period=data['period']
    mes=data['mes']
    passed=passes(data,s)
    prv.plotGrid(period,mes,passed,climRange)
    plt.title(title,fontsize=11)
    

def estimateReliability(invdata,opsdata,s=(0.0,1.0)):
    """
    Given an inversion run cut down to the population you want to consider.
    Return a reliability number.
    1-Reliability = Nfp/Npc * (1-E)/E
    E is the effectiveness as measured by Inversion.
    """
    
    invwant=((invdata['disp'] == 'PC') & (invdata['score']>=s[0])) |\
            ((invdata['disp']=='FP') & (invdata['score']>=s[1]))
            
    opsPC=((opsdata['disp'] == 'PC') & (opsdata['score']>=s[0])) |\
            ((opsdata['disp']=='FP') & (opsdata['score']>=s[1]))   
    
    if len(invwant)>0:
        E=len(invwant[~invwant])/float(len(invwant))
        if (E !=0) & (len(opsPC[opsPC])!=0):
            R= 1-len(opsPC[~opsPC])/np.float(len(opsPC[opsPC])) * (1.0-E)/E
        else:   
            R=-1
    else:
        E=-1
        R=-1
        
    print len(invwant[invwant])
    print len(opsPC[opsPC])
    return R,E


def estimateReliabilityErr(invdata,opsdata,s=(0.0,1.0)):
    """
    Given an inversion run cut down to the population you want to consider.
    Return a reliability number.
    1-Reliability = Nfp/Npc * (1-E)/E
    E is the effectiveness as measured by Inversion.
    """
    
    invwant=((invdata['disp'] == 'PC') & (invdata['score']>=s[0])) |\
            ((invdata['disp']=='FP') & (invdata['score']>=s[1]))
            
    opsPC=((opsdata['disp'] == 'PC') & (opsdata['score']>=s[0])) |\
            ((opsdata['disp']=='FP') & (opsdata['score']>=s[1]))   
    
    if len(invwant)>0:
        E=len(invwant[~invwant])/float(len(invwant))
        if (E !=0) & (len(opsPC[opsPC])!=0):
            R= 1-len(opsPC[~opsPC])/len(opsPC[opsPC]) * (1.0-E)/E
        else:   
            R=-1
    else:
        E=-1
        R=-1
        
    if (R !=-1) & (E!=-1):
        errR,errE=reliabilityErrorPoisson(R,E,len(opsPC[~opsPC]),len(opsPC[opsPC]),len(invwant),len(invwant[~invwant]))        
    else:
        errR=0
        
    return R,errR,E,errE
    
def reliabilityErrorPoisson(R,E,Nfp,Npc,Nfa,Nffa):
    """
    return the poisson derived error on reliability
    R = measured reliability
    E = effectiveness
    Nfp number of ops false positives
    Npc number of ops planet candidates
    Nfa number of false alarms (inversion, ss1)
    Nffa number of false alarms determined to be false alarms
    """
    print Nfp,Npc,Nfa,Nffa
    
    if E == 1.0:
        E = 0.9999
    else:
        E=np.float(E)
        
    errE=np.sqrt(Nffa)/Nfa
    
    #Check this algebra before you use. It should be for the full error term wtih poisson errors on Nfp and Nffa
    #errRoverRsq = (R**2)*(np.float(Nfp)/np.float(Npc)) * errE/(E)**2 + 1.0/np.float(Nfp) + 1.0/np.float(Nffa)   
    
    #errR = np.float(R) * np.sqrt(errRoverRsq)
    
    errR = (np.float(Nfp)/np.float(Npc)) * errE/(E)**2  #(approximation considering only the one dominate term)
    #Full thing
    errR = np.sqrt(np.float(Nfp)/np.float(Npc) * errE/(E)**2 + R/np.float(Nfp) + R/np.float(Npc) )
    
    
    return errR,errE
    
def reliabiltiyErrorBinomial(R,E,Nfp,Npc,Nfa,Npcfa):
    """
    return the poisson derived error on reliability
    R = measured reliability
    E = effectiveness
    Nfp number of ops false positives
    Npc number of ops planet candidates
    Nfa number of total false alarms (inversion, ss1)
    Npcfa number of false alarms determined to be PCs
    """
    
    Ntot=Nfp+Npc
    sig_Npc=np.sqrt(Npc - (Npc**2)/np.float(Ntot))    
    sig_E=np.sqrt(Npcfa)/np.float(Nfa)    
    
    part1= (Ntot*sig_Npc/np.float(Nfp*Npc))
    part2 = (sig_E/(E*(1-E)))
    
    sig_R = (1-R) * np.sqrt(part1**2 + part2**2)
    
    return sig_R,sig_E
    


def getSvnNum(fname):
    """
    Get the first line of file that contains the id string
    """
    fid=open(fname,'r')
    line=fid.readline()
    line=line.replace('$',' ')
    
    return line
    
def createScoreCard(opsdata,opsfgkData,invdata,invfgkData,injdata,injfgkData,Srange,pRadius,\
                        hzpercent,hztotal,hzpcs,climrange,figureTitle,outname):
    """
    Create 9 plot scorecard of Inj,INV and OPS.
    Inlude HZ estimate
    """

    fig=plt.figure(figsize=(8.5,11))
    gs=gridspec.GridSpec(7,2)
    plt.subplots_adjust(hspace=0.75, wspace=0.75)

    ax=plt.subplot(gs[0:2,0])
    plotMPGrid(opsdata,'All OPS PC Rate',climrange)
   
    ax=plt.subplot(gs[0:2,1])
    plotMPGrid(opsfgkData,'FGK Dwarf OPS PC Rate',climrange)
    
    ax=plt.subplot(gs[2:4,0])
    plotMPGrid(invdata,'All RelType PC Rate',(climrange[0],climrange[1]/2))
    
    ax=plt.subplot(gs[2:4,1])
    plotMPGrid(invfgkData,'FGK Dwarf RelType PC Rate',(climrange[0],climrange[1]/2))
    ax.plot([-2,1.1],[1.15,1.15],transform=ax.transAxes, clip_on=False,color='k')

    plt.subplot(gs[4:6,0])
    plotMPGrid(injdata,'All INJ PC Rate',(climrange[1]/2,climrange[1]))
    
    ax=plt.subplot(gs[4:6,1])
    plotMPGrid(injfgkData,'FGK Dwarf INJ PC Rate',(climrange[1]/2,climrange[1]))
    ax.plot([-2,1.1],[1.15,1.15],transform=ax.transAxes, clip_on=False,color='k')
    
       
    #---Habitable Zone --- #
    hzstr=np.array(["%.2f%%" % x for x in hzpercent.reshape(hzpercent.size)])
    hzstr = hzstr.reshape(hzpercent.shape)
    
    hztotstr=np.array(["%u tot" % x for x in hztotal.reshape(hztotal.size)])
    hztotstr = hztotstr.reshape(hztotal.shape)

    hzpcstr=np.array(["%u pcs" % x for x in hzpcs.reshape(hzpcs.size)])
    hzpcstr = hzpcstr.reshape(hzpcs.shape)
    
    tabdata=[hzstr[:,0], hztotstr[:,0], hzpcstr[:,0]]
    ax = fig.add_subplot(gs[6,0])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels=("Ops","Inv","Inj")
    
    ax.table(cellText=tabdata,colLabels=colLabels,loc='center')
    atitle="All HZ PCs S=%s Pr=%s" % (str(Srange), str(pRadius))
    plt.title(atitle,fontsize=11)
    
    tabdata=[hzstr[:,1], hztotstr[:,1], hzpcstr[:,1]]
    ax = fig.add_subplot(gs[6,1])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels=("Ops","Inv","Inj")

    ax.table(cellText=tabdata,colLabels=colLabels,loc='center', fontsize=12)
    atitle="FGK HZ PCs S=%s Pr=%s" % (str(Srange), str(pRadius))
    plt.title(atitle,fontsize=11)
    
    plt.annotate(figureTitle,xy=(0.20,0.94),xycoords='figure fraction',  fontsize=12)
    plt.savefig(outname)
    plt.close()

def createScoreCardScore(opsdata,opsfgkData,invdata,invfgkData,injdata,injfgkData,Srange,pRadius,\
                        hzpercent,hztotal,hzpcs,hzinvwant,hzopswant,climrange,figureTitle,outname,scores=[0.0,1.0]):
    """
    Create 9 plot scorecard of Inj,INV and OPS.
    Inlude HZ estimate
    In this case make a cut on Score s1 refers to the PC score and s2 the FP score.
    """

    fig=plt.figure(figsize=(8.5,11))
    gs=gridspec.GridSpec(7,2)
    plt.subplots_adjust(hspace=0.75, wspace=0.75)

    ax=plt.subplot(gs[0:2,0])
    plotMPGridS(opsdata,'All OPS PC Rate',climrange,s=scores)
   
    ax=plt.subplot(gs[0:2,1])
    plotMPGridS(opsfgkData,'FGK Dwarf OPS PC Rate',climrange,s=scores)
    
    ax=plt.subplot(gs[2:4,0])
    plotMPGridS(invdata,'All False PC Rate',(climrange[0],climrange[1]/2),s=scores)
    
    ax=plt.subplot(gs[2:4,1])
    plotMPGridS(invfgkData,'FGK Dwarf False PC Rate',(climrange[0],climrange[1]/2),s=scores)
    ax.plot([-2,1.1],[1.15,1.15],transform=ax.transAxes, clip_on=False,color='k')

    plt.subplot(gs[4:6,0])
    plotMPGridS(injdata,'All INJ PC Rate',(climrange[1]/2,climrange[1]),s=scores)
    
    ax=plt.subplot(gs[4:6,1])
    plotMPGridS(injfgkData,'FGK Dwarf INJ PC Rate',(climrange[1]/2,climrange[1]),s=scores)
    ax.plot([-2,1.1],[1.15,1.15],transform=ax.transAxes, clip_on=False,color='k')
    
       
    #---Habitable Zone --- #
    hzstr=np.array(["%.2f%%" % x for x in hzpercent.reshape(hzpercent.size)])
    hzstr = hzstr.reshape(hzpercent.shape)
    
    hztotstr=np.array(["%u tot" % x for x in hztotal.reshape(hztotal.size)])
    hztotstr = hztotstr.reshape(hztotal.shape)

    hzpcstr=np.array(["%u pcs" % x for x in hzpcs.reshape(hzpcs.size)])
    hzpcstr = hzpcstr.reshape(hzpcs.shape)
    
    tabdata=[hzstr[:,0], hztotstr[:,0], hzpcstr[:,0]]
    ax = fig.add_subplot(gs[6,0])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels=("Ops","Inv","Inj")
    
    ax.table(cellText=tabdata,colLabels=colLabels,loc='center')
    atitle="All HZ PCs S=%s Pr=%s" % (str(Srange), str(pRadius))
    plt.title(atitle,fontsize=11)
    
    tabdata=[hzstr[:,1], hztotstr[:,1], hzpcstr[:,1]]
    ax = fig.add_subplot(gs[6,1])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels=("Ops","Inv","Inj")

    ax.table(cellText=tabdata,colLabels=colLabels,loc='center', fontsize=12)
    atitle="FGK HZ PCs S=%s Pr=%s" % (str(Srange), str(pRadius))
    plt.title(atitle,fontsize=11)
    
    R,E = estimateReliability(invdata[hzinvwant],opsdata[hzopswant],s=scores)    
    hzRelStr="HZ (%4.2f, %4.2f): E= %5.2f%%  R=%5.2f%%" % (Srange[0],Srange[1],E*100,R*100)
    plt.annotate(hzRelStr,xy=(0.11,0.06),xycoords='figure fraction', fontsize=11)
    
    plt.annotate(figureTitle,xy=(0.20,0.94),xycoords='figure fraction',  fontsize=12)
    plt.savefig(outname)
    plt.close()


def createOnePageStats(opsdata,injdata,invdata,scores=[0.0,1.0],figureTitle='onepage',hz=(0.25,1.7)):
    """
    Create a page with reliability and completeness numbers for large bins in a table.
    """
    injpc=passes(injdata,s=scores)
    opspc=passes(opsdata,s=scores)
    invpc=passes(invdata,s=scores)    
    
    injwant1=injdata.period>0
    opswant1=opsdata.period>0
    invwant1=invdata.period>0
    
    injwant2=(injdata.period>=200) & (injdata.period<=500)
    opswant2=(opsdata.period>=200) & (opsdata.period<=500)
    invwant2=(invdata.period>=200) & (invdata.period<=500)
    
    injwant3=(injdata.mes<=10)
    opswant3=(opsdata.mes<=10)
    invwant3=(invdata.mes<=10)
    
    injwanthz=getHZarray(injdata,Srange=hz,pRadius=hz,s=scores)
    opswanthz=getHZarray(opsdata,Srange=hz,pRadius=hz,s=scores)
    invwanthz=getHZarray(invdata,Srange=hz,pRadius=hz,s=scores)
    
    #initialize arrays
    C=np.zeros(4)
    E=np.zeros(4)
    R=np.zeros(4)
    N=np.zeros(4)
    SNR=np.zeros(4)
    
    
    #Calculate completeness for each desired range.
    C[0]=len(injdata[injwant1 & injpc])/np.float(len(injwant1[injwant1]))
    C[1]=len(injdata[injwant2 & injpc])/np.float(len(injwant2[injwant2]))
    C[2]=len(injdata[injwant3 & injpc])/np.float(len(injwant3[injwant3]))
    C[3]=len(injdata[injwanthz & injpc])/np.float(len(injwanthz[injwanthz]))
    
    
    N[0]=len(opsdata[opswant1 & opspc])
    N[1]=len(opsdata[opswant2 & opspc])
    N[2]=len(opsdata[opswant3 & opspc])
    N[3]=len(opsdata[opswanthz & opspc])
    #Calculate Effectiveness
    #E[0]=len(invdata[invwant1 & invpc])/len(invwant1[invwant1])
    #E[1]=len(invdata[invwant2 & invpc])/len(invwant2[invwant2])
    #E[2]=len(invdata[invwant3 & invpc])/len(invwant3[invwant3])
    
    #Calculate Reliability
    R[0],E[0]=estimateReliability(invdata[invwant1],opsdata[opswant1],s=scores)
    R[1],E[1]=estimateReliability(invdata[invwant2],opsdata[opswant2],s=scores)
    R[2],E[2]=estimateReliability(invdata[invwant3],opsdata[opswant3],s=scores)
    R[3],E[3]=estimateReliability(invdata[invwanthz],opsdata[opswanthz],s=scores)
    
    #Calculate SNR
    SNR[0]=RelSNR(R[0],len(opsdata[opspc & opswant1]))
    SNR[1]=RelSNR(R[1],len(opsdata[opspc & opswant2]))
    SNR[2]=RelSNR(R[2],len(opsdata[opspc & opswant3]))
    SNR[3]=RelSNR(R[3],len(opsdata[opspc & opswanthz]))

    Cstr=np.array(["%.4f%%" % x for x in 100*C.reshape(C.size)])
    Cstr=Cstr.reshape(C.shape)
    
    Nstr=np.array(["%.1f" % x for x in N.reshape(N.size)])
    Nstr=Nstr.reshape(N.shape)

    Estr=np.array(["%.4f%%" % x for x in 100*E.reshape(E.size)])
    Estr=Estr.reshape(E.shape)

    Rstr=np.array(["%.4f" % x for x in R.reshape(R.size)])
    Rstr=Rstr.reshape(R.shape)
    
    SNRstr=np.array(["%.4f" % x for x in SNR.reshape(SNR.size)])
    SNRstr=SNRstr.reshape(SNR.shape) 
    
    tabdata=[Cstr[:], Estr[:], Rstr[:], SNRstr[:], Nstr[:]]
    
    plt.figure(figsize=(8.5,11))    
    ax = plt.subplot(311)
    the_table=ax.table(cellText=tabdata,colLabels=['All','P:200-500','MES<10','HZ'],rowLabels=['C','E','R','SNR','Npc'],loc='center', fontsize=12)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    
    #Change height of the cells
    table_props = the_table.properties()
    table_cells = table_props['child_artists']
    for cell in table_cells: 
        cell.set_height(0.1)
    
    plt.annotate(figureTitle,xy=(0.20,0.94),xycoords='figure fraction',  fontsize=11) 
    
def marginalizedRunPlots(opsdata,invdata,scores=[0.0,1.0],figureTitle='onepage',hz=(0.25,1.7),Effect=False):
    plt.figure(figsize=(8.5,11))
    plt.subplot(311)
    plotRunningRel(opsdata[opsdata.mes<20],invdata[invdata.mes<25],binbeg=np.arange(20,450,10),\
            metric='period',scores=scores,binwidth=60,label='mes<20, bin=60d',Effect=Effect)
    
    plt.subplot(312)
    plotRunningRel(opsdata[opsdata.period>50],invdata[invdata.period>70],binbeg=np.arange(7.0,20,.5),\
            metric='mes',scores=scores,binwidth=1.5,label='period>70d, bin=1.5',Effect=Effect)
    
    plt.subplot(313)
    plotRunningRel(opsdata[opsdata.mes<20],invdata[invdata.mes<25],binbeg=np.arange(2.5,15.5,1),\
            metric='ngood',scores=scores,binwidth=0.95,label='mes<20, bin=1',Effect=Effect)
    
    plt.annotate(figureTitle,xy=(0.20,0.94),xycoords='figure fraction',  fontsize=11)    
    
def getHZarray(data,Srange=(0.75,1.25),pRadius=(0.75,1.25),s=(0.0,1.0)):
    """
    Return the array of hz objects
    Use the input ranges.
    """
    
    hzwant=(data['srad'] >=Srange[0]) & (data['srad']<=Srange[1]) & \
           (data['rplanet'] >=pRadius[0]) & (data['rplanet']<=pRadius[1])

    return hzwant


def RelSNR(R,npc):
    """
    """
    snr=(R * npc**(0.5)) / (2-R)**(0.5)   
    
    return snr

def plotRunningRel(opsdata,invdata,binbeg=[0,100,200,300],metric='period',scores=[0.0,1.0],binwidth=50,label='',Effect=False):
    """
    Reliability and Completeness as a running bin.
    """
    R=np.zeros(len(binbeg))
    E=np.zeros(len(binbeg))
    Rs=np.zeros(len(binbeg))
    Es=np.zeros(len(binbeg))
    for (i,b) in enumerate(binbeg):
        opswant=(opsdata[metric] >= b) & (opsdata[metric] <= b+binwidth)
        invwant=(invdata[metric] >= b) & (invdata[metric] <= b+binwidth)
        R[i],E[i]=estimateReliability(invdata[invwant],opsdata[opswant],s=[0.0,1.0])
        Rs[i],Es[i]=estimateReliability(invdata[invwant],opsdata[opswant],s=scores)
    
    if Effect:
        R=E
        Rs=Es
        ylab='Effectiveness'
    else:
        ylab='Reliability'

    plt.plot(binbeg+(binwidth/2.0),R,'.b-',label=label)
    plt.plot(binbeg+(binwidth/2.0),Rs,'.r--',label='score=%s'% str(scores))
    plt.ylabel(ylab)
    plt.xlabel(metric)
    plt.legend(loc='lower right',fontsize='small')
    
    

def fineGridRelCompleteness(opsdata,invdata,injdata,opsFgkWant,invFgkWant,injFgkWant,climrange,outname,figTitle,s=[0.0,1.0], \
        xbins=[0,100,200,300,400,500],ybins=[7,9,11,13,15,20] ):
    plt.figure(figsize=(18,13))
    plt.subplots_adjust(hspace=0.15, wspace=0.15)    
    plt.subplot(231)
    #xbins=[0,50,100,150,200,250,350,400,500]
    #ybins=[7,8,9,10,12,15,20]
    #xbins=[0,100,200,300,400,500]
    #ybins=[7,9,11,13,15,20]
    #Plot reliability
    prv.plotOnlyReliability(invdata,opsdata,xBins=xbins,yBins=ybins,atitle="All Reliability\n#PCs/#FPs",s=s)
    
    plt.subplot(232)
    plotMPGrid2(invdata,'InEffectiveness\nAll INV PC Rate',(climrange[0],climrange[1]/4),xBins=xbins,yBins=ybins,s=s)
    
    plt.subplot(233)
    #Plot completeness
    plotMPGrid2(injdata,'Completeness\nAll INJ PC Rate',(climrange[1]/2,climrange[1]),xBins=xbins,yBins=ybins,s=s)
    
    plt.subplot(234)
    #FGK reliability
    prv.plotOnlyReliability(invdata[invFgkWant],opsdata[opsFgkWant],xBins=xbins,yBins=ybins,atitle="FGK Reliability\n#PCs/#FPs",s=s)

    plt.subplot(235)
    plotMPGrid2(invdata[invFgkWant],'InEffectiveness\nFGK INV PC Rate',(climrange[0],climrange[1]/4),xBins=xbins,yBins=ybins,s=s)

    plt.subplot(236)    
    #FGK Completeness
    plotMPGrid2(injdata[injFgkWant],'FGK Completeness\nINJ PC Rate',(climrange[1]/2,climrange[1]),xBins=xbins,yBins=ybins,s=s)
    
    plt.annotate(figTitle,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10) 
    
    plt.savefig(outname)
    plt.close()

def completenessAndReliability(opsdata,invdata,injdata,opsFgkWant,invFgkWant,injFgkWant,climrange,outname,figTitle,s=[0.0,1.0]):
    plt.figure(figsize=(8.5,6))
    plt.subplots_adjust(hspace=0.25, wspace=0.45)    
    plt.subplot(221)
    #Plot reliability
    prv.plotOnlyReliability(invdata,opsdata,atitle="Reliability, All Stars \n#PCs/#FPs",s=s,xlabel='Period (days)',ylabel='MES')
    
    plt.subplot(222)
    #Plot completeness
    plotMPGridS(injdata,'Completeness, All Stars\n#INJ-TCE',(climrange[1]/2,climrange[1]),s=s)
    
    plt.subplot(223)
    #FGK reliability
    prv.plotOnlyReliability(invdata[invFgkWant],opsdata[opsFgkWant],atitle="Reliability, FGK Dwarfs\n#PCs/#FPs",s=s,xlabel='Period (days)',ylabel='MES')

    plt.subplot(224)    
    #FGK Completeness
    plotMPGridS(injdata[injFgkWant],'Completeness, FGK Dwarfs\n#INJ-TCE',(climrange[1]/2,climrange[1]),s=s)
    
    plt.annotate(figTitle,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10) 
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()
    
def plotPeriodRadius(data,outname, figTitle,label="",colorkey="score",colorlim=[0,1],s=[0.0,1.0]):
    """
    Plot the period vs radius plot, colored by the specified key
    passes according to disposition and score
    FPs that end up passing are given a different color.
    """
    
    pcs=passes(data,s) 
    fppass=(data['disp']=='FP') & (data['score']>s[1])
    plotsortdata=data.sort_values(by=colorkey,ascending=False)
    period=plotsortdata.loc[pcs,'period']
    color=plotsortdata.loc[pcs,colorkey]
    #print colorkey,len(color), len(pcs[pcs]),color[0],color[2]
    prad=plotsortdata.loc[pcs,'rplanet']
    
    fpperiod=plotsortdata.loc[fppass,'period']
    #fpcolor=plotsortdata.loc[fppass,colorkey]
    fpprad=plotsortdata.loc[fppass,'rplanet']

    plt.figure()
    plt.scatter(np.log10(period),np.log10(prad),c=color,marker='o',linewidth=0)

    
    plt.set_cmap('plasma')
    cbar=plt.colorbar()
    cbar.set_label(colorkey)
    lab=np.array([.4,1,4,10,40])
    locat=np.log10(lab)
    plt.yticks(locat,lab.astype(str))
    plt.ylabel(r"Planet Radius (Earth radii)")

    lab=np.array([0.3,1,3,10,30,100,300])
    locat=np.log10(lab)
    plt.xticks(locat,lab.astype(str))
    plt.xlabel('Period (days)')
    
    plt.plot(np.log10(fpperiod),np.log10(fpprad),'g.',ms=3)
    plt.xlim(-0.6,3.0)
    plt.ylim(-.5,1.9)
    plt.clim(colorlim)
    plt.title('%s Planet Candidates\n low value plotted on top' % (label))
    plt.annotate(figTitle,xy=(0.55,0.001),xycoords='figure fraction',fontsize=6)
    plt.savefig(outname)
    plt.close()

def plotInsolationRadius(data,outname, figTitle,label="",colorkey="score",colorlim=[0,1],s=[0.0,1.0]):
    """
    Plot the period vs radius plot, colored by the specified key
    passes according to disposition and score
    FPs that end up passing are given a different color.
    """
    
    pcs=passes(data,s) 
    fppass=(data['disp']=='FP') & (data['score']>s[1])
    plotsortdata=data.sort_values(by=colorkey,ascending=False)
    insol=plotsortdata.loc[pcs,'srad']
    color=plotsortdata.loc[pcs,colorkey]
    print colorkey,len(color), len(pcs[pcs]),color[0],color[2]
    prad=plotsortdata.loc[pcs,'rplanet']
    
    fpinsol=plotsortdata.loc[fppass,'srad']
    #fpcolor=plotsortdata.loc[fppass,colorkey]
    fpprad=plotsortdata.loc[fppass,'rplanet']

    plt.figure()
    plt.scatter(np.log10(insol),np.log10(prad),c=color,marker='o',linewidth=0)

    
    plt.set_cmap('viridis')
    cbar=plt.colorbar()
    cbar.set_label(colorkey)
    lab=np.array([.4,1,4,10,40])
    locat=np.log10(lab)
    plt.yticks(locat,lab.astype(str))
    plt.ylabel(r"Planet Radius (Earth radii)")

    lab=np.array([.3,1,3,10,100,1000,10000])
    locat=np.log10(lab)
    plt.xticks(locat,lab.astype(str))
    plt.xlabel('Insolation Flux (Earth)')
    
    plt.plot(np.log10(fpinsol),np.log10(fpprad),'r.',ms=2)
    plt.plot(0,0,'y*',ms=8)
    plt.xlim(-1,5.5)
    plt.ylim(-.5,1.9)
    plt.clim(colorlim)
    plt.gca().invert_xaxis()
    plt.title('%s Planet Candidates\n low value plotted on top' % (label))
    plt.annotate(figTitle,xy=(0.55,0.001),xycoords='figure fraction',fontsize=6)
    plt.savefig(outname)
    plt.close()
    
def ebEfficiency(data,outname,figTitle):
    """
    Efficiency of the EBs
    """
    period=data['period']
    mes=data['mes']
    pcIdx=prv.passes(data,s=[0.0,1.0])
    plt.figure(figsize=(10,5))
    plt.subplot(1,2,1)
    prv.plotGrid(period, mes, pcIdx, drange=(0,100),xBins = [0, 10, 200, 500],yBins = [7, 10, 20, 2000])
    plt.title('INJ-EBs Fraction called PC')
    plt.subplot(1,2,2)
    pcIdx=(data['S']==1) & (data['N']==0)
    prv.plotGrid(period, mes, pcIdx, drange=(0,100),xBins = [0, 10, 200, 500],yBins = [7, 10, 20, 2000])
    plt.title('INJ-EBs Fraction with S Flag\n (and not N)')
    plt.annotate(figTitle,xy=(0.55,0.001),xycoords='figure fraction',fontsize=6)
    plt.savefig(outname)
    plt.close()

def plotTceDist(ops,inv,bins=200,metrics=('period','mes'),labels=('OPS','INV')):
    """
    Compare the histograms along two metrics for the two input metrics
    """
    plt.figure(figsize=(8,10))
    plt.subplot(211)
    plt.hist(np.log10(ops[metrics[0]]),bins=bins,color='red',histtype='step',label=labels[0])
    plt.hist(np.log10(inv[metrics[0]]),bins=bins,color='blue',histtype='step',label=labels[1])
    plt.xlabel(metrics[0])
    plt.legend(loc='upper center')
    
    plt.subplot(212)
    plt.hist(np.log10(ops[metrics[1]]),bins=bins,color='red',histtype='step',label=labels[0])
    plt.hist(np.log10(inv[metrics[1]]),bins=bins,color='blue',histtype='step',label=labels[1])
    plt.xlabel(metrics[1])
    plt.legend(loc='upper center')

def plotTceDist3(ops,inv,inj,bins=200,metric=('period'),labels=('OPS','INV','INJ'),ylim=(0,1200)):
    """
    Plot three distributions in a columns
    """
    
    plt.figure(figsize=(7,10))
    plt.subplot(311)
    plt.hist(np.log10(ops[metric]),bins=bins,color='black',histtype='step',label=labels[0])
    plt.xlabel(metric)
    plt.ylim(ylim)
    plt.legend()
    plt.legend(loc='upper left')
    plt.subplot(312)
    plt.hist(np.log10(inv[metric]),bins=bins,color='blue',histtype='step',label=labels[1])
    plt.xlabel(metric)
    plt.ylim(ylim)
    plt.legend()
    plt.legend(loc='upper left')
    plt.subplot(313)
    plt.hist(np.log10(inj[metric]),bins=bins,color='red',histtype='step',label=labels[2])
    plt.xlabel(metric)
    plt.ylim(ylim)
    plt.legend()
    plt.legend(loc='upper left')

def plotPCs(ops,inv,metricx,metricy,score=[0.0,1.0],xlim=None,ylim=None,figTitle=None):
    """
    plot metrixy vs metricx for those that pass according to the score cut.
    """
    
    opspc=passes(ops,s=score)
    invpc=passes(inv,s=score)
    
    plt.figure(figsize=(9,9))
    
    plt.plot(ops.loc[opspc][metricx],ops.loc[opspc][metricy],'ro',alpha=0.6,markeredgewidth=0,markersize=14)
    plt.plot(inv.loc[invpc][metricx],inv.loc[invpc][metricy],'bs',alpha=0.6,markeredgewidth=0,markersize=13)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.xlabel(metricx)
    plt.ylabel(metricy)
    plt.title('PCs for OPS (red) and INV (blue)')
    plt.annotate(figTitle,xy=(0.1,0.91),xycoords='figure fraction',fontsize=10)
    