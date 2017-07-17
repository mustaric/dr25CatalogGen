# -*- coding: utf-8 -*-
"""
Created on Thu May 19 13:54:18 2016

@author: smullall

Code to look at individual metrics.
Contains code to look at the major flags provided by Jeff's Robovetter.

"""

import rvIO as io
import numpy as np

#from multiprocessing import Pool

def uniqueFlags(dataf,fieldname='flags'):
    """
    Given a data frame with the column flags
    Return all the unique reasons objects failed.
    And count the number of times that flag is used.
    """
    
    all=[]
    
    fl=dataf[fieldname]
    
    for f in fl:
        if type(f)==float:
            f=""
            
        all.extend(f.split('---'))
    
    [uniqueFlags,flagCount]=np.unique(all,return_counts=True)
   
    return uniqueFlags,flagCount

def flagCount(dataf,flags):
    """
    Return array of the number that are true for each flag
    """
    numTrue=np.zeros(len(flags))
    
    for i,f in enumerate(flags):
        try:
            numTrue[i]=sum(dataf[f])
        except KeyError:
            numTrue[i]=1

    print len(flags),len(numTrue)
    return numTrue

def groupFlagsPost(dataf):
    """
    Create new columns of true/false based on desired logic
    """
    
    dataf['g_EPHEM_MATCH']=np.array(map(lambda x: x.find("PARENT_IS")>=0, dataf['flags']),dtype=bool)
    dataf['g_RUBBLE'] = np.array(map(lambda x: x.find("RUBBLE")>=0, dataf['flags']),dtype=bool)
    dataf['g_CHASES'] = np.array(map(lambda x: x.find("CHASES")>=0, dataf['flags']),dtype=bool)       
    dataf['g_MARSHALL'] = np.array(map(lambda x: x.find("MARSHALL")>=0, dataf['flags']),dtype=bool)
    dataf['g_ROCKY'] = np.array(map(lambda x: x.find("ROCKY")>=0, dataf['flags']),dtype=bool)
    dataf['g_SKYE'] = np.array(map(lambda x: x.find("SKYE")>=0, dataf['flags']),dtype=bool)
    dataf['g_ODD_EVEN'] = np.array(map(lambda x: x.find("ODD_EVEN")>=0, dataf['flags']),dtype=bool)
    dataf['g_ITRANS_ALL'] = np.array(map(lambda x: x.find("ITRANS")>=0, dataf['flags']),dtype=bool)
    try:  
        dataf['g_LPP_HIGH']= dataf['LPP_ALT_TOO_HIGH'] | dataf['LPP_DV_TOO_HIGH']
    except:
        pass
    try:
        dataf['g_NO_CENTROID']=dataf['CROWDED_DIFF'] | dataf['EYEBALL'] | dataf['CENTROID_SIGNIF_UNCERTAIN'] | dataf['TOO_FEW_CENTROIDS'] | dataf['TOO_FEW_QUARTERS'] | dataf['KIC_OFFSET']
    except:
        pass    
    try:
        dataf['g_SEASONAL']=dataf['SEASONAL_DEPTH_DIFFS_IN_ALT'] | dataf['SEASONAL_DEPTH_DIFFS_IN_DV']
    except:
        pass    
    try:
        dataf['g_SIG_PRI_OVER_FRED']=dataf['ALT_SIG_PRI_OVER_FRED_TOO_LOW'] | dataf['DV_SIG_PRI_OVER_FRED_TOO_LOW']
    except:
        pass    
    try:    
        dataf['g_TAT']=dataf['MS_ALT_TAT_FAIL'] | dataf['MS_DV_TAT_FAIL']
    except:
        pass    
    try:    
        dataf['g_MS_SHAPE']=dataf['MS_DV_SHAPE_FAIL'] | dataf['MS_ALT_SHAPE_FAIL']
    except:
        pass
    try:
        dataf['g_DMM']=dataf['MS_DV_DMM_FAIL'] | dataf['MS_ALT_DMM_FAIL']
    except:
        pass
    try:    
        dataf['g_ITRANS_ALONE']=dataf['g_ITRANS_ALL'] & ~dataf['g_LPP_HIGH'] & ~dataf['g_SIG_PRI_OVER_FRED'] & ~dataf['g_DMM'] & ~dataf['g_TAT'] & ~dataf['g_ODD_EVEN']
    except:
        pass
    
    newflags=['g_EPHEM_MATCH','g_RUBBLE','g_CHASES','g_MARSHALL','g_ROCKY','g_ODD_EVEN','g_ITRANS_ALL','g_LPP_HIGH',
              'g_NO_CENTROID','g_SEASONAL','g_SIG_PRI_OVER_FRED','g_TAT','g_MS_SHAPE','g_DMM','g_ITRANS_ALONE','g_SKYE']
    
    return dataf,newflags
    
#import pdb
def accountFlag(rvdata,uniFlags):
    """
    Create True false columns out of the minor flags
    """ 

    for uf in uniFlags:
        if (uf.find("PARENT_IS")<0):
            rvdata[uf] = np.array(map(lambda x: x.find(uf)>=0, rvdata['flags']), dtype=bool)

    return rvdata    
    
def reportFlags(datafile,tcefile,ops=True):
    """Run uniqueFlags and create an output file
    """
    if ops:
        fedfile='/soc/nfs/so-nfs/DR25/other/koimatch_DR25_03032016-DMedit.txt'
        cumfile='/home/smullall/Kepler/cumulativeTable/q1q17/q1q17-status.csv'
        rvdata=io.createAllResults(datafile,tcefile,fedfile,cumfile)
    else:
        rvdata=io.createRVResults(datafile,tcefile)
    rvdata['flags']=rvdata['flags'].fillna('ABSENT_FLAG')
    #pdb.set_trace()
    #rvdata=groupFlagsPre(rvdata)
    #want=(rvdata['period']<percut[1]) & (rvdata['period']>percut[0]) & (rvdata['mes']<mescut[1]) & (rvdata['mes']>mescut[0])
    [uniFlags, fc]=uniqueFlags(rvdata)
    rvdata=accountFlag(rvdata,uniFlags)
    rvdata,newflags = groupFlagsPost(rvdata)
    extFlags=list(uniFlags)    
    extFlags.extend(newflags)
    arrFlags=np.array(extFlags)
    flagCnt=flagCount(rvdata,extFlags)
    data=[arrFlags,flagCnt]
    
    return (rvdata,data)
#%%
def reportFlagDf(df):
    """
    Given a dataframe, report the unique flags and their frequency
    """
    [uniFlags, fc]=uniqueFlags(df)
    flagCnt=flagCount(df,uniFlags)
    data=[uniFlags,flagCnt]  
    
    return (df,data)
#%%    
import matplotlib.pyplot as plt
plt.switch_backend('Agg') 
def plotUniqueFlags(data,x,y,flags,markSize=1):
    """
    For a parameter space of x vs. y, plot each flag name in a subplot. Good for seeing why thigns failed.
    """
    L=len(flags)+1
    rows=np.ceil(L/3.0)
    for j in np.arange(0,int(rows),1):
        for i in np.arange(0,3,1):
            
            it=int(j*3+i)

            if (it) < len(flags):

                f=flags[it]
                plt.subplot(rows,3,it+1)
                
                try:
                    want=data[f]
                    n=len(want[want])                    
                    plt.plot(np.log10(data[want][x]),np.log10(data[want][y]),'k.',ms=markSize)
                    plt.title("%s %u" % (f,n),fontsize=10)
                    plt.xlim((-.3,2.8))
                    plt.ylim((.7,3))
                except:
                    plt.title("%s %u" % (f,0),fontsize=10)
                    pass
    
    plt.subplot(rows,3,int(rows*3))
    plt.plot(np.log10(data[x]),np.log10(data[y]),'r.',ms=2)
    plt.title('All Population')
    plt.xlim((-.3,2.8))
    plt.ylim((.7,3))
#%%
import datetime as dt
def printFlags(df,uniFlags,filename,svnId):
    """
    Print to a text file the TCE name, N,S,C,E and then all the flags associated with that target.
    """
    date=dt.date.today()
    header="#minorFlagsReport\n#Created %s\n#Revision %s\n#\n" % (date,svnId)
    fid=open(filename,'w')
    fid.write(header)
    
    columns=['N','S','C','E','disp']
    columns.extend(uniFlags)
    
    #Recast all the T/F to 1 and 0
    intdf=df*1
    
    astring=intdf.to_csv(sep=',',columns=columns,header=True,index=True,float_format='%1',na_rep='0')
    fid.write(astring)
    fid.close()
    
    

#Change the next line's revision number.
id='r62353'
injfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-INJ-PlanetOn.txt'
inj2file='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-INJ-PlanetOff.txt'
opsfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-OPS.txt'
invfile='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-INV.txt'
ss1file='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/RoboVetterOut-SS1.txt'

tcefile='/soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt'
injtfile='/soc/nfs/so-nfs/DR25/INJ/DATA/TCEs.txt'
invtfile='/soc/nfs/so-nfs/DR25/INV/DATA/TCEs.txt'
ss1tfile='/soc/nfs/so-nfs/DR25/SS1/DATA/TCEs.txt'

outinj='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/uniqueFlags-Inj.txt'
outinj2='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/uniqueFlags-Inj.txt'
outinv='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/uniqueFlags-Inv.txt'
outops='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/uniqueFlags-Ops.txt'
outss1='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/uniqueFlags-Ss1.txt'

(rvdata,dataops)=reportFlags(opsfile,tcefile)
np.savetxt(outops,np.transpose(dataops),fmt=['%s','\t\t%s'])
uniflags=dataops[0]
flagswant=np.array(map(lambda x: x.find("PARENT_IS")<0,uniflags),dtype=bool)


(injdata,datainj)=reportFlags(injfile,injtfile,ops=False)
np.savetxt(outinj,np.transpose(datainj),fmt=['%s','\t\t%s'])

(inj2data,datainj2)=reportFlags(inj2file,injtfile,ops=False)
np.savetxt(outinj,np.transpose(datainj2),fmt=['%s','\t\t%s'])

(invdata,datainv)=reportFlags(invfile,invtfile,ops=False)
np.savetxt(outinv,np.transpose(datainv),fmt=['%s','\t\t%s'])

(ss1data,datass1)=reportFlags(ss1file,ss1tfile,ops=False)
np.savetxt(outss1,np.transpose(datass1),fmt=['%s','\t\t%s'])
#%%
flags=[ u'g_LPP_HIGH', u'g_EPHEM_MATCH',u'GHOST_HALO_TO_CORE_TOO_HIGH', u'g_ITRANS_ALL',
       u'g_RUBBLE',u'g_CHASES',u'g_MARSHALL',u'g_SKYE', u'g_ITRANS_ALONE',
       u'g_ROCKY',u'g_DMM',u'g_SIG_PRI_OVER_FRED',u'SWEET_NTL',u'SWEET_EB',
       u'FERGAL_SHAPE_FAIL',u'CLEAR_APO',u'g_ODD_EVEN']

plt.figure(figsize=(15,15))
plt.subplots_adjust(hspace=0.45, wspace=0.3) 
plotUniqueFlags(rvdata,'period','mes',flags)
plt.annotate('Reasons TCEs are Labeled FP  '+id,xy=(0.2,0.94),xycoords='figure fraction')
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/' + id+'/UniqueFlagOpsFPs-'+id+'.png')


plt.figure(figsize=(15,15))
plt.subplots_adjust(hspace=0.45, wspace=0.3)
rvdata['iskoi']=rvdata['koi'].notnull()
koisnomore=(rvdata['iskoi']) & (rvdata['N']==1)
plotUniqueFlags(rvdata[koisnomore],'period','mes',flags,markSize=2)
plt.annotate('Reason KOIs are now listed as not-transit-like'+ id,xy=(0.2,0.94),xycoords='figure fraction')
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/UniqueFlagKOI-NFPs-'+id+'.png')
#%
plt.figure(figsize=(15,15))
plt.subplots_adjust(hspace=0.45, wspace=0.3) 
plotUniqueFlags(injdata,'period','mes',flags)
plt.annotate('Reason on-Target Injections are FPs ' + id,xy=(0.2,0.94),xycoords='figure fraction')
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/UniqueFlagInjFPs-'+id+'.png')

plt.figure(figsize=(15,15))
plt.subplots_adjust(hspace=0.45, wspace=0.3) 
plotUniqueFlags(invdata,'period','mes',flags)
plt.annotate('Reason Inversions are labeled FPs ' + id,xy=(0.2,0.94),xycoords='figure fraction')
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/UniqueFlagInvFPs-'+id+'.png')

plt.figure(figsize=(15,15))
plt.subplots_adjust(hspace=0.45, wspace=0.3) 
plotUniqueFlags(ss1data,'period','mes',flags)
plt.annotate('Reason SS1 are labeled FPs ' + id,xy=(0.2,0.94),xycoords='figure fraction')
plt.savefig('/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/UniqueFlagSs1FPs-'+id+'.png')
#%%
svnId=id
filename='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsOPS-' + id + '.csv'
printFlags(rvdata,uniflags[flagswant],filename,svnId)
filename='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsINV-'+ id + '.csv'
printFlags(invdata,uniflags[flagswant],filename,svnId)
filename='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsSS1-'+ id + '.csv'
printFlags(ss1data,uniflags[flagswant],filename,svnId)
filename='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsINJ-PlanetOn-'+ id + '.csv'
printFlags(injdata,uniflags[flagswant],filename,svnId)
filename='/home/smullall/Kepler/RoboVetter/DR25/stats/Versions/'+id+'/IndividualFlagsINJ-PlanetOff-'+ id + '.csv'
printFlags(inj2data,uniflags[flagswant],filename,svnId)

#if __name__ == "__main__":
#    main()