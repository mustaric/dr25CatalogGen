# -*- coding: utf-8 -*-
"""
Created on Fri May  6 16:41:52 2016

Tools for dealing with DR25

@author: fmullall
"""

__version__ = "$Id: dr25.py 2358 2016-06-24 19:38:43Z fmullall $"
__URL__ = "$URL: svn+ssh://flux.amn.nasa.gov/home/fmullall/svn/kepler/py/dr25.py $"

from scipy.io import matlab
import pandas as pd
import numpy as np
import tools
import re
import os

#Define some paths
opsDvPath='/soc/nfs/production-nfs4/mq-q1-q17/pipeline_results/dv-v4'
opsTpsPath='/soc/nfs/production-nfs3/mq-q1-q17/pipeline_results/tps-v4'
opsRbvFile='/soc/nfs/so-nfs/DR25/OPS/RoboVet/RoboVetterOut-OPS.txt'
opsRbv24File='/soc/nfs/so-nfs/DR25/OPS/RoboVet/RoboVetterOut-OPS-vDR24.txt'


injDvPath1='/soc/nfs/test-nfs3/ksop-2540-dv-batch1'
injDvPath2='/soc/nfs/test-nfs3/ksop-2540-dv-batch2'
injDvPath3='/soc/nfs/test-nfs4/ksop-2540-dv-batch3'
injDvPath4='/soc/nfs/test-nfs1/ksop-2540-dv-batch4'
injTpsPath='/soc/nfs/test-nfs1/ksop-2539-tps-q1-q17-transit-injection-reRun'

injRbvPath='/net/spmdb-fs/so-nfs/DR25/INJ/RoboVet/'
injRbv24File='/net/spmdb-fs/so-nfs/DR25/INJ/RoboVet/RoboVetterOut-INJ-All-vDR24.txt'
injRbvFile='/net/spmdb-fs/so-nfs/DR25/INJ/RoboVet/RoboVetterOut-INJ-All.txt'
injFedFile='/net/spmdb-fs/so-nfs/DR25/INJ/injmatch_DR25_03182016.txt'


invDvPath='/soc/nfs/test-nfs4/ksop-2543-transit-inversion-dv-rerun'
invTpsPath='/soc/nfs/test-nfs5/ksop-2542-transit-inversion-tps-rerun'
invRbvFile='/net/spmdb-fs/so-nfs/DR25/INV/RoboVet/RoboVetterOut-INV.txt'
invRbv24File='/net/spmdb-fs/so-nfs/DR25/INV/RoboVet/RoboVetterOut-INV-vDR24.txt'
invFedFile='/net/spmdb-fs/so-nfs/DR25/INV/koimatch_DR25INV_03172016.txt'

def searchColumns(df, key):
    """Returns a list of column names in a dataframe that match key.

    Example:
    -------------
    searchColumns(df, 'Radius')
    allTransitFits_planetRadius_value
    allTransitFits_planetRadius_uncertainty
    ....
    """

    f = lambda x: re.search(key, x, re.IGNORECASE)
    wh = np.where(map(f, df.columns))
    return df.columns[wh]



def loadFederation(fedFile):
    """
    Return a list of TCEs that are federated according to input file

    Inputs:
    ----------
    fedFile  (string) Path to file containing federation.

    Returns:
    ------------
    A DataFrame of the koi names of the federated TCEs. The index is the
    tceId


    Notes:
    -------
    Chris' standard federation format is understood by this function.
    In general, function expects the following quanities in the following
    columns
    0: koi number (eg 5.01)
    1: kepid
    2: tce
    -1: 1 if koi federated, zero otherwise. -1 means the last
       column of the file.

    """
    data = np.loadtxt(fedFile)
    idx = data[:,-1] == 1
    data = data[idx]

    tceId = tools.createTceId(data[:,1], data[:,2]).astype(np.int64)
    koiName = map(tools.fixKoiString, data[:,0].astype(str))
    return pd.DataFrame(koiName, index=tceId, columns=["koiName"])


def loadRobovetter(rbvFile):
    columns="""tceStr Score Disposition nFlag sFlag cFlag eFlag""".split()

    rbv = pd.read_csv(rbvFile, sep=' ', usecols=range(len(columns)),
                      names=columns, comment='#')

    tmp = rbv['tceStr'].as_matrix()
    kepid = np.array( map(lambda x: int(x[:9]), tmp) )
    tce = np.array( map(lambda x: int(x[-2:]), tmp) )
    tceId = tools.createTceId(kepid, tce)
    rbv.index = tceId

    return rbv


def loadDVForInjection(doFederation=True):
    """
    Loading DV params for injection.

    Loading Dv for injection is harder because it's split over
    4 files. This function implements this frequent, and slightly
    annoying task

    Optional Inputs:
    ----------------
    federate
        Cull to only TCEs that federate (i.e were injected
        and recovered)

    Returns:
    ---------------
    A dataframe.
    """

    dvPath = [injDvPath1, injDvPath2, injDvPath3, injDvPath4]

    data = []
    for p in dvPath:
        f = os.path.join(p, "dvOutputMatrix.mat")
        dv = loadDvAsDataFrame(f)
        data.append(dv)

    allData = pd.concat(data)

    if doFederation:
        fed = loadFederation(injFedFile)
        allData = federate(allData, fed)

    return allData



def federate(data, fed):
    """
    Federate data

    Return only the rows of data that are included in fed.

    Matching is done on index. If fed is created by loadFederation()
    the index of data should be tceId. If data is loaded by, eg.
    loadDvAsDataFrame, everything will work smoothly.

    Inputs:
    ---------------
    data
        (DataFrame) Data to federate
    fed
        (DataFrame) Dataframe as returned by loadFederation()

    Returns:
    -------------
    A dataframe.
    """

    return pd.merge(fed, data, left_index=True, right_index=True)


def antiFederate(data, fed):
    """
    Return the rows of data that are **not** in fed

    See federate for more details
    """

    idx = data.index.difference( fed.index)
    return data.loc[idx, :]


def tceStrToId(tceStr):
    """Convert Jeff's TCE string to a tceId float

    Jeff uses a string as the unique identifier of a TCE, where as I
    use an integer. This function converts from the string form to the
    integer form.

    Inputs:
    ------------
    tceStr  (string) A string consisting of "%s-%s" storing the kepid and
            the planet number. Lists (or arrays, or Series) of tceStrs can
            also be passed

    Returns:
    ------------
    tceId. An integer corresponding to 100*kepid + planNum
    """

    isSingle = False
    if isinstance(tceStr, str):
        isSingle = True
        tceStr = [tceStr]

    kepid   = np.array( map(lambda x: float(x[:-3]), tceStr) )
    planNum = np.array( map(lambda x: float(x[-2:]), tceStr) )
    tceId = tools.createTceId(kepid, planNum).astype(np.int32)

    if isSingle:
        tceStr = tceStr[0]

    return tceId


def tceIdToStr(tceId, sep="-"):
    """Convert a tceId (which is a number) to a tce string (with a defined format).

    A tce string is of the format "%08i-%02i" %(kepid, planetNumber)

    Inputs:
    ----------
    tceId
        (int, float, or iterable) If an int or float, this is the tceId to
        be converted. If it's an iterable (e.g a list) it's a list of tceIds.

    Optional Inputs:
    ----------------
    sep
        (string) String separator of kepid and planet number.

    Returns:
    -----------
    A string, or a list of strings. If the input is an iterable, the
    return value is always a list (not an array, or a tuple)
    """

    isSingle = False
    if not hasattr(tceId, "__len__"):
        isSingle = True
        tceId = np.array([tceId])

    kepid   = np.array( map(lambda x: int(np.floor(x/100.)), tceId) )
    planNum = np.array( map(lambda x, k: int(x - 100*k), tceId, kepid) )

    tceStr = map(lambda k,p: "%09i%s%02i" %(k, sep, p), kepid, planNum)

    if isSingle:
        tceStr = tceStr[0]

    return tceStr


def getRbvVersion(fn):
    """Return the version string in a robovetter file

    For the DR24 versions of files, it just returns "DR24"
    """
    fp = open(fn)

    ver = fp.readline()
    if ver.find("$Id") < 0:
        ver = "DR24"
    else:
        ver = ver[3:-1]
    return ver


def loadDvAsDataFrame(fn):
    mat = matlab.loadmat(fn, struct_as_record=False, squeeze_me=True)

    pars = mat['dvOutputMatrix']
    cols = mat['dvOutputMatrixColumns']

    df = pd.DataFrame(pars, columns=cols)

    kepid = df.loc[:, 'keplerId'].astype(int)
    tce = df.loc[:, 'planetIndexNumber'].astype(int)
    tceId = tools.createTceId(kepid, tce).as_matrix()
    tceId = np.round(tceId).astype(int)   #Round off error is a problem here!
    assert(len(tceId) == len(np.unique(tceId)))
    df.index = tceId

    return df