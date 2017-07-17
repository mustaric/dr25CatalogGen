
"""Tools for interacting with DV ouputs"""

from scipy.io import matlab
import pandas as pd
import numpy as np
import tools

__version__ = "$Id: dv.py 2291 2016-03-09 23:57:52Z fmullall $"
__URL__ = "$URL: svn+ssh://flux.amn.nasa.gov/home/fmullall/svn/kepler/py/dv.py $"



def loadDvOutputMatrix(fn):
    mat = matlab.loadmat(fn, struct_as_record=False, squeeze_me=True)

    pars = mat['dvOutputMatrix']
    cols = mat['dvOutputMatrixColumns']

    out = nca.Nca(pars)
    out.setLookup(1, cols)
    return out



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