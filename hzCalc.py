import numpy as np

def hzBorder(teff):
    """
    # The coefficients are as follows. The columns, i, are arranged according to
    # the HZ limits given in the paper.
    #
    # i = 0 --> Recent Venus
    # i = 1 --> Runaway Greenhouse
    # i = 2 --> Maximum Greenhouse
    # i = 3 --> Early Mars
    # i = 4 --> Runaway Greenhouse for 5 ME
    # i = 5 --> Runaway Greenhouse for 0.1 ME
    # First row: S_effSun(i) 
    # Second row: a(i)
    # Third row:  b(i)
    # Fourth row: c(i)
    # Fifth row:  d(i)
    """
    seffsun=np.zeros(6)
    seffsun[0] = 1.77600E+00;
    seffsun[1] = 1.10700E+00;
    seffsun[2] = 3.56000E-01;   
    seffsun[3] = 3.20000E-01;   
    seffsun[4] = 1.18800E+00;   
    seffsun[5] = 9.90000E-01;
    
    a=np.zeros(6)
    a[0] = 2.13600E-04;
    a[1] = 1.33200E-04;
    a[2] = 6.17100E-05;
    a[3] = 5.54700E-05;
    a[4] = 1.43300E-04;
    a[5] = 1.20900E-04;
    
    b=np.zeros(6)
    b[0] = 2.53300E-08;
    b[1] = 1.58000E-08;
    b[2] = 1.69800E-09;
    b[3] = 1.52600E-09;   
    b[4] = 1.70700E-08;
    b[5] = 1.40400E-08;
    
    c=np.zeros(6)
    c[0] = -1.33200E-11;
    c[1] = -8.30800E-12;
    c[2] = -3.19800E-12;
    c[3] = -2.87400E-12;
    c[4] = -8.96800E-12;
    c[5] = -7.41800E-12;
    
    d=np.zeros(6)
    d[0] = -3.09700E-15;
    d[1] = -1.93100E-15;
    d[2] = -5.57500E-16;
    d[3] = -5.01100E-16;
    d[4] = -2.08400E-15;
    d[5] = -1.71300E-15;
    
    seff=np.zeros((len(teff),6))
    
    for i,v in enumerate(teff):
        for j in np.arange(0,5,1):
            tstar=teff[i]-5777.0
            seff[i,j]=seffsun[j]+ a[j]*tstar+ b[j]*tstar**2.0+ c[j]*tstar**3.0+ d[j]*tstar**4.0 

    return seff

def inHZ(teff,S,which=[0,3]):
    """
    Given a stellar temerpature in K and an insolation flux
    Return whether it lies in the HZ
    Which specifies the upper and lower boundaries of the HZ.
    returns a T/F
    """
    
    hzb=hzBorder([teff])    
    inside= (S<=hzb[0,which[0]]) & (S>=hzb[0,which[1]])
    
    return inside