import pandas as pd
import numpy as np
import glob
import os
from scipy.interpolate import griddata
import pkg_resources



def _get_MTR(filename):
    #print(filename)
    Rsun = 695800000.0

    data = np.loadtxt(filename)

    log10T = data[:,0]
    R = data[:,4]
    M = float(os.path.basename(filename).replace('Msun.trk','')) * np.ones_like(R)
    
    return np.c_[M,log10T,R]



def _get_MTage(filename):
    #print(filename)
    Rsun = 695800000.0

    data = np.loadtxt(filename)

    log10T = data[:,0]
    R = data[:,4]
    M = float(os.path.basename(filename)[:3]) * np.ones_like(R)
    age = data[:,3]    

    return np.c_[M,log10T,age]



def _make_MTR(ddir=None):
    if ddir is None:
        ddir = pkg_resources.resource_filename('WDfunctions', 'data/')
    files = glob.glob(ddir+'WDtracks/*.trk')
    files.sort()
    alldata = np.vstack(np.array([_get_MTR(f) for f in files],dtype=object))

    def WD_MTR(M,logT):
        """ 
        input:
        M : float or array-like
            the white dwarf mass in solar units
        logT : float or array-like
            the white dwarf temperature in log10(Kelvin)

        output:
        R : float or array-like
            the white dwarf radius in solar units
        """
        if (isinstance(M,float) or isinstance(M,int)) and (isinstance(logT,float) or isinstance(logT,int)):
            input_values = np.array([M,logT])
        else:
            # assuming they are arrays
            M = np.array(M,dtype=float)
            logT = np.array(logT,dtype=float)
            input_values = np.c_[M,logT]
        return griddata(alldata[:,:2],alldata[:,2],input_values,method='cubic')
    return WD_MTR



def _make_MTage(ddir=None):
    if ddir is None:
        ddir = pkg_resources.resource_filename('WDfunctions', 'data/')
    files = glob.glob(ddir+'WDtracks/*.trk')
    files.sort()
    alldata = np.vstack(np.array([_get_MTage(f) for f in files],dtype=object))

    def WD_MTage(M,logT):
        """ 
        input:
        M : float or array-like
            the white dwarf mass in solar units
        logT : float or array-like
            the white dwarf temperature in log10(Kelvin)

        output:
        age : float or array-like
            the white dwarf age in Gyr
        """
        if (isinstance(M,float) or isinstance(M,int)) and (isinstance(logT,float) or isinstance(logT,int)):
            input_values = np.array([M,logT])
        else:
            # assuming they are arrays
            M = np.array(M,dtype=float)
            logT = np.array(logT,dtype=float)
            input_values = np.c_[M,logT]
        return griddata(alldata[:,:2],alldata[:,2],input_values,method='cubic')
    return WD_MTage

    f = lambda pos : griddata(alldata[:,:2],alldata[:,2],pos,method='cubic')

    return f



def WD_MR_Eggleton(M):
    '''calculate the radius of a zero temperature wd given a mass

    input :
    M : float or array-like
        the mass of the wd in solar units
    output : float or array-like
        the radius in solar units
    '''

    # constants in SI units
    G = 6.673*10**-11
    M_sun = 1.988*10**30
    R_sun = 6.955*10**8

    M_ch = 1.44
    M_p = 0.00057
    R = 0.0114 * (( M/M_ch )**(-2./3.)-(M/M_ch)**(2./3.))**(0.5) * (1+3.5*(M/M_p)**(-2./3.)+(M/M_p)**(-1))**(-2./3.)
    return R



def WD_MR_Nauenberg(M,mu=2.001833):
    """
    A mass radius function for white dwarfs
    from https://articles.adsabs.harvard.edu//full/1972ApJ...175..417N/0000420.000.html

    
    mu is average molecular weight per electron!
    pure He: 4.0026/2 = 2.0013
    pure C: 12.011/6 = 2.001833
    pure O: 15.999/8 = 1.999875
    pure Ne: 20.1797/10 = 2.01797

    input :
    M : float or array-like
        mass of the wd in solar units
    mu : float
        average molecular weight per electron, in proton mass
    output : float or array-like
        the radius in solar units
    """

    # calculate M3
    M3 = 5.816 * mu**(-2.) # in solar units
    
    # calculate the radius
    R = 0.0225 / mu 
    R *= np.sqrt(1-(M/M3)**(4./3.))
    R /= (M/M3)**(1./3.)

    return R


    
# run some functions to make them on import
WD_MTR = _make_MTR()
WD_MTage = _make_MTage()
