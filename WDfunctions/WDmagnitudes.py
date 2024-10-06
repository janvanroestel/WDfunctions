import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator as RGI
import copy
import pkg_resources

def calc_logg(M,R):
    '''calculate logg [cgs] from solarmass and solarradius

    input :
    M : float or array-like
        mass [M_sun]
    R : float or array-like
        Radius [R_sun]

    output : float or array-like
        logg [cgs]
    '''

    G = 6.673*10**-11
    M_sun = 1.988*10**30
    R_sun = 6.955*10**8

    return np.log10(G*M*M_sun*(R*R_sun)**-2)+2

def calc_R(M,logg):
    '''calculate R from M and logg [cgs]

    input :
    M : float or array-like
        mass [M_sun]
    output : float or array-like
        logg [cgs]
        
    output:
    R : float or array-like
        Radius [R_sun]
    '''

    G = 6.673*10**-11
    M_sun = 1.988*10**30
    R_sun = 6.955*10**8

    #np.log10(G*M*M_sun*(R*R_sun)**-2)+2

    return np.sqrt(G*M*M_sun / (10**(logg-2)) )/R_sun


class WDSED():
    "Class that returns WD magnitudes"

    def __init__(self,Table):
        """ Load a table from 
            http://www.astro.umontreal.ca/~bergeron/CoolingModels/ as input
        """

        # replace 'log g' with 'logg' in header

        # load the data
        table_path = ddir = pkg_resources.resource_filename('WDfunctions', 'data/Table_DA.txt')
        self.data = np.genfromtxt(table_path,skip_header=1,names=True)
    
        # make the interpolation functions for each filter
        cols = {}
        cols['FUV'] = 'FUV'
        cols['NUV'] = 'NUV'
        cols['SDSSu'] = 'u'
        cols['SDSSg'] = 'g'
        cols['SDSSr'] = 'r'
        cols['SDSSi'] = 'i'
        cols['SDSSz'] = 'z'
        cols['PSg'] = 'g_1'
        cols['PSr'] = 'r_1'
        cols['PSi'] = 'i_1'
        cols['PSz'] = 'z_1'
        cols['PSy'] = 'y'
        cols['J'] = 'J'
        cols['Y'] = 'Y'
        cols['W1'] = 'W1'
        cols['W2'] = 'W2'
        cols['R'] = 'R'

        self.cols = cols

        # make interpolators
        NT = np.size(np.unique(self.data['Teff']))
        points = [self.data['Teff'].reshape(5,NT)[0],
                  self.data['logg'].reshape(5,NT)[:,0]]

        interpolators = dict()
        for key, value in cols.items():
            interpolators[key] = RGI(points,self.data[value].reshape(5,NT).T)
        self._interpolators = interpolators

        self.model_mass = RGI(points,self.data["MMo"].reshape(5,NT).T)



    # for testing, this method is slow... (relativly)
    def _interpolator(self,Teff,logg,colname):
        output = griddata(points=(self.data['Teff'],self.data['logg']),
                        values=self.data[self.cols[colname]],xi=np.c_[Teff,logg],
                        method='linear')
        return output

    def get_magnitudes(self,filters,T,M=None,R=None,logg=None,distance=10.):
    
        if logg is None and M is not None and R is not None:
            logg = calc_logg(M,R)
        elif logg is not None and M is None and R is None:
            logg = logg
        else:
            print('Error; provide only logg OR M and R')

        output = [self._interpolators[f]([T,logg]) for f in filters]
        output = np.atleast_2d(np.array(output))

        # convert to right radius
        if R is not None:
            M_model = self.model_mass([T,logg])
            R_model = calc_R(M_model,logg)
            radius_mag = -2.5*np.log10((R/R_model)**2)
            print(radius_mag)
            output += radius_mag    
        
        # set to the right distance
        mu = 5*np.log10(distance)-5 
        output += mu

        return output
        

