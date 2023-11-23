import numpy as np


class Parameters:
    def __int__(self):
        pass
    @staticmethod
    def parameters():

        R          = 8.314 #----------------------------------------------------------> J/K*mol
        sigma      = np.array([3.0007, 2.7852])  #------------------------------------> Diameter of molecule (A°)
        Ek         = np.array([366.51, 169.21]) #------------------------------------->
        mi         = np.array([1.0656, 2.0729]) #------------------------------------->
        Eabk       = np.array([2500.7, 0]) #------------------------------------------>
        Kab        = np.array([0.034868, 0]) #---------------------------------------->
        Kij        = np.array([0.0]) #------------------------------------------------> Binary interaction parameters
        x          = np.array([0.9882, 0.0118]) #-------------------------------------> Molar fraction
        A          = np.array([3.470, 5.457]) #---------------------------------------> Heat capacity coefficient
        B          = np.array([1.450e-3, 1.045e-3]) #---------------------------------> Heat capacity coefficient
        C          = np.array([0, 0]) #-----------------------------------------------> Heat capacity coefficient
        D          = np.array([0.121e5, -1.157e5]) #----------------------------------> Heat capacity coefficient
        molar_mass = np.array([18.01, 44.01]) #---------------------------------------> Molar mass
        density    = np.array([55525, 473.43]) #--------------------------------------> Molar density (mol/m³)
        code       = [0,1,1,1,1,0] #--------------------------------------------------> Code of sites
        num_sites  = np.array([4, 0]) #-----------------------------------------------> Number of sites for each compound

        return R, sigma, Ek, mi, Eabk, Kab, Kij, x, A, B, C, D, molar_mass, density, code, num_sites
