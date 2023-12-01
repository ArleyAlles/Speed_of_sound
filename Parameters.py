import numpy as np


class Parameters:
    """
        * Parameters of PC-SAFT EoS.
        In this class the following methods are present:
            - parameters: All parameters for PC-SAFT EoS
    """
    def __int__(self):
        pass
    @staticmethod
    def parameters(dataset):
        """
        All parameters for PC-SAFT.
        In this case, cross-association between H2O and CO2 was not considered. Only H2O association is available
            :param dataset: Which dataset will be used
        """
        #The parameters of PC-SAFT for H20 were obtained in the paper: Comparison of two Association Models
        #(Elliott-Suresh-Donohue and Simplified PC-SAFT) for complex Phase Equilibria of Hydrocarbon-Water
        #and Amine-Containing Mixtures (Ref. 17; Table 2.)
        R          = 8.314 #-----------------------------------------------------------> J/K*mol
        #sigma      = np.array([3.0007, 2.7852])  #------------------------------------> Diameter of molecule (A°)
        #Ek         = np.array([207.84, 169.21]) #------------------------------------->
        #mi         = np.array([2.0, 2.0729]) #---------------------------------------->
        #Eabk       = np.array([1506.4, 0]) #------------------------------------------>
        #Kab        = np.array([0.1550, 0]) #------------------------------------------>
        sigma      = np.array([2.229, 2.731])  #---------------------------------------> Temperature-independent segment diameter (A°)
        Ek         = np.array([141.66, 157.25]) #--------------------------------------> Depth of pair potential/Boltzmann constant (K)
        mi         = np.array([2.1945, 2.2280]) #--------------------------------------> Number of segments (dimensionless)
        Eabk       = np.array([1804.17, 307.41]) #-------------------------------------> Association energy (K)
        Kab        = np.array([0.2039, 0.0287]) #--------------------------------------> Association volume (dimensionless)
        Kij        = np.array([0.0]) #-------------------------------------------------> Binary interaction parameters
        A          = np.array([3.470, 5.457]) #----------------------------------------> Heat capacity coefficient
        B          = np.array([1.450e-3, 1.045e-3]) #----------------------------------> Heat capacity coefficient
        C          = np.array([0, 0]) #------------------------------------------------> Heat capacity coefficient
        D          = np.array([0.121e5, -1.157e5]) #-----------------------------------> Heat capacity coefficient
        molar_mass = np.array([18.01, 44.01]) #----------------------------------------> Molar mass
        code       = [0,1,1,1,1,0] #---------------------------------------------------> Code of sites
        num_sites  = np.array([4, 0])  # ----------------------------------------------> Number of sites for each compound

        if dataset == 0:
            x = np.array([0.9882, 0.0118])  # ------------------------------------------> Molar fraction
        elif dataset == 1:
            x = np.array([0.9934, 0.0066])  # -----------------------------------------> Molar fraction
        elif dataset == 2:
            x = np.array([0.9985, 0.0015])  # -----------------------------------------> Molar fraction
        else:
            raise TypeError("There is no such experimental points in this composition. The options are:\n"
                            "x (H2O): 0.9882 ----------> dataset=0\n"
                            "x (H2O): 0.9934 ----------> dataset=1\n"
                            "x (H2O): 0.9985 ----------> dataset=2")

        return R, sigma, Ek, mi, Eabk, Kab, Kij, x, A, B, C, D, molar_mass, code, num_sites
