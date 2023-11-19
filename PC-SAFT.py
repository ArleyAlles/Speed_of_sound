import numpy as np
import math
from sympy import *

class PCSAFT:
    """
    * Calculation of properties using PC-SAFT EoS.
    In this class the following methods are present:

        - a_hc: Hard-Chain contribution
        - a_disp: Dispersion contribution
        - helmholtz: Helmholtz free energy
        - compressibility: Compressibility factor
        - pressure: Pressure of system
        - first_derivative: First derivative of Helmholtz free energy for reduced density (eta)
        - second_derivative: Second derivative of Helmholtz free energy for reduced density (eta)
        - roots: Calculation of root (eta) using Newton's method plus bisection method
        - coef_fugacity: Fugacity coefficient in a specific phase
        - enthalpy: Residual enthalpy
    """

    def __init__(self, R, T, sigma, Ek, mi, eabk, kab, Kij, A, B, C, D):
        self.R     = R #---------------------> Gas constant (J.mol-1.K-1)
        self.T     = T #---------------------> Temperature (K)
        self.sigma = sigma #-----------------> Temperature-independent segment diameter (A°)
        self.Ek    = Ek #--------------------> Depth of pair potential/Boltzmann constant (K)
        self.mi    = mi #--------------------> Number of segments (dimensionless)
        self.eabk  = eabk #------------------> Association energy (K)
        self.kab   = kab #-------------------> Association volume (dimensionless)
        self.Kij   = Kij #-------------------> Binary interaction parameters (dimensionless). Order for 3 comp. = [K01, K02, K12]
        self.k     = 1.380649e-23  # --------> Boltzmann's constant (J.K-1)
        self.Na    = 6.0221e+23  # ----------> Avogadro's number (mol-1)

        # Universal constants a = [[a00, a01, a02 ...], [a10, a11, a12 ...], [a20, a21, a22 ...]]
        self.a = np.array([[0.910563, 0.636128, 2.686135, -26.547362, 97.759209, -159.591541, 91.297774],
                           [-0.308402, 0.186053, -2.503005, 21.419794, -65.255885, 83.318680, -33.746923],
                           [-0.090615, 0.452784, 0.596270, -1.724183, -4.130211, 13.776632, -8.672847]])

        # Universal constants b = [[b00, b01, b02 ...], [b10, b11, b12 ...], [b20, b21, b22 ...]]
        self.b = np.array([[0.724095, 2.238279, -4.002585, -21.003577, 26.855641, 206.551338, -355.602356],
                           [-0.575550, 0.699510, 3.892567, -17.215472, 192.672264, -161.826462, -165.207693],
                           [0.097688, -0.255757, -9.155856, 20.642076, -38.804430, 93.626774, -29.666906]])

        self.A = A #Heat capacity coefficient
        self.B = B #Heat capacity coefficient
        self.C = C #Heat capacity coefficient
        self.D = D #Heat capacity coefficient

    def kij(self):
        cont   = 0
        row    = len(self.mi)
        column = len(self.mi)
        matrix = np.zeros((row, column))
        Kij    = []

        for i in range(row):
            for j in range(column):
                if j>i:
                    matrix[i,j] = self.Kij[cont]
                    cont += 1
                else:
                    matrix[i,j] = 0

        matrix_T   = np.transpose(matrix)
        matrix_kij = matrix+matrix_T

        for i in range(row):
            for j in range(column):
                Kij.append(matrix_kij[i,j])

        return Kij

    def a_hc(self, eta, T, x):
        """
        * Hard-chain contribution
            :param eta: Reduced density
            :param x: Molar fraction
        """

        di      = self.sigma*(1-0.12*np.exp(-3*(self.Ek/T)))
        ro      = (6/math.pi)*eta*(sum(x*self.mi*(di**3)))**-1
        epsilon = np.zeros(4)

        for i in range(4):
            epsilon[i] = (math.pi/6)*ro*sum(x*self.mi*(di**i))

        Pa = 1/(1-epsilon[3])
        Pb = (di/2)*(3*epsilon[2])/((1-epsilon[3])**2)
        Pc = ((di/2)**2)*2*(epsilon[2]**2)/((1-epsilon[3])**3)
        gij = Pa+Pb+Pc
        Pa = (3*epsilon[1]*epsilon[2])/(1-epsilon[3])
        Pb = (epsilon[2]**3)/(epsilon[3]*((1-epsilon[3])**2))
        Pc = (((epsilon[2]**3)/(epsilon[3]**2))-epsilon[0])*np.log(1-epsilon[3])
        ahs = (1/epsilon[0])*(Pa+Pb+Pc)
        m = sum(x*self.mi)
        ahc = m*ahs-sum(x*(self.mi-1)*np.log(gij))

        return ahc

    def a_disp(self, eta, T, x):
        """
        * Dispersion contribution
            :param ro_molar: molar density of mixture
            :param T: temperature
            :param x: Molar fraction
        """
        di = self.sigma*(1-0.12*np.exp(-3*(self.Ek/T)))
        ro = (6/math.pi)*eta*(sum(x*self.mi*(di**3)))**-1
        m  = sum(x*self.mi)
        a  = self.a[0] + (((m-1)/m)*self.a[1]) + (((m-1)/m)*((m-2)/m)*self.a[2])
        b  = self.b[0] + (((m-1)/m)*self.b[1]) + (((m-1)/m)*((m-2)/m)*self.b[2])

        I1 = 0
        I2 = 0

        for i in range(7):
            I1 += a[i]*(eta**i)
            I2 += b[i]*(eta**i)


        Pa = m*((8*eta)-(2*(eta**2)))/((1-eta)**4)
        Pb = ((20*eta)-(27*eta**2)+(12*eta**3)-(2*eta**4))/(((1-eta)*(2-eta))**2)
        C1 = 1/(1+Pa+(1-m)*Pb)
        cont = 0
        m2esig3 = 0
        m2esi2g3 = 0
        Kij = self.kij()

        for l in range(len(x)):
            for n in range(len(x)):
                eij       = (1-Kij[cont])*np.sqrt(self.Ek[l]*self.Ek[n])
                sigma_ij  = 0.5*(self.sigma[l]+self.sigma[n])
                m2esig3  += x[l]*x[n]*self.mi[l]*self.mi[n]*(eij/T)*(sigma_ij**3)
                m2esi2g3 += x[l]*x[n]*self.mi[l]*self.mi[n]*((eij/T)**2)*(sigma_ij**3)
                cont += 1

        adisp = (-2*math.pi*ro*I1*m2esig3) - (math.pi*ro*m*C1*I2*m2esi2g3)
        return adisp

    def association_parameters(self, eta, T, x, code):
        """
        Generation of all pair combinations of sites (without pairs with same site,
        like AA, BB and so on). Such procedure is accomplished with the code input,
        which is the label of all possible sites combination. Each number will be
        responsible for a specific value of association strength. For example:

        For a system of Methanol(i) and ethanol(j), which has 3 sites of association each one,
        is given the code below:
            code = [0,1,0,0,2,1,0,0,2,2,2,0,0,3,3]
            0 ---> Association between sites with same characteristic (donors OR receptors)
            1 ---> Association between sites of the same component (i) but with different
                   characteristics (donor AND receptor)
            2 ---> Association between sites of different components (i with j or vice-versa)
                   but with different characteristics (donor AND receptor)
            3 ---> Association between sites of the same component (j) but with different
                   characteristics (donor AND receptor)

        Which is the same of the combinations:
            combinations = [AiBi, AiCi, AiAj, AiBj, AiCj, BiCi, BiAj, BiBj, BiCj,
                            CiAj, CiBj, CiCj, AjBj, AjCj, BjCj]

            A, B, C ---> Sites of molecules

        OBS: For a given code, the order of label MUST follow the order:
            For a system of 2 molecules:
                * Label 1 ---> molecule (i) with molecule (i)
                * Label 2 ---> molecule (i) with molecule (j)
                * Label 3 ---> molecule (j) with molecule (j)

            For a system of 3 molecules:
                * Label 1 ---> molecule (i) with molecule (i)
                * Label 2 ---> molecule (i) with molecule (j)
                * Label 3 ---> molecule (i) with molecule (k)
                * Label 4 ---> molecule (j) with molecule (j)
                * Label 5 ---> molecule (j) with molecule (k)
                * Label 6 ---> molecule (k) with molecule (k)

            And so on ...

        :param code: list of labels for each pairwise association
        :return: association strength for each combination
        """

        num_comp = len(self.mi)
        di = self.sigma * (1 - 0.12 * np.exp(-3 * (self.Ek / T)))
        ro = (6/math.pi)*eta*(sum(x*self.mi*(di**3)))**-1
        epsilon  = np.zeros(4)

        for i in range(4):
            epsilon[i] = (math.pi/6)*ro*sum(x*self.mi*(di**i))

        #Energy and volume of association
        delta = []
        for i in range(num_comp):
            for j in range (num_comp):
                if j>=i:
                    dij     = (di[i]*di[j])/(di[i]+di[j])
                    Pa      = 1/(1-epsilon[3])
                    Pb      = 3*dij
                    Pc      = epsilon[2]/((1-epsilon[3])**2)
                    Pd      = 2*dij**2
                    Pe      = (epsilon[2]**2)/((1-epsilon[3])**3)
                    gij_seg = Pa+(Pb*Pc)+(Pd*Pe)
                    eaibjk  = (self.eabk[i]+self.eabk[j])*0.5
                    Pa      = math.sqrt(self.kab[i]*self.kab[j])
                    Pb      = math.sqrt(self.sigma[i]*self.sigma[j])
                    Pc      = 0.5*(self.sigma[i]+self.sigma[j])
                    kaibj   = Pa*(Pb/Pc)**3
                    Pa      = np.exp(eaibjk/T)-1
                    delta.append(gij_seg*(dij**3)*Pa*kaibj)

        #Creating vectors for combination sites
        Delta = np.zeros(len(code))
        for i, num in enumerate(code):
            if num == 0:
                Delta[i] = 0
            else:
                Delta[i] = delta[num-1]

        return Delta

    def a_assoc(self, num_sites, eta, T, x, code):
        """
        * Calculation of association Helmholtz free energy term, using successive substitution method.
            :param num_sites: number of sites for each component
            :param ro_molar: molar density of mixture
            :param T: temperature
            :param x: molar fraction of each component in the mixture
            :param code: list of labels for each pairwise association
        """

        #Creating the energy of association matrix
        delta  = self.association_parameters(eta, T, x, code)
        cont   = 0
        row    = sum(num_sites)
        column = sum(num_sites)
        matrix = np.zeros((row, column))

        for i in range(row):
            for j in range(column):
                if j > i:
                    matrix[i, j] = delta[cont]
                    cont += 1
                else:
                    matrix[i, j] = 0

        matrix_T = np.transpose(matrix)
        Delta    = matrix + matrix_T

        #Initial guess for X
        X          = 0.1*np.ones(sum(num_sites))

        # Calculation of molar density of molecules
        di         = self.sigma*(1-0.12*np.exp(-3*(self.Ek/T)))
        ro         = (6/math.pi)*eta*(sum(x*self.mi*(di**3)))**-1
        ro_comp    = ro*x
        Ro         = np.array([])
        for i, n_sites in enumerate(num_sites):
            ro     = np.array([ro_comp[i]]*n_sites)
            Ro     = np.concatenate((Ro,ro), axis=0)


        #Newton-Raphson method
        error = 1
        tol = 1e-10
        f = []

        while error > tol:

            #Jacobian matrix
            jacobian = np.zeros((row, column))
            for i in range(row):
                for j in range(column):
                    if i == j:
                        jacobian[i,j] = 1+sum(Ro*X*Delta[i,:])
                    else:
                        jacobian[i,j] = Ro[j]*X[i]*Delta[i,j]

                f.append(X[i]+X[i]*sum(Ro*X*Delta[i,:])-1)


            #Pseudoinverse
            J_pseudo = np.linalg.inv(np.transpose(jacobian)*jacobian)*np.transpose(jacobian)

            X_mod = np.reshape(X,(row,1))

            #Vector of functions
            F = np.reshape(np.array(f),(row,1))
            f.clear()

            #Actualization of X
            X_new = X_mod - np.dot(J_pseudo,F)
            error = sum((X_new-X_mod)**2)
            X     = np.reshape(X_new,(row,))

        #Association Helmholtz free energy
        a_association = 0
        for i, n_sites in enumerate(num_sites):
            Pa = 0
            for j in range(n_sites):
                Pa   += np.log(X[j])-(0.5*X[j])
            a_association += x[i]*(Pa+(0.5*n_sites))
        return a_association

    def a_res(self, eta, T, x, num_sites, code, association):
        """
        * Helmholtz free energy for all contributions
            :param ro_molar: molar density of mixture
            :param T: temperature
            :param x: Molar fraction
        :return:
        """
        a_hc   = self.a_hc(eta, T, x)
        a_disp = self.a_disp(eta, T, x)

        if association == 0:
            a_res = a_hc+a_disp
        elif association == 1:
            a_assoc = self.a_assoc(num_sites, eta, T, x, code)
            #a_assoc = self.association(num_sites, eta, T, x, code)
            a_res = a_hc+a_disp+a_assoc

        return a_res

    def compressibility(self, eta, T, x, num_sites, code, association):
        """
        Compressibility factor
        :param ro: Total number density of spherical segments (for mixture)
        :param x: Molar fraction
        """
        h = eta*1e-6
        am_2 = self.a_res(eta-(2*h), T, x, num_sites, code, association)
        am2  = self.a_res(eta+(2*h), T, x, num_sites, code, association)
        am_1 = self.a_res(eta-h, T, x, num_sites, code, association)
        am1  = self.a_res(eta+h, T, x, num_sites, code, association)
        da_m = (am_2-8*am_1+8*am1-am2)/(12*h)
        Z    = 1+(eta*da_m)

        return Z

    def pressure(self, eta, T, x, num_sites, code, association):
        """
        * Pressure of system
            :param eta: Reduced density
            :param x: Molar fraction
        """

        Z = self.compressibility(eta, T, x, num_sites, code, association)
        di = self.sigma * (1 - 0.12 * np.exp(-3 * (self.Ek / T)))
        ro = (6 / math.pi) * eta * (sum(x * self.mi * (di ** 3))) ** -1
        P = Z*self.k*T*ro*1e30 #Pascal

        return P

    def dpdro(self, eta, T, x, num_sites, code, association):
        """
        * First derivative of pressure by density
        :param ro_molar: molar density
        :param T: temperature
        :param x: molar fraction
        """
        h = eta*1e-7
        am_2     = self.pressure(eta-(2*h), T, x, num_sites, code, association)
        am2      = self.pressure(eta+(2*h), T, x, num_sites, code, association)
        am_1     = self.pressure(eta-h, T, x, num_sites, code, association)
        am1      = self.pressure(eta+h, T, x, num_sites, code, association)
        dPdeta = (am_2 - 8 * am_1 + 8 * am1 - am2) / (12 * h)

        di = self.sigma * (1 - 0.12 * np.exp(-3 * (self.Ek / T)))
        detadro  = ((math.pi/6)*self.Na*(sum(x*self.mi*(di**3))))*1e-30
        dpdro    = dPdeta*detadro
        return dpdro

    def dpdt(self, eta, T, x, num_sites, code, association):
        """
        * First derivative of pressure by temperature
        :param ro_molar: molar density
        :param T: temperature
        :param x: molar fraction
        """

        h = T*1e-7
        am_2 = self.compressibility(eta, T-(2*h), x, num_sites, code, association)
        am2  = self.compressibility(eta, T+(2*h), x, num_sites, code, association)
        am_1 = self.compressibility(eta, T-h, x, num_sites, code, association)
        am1  = self.compressibility(eta, T+h, x, num_sites, code, association)
        dZdt = (am_2-8*am_1+8*am1-am2)/(12*h)

        di = self.sigma * (1 - 0.12 * np.exp(-3 * (self.Ek / T)))
        ro = (6 / math.pi) * eta * (sum(x * self.mi * (di ** 3))) ** -1
        Z = self.compressibility(eta, T, x, num_sites, code, association)
        dpdt = self.k*T*ro*1e30*(dZdt+(Z/T))
        return dpdt

    def dadt(self, eta, T, x, num_sites, code, association):
        """
        * First derivative of pressure by temperature
        :param ro_molar: molar density
        :param T: temperature
        :param x: molar fraction
        """

        h = T*1e-7
        am_2 = self.a_res(eta, T-(2*h), x, num_sites, code, association)
        am2  = self.a_res(eta, T+(2*h), x, num_sites, code, association)
        am_1 = self.a_res(eta, T-h, x, num_sites, code, association)
        am1  = self.a_res(eta, T+h, x, num_sites, code, association)
        dadt = (am_2-8*am_1+8*am1-am2)/(12*h)

        return dadt

    def d2adt2(self, eta, T, x, num_sites, code, association):
        """
        * Second derivative of Hemholtz free energy by temperature
        :param ro_molar: molar density
        :param T: temperature
        :param x: molar fraction
        """

        h      = T*1e-7
        am_2   = self.dadt(eta, T-(2*h), x, num_sites, code, association)
        am2    = self.dadt(eta, T+(2*h), x, num_sites, code, association)
        am_1   = self.dadt(eta, T-h, x, num_sites, code, association)
        am1    = self.dadt(eta, T+h, x, num_sites, code, association)
        am     = self.dadt(eta, T, x, num_sites, code, association)
        d2adt2 = (am_2-8*am_1+8*am1-am2)/(12*h)

        return d2adt2

    def cv(self, eta, T, x, num_sites, code, association):
        """
        Heat capacity at constant volume
        :param ro_molar: molar density
        :param T: temperature
        :param x: molar fraction
        """
        dadt     = self.dadt(eta, T, x, num_sites, code, association)
        d2adt2   = self.d2adt2(eta, T, x, num_sites, code, association)
        cv_res   = (-self.R*d2adt2*T**2)+(-2*self.R*T*dadt)
        cp_ideal = sum(((self.A+(self.B*T)+(self.C*T**2)+(self.D/(T**2)))*self.R))
        cv_ideal = cp_ideal-self.R
        cv       = cv_ideal+cv_res
        return cv

    def cp(self, eta, T, x, num_sites, code, association):
        """
        * Heat capacity at constant pressure
        :param ro_molar: molar density
        :param T: temperature
        :param x: molar fraction
        """
        di       = self.sigma * (1-0.12*np.exp(-3 * (self.Ek / T)))
        ro       = (6 / math.pi) * eta * (sum(x * self.mi * (di ** 3))) ** -1
        density  = ro/(self.Na*1e-30)
        dpdt     = self.dpdt(eta, T, x, num_sites, code, association)
        dpdro    = self.dpdro(eta, T, x, num_sites, code, association)
        cp_ideal = sum((self.A+(self.B*T)+(self.C*T**2)+(self.D/(T**2)))*self.R)
        Pa       = T/(density**2)
        cv       = self.cv(eta,T,x,num_sites,code,association)
        cp       = cv+(Pa*((dpdt**2)/dpdro))

        return cp

    def speed_of_sound(self, molar_mass, eta, T, x, num_sites, code, association):
        """
        Speed of sound
        :param molar_mass: molar mass
        :param ro_molar: molar density
        :param T: temperature
        :param x: molar fraction
        """
        Mw = sum(x*molar_mass)
        dpdro            = self.dpdro(eta, T, x, num_sites, code, association)
        cp               = self.cp(eta, T, x, num_sites, code, association)
        cv               = self.cv(eta, T, x, num_sites, code, association)
        gamma            = cp/cv
        u                = math.sqrt((gamma/Mw)*dpdro*1e3)

        print(f"Isocoric heat: {cv:.2f} J/mol*K")
        print(f"Isobaric heat: {cp:.2f} J/mol*K")
        print(f"Speed of sound: {u:.2f} m/s")

    def d2pdro2(self, eta, T, x, num_sites, code, association):

        h = 0.000001 * eta
        am_2 = self.dpdro(eta - (2 * h), T, x, num_sites, code, association)
        am2 = self.dpdro(eta + (2 * h), T, x, num_sites, code, association)
        am_1 = self.dpdro(eta - h, T, x, num_sites, code, association)
        am1 = self.dpdro(eta + h, T, x, num_sites, code, association)
        d2pdro2 = (am_2 - 8 * am_1 + 8 * am1 - am2) / (12 * h)
        return d2pdro2


    def roots(self, Po, T, x, num_sites, code, association):
        """
        * Calculation of root (eta) using Newton's method plus bisection method.
            Initially, Newton's optimization method is used for finding stationary points.
            Such extrema points are used for establishing regions of roots. All roots are
            obtained by bisection method and Gibbs free energy was calculated for all roots
            for finding the corrected roots (that one which minimizes Gibbs).

            :param Po: Pressure of system
            :param x: Molar fraction
                """
        n     = 60
        tol   = 1e-12
        E     = 1e-15
        eta   = np.zeros(n)
        P     = np.zeros(n+2)
        Eta0  = np.linspace(E,(1-E),num=n)
        error = 1
        eta2  = np.zeros(n+2)
        intervals_P   = []
        intervals_eta = []

        #Newton's method
        while error>tol:
            for i, eta0 in enumerate(Eta0):
                df     = self.dpdro(eta0, T, x, num_sites, code, association)
                d2f    = self.d2pdro2(eta0, T, x, num_sites, code, association)
                eta[i] = eta0-(df/d2f)

            error = sum(abs(eta-Eta0))
            Eta0  = eta

        #P and eta extension
        for i in range(len(eta)):
            if i != n+1:
                P[i+1] = self.pressure(eta[i], T, x, num_sites, code, association)
            if i == n-1:
                P[n+1] = (self.pressure(eta[i], T, x, num_sites, code, association))**4

        eta2[1:n+1] = eta
        eta2[n+1] = 1

        #Finding intervals
        for i in range(len(P)):
            if i != 0:
                if (P[i-1]<Po<P[i]) or (P[i-1]>Po>P[i]):
                    intervals_P.append([P[i-1], P[i]])
                    intervals_eta.append([eta2[i-1], eta2[i]])


        #Finding roots (bisection)
        eta = []

        for _, int_eta in enumerate(intervals_eta):
            error = 1
            tol  = 1e-10
            xa   = int_eta[0]
            xb   = int_eta[1]
            xro  = 120 #guess

            while error > tol:
                xr = (xa+xb)/2
                fa = self.pressure(xa, T, x, num_sites, code, association)-Po
                fr = self.pressure(xr, T, x, num_sites, code, association)-Po
                if fa*fr<0:
                    xb = xr
                else:
                    xa = xr
                error = abs(xr-xro)
                xro = xr

            eta.append(xr)

        # Removing negative values
        eta = list(filter(lambda x: x > 0, eta))

        #Choosing the correct root
        if len(eta) == 1:
            return eta[0]
        else:
            G_res = []

            for _, e in enumerate(eta):
                a_res = self.a_res(e, T, x, num_sites, code, association)
                Z = self.compressibility(e, T, x, num_sites, code, association)
                G_res.append(a_res+(Z-1)-np.log(Z))

            # Removing 'nan' values
            G_res_new = [x for x in G_res if math.isnan(x) == False]
            min_G_res = min(G_res_new)
            Pos_min = G_res.index(min_G_res)

            return eta[Pos_min]


    def density_to_eta(self, x, density):

        mix_density = sum(x*density)
        di          = sigma*(1-0.12*np.exp(-3*(self.Ek/self.T)))
        ro          = mix_density*6.0221e-7
        eta         = (math.pi/6)*ro*sum(x*mi*di**3)

        return eta



class PCSAFT2:

    def __init__(self, R, T, sigma, Ek, mi, eabk, kab, Kij, A, B, C, D):
        self.R     = R #---------------------> Gas constant (J.mol-1.K-1)
        self.T     = T #---------------------> Temperature (K)
        self.sigma = sigma #-----------------> Temperature-independent segment diameter (A°)
        self.Ek    = Ek #--------------------> Depth of pair potential/Boltzmann constant (K)
        self.mi    = mi #--------------------> Number of segments (dimensionless)
        self.eabk  = eabk #------------------> Association energy (K)
        self.kab   = kab #-------------------> Association volume (dimensionless)
        self.Kij   = Kij #-------------------> Binary interaction parameters (dimensionless). Order for 3 comp. = [K01, K02, K12]
        self.k     = 1.380649e-23  # --------> Boltzmann's constant (J.K-1)
        self.Na    = 6.0221e+23  # ----------> Avogadro's number (mol-1)

        # Universal constants a = [[a00, a01, a02 ...], [a10, a11, a12 ...], [a20, a21, a22 ...]]
        self.a = np.array([[0.910563, 0.636128, 2.686135, -26.547362, 97.759209, -159.591541, 91.297774],
                           [-0.308402, 0.186053, -2.503005, 21.419794, -65.255885, 83.318680, -33.746923],
                           [-0.090615, 0.452784, 0.596270, -1.724183, -4.130211, 13.776632, -8.672847]])

        # Universal constants b = [[b00, b01, b02 ...], [b10, b11, b12 ...], [b20, b21, b22 ...]]
        self.b = np.array([[0.724095, 2.238279, -4.002585, -21.003577, 26.855641, 206.551338, -355.602356],
                           [-0.575550, 0.699510, 3.892567, -17.215472, 192.672264, -161.826462, -165.207693],
                           [0.097688, -0.255757, -9.155856, 20.642076, -38.804430, 93.626774, -29.666906]])

        self.A = A #Heat capacity coefficient
        self.B = B #Heat capacity coefficient
        self.C = C #Heat capacity coefficient
        self.D = D #Heat capacity coefficient

    def density_to_eta(self, x, density):

        mix_density = sum(x*density)
        di          = sigma*(1-0.12*np.exp(-3*(self.Ek/self.T)))
        ro          = mix_density*6.0221e-7
        eta         = (math.pi/6)*ro*sum(x*mi*di**3)

        return eta

    def kij(self):
        cont   = 0
        row    = len(self.mi)
        column = len(self.mi)
        matrix = np.zeros((row, column))
        Kij    = []

        for i in range(row):
            for j in range(column):
                if j>i:
                    matrix[i,j] = self.Kij[cont]
                    cont += 1
                else:
                    matrix[i,j] = 0

        matrix_T   = np.transpose(matrix)
        matrix_kij = matrix+matrix_T

        for i in range(row):
            for j in range(column):
                Kij.append(matrix_kij[i,j])

        return Kij

    def a_hc(self, x):
        """
        * Hard-chain contribution
            :param eta: Reduced density
            :param x: Molar fraction
        """
        eta, T = symbols('eta T ')
        di     = []
        gij    = []
        soma_0 = 0
        soma_1 = 0
        soma_2 = 0
        soma_3 = 0
        for i in range(len(self.Ek)):
            di.append(self.sigma[i]*(1-0.12*exp(-3*(self.Ek[i]/T))))
            soma_0 += x[i]*self.mi[i]*(di[i]**0)
            soma_1 += x[i]*self.mi[i]*(di[i]**1)
            soma_2 += x[i]*self.mi[i]*(di[i]**2)
            soma_3 += x[i]*self.mi[i]*(di[i]**3)

        for i in range(len(x)):
            ro       = (6/math.pi)*eta*soma_3**-1
            epsilon0 = (math.pi/6)*ro*soma_0
            epsilon1 = (math.pi/6)*ro*soma_1
            epsilon2 = (math.pi/6)*ro*soma_2
            epsilon3 = (math.pi/6)*ro*soma_3

            Pa  = 1/(1-epsilon3)
            Pb  = (di[i]/2)*(3*epsilon2)/((1-epsilon3)**2)
            Pc1 = (di[i]**2/4)
            Pc2 = 2*(Pow(epsilon2,2))/(Pow((1-epsilon3),3))
            gij.append(Pa+Pb+Pc1*Pc2)

        Pa  = (3*epsilon1*epsilon2)/(1-epsilon3)
        Pb  = (epsilon2**3)/(epsilon3*((1-epsilon3)**2))
        Pc  = (((epsilon2**3)/(epsilon3**2))-epsilon0)*log(1-epsilon3)
        ahs = (1/epsilon0)*(Pa+Pb+Pc)

        m     = 0
        soma2 = 0
        for i in range(len(x)):
            m     += x[i]*self.mi[i]
            soma2 += x[i]*(self.mi[i]-1)*log(gij[i])

        ahc = m*ahs-soma2
        return ahc

    def a_disp(self, x):

        eta, T = symbols('eta T ')
        di   = []
        m    = 0
        I1   = 0
        I2   = 0
        soma = 0
        for i in range(len(self.Ek)):
            di.append(self.sigma[i]*(1-0.12*exp(-3*(self.Ek[i]/T))))
            soma += x[i]*self.mi[i]*(di[i]**3)

        ro = (6/math.pi)*eta*soma**-1

        for i in range(len(x)):
            m += x[i]*self.mi[i]

        for i in range(7):
            a   = self.a[0,i]+(((m-1)/m)*self.a[1,i])+(((m-1)/m)*((m-2)/m)*self.a[2,i])
            b   = self.b[0,i]+(((m-1)/m)*self.b[1,i])+(((m-1)/m)*((m-2)/m)*self.b[2,i])
            I1 += a*eta**i
            I2 += b*eta**i

        Pa = m*((8*eta)-(2*(eta**2)))/((1-eta)**4)
        Pb = ((20*eta)-(27*(eta**2))+(12*(eta**3))-(2*(eta**4)))/(((1-eta)*(2-eta))**2)
        C1 = 1/(1+Pa+((1-m)*Pb))

        cont = 0
        m2esig3 = 0
        m2esi2g3 = 0
        Kij = self.kij()
        for l in range(len(x)):
            for n in range(len(x)):
                eij       = (1-Kij[cont])*np.sqrt(self.Ek[l]*self.Ek[n])
                sigma_ij  = 0.5*(self.sigma[l]+self.sigma[n])
                m2esig3  += x[l]*x[n]*self.mi[l]*self.mi[n]*(eij/T)*(sigma_ij**3)
                m2esi2g3 += x[l]*x[n]*self.mi[l]*self.mi[n]*((eij/T)**2)*(sigma_ij**3)
                cont     += 1

        adisp = (-2*math.pi*ro*I1*m2esig3)-(math.pi*ro*m*C1*I2*m2esi2g3)

        return adisp

    def association_parameters(self, eta, T, x, code):
        """
        Generation of all pair combinations of sites (without pairs with same site,
        like AA, BB and so on). Such procedure is accomplished with the code input,
        which is the label of all possible sites combination. Each number will be
        responsible for a specific value of association strength. For example:

        For a system of Methanol(i) and ethanol(j), which has 3 sites of association each one,
        is given the code below:
            code = [0,1,0,0,2,1,0,0,2,2,2,0,0,3,3]
            0 ---> Association between sites with same characteristic (donors OR receptors)
            1 ---> Association between sites of the same component (i) but with different
                   characteristics (donor AND receptor)
            2 ---> Association between sites of different components (i with j or vice-versa)
                   but with different characteristics (donor AND receptor)
            3 ---> Association between sites of the same component (j) but with different
                   characteristics (donor AND receptor)

        Which is the same of the combinations:
            combinations = [AiBi, AiCi, AiAj, AiBj, AiCj, BiCi, BiAj, BiBj, BiCj,
                            CiAj, CiBj, CiCj, AjBj, AjCj, BjCj]

            A, B, C ---> Sites of molecules

        OBS: For a given code, the order of label MUST follow the order:
            For a system of 2 molecules:
                * Label 1 ---> molecule (i) with molecule (i)
                * Label 2 ---> molecule (i) with molecule (j)
                * Label 3 ---> molecule (j) with molecule (j)

            For a system of 3 molecules:
                * Label 1 ---> molecule (i) with molecule (i)
                * Label 2 ---> molecule (i) with molecule (j)
                * Label 3 ---> molecule (i) with molecule (k)
                * Label 4 ---> molecule (j) with molecule (j)
                * Label 5 ---> molecule (j) with molecule (k)
                * Label 6 ---> molecule (k) with molecule (k)

            And so on ...

        :param code: list of labels for each pairwise association
        :return: association strength for each combination
        """

        num_comp = len(self.mi)
        di = self.sigma * (1 - 0.12 * np.exp(-3 * (self.Ek / T)))
        ro = (6/math.pi)*eta*(sum(x*self.mi*(di**3)))**-1
        epsilon  = np.zeros(4)

        for i in range(4):
            epsilon[i] = (math.pi/6)*ro*sum(x*self.mi*(di**i))

        #Energy and volume of association
        delta = []
        for i in range(num_comp):
            for j in range (num_comp):
                if j>=i:
                    dij     = (di[i]*di[j])/(di[i]+di[j])
                    Pa      = 1/(1-epsilon[3])
                    Pb      = 3*dij
                    Pc      = epsilon[2]/((1-epsilon[3])**2)
                    Pd      = 2*dij**2
                    Pe      = (epsilon[2]**2)/((1-epsilon[3])**3)
                    gij_seg = Pa+(Pb*Pc)+(Pd*Pe)
                    eaibjk  = (self.eabk[i]+self.eabk[j])*0.5
                    Pa      = math.sqrt(self.kab[i]*self.kab[j])
                    Pb      = math.sqrt(self.sigma[i]*self.sigma[j])
                    Pc      = 0.5*(self.sigma[i]+self.sigma[j])
                    kaibj   = Pa*(Pb/Pc)**3
                    Pa      = np.exp(eaibjk/T)-1
                    delta.append(gij_seg*(dij**3)*Pa*kaibj)

        #Creating vectors for combination sites
        Delta = np.zeros(len(code))
        for i, num in enumerate(code):
            if num == 0:
                Delta[i] = 0
            else:
                Delta[i] = delta[num-1]

        return Delta

    def X_assoc(self, num_sites, eta, T, x, code):
        """
        * Calculation of association Helmholtz free energy term, using successive substitution method.
            :param num_sites: number of sites for each component
            :param ro_molar: molar density of mixture
            :param T: temperature
            :param x: molar fraction of each component in the mixture
            :param code: list of labels for each pairwise association
        """

        #Creating the energy of association matrix
        delta  = self.association_parameters(eta, T, x, code)
        cont   = 0
        row    = sum(num_sites)
        column = sum(num_sites)
        matrix = np.zeros((row, column))

        for i in range(row):
            for j in range(column):
                if j > i:
                    matrix[i, j] = delta[cont]
                    cont += 1
                else:
                    matrix[i, j] = 0

        matrix_T = np.transpose(matrix)
        Delta    = matrix + matrix_T

        #Initial guess for X
        X          = 0.1*np.ones(sum(num_sites))

        # Calculation of molar density of molecules
        di         = self.sigma*(1-0.12*np.exp(-3*(self.Ek/T)))
        ro         = (6/math.pi)*eta*(sum(x*self.mi*(di**3)))**-1
        ro_comp    = ro*x
        Ro         = np.array([])
        for i, n_sites in enumerate(num_sites):
            ro     = np.array([ro_comp[i]]*n_sites)
            Ro     = np.concatenate((Ro,ro), axis=0)


        #Newton-Raphson method
        error = 1
        tol = 1e-10
        f = []

        while error > tol:

            #Jacobian matrix
            jacobian = np.zeros((row, column))
            for i in range(row):
                for j in range(column):
                    if i == j:
                        jacobian[i,j] = 1+sum(Ro*X*Delta[i,:])
                    else:
                        jacobian[i,j] = Ro[j]*X[i]*Delta[i,j]

                f.append(X[i]+X[i]*sum(Ro*X*Delta[i,:])-1)


            #Pseudoinverse
            J_pseudo = np.linalg.inv(np.transpose(jacobian)*jacobian)*np.transpose(jacobian)

            X_mod = np.reshape(X,(row,1))

            #Vector of functions
            F = np.reshape(np.array(f),(row,1))
            f.clear()

            #Actualization of X
            X_new = X_mod - np.dot(J_pseudo,F)
            error = sum((X_new-X_mod)**2)
            X     = np.reshape(X_new,(row,))

        return X

    def ass_parameters(self, x, code):

        eta, T   = symbols('eta T ')
        num_comp = len(self.mi)
        di       = []
        gij      = []
        soma_0   = 0
        soma_1   = 0
        soma_2   = 0
        soma_3   = 0
        for i in range(len(self.Ek)):
            di.append(self.sigma[i] * (1 - 0.12 * exp(-3 * (self.Ek[i] / T))))
            soma_0 += x[i] * self.mi[i] * (di[i] ** 0)
            soma_1 += x[i] * self.mi[i] * (di[i] ** 1)
            soma_2 += x[i] * self.mi[i] * (di[i] ** 2)
            soma_3 += x[i] * self.mi[i] * (di[i] ** 3)

        for i in range(len(x)):
            ro       = (6/math.pi)*eta*soma_3**-1
            epsilon2 = (math.pi / 6) * ro * soma_2
            epsilon3 = (math.pi / 6) * ro * soma_3

            # Energy and volume of association
            delta = []
            for i in range(num_comp):
                for j in range(num_comp):
                    if j >= i:
                        dij     = (di[i]*di[j])/(di[i]+di[j])
                        Pa      = 1/(1-epsilon3)
                        Pb      = 3*dij
                        Pc      = epsilon2/((1-epsilon3)**2)
                        Pd      = 2*dij**2
                        Pe      = (epsilon2**2)/((1-epsilon3)**3)
                        gij_seg = Pa+(Pb*Pc)+(Pd*Pe)
                        eaibjk  = (self.eabk[i]+self.eabk[j])*0.5
                        Pa      = math.sqrt(self.kab[i]*self.kab[j])
                        Pb      = math.sqrt(self.sigma[i]*self.sigma[j])
                        Pc      = 0.5 * (self.sigma[i]+self.sigma[j])
                        kaibj   = Pa*(Pb/Pc)**3
                        Pa      = exp(eaibjk/T)-1
                        delta.append(gij_seg*(dij**3)*Pa*kaibj)

            # Creating vectors for combination sites
            Delta = [0]*len(code)
            for i, num in enumerate(code):
                if num == 0:
                    Delta[i] = 0
                else:
                    Delta[i] = delta[num-1]

            return Delta

    def a_association(self, num_sites, eta1, T1, x, code2, code1):
        """

        :param num_sites:
        :param eta1:
        :param T1:
        :param x:
        :param code2: Complete association combinations
        :param code1: Not complete association combinations
        :return:
        """
        eta, T, = symbols('eta T ')

        delta = self.ass_parameters(x, code2)

        # Calculation of molar density of molecules
        soma = 0
        di   = []
        for i in range(len(self.Ek)):
            di.append(self.sigma[i]*(1-0.12*exp(-3*(self.Ek[i]/T))))
            soma += x[i] * self.mi[i] * (di[i] ** 3)

        ro = (6/math.pi)*eta*soma**-1
        ro_comp = ro*x
        Ro      = np.array([])
        X       = self.X_assoc(num_sites, eta1, T1, x, code1)

        for i, n_sites in enumerate(num_sites):
            ro  = np.array([ro_comp[i]]*n_sites)
            Ro  = np.concatenate((Ro, ro), axis=0)

        row  = sum(num_sites)
        resp = []
        cont = 0
        for j in range(row):
            soma = 0
            for i in range(row):
                soma += Ro[i]*X[i]*delta[cont]
                cont += 1
            resp.append((1+soma)**-1)

        # Association Helmholtz free energy
        a_association = 0
        for i, n_sites in enumerate(num_sites):
            Pa = 0
            for j in range(n_sites):
                Pa += log(resp[j]) - (0.5 * resp[j])
            a_association += x[i]*(Pa+(0.5*n_sites))

        return a_association


    def a_res(self, num_sites, eta1, T1, x, code2, code1, association):

        a_hc   = self.a_hc(x)
        a_disp = self.a_disp(x)

        if association == 0:
            a_res  = a_hc+a_disp

        elif association == 1:
            a_association = self.a_association(num_sites, eta1, T1, x, code2, code1)
            a_res = a_hc+a_disp+a_association

        return a_res


    def compressibility(self, num_sites, eta1, T1, x, code2, code1, association):
        eta, T = symbols('eta T')
        di = []
        soma = 0
        a_res = self.a_res(num_sites, eta1, T1, x, code2, code1, association)
        dadn = diff(a_res, eta)
        Z = 1 + (eta * dadn)
        resp = Z.subs([(eta, eta1), (T, T1)])
        return resp

    def pressure(self, num_sites, eta1, T1, x, code2, code1, association):

        eta, T = symbols('eta T')
        di     = []
        soma   = 0
        a_res  = self.a_res(num_sites, eta1, T1, x, code2, code1, association)
        dadn   = diff(a_res, eta)
        Z      = 1+(eta*dadn)
        for i in range(len(self.Ek)):
            di.append(self.sigma[i]*(1-0.12*exp(-3*(self.Ek[i]/T))))
            soma += x[i]*self.mi[i]*(di[i]**3)

        for i in range(len(x)):
            ro = (6/math.pi)*eta*soma**-1

        P = Z*self.k*T*ro*1e30  #Pascal
        dpdeta = diff(P, eta, 2).subs([(eta, eta1), (T, T1)])
        print(dpdeta)
        return P

    def cv(self, num_sites, eta1, T1, x, code2, code1, association):

        eta, T = symbols('eta T')
        a_res    = self.a_res(num_sites, eta1, T1, x, code2, code1, association)
        dadt     = diff(a_res,T).subs([(eta,eta1),(T,T1)])
        d2adt2   = diff(a_res,T,2).subs([(eta,eta1),(T,T1)])
        cv_res   = (-self.R*d2adt2*T1**2)+(-2*self.R*T1*dadt)
        cp_ideal = sum(((self.A+(self.B*T1)+(self.C*T1**2)+(self.D/(T1**2)))*self.R))
        cv_ideal = cp_ideal-self.R
        cv       = cv_ideal+cv_res
        return cv

    def cp(self, num_sites, eta1, T1, x, code2, code1, association):

        eta, T   = symbols('eta T')
        P        = self.pressure(num_sites, eta1, T1, x, code2, code1, association)
        dpdt     = diff(P, T).subs([(eta, eta1), (T, T1)])
        dpdeta   = diff(P, eta).subs([(eta, eta1), (T, T1)])
        di       = self.sigma*(1-0.12*np.exp(-3*(self.Ek/T1)))
        detadro  = ((math.pi/6)*self.Na*(sum(x*self.mi*(di**3))))*1e-30
        dpdro    = dpdeta * detadro

        ro       = (6/math.pi)*eta1*(sum(x*self.mi*(di**3)))**-1
        density  = ro/(self.Na*1e-30)
        Pa       = T1/(density**2)
        cv       = self.cv(num_sites, eta1, T1, x, code2, code1, association)
        cp       = cv+(Pa*((dpdt**2)/dpdro))

        return cp


    def speed_of_sound(self, molar_mass, num_sites, eta1, T1, x, code2, code1, association):

        eta, T  = symbols('eta T')

        P       = self.pressure(num_sites, eta1, T1, x, code2, code1, association)
        dpdeta  = diff(P,eta).subs([(eta, eta1), (T, T1)])
        di      = self.sigma*(1-0.12*np.exp(-3*(self.Ek/T1)))
        detadro = ((math.pi/6)*self.Na*(sum(x*self.mi*(di**3))))*1e-30
        dpdro   = dpdeta*detadro

        Mw      = sum(x*molar_mass)
        cp      = self.cp(num_sites, eta1, T1, x, code2, code1, association)
        cv      = self.cv(num_sites, eta1, T1, x, code2, code1, association)
        gamma   = cp/cv
        u       = math.sqrt((gamma/Mw)*dpdro*1e3)

        print(f"Isocoric heat: {cv:.2f} J/mol*K")
        print(f"Isobaric heat: {cp:.2f} J/mol*K")
        print(f"Speed of sound: {u:.2f} m/s")



"""
################################ METANO/ETANO/PROPANO ###########################
R          = 8.314 #----------------------------------------------------> J/K*mol
T          = 233.15 #---------------------------------------------------> K
sigma      = np.array([3.7039, 3.5206, 3.6184])  #----------------------> A°
Ek         = np.array([150.03, 191.42, 208.11]) #----------------------->
mi         = np.array([1.0, 1.6069, 2.0020]) #-------------------------->
Kij        = np.array([3.0e-4, 1.15e-2, 5.10e-3]) #--------------------------------------------> binary interaction parameters
x          = np.array([0.1, 0.3, 0.6]) #--------------------------------> molar fraction
eta        = 0.402 #---------------------------------------------------->
A          = np.array([5.457]) #----------------------------------------> Heat capacity coefficient
B          = np.array([1.045e-3]) #-------------------------------------> Heat capacity coefficient
C          = np.array([0]) #--------------------------------------------> Heat capacity coefficient
D          = np.array([-1.157e5]) #-------------------------------------> Heat capacity coefficient
P = 5.022e7

#resp = PCSAFT(R, T, sigma, Ek, mi, 0, 0, Kij, A, B, C, D).coef_fugacity(P, T, x, 0, 0, 0)
resp = PCSAFT(R, T, sigma, Ek, mi, 0, 0, Kij, A, B, C, D).a_res(eta, T, x, 0, 0, 0)
print(resp)"""



"""
########################################## CO2 ###############################################
R          = 8.314 #----------------------------------------------------> J/K*mol
T          = 300 #------------------------------------------------------> K
sigma      = np.array([2.7852])  #---------------------------------------> A°
Ek         = np.array([169.21]) #--------------------------------------->
mi         = np.array([2.0729]) #--------------------------------------->
Kij        = np.array([0]) #--------------------------------------------> binary interaction parameters
x          = np.array([1]) #--------------------------------------------> molar fraction
A          = np.array([5.457]) #----------------------------------------> Heat capacity coefficient
B          = np.array([1.045e-3]) #-------------------------------------> Heat capacity coefficient
C          = np.array([0]) #--------------------------------------------> Heat capacity coefficient
D          = np.array([-1.157e5]) #-------------------------------------> Heat capacity coefficient
molar_mass = np.array([44.01]) #----------------------------------------> Molar mass
density    = 651.53 #----------------------------------------> Molar density (mol/m³)
di         = sigma*(1-0.12*np.exp(-3*(Ek/T)))
ro         = density*6.0221e-7
eta        = (math.pi/6)*ro*(sum(x*mi*(di**3)))
#resp       = PCSAFT(R, T, sigma, Ek, mi, 0, 0, Kij, A, B, C, D).speed_of_sound(molar_mass,eta,T,x,0,0,0)
resp       = PCSAFT2(R, T, sigma, Ek, mi, 0, 0, Kij, A, B, C, D).speed_of_sound(molar_mass, eta, T, x)"""


"""
#################################### METANO #######################################
R          = 8.314 #----------------------------------------------------> J/K*mol
T          = 92 #------------------------------------------------------> K
sigma      = np.array([3.7039])  #---------------------------------------> A°
Ek         = np.array([150.03]) #--------------------------------------->
mi         = np.array([1.0]) #--------------------------------------->
Kij        = np.array([0]) #--------------------------------------------> binary interaction parameters
x          = np.array([1]) #--------------------------------------------> molar fraction
A          = np.array([5.457]) #----------------------------------------> Heat capacity coefficient
B          = np.array([1.045e-3]) #-------------------------------------> Heat capacity coefficient
C          = np.array([0]) #--------------------------------------------> Heat capacity coefficient
D          = np.array([-1.157e5]) #-------------------------------------> Heat capacity coefficient
molar_mass = np.array([44.01]) #----------------------------------------> Molar mass
density    = 28074
di         = sigma*(1-0.12*np.exp(-3*(Ek/T)))
ro         = density*6.0221e-7
eta        = (math.pi/6)*ro*(sum(x*mi*(di**3)))
resp = PCSAFT2(R, T, sigma, Ek, mi, 0, 0, Kij, A, B, C, D).pressure(0, eta, T, x, 0, 0, association=0)
print(resp)"""



"""
############################################# ÁGUA ################################################
R          = 8.314 #----------------------------------------------------------> J/K*mol
T          = 800 #------------------------------------------------------------> Temperature (K)
sigma      = np.array([3.0007])  #--------------------------------------------> Diameter of molecule (A°)
Ek         = np.array([366.51]) #--------------------------------------------->
mi         = np.array([1.0656]) #--------------------------------------------->
Eabk       = np.array([2500.7]) #--------------------------------------------->
Kab        = np.array([0.034868]) #------------------------------------------->
Kij        = np.array([0]) #--------------------------------------------------> Binary interaction parameters
x          = np.array([1]) #--------------------------------------------------> Molar fraction
A          = np.array([3.470]) #----------------------------------------------> Heat capacity coefficient
B          = np.array([1.450e-3]) #-------------------------------------------> Heat capacity coefficient
C          = np.array([0]) #--------------------------------------------------> Heat capacity coefficient
D          = np.array([0.121e5]) #--------------------------------------------> Heat capacity coefficient
molar_mass = np.array([18.01]) #----------------------------------------------> Molar mass
density    = np.array([12105]) #----------------------------------------------------------> Molar density (mol/m³)
code       = [0,1,1,1,1,0] #--------------------------------------------------> Code of sites
num_sites  = np.array([4]) #--------------------------------------------------> Number of sites for each compound
SAFT       = PCSAFT(R, T, sigma, Ek, mi, Eabk, Kab, Kij, A, B, C, D)
eta        = SAFT.density_to_eta(x,density)
resp       = SAFT.speed_of_sound(molar_mass, eta, T, x, num_sites, code, association=1)
"""





############################################# ÁGUA+CO2 ################################################
#COMPOSIÇÕES BASEADAS DO TRABALHO: "Speeds of Sound in Binary Mixtures of Water and Carbon Dioxide at Temperatures from
# 273 K to 313 K and at Pressures up to 50 MPa"

R          = 8.314 #----------------------------------------------------------> J/K*mol
T          = 273.16 #-----------------------------------------------------------> Temperature (K)
sigma      = np.array([3.0007, 2.7852])  #--------------------------------------------> Diameter of molecule (A°)
Ek         = np.array([366.51, 169.21]) #--------------------------------------------->
mi         = np.array([1.0656, 2.0729]) #--------------------------------------------->
Eabk       = np.array([2500.7, 0]) #--------------------------------------------->
Kab        = np.array([0.034868, 0]) #------------------------------------------->
Kij        = np.array([0.0]) #--------------------------------------------------> Binary interaction parameters
x          = np.array([0.9882, 0.0118]) #-------------------------------------> Molar fraction
A          = np.array([3.470, 5.457]) #----------------------------------------------> Heat capacity coefficient
B          = np.array([1.450e-3, 1.045e-3]) #-------------------------------------------> Heat capacity coefficient
C          = np.array([0, 0]) #--------------------------------------------------> Heat capacity coefficient
D          = np.array([0.121e5, -1.157e5]) #--------------------------------------------> Heat capacity coefficient
molar_mass = np.array([18.01, 44.01]) #----------------------------------------------> Molar mass
density    = np.array([55525, 473.43]) #----------------------------------------------> Molar density (mol/m³)
code       = [0,1,1,1,1,0] #--------------------------------------------------> Code of sites
num_sites  = np.array([4, 0]) #--------------------------------------------------> Number of sites for each compound
SAFT       = PCSAFT(R, T, sigma, Ek, mi, Eabk, Kab, Kij, A, B, C, D)
eta        = SAFT.density_to_eta(x,density)
resp       = SAFT.speed_of_sound(molar_mass, eta, T, x, num_sites, code, association=1)


"""
Components = np.array(['HFC 32', 'HFC 152a']) #---------> Components
R          = 8.314  #--------------------------------------------> Gas constant (J*mol-1*K-1)
mi         = np.array([1.730, 1.699])
sigma      = np.array([3.231, 3.687])
Ek         = np.array([149.988, 193.413])
T          = 323.15
eta        = 0.999999999998866
x          = np.array([0.4, 1-0.4])
A          = np.array([0, 0]) #----------------------------------------------> Heat capacity coefficient
B          = np.array([0, 0]) #-------------------------------------------> Heat capacity coefficient
C          = np.array([0, 0]) #--------------------------------------------------> Heat capacity coefficient
D          = np.array([0, 0])
Kij        = np.array([0])
SAFT       = PCSAFT2(R, T, sigma, Ek, mi, 0, 0, Kij, A, B, C, D).compressibility(0, eta, T, x, 0, 0, association=0)
print(SAFT)"""

