from EoS import PCSAFT
from Parameters import Parameters
import numpy as np


############################################# ÁGUA+CO2 ################################################
#COMPOSIÇÕES BASEADAS DO TRABALHO: "Speeds of Sound in Binary Mixtures of Water and Carbon Dioxide at
# Temperatures from 273 K to 313 K and at Pressures up to 50 MPa"

class Data:

    def __int__(self):
        pass

    @staticmethod
    def data():

        T = np.array([273.17, 273.19, 273.20, 273.23, 273.14, 273.13, 273.13,
                      273.13, 273.13, 273.13, 273.13, 273.12, 273.13, 273.13,
                      273.13, 278.12, 278.19, 278.22, 278.19, 278.18, 278.14,
                      278.14, 278.17, 278.19, 278.18, 278.17, 278.17, 278.16,
                      278.16, 278.16, 278.15, 278.15, 278.15, 283.18, 283.18,
                      283.18, 283.20, 283.19, 283.19, 283.17, 283.18, 283.19,
                      283.19])

        P = 1e6*np.array([3.97, 5.97, 7.97, 9.97, 9.97, 11.95, 14.94, 17.93,
                          19.93, 24.92, 29.90, 34.89, 39.87, 44.85, 49.82,
                          3.63, 3.99, 5.99, 7.99, 9.98, 9.98, 11.98, 11.96,
                          14.93, 17.93, 19.92, 24.91, 27.90, 29.90, 34.88,
                          39.86, 44.84, 49.82, 3.99, 4.99, 7.99, 9.99, 11.95,
                          14.96, 17.93, 19.94, 24.90, 29.89])

        u_exp = np.array([1433.28, 1436.05, 1438.84, 1441.50, 1441.25, 1444.05,
                          1448.15, 1452.38, 1455.33, 1462.65, 1470.03, 1477.73,
                          1485.39, 1493.30, 1501.26, 1452.58, 1453.50, 1456.24,
                          1459.12, 1462.00, 1462.00, 1464.75, 1464.84, 1469.35,
                          1473.80, 1476.68, 1484.21, 1488.81, 1491.85, 1499.39,
                          1507.17, 1515.13, 1523.13, 1472.01, 1473.63, 1477.97,
                          1480.84, 1483.93, 1488.38, 1492.84, 1495.84, 1503.55,
                          1511.23])

        return T, P, u_exp
    @classmethod
    def ARD(cls):

        R, sigma, Ek, mi, Eabk, Kab, Kij, x, A, B, C, D, molar_mass, density, code, num_sites = Parameters().parameters()
        T, P, u_exp = cls.data()
        u_calc      = np.zeros(len(T))
        SAFT        = PCSAFT(R, sigma, Ek, mi, Eabk, Kab, Kij, A, B, C, D)

        for i in range(len(T)):
            u_calc[i] = SAFT.speed_of_sound(molar_mass, P[i], T[i], x, num_sites, code, association=0, specification=0)

        ARD = (100/len(T))*sum(abs((u_calc-u_exp)/u_exp))

        print(f'\nAverage relative deviation (ARD%): {ARD:.4f}')
