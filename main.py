"""
Authors: Arley A. Cruz

    Some basic information:

    * The aim of this program is calculate speed of sound for H2O+CO2 using PC-SAFT with association

    * The user must insert some input parameters in order to obtain the results. Such parameters are:
          - dataset: Which dataset will be used
          - association: Use or not association contribution (0 --> No; 1 --> Yes)

          OBS: The experimental data used in this code are found in the paper: "Speeds of Sound in Binary Mixtures of Water
          and Carbon Dioxide at Temperatures from 273 K to 313 K and at Pressures up to 50 MPa". The conditions
          are:
                - x(H20) --------> 0.9882; 0.9934; 0.9985
                - Temperature ---> 273 K to 313 K
                - Pressure ------> 4 MPa to 50 MPa

    * The output of such code is the calculated speed of sound for each experimental point (P,T) and the
    Average Relative Deviation (ARD%)
"""

from Data import Data
from Parameters import Parameters

dataset     = 0
association = 1
Parameters().parameters(dataset)
Data().ARD(dataset=dataset, association=association)
