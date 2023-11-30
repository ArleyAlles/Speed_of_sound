from Data import Data
from Parameters import Parameters


set           = 0
association   = 1
specification = 0
Parameters().parameters(set)
Data().ARD(set=set, association=association, specification=specification)

