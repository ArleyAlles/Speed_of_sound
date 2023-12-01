The present program has the goal to calculate speed of sound for H2O+CO2 mixtures. The results are compared with experimental dataset of a paper of DHAKAL, et al. 2023

The structure of code is divided into:
	
	  - main.py
	  - EoS.py
	  - Data.py
	  - Parameters.py
	  - requirements.txt

--> The file **main.py** is the user-executable part of the program. The input parameters necessary for execution of program are:

	  - dataset: Which dataset will be used;
          - association: Use or not association contribution (0 --> No; 1 --> Yes).

The output of **main.py** code is:
	  - Calculated speed of sound for each experimental point (P,T) (given in m/s²);
	  - Average Relative Deviation (ARD%) between calculated and experimental speed of sound.


--> **Equation_of_state.py** is the file which Pc_SAFT(with association contibution) EoS remains. The user doesn't need to change this part of program.

--> **Data.py** contains all all experimental data (Temperature, Pressure , experimental speed of sound). If user wants to add a new
dataset, the inclusions MUST be done in the method **data** as the same way of data already entered. 

--> **Parameters.py** gathers all parameters for running PC-SAFT. If a different mixture (or compound) need to be evaluated, the user MUST to change the parameters. OBS:The *code* parameter treats about all the possible sites combinations of system for association contribution. An example for how to add such parameter is in the docstring of **association_parameters** at **EoS.py** file. 

--> At the end, *requirements.txt* is the list of all python extensions needed for running the code. To install such extensions use the command: pip install -r requirements.txt


REFERENCES:
	* [dataset] - Speeds of Sound in Binary Mixtures of Water and Carbon Dioxide at Temperatures from 273 K to 313 K and at Pressures up to 50 MPa;
	* [Parameters] - Modelling the phase equilibria of a H2O-CO2 mixture with PC-SAFT and tPC-PSAFT equations of state