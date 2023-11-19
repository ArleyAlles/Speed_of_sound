import numpy as np

density     = 100
n_particles = 864
molar_mass  = 16.0427
avogadro    = 6.0221409e23
rx          = np.zeros(n_particles+1)
ry          = np.zeros(n_particles+1)
rz          = np.zeros(n_particles+1)

#Converting molar mass from g/mol to Kg/mol
molar_mass = molar_mass*1e-3

#Size of unit cell (Angstroms)
unit_cell_size = ((n_particles*molar_mass/density/avogadro)**(1/3))*1e10

num = int(round((n_particles)/4,0))

unit_cell = int(round(num**(1/3),0))

particle = 1
for i in range(unit_cell):
    for j in range(unit_cell):
        for k in range(unit_cell):

            #(0,0,0)
            rx[particle] = unit_cell_size*i
            ry[particle] = unit_cell_size*j
            rz[particle] = unit_cell_size*k
            particle += 1

            # (1,1,0)
            rx[particle] = unit_cell_size * (i+0.5)
            ry[particle] = unit_cell_size * (j+0.5)
            rz[particle] = unit_cell_size * k
            particle += 1

            # (1,0,1)
            rx[particle] = unit_cell_size * (i+0.5)
            ry[particle] = unit_cell_size * j
            rz[particle] = unit_cell_size * (k+0.5)
            particle += 1

            # (0,1,1)
            rx[particle] = unit_cell_size * i
            ry[particle] = unit_cell_size * (j+0.5)
            rz[particle] = unit_cell_size * (k+0.5)
            particle += 1
