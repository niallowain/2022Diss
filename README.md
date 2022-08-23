# 2022 Dissertation Work

## Required Packages

`numpy`

## Brief Functions Overview

### adjmatgen.py

`adjmatgen.py` contains `adj_mat_gen` which generates a random adjacency matrix, for a given number of dimensions/Nodes

This is developed on with `temp_adj_mats` which develops a unique adjacenecy matrix for each time step

### montecarlosir.py

This file gives `montecarlosir` which returns the average probability for each node being S,I,R over a number of simulations


### models.py

Contains `models` which hopes to return for given nodes the trajectories of Temporal Pair Based and Contact Based Modelling

### updatedcp.py

Hopefully this will not need to be used, but contains the updated code for Koher et al (2017) which modelled Contact Based Model with an edge-based approach, it may be quicker but would need different inputs to other scripts provided here

Also, requires networkx for execution

