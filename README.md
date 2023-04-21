# 2D-ShwarzSolver-for-Poisson
This is the AD for the paper: Physics-informed Quantum Deep Neural Network for Solving PDEs.
The main usage of this code is to generate training data for the PIQDNN

The code employs Alternating FEM scheme to solve 2D Poisson problem.

Download the code and use make to compile executable.

To generate data there are two ways:
    1. mpirun -np 2 ./main Mesh_size Overlapping
    The overlaping can be set as 2. The code will generate a csv file
    2. bash ./run_DifferentMesh.sh
    The script will generate automatically up to 45 data file stored under ./outputfile
