# Running PRISMS-plasticity simulations from EBSD: an example
The MATLAB script named [Copper.m](Copper.m) illustrates the steps by step procedure used for converting EBSD into data suitable for PRISMS-Plasticity.

## Generate input data for PRISMS-Plasticity
All you have to do is running this script file to generate the files needed for running the crystal plasticity simulation, namely the mesh file and the orienation file.

## Running PRISMS-Plasticity
The PRISMS-Plasticity parameters are already set in [``prm.prm``](prm.prm) for simulating a tensile test along the x direction, up to 3% elongation. Just run
   
    path/to/prisms_binary prm.prm
 
 It is higly advised to run this simulation on multiple threads in order to speed it up. This can be done with MPI, e.g.:
 
     mpirun -np 8 path/to/prisms_binary prm.prm
     
 ---
 **NOTE**
 
This simulation takes about 5 hours on 32 threads on Intel Xeon Gold 6242 CPU @ 2.80GHz.
 
 ---
