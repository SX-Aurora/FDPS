# Samples and Makefile
SPH sample and N-body sample are available.

Executable of SHP is "sph.out" and the one of N-body is "nbody.out".
Additionally, shared library file "libvesph.so" or "libvenbody.so" is necessary for VEOffload.

Each sample includes the following Makefiles.

```
Makefile_g++: For C++ sample with g++ compiler.
Makefile_mpic++: For C++ sample with gcc MPI compiler.
Makefile_icpc: For C++ sample with intel compiler.
Makefile_gfortran: For Fortran sample with gfortran compiler.
Makefile_mpifort: For Fortran sample with gcc MPI compiler.
Makefile_ifort: For Fortran sample with intel compiler.
Makefile_nfort: For shared library run on VE for Fortran sample with NEC compiler.
```

Invoke make command like,

```
$ make -f Makefile_g++
```

Made executable "sph.out" or "nbody.out" is thread parallelized with OpenMP.

To enable MPI parallelization, enable "PARTICLE_SIMULATOR_MPI_PARALLEL" and make with MPI compiler.
Executables "sph_mpi.out" and "nbody_mpi.out" made with "Makefile_mpic++" and "Makefile_mpifort" are MPI+OpenMP parallelized.

# OpenMP parallelization
To execute FDPS application with VEOffload, specify number fo threads with environment variable "OMP_NUM_THREADS" and "VE_OMP_NUM_THREADS=1".
Furthermore, specify VE number with environment variable "VE_NODE_NUMBER".

For example, to execute SPH sample with 8 threads on VE number 0.

```
$ export OMP_NUM_THREADS=8
$ export VE_OMP_NUM_THREADS=1
$ export VE_NODE_NUMBER=0
$ ./sph.out
```

To execute SPH sample with 16 threads on VE number 0 and 1.

```
$ export OMP_NUM_THREADS=16
$ export VE_OMP_NUM_THREADS=1
$ export VE_NODE_NUMBER=0-1
$ ./sph.out
```

# MPI+OpenMP parallelization
Environment variable specifications are same as OpenMP parallelization.
Utilized VE number has to be specified for each MPI process.

For example, to execute SPH sample made with gcc MPI with 8 threads x 2 MPI processes on VE number 0 and 1.

```
$ export OMP_NUM_THREADS=8
$ export VE_OMP_NUM_THREADS=1
$ export VE_NODE_NUMBER=0-1
$ /usr/mpi/gcc/openmpi-4.0.3rc4/bin/mpirun -np 2 sph_mpi.out
```

Discontinuous VE numbers can be specified like,

```
$ export VE_NODE_NUMBER=0,4
```

Use appropriate "mpirun" command depending on your environment.

# "n_group_limit" parameter
On N-body sample, number of "EPI", "EPJ" particles processed on each user defined interaction function call can be specified with "n_group_limit" parameter.

Call count is decreased and overhead of VE-VH data copy is mitigated if larger number is set on "n_group_limit". (i.e. "n_group_limit=512")

# How to offload with VEOffload
In default, 1 process is created on 1 VE and each process runs with 8 threads to process user defined interaction functions.

If you make with "-DVEO_METHOD2" option, 8 processes are created on 1 VE and each process run with 1 thread. In some cases this is faster.
