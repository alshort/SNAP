# SNAP: SN (Discrete Ordinates) Application Proxy

snap-c: C port of SNAP 1.0 produced by Intel Corporation


## Installation

* sudo apt-get install mpich
* sudo apt-get install libomp-dev
* gcc (8.3)

* Downloaded and install score-p
  * Followed instructions in README
  * ../configure --with-mpi=mpich3

* `sudo ln -s /usr/local/cuda/lib64/libcudart.so /usr/lib/libcudart.so`


## Running

## Locally
* `. intel.sh`
* `make clean all`


## Tinis
* module load intel
* module load intel impi


## Description

SNAP serves as a proxy application to model the performance of a modern discrete ordinates neutral particle transport application. SNAP may be considered an update to [Sweep3D](http://www.ccs3.lanl.gov/PAL/software.shtml), intended for hybrid computing architectures. It is modeled off the Los Alamos National Laboratory code PARTISN. PARTISN solves the linear Boltzmann transport equation (TE), a governing equation for determining the number of neutral particles (e.g., neutrons and gamma rays) in a multi-dimensional phase space. SNAP itself is not a particle transport application; SNAP incorporates no actual physics in its available data, nor does it use numerical operators specifically designed for particle transport. Rather, SNAP mimics the computational workload, memory requirements, and communication patterns of PARTISN. The equation it solves has been composed to use the same number of operations, use the same data layout, and load elements of the arrays in approximately the same order. Although the equation SNAP solves looks similar to the TE, it has no real world relevance.

The solution to the time-dependent TE is a "flux" function of seven independent variables: three spatial (3-D spatial mesh), two angular (set of discrete ordinates, directions in which particles travel), one energy (particle speeds binned into "groups"), and one temporal. PARTISN, and therefore SNAP, uses domain decomposition over these dimensions to coherently distribute the data and the tasks associated with solving the equation. The parallelization strategy is expected to be the most efficient compromise between computing resources and the iterative strategy necessary to converge the flux.

The iterative strategy is comprised of a set of two nested loops. These nested loops are performed for each step of a time-dependent calculation, wherein any particular time step requires information from the preceding one. No parallelization is performed over the temporal domain. However, for time-dependent calculations two copies of the unknown flux must be stored, each copy an array of the six remaining dimensions. The outer iterative loop involves solving for the flux over the energy domain with updated information about coupling among the energy groups. Typical calculations require tens to hundreds of groups, making the energy domain suitable for threading with the node's (or nodes') provided accelerator. The inner loop involves sweeping across the entire spatial mesh along each discrete direction of the angular domain. The spatial mesh may be immensely large. Therefore, SNAP spatially decomposes the problem across nodes and communicates needed information according to the KBA method. KBA is a transport-specific application of general parallel wavefront methods. Nested threads, spawned by the energy group threads, are available to use in one of two ways. Per one approach, nested threads may be used to further parallelize the work to sweep different energy groups assigned to a main-level thread. This option is still experimental and has only been implemented to work in the case of using a single MPI process. Alternatively, nested threads are used to perform "mini KBA" sweeps by concurrently operating on cells lying on the same diagonal of spatial sub-domains already decomposed across the distributed memory architecture (i.e., different MPI ranks). Lastly, although KBA efficiency is improved by pipelining operations according to the angle, current chipsets operate best with vectorized operations. During a mesh sweep, SNAP operations are vectorized over angles to take advantage of the modern hardware.

SNAP should be tested with problem sizes that accurately reflect the types of calculations PARTISN frequently handles. The spatial domain shall be decomposed to 2,000--4,000 cells per node (MPI rank). Each node will own all the energy groups and angles for that group of cells; typical calculations feature 10--100 energy groups and as few as 100 to as many as 2,000 angles. Moreover, sufficient memory must be provided to store two full copies of the solution vector for time-dependent calculations. The preceding parameters assume current trends in available per core memory. Significant advances or detriments affecting this assumption shall require reconsideration of appropriate parameters per compute node.

## Compilation

SNAP has been written to the Fortran 90/95 standard primarily. The retrieval of command line arguments, which contain file names, is handled with a standard Fortran 2003 intrinsic subroutine. It has been successfully built with, but not necessarily limited to, gfortran and ifort. Moreover, the code has been built with the profiling tool [Byfl](https://github.com/losalamos/byfl). The accompanying Makefile provides sample build options for gfortran and ifort. The build system depends on the availability of MPI. Both example builds assume the usage of mpif90 from an MPI installation. Builds may be selected by switching the COMPILER option in the Makefile or choosing one with the "make COMPILER=[]" command. The builds also assume the availability of OpenMP. Compiling SNAP without MPI or OpenMP will require modification to the source code to remove related subroutine calls and directives.

MPI implementations typically suggest using a "wrapper" compiler to compile the code. SNAP has been built and tested with OpenMPI and MPICH. OpenMPI allows one to set the underlying Fortran compiler with the environment variable OMPI_FC, where the variable is set to the (path and) compiler of choice, e.g., ifort, gfortran, etc.

The makefile currently is set up for several build options using different MPI wrappers and Fortran compilers. One example uses:

    FORTRAN = mpif90
    COMPILER = ifort

and testing has been performed with

    OMPI_FC = [path]/ifort

Fortran compilation flags can be set according to the underlying compiler. The current flags are set for the ifort compiler and using OpenMP for parallel threading.

    TARGET = isnap
    FFLAGS = -03 -[q]openmp -align array32byte -fp-model fast -fp-speculation fast -xHost
    FFLAG2 =

where `FFLAG2` is reserved for additional flags that may need applied differently, depending on the compiler. To make SNAP with these default settings, simply type

    make

on the command line within the SNAP directory.

A debugging version of SNAP can be built by typing

    make OPT=no

on the command line. The unoptimized, debugging version of SNAP features bounds checking, back-tracing an error, and the necessary debug compiler flags. With ifort, these flags appear as:

    FFLAGS = -O0 -[q]openmp -g -check bounds -traceback -warn unused
    FFLAG2 =

The values for these compilation variables have been modified for various Fortran compilers and the Makefile provides details of what has been used previously. These lines are commented out for clarity at this time and to ensure that changes to the build system are made carefully before attempting to rebuild with a different compiler.

The SNAP directory can be cleaned up of its module and object files if the user desires with:

    make clean

This removes all the `*.mod` and `*.o` files, as well as `*.bc` files from Byfl builds. Moreover, it will enforce complete recompilation of all files upon the next instance of `make` or `make OPT=no`. Currently, there is no separate directory for the compilation files of separate optimized and unoptimized builds. The user must do a `make clean` before building the code if the previous build used the opposite command.

Pre-processing has been added for the inclusion/exclusion of MPI and OpenMP. To build without MPI, OpenMP, or both, use the command lines, respectively:

    make MPI=no
    make OPENMP=no
    MAKE MPI=no OPENMP=no

Default make settings will build with MPI and OpenMP included. These options are further available with unpotimized settings, OPT=no.

Lastly, a line count report is generated with:

    make count

The line count report excludes blank lines and comments. It counts the number of code lines in all `*.f90` files and sums the results. The information is printed to the the `Lines` file.

## Usage

When SNAP is built with MPI, to execute SNAP, use the following command:

    mpirun -np [#] [path]/snap --fi [infile] --fo [outfile]

This command will automatically run with the number of threads specified by the input file, which is used to set the number of OpenMP threads, overwriting any environment variable to set the number of threads. Testing has shown that to ensure proper concurrency of work, the above command can be modified to

    mpirun -cpus-per-proc [#threads] -np [#procs] [path]/snap [infile] [outfile]

The command line is read for the input/output file names. If one of the names is missing, the code will not execute. Moreover, the output file overwrites any pre-existing files of the same name.

The specific command to invoke a run with MPI and the corresponding options may be dependent on the specific machine used to execute. Most notably, the "aprun" command is used on Cray systems.