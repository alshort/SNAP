#!/bin/sh
. ~/intel/bin/compilervars.sh -arch intel64 -platform linux
. ~/intel/mkl/bin/mklvars.sh intel64
. ~/intel/bin/iccvars.sh -arch intel64 -platform linux
. ~/intel/compilers_and_libraries/linux/mpi/intel64/bin/mpivars.sh 
