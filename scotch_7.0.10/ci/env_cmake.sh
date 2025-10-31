#!/bin/bash

case "${CC}" in
     icx )
       source /disk1/builds/intel/setvars.sh intel64 && export MPI_HOME=/disk1/builds/intel/mpi/latest/ ;;
     pgcc )
       source /disk1/builds/nvidia/hpc_sdk/load_env ;;
esac
