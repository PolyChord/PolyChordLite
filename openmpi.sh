#!/bin/bash

mpiVersion="2.1.5"
mpiFile="openmpi-${mpiVersion}"

# Download
curl http://download.open-mpi.org/release/open-mpi/v2.1/${mpiFile}.tar.gz --output "/${mpiFile}.tar.gz"
# Unzip
tar -xf "/${mpiFile}.tar.gz"

cd "/${mpiFile}" && ./configure --prefix=/usr/local && make all install

/io/build.sh MPI

# Set path variables
# export PATH=${mpiInstall}/mpi/bin:$PATH
# export LD_LIBRARY_PATH=/usr/lib/mpi/lib:$LD_LIBRARY_PATH