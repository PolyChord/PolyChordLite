#!/bin/bash
set -e -x

# --------------------------------------------------
# Run these commands to downloads images :
# sudo docker pull quay.io/pypa/manylinux1_x86_64
# sudo docker pull quay.io/pypa/manylinux1_i686
#
# Run these commands to build the wheels :
# sudo docker run --rm -v `pwd`:/io quay.io/pypa/manylinux1_x86_64 /io/build.sh
# sudo docker run --rm -v `pwd`:/io quay.io/pypa/manylinux1_i686 /io/build.sh
# --------------------------------------------------

if ! hash mpicxx 2>/dev/null; then
    # Install openmpi
    mpiVersion="2.1.5"
    mpiFile="openmpi-${mpiVersion}"

    # Download
    curl http://download.open-mpi.org/release/open-mpi/v2.1/${mpiFile}.tar.gz --output "/${mpiFile}.tar.gz"
    # Unzip
    tar -xf "/${mpiFile}.tar.gz"

    cd "/${mpiFile}" && ./configure --prefix=/usr/local && make all install
fi

MPI=0
wheelsPath="/io/wheelhouse"
distPath="/io/dist"

for i in `seq 1 2`; do
    manylinuxWheelsPath="${wheelsPath}/manylinux"
    outputPath="${wheelsPath}/output"

    mkdir -p $distPath $wheelsPath $manylinuxWheelsPath $outputPath

    # Cleans and makes the shared libraries
    cd /io && make veryclean && make MPI=$MPI 

    # Build a wheel for each version of Python
    for PYBIN in /opt/python/*/bin; do 
        "${PYBIN}/pip" install -r /io/requirements.txt
        if [[ $MPI == 1 ]];
        then 
            "${PYBIN}/pip" install mpi4py
            # wheels are stored in the '/io/wheelhouse' directory
            cd /io && CC=mpicc CXX=mpicxx "${PYBIN}/python" setup.py bdist_wheel -d $wheelsPath mpi
        else
            # wheels are stored in the '/io/wheelhouse' directory
            cd /io && "${PYBIN}/python" setup.py bdist_wheel -d $wheelsPath
        fi
    done

    # Bundle external shared libraries into the wheels
    for whl in $wheelsPath/py*$(uname -m).whl; do
        auditwheel repair "$whl" -w $manylinuxWheelsPath -L pypolychord/libs
    done

    cp $manylinuxWheelsPath/*.whl $distPath

    # cd /io && make veryclean

    # Test all the wheels

    # for PYBIN in /opt/python/*/bin/; do
    #     # Split the path to get the python version ($PYPATH[3])
    #     IFS='/' read -ra PYPATH <<< "$PYBIN"

    #     # Install and run all the weels with run_pypolychord
    #     "${PYBIN}/pip" install pypolychord --no-index -f $manylinuxWheelsPath
    #     cd "$HOME" && "${PYBIN}/python" /io/run_pypolychord.py > "${outputPath}/${PYPATH[3]}"_$(uname -m)_output.txt
    # done
    
    MPI=1
    wheelsPath="/io/wheelhouse/mpi"
done


echo "Wheels builded successfully"
