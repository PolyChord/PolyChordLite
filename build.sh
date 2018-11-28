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

wheelsPath="/io/wheelhouse"
manylinuxWheelsPath="${wheelsPath}/manylinux"
outputPath="${wheelsPath}/output"

mkdir -p $wheelsPath $manylinuxWheelsPath $outputPath

# Cleans and makes the shared libraries 
cd /io && make veryclean && make

# Build a wheel for each version of Python
for PYBIN in /opt/python/*/bin; do 
    "${PYBIN}/pip" install -r /io/requirements.txt
    # wheels are stored in the '/io/wheelhouse' directory
    cd /io && "${PYBIN}/python" setup.py bdist_wheel -d $wheelsPath
done

# Bundle external shared libraries into the wheels
for whl in $wheelsPath/py*$(uname -m).whl; do
    auditwheel repair "$whl" -w $manylinuxWheelsPath
done

cd /io && make veryclean

# Test all the wheels

for PYBIN in /opt/python/*/bin/; do
    # Split the path to get the python version ($PYPATH[3])
    IFS='/' read -ra PYPATH <<< "$PYBIN"

    # Install and run all the weels with run_pypolychord
    "${PYBIN}/pip" install pypolychord --no-index -f $manylinuxWheelsPath
    cd "$HOME" && "${PYBIN}/python" /io/run_pypolychord.py > "${outputPath}/${PYPATH[3]}"_output_$(uname -m).txt
done
