#!/bin/bash
set -e -x

# --------------------------------------------------
# Run this commads to build wheels :
# sudo docker pull quay.io/pypa/manylinux1_x86_64
# sudo docker run --rm -v `pwd`:/io quay.io/pypa/manylinux1_x86_64 /io/build.sh
# --------------------------------------------------

wheelsPath="/io/wheelhouse"
manylinuxWheelsPath="${wheelsPath}/manylinux"

# Compile wheels

cd /io && make veryclean && make

for PYBIN in /opt/python/*/bin; do 
    "${PYBIN}/pip" install -r /io/requirements.txt
    # wheels are stored in the '/io/wheelhouse' directory
    cd /io && "${PYBIN}/python" setup.py bdist_wheel -d $wheelsPath
done

# Bundle external shared libraries into the wheels
for whl in $wheelsPath/py*.whl; do
    auditwheel repair "$whl" -w $manylinuxWheelsPath
done

# Test of all the wheels
outputPath="/io/output"

mkdir -p $outputPath
for PYBIN in /opt/python/*/bin/; do
    # Split the path
    IFS='/' read -ra PYPATH <<< "$PYBIN"

    # Install and run all the weels with run_pypolychord
    "${PYBIN}/pip" install pypolychord --no-index -f $manylinuxWheelsPath
    cd "$HOME" && "${PYBIN}/python" /io/run_pypolychord.py > "${outputPath}/${PYPATH[3]}"_output.txt
done
