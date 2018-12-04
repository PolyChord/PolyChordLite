PolyChord Wheel Documentation
=============================

Please follow the documentation on this site to
install docker Community Edition (Docker CE):
https://docs.docker.com/install/linux/docker-ce/fedora/

    .. code:: bash

        $ sudo dnf -y install dnf-plugins-core
        $ sudo dnf config-manager --add-repo
            https://download.docker.com/linux/fedora/docker-ce.repo
        $ sudo dnf install docker-ce

Creating manylinux tags for Python 2.7 and 3.4+
===============================================
The following commands create the CentOS 5 64 & 32 bit docker images:

    .. code:: bash

        $ sudo docker pull quay.io/pypa/manylinux1_x86_64
        $ sudo docker pull quay.io/pypa/manylinux1_i686

**Note:** If you want to upload the generated wheels,
don't forget to update PolyChordLite's version number.

The following commands make the wheels,
add manlinux tags using auditwheel,
and test them.

    .. code:: bash

        $ cd <PolyChordLite directory>
        $ sudo docker run --rm -v `pwd`:/io quay.io/pypa/manylinux1_x86_64 /io/build.sh
        $ sudo docker run --rm -v `pwd`:/io quay.io/pypa/manylinux1_i686 /io/build.sh

The output texts are stored in wheelhouse/output and the
manylinux wheels are stored in wheelhouse/manylinux.

Manually check if wheels work:
------------------------------

You should first clean the cache and sites-package
of your Python installation using this script:

    .. code:: bash

        #!/usr/bin/env bash

        for PYLIB in /usr/lib64/python*; do
            rm -rf "${PYLIB}"/site-packages/pypolychord*
            rm -rf "${PYLIB}"/site-packages/_pypolychord*
        done

        rm -rf ~/.cache/Python-Eggs/*

        for PYLIB in ~/.local/lib/python*; do
            rm -rf "${PYLIB}"/site-packages/pypolychord*
            rm -rf "${PYLIB}"/site-packages/_pypolychord*
        done

        rm -rf .local/lib/python(version)/site-packages/.libs_pypolychord
        rm -rf.cache/python-eggs

We recommend using virtualenv to test the different python version.
Instructions are in here:
https://docs.python.org/3/tutorial/venv.html

Uploading to PyPi:
==================
Make sure twine is installed:

    .. code:: bash

        $ pip install twine

Upload your wheels with the following command:

    .. code:: bash

        $ twine upload wheelhouse/manylinux/*

Downloading from PyPi:
======================

    .. code:: bash

        $ pip install pypolychord --user
        $ python run_pypolychord.py


