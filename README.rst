.. image:: https://github.com/PolyChord/PolyChordLite/workflows/CI/badge.svg?branch=master
   :target: https://github.com/PolyChord/PolyChordLite/actions?query=workflow%3ACI+branch%3Amaster
   :alt: Build Status
.. image:: https://img.shields.io/badge/arXiv-1506.00171-b31b1b.svg
   :target: https://arxiv.org/abs/1506.00171
   :alt: Open-access paper

PolyChord v 1.20.1

Will Handley, Mike Hobson & Anthony Lasenby

wh260@mrao.cam.ac.uk

arXiv:1502.01856

arXiv:1506.00171

Latest version Released Apr 2020


PolyChord Licence
=================

Users are required to accept the licence agreement given in LICENCE
file. PolyChord is free for academic usage

Users are also required to cite the PolyChord papers: 

- arXiv:1502.01856
- arXiv:1506.00171

in their publications.

Python quickstart
=================

For Python users in a hurry:

.. code:: bash

    pip install git+https://github.com/PolyChord/PolyChordLite@master
    wget https://raw.githubusercontent.com/PolyChord/PolyChordLite/master/run_pypolychord.py
    python run_pypolychord.py

You can then modify the file run_pypolychord.py to your needs. If you have mpi compilers available, this version can be run in parallel with mpi.

You should make sure that you have gfortran (or equivalent) fortran compilers installed. 

If any of the above steps fail (this can in general happen for certain Mac OSX versions), then try installing without pip:

.. code:: bash

    git clone https://github.com/PolyChord/PolyChordLite.git
    cd PolyChordLite
    python setup.py install

If you do not have sudo access/virtual environments/anaconda, then appending `--user` to the install command may be necessary.

Post Processing
===============

We recommend the tool `anesthetic <https://github.com/williamjameshandley/anesthetic>`_ for post-processing your nested sampling runs. A plot gallery can be found `here <http://htmlpreview.github.io/?https://github.com/williamjameshandley/cosmo_example/blob/master/demos/demo.html>`_


https://github.com/williamjameshandley/anesthetic

MPI Support
===========

The code is MPI compatible with openMPI. To disable the MPI parallelization, 
set MPI=0 in ./Makefile, or compile with

.. code::

    make <target>  MPI=0

Additional Libraries  
====================

PolyChord requires no additional libraries to run in linear mode
To run with MPI it requires the openMPI library


Compilers
=========

PolyChord compiles with both gfortran and intel compilers. 

Compiler type is chosen in the Makefile with the COMPILER_TYPE flag;

set
COMPILER_TYPE = gnu
for gfortran compilers (free)

set
COMPILER_TYPE = intel
for intel compilers (proprietary, much faster)


Running PolyChord
=================

Examples
--------
First, try a couple of quick examples:

1) 20 dimensional Gaussian

Run the commands:

.. code::

    $  make gaussian
    $  ./bin/gaussian ini/gaussian.ini

2) Rastrigin

Run the commands:

.. code::

    $ make rastrigin
    $ ./bin/rastrigin ini/rastrigin.ini

This runs the rastrigin 'bunch of grapes' loglikelihood.

In general, binary executables are stored in the directory ./bin, and ini files are
stored in the directory ./ini.

You can create new likelihoods by modelling them on the ones in
likelihoods/examples, and triggering them with their own ini files

Alternatively you can take a more "MultiNest" like approach, and manually
generate the prior transformations. PolyChord's settings are then modified in
the driver files src/drivers.


Fortran likelihoods
-------------------
You should place your likelihood code in the function loglikelihood and your
prior code in the function prior, contained in:

./likelihoods/fortran/likelihood.f90 

Any setup required (such as reading in input files) should be conducted in the
function setup_loglikelihood. In most cases, this will likely just be a call
to your own pre-written library.

You should then alter the polychord run-time settings within the driver file:

./src/drivers/polychord_fortran.f90

Your code can be compiled and run with the commands:

.. code::

    $  make polychord_fortran
    $  ./bin/polychord_fortran



C++/C likelihoods
-----------------
You should place your likelihood code in the function loglikelihood,
contained in 

./likelihoods/CC/CC_likelihood.cpp

Any setup required (such as reading in input files) should be conducted in the
function setup_loglikelihood.  In most cases, this will likely just be a call
to your own pre-written library.

You should then alter the polychord run-time settings within the driver file:

./src/drivers/polychord_CC.cpp

or use the ini file version:

./likelihoods/CC_ini/CC_ini_likelihood.cpp
./src/drivers/polychord_CC_ini.cpp

Your code can be compiled and run with the commands:

.. code::

    $  make polychord_CC
    $  ./bin/polychord_CC 

or

.. code::

    $  make polychord_CC_ini
    $  ./bin/polychord_CC_ini ini/gaussian_CC.ini

If you have an additional suggestions to make the c++ wrapper more easy to use, 
please email Will (wh260@mrao.cam.ac.uk).



Python likelihoods (pypolychord)
--------------------------------
Being python, this interface is the most self-explanatory. 
You can install direct from the git repository using:

.. code:: bash

    pip install https://github.com/PolyChord/PolyChordLite/archive/master.zip

(N.B. PyPi coming soon)
or you can install locally with the command:

.. code:: bash

   git clone https://github.com/PolyChord/PolyChordLite.git
   cd PolyChordLite
   pip install . --user

This has the advantage of using intel compilers if you have them (e.g. on a HPC machine). You may wish to consider installing pypolychord in a `virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments>`, in which case you don't need the --user argument.

Once installed, you can then import pypolychord from anywhere with the lines:

.. code:: python

   import pypolychord

and check that it's working by running:

.. code:: bash

    $  python run_pypolychord.py

or in MPI:

.. code:: bash

    $  mpirun -np 4 python run_pypolychord.py

If so, the rest of the interface is relatively painless. Follow the example in
run_pypolychord.py, and consult the docstring if you need help:

.. code:: python

    import pypolychord
    from pypolychord.settings import PolyChordSettings

    help(pypolychord.run_polychord)
    help(PolyChordSettings)

There is also a demo `python notebook <https://github.com/PolyChord/PolyChordLite/blob/master/run_pypolychord.ipynb>`_.

Output files 
=============
PolyChord produces several output files depending on which settings
are chosen


[root].stats
------------
Run time statistics

[root].resume
-------------
Files for resuming a stopped run. Semi-human readable.
This is produced if settings%write_resume=.true.
This is used if settings%read_resume=.true.

[root].txt
----------
File containing weighted posterior samples. Compatable with the format
required by getdist package which is part of the CosmoMC package.
Contains ndims+nderived+2 columns:

.. code::

    weight -2*loglike <params> <derived params>

Refer to the following website in order to download or get more
information about getdist:
http://cosmologist.info/cosmomc/readme.html#Analysing

If settings%cluster_posteriors=.true. there are additional cluster files in
clusters/[root]_<integer>.txt 

[root]_equal_weights.txt
------------------------
As above, but the posterior points are equally weighted. This is
better for 'eyeballing' the posterior, and provides a natural ~4 fold
compression of the .txt file. 


[root]_phys_live.txt
--------------------
Live points in the physical space. This is produced if
settings%write_phys_live=.true.
This file contains ndims+nderived+1 columns, indicating the physical
parameters, derived parameters and the log-likelihood. This is useful
for monitoring a run as it progresses. 

[root]_dead.txt
---------------
Points that have been killed off. This is produced if
settings%write_dead=.true.
This file contains ndims+nderived+1 columns, indicating the loglikelihood,
physical parameters, derived parameters and the log-likelihood. This is useful
for monitoring a run as it progresses, and for performing alternative
calculations and checks on evidence and posterior computations

[root].paramnames
-----------------
Parameter names file for compatibility with getdist


[root]phys_live-birth.txt & [root]dead-birth.txt 
------------------------------------------------

These can be used to reconstruct a full nested sampling run, as well as
simulate dynamic nested sampling.  The format & contents of these two files
are as follows: They have has ndims+nderived+2 columns. The first
ndims+nderived columns are the ndim parameter values along with the nderived
additional parameters that are being passed by the likelihood routine for
PolyChord to save along with the ndims parameters. The ndims+nderived+1 column
is the log-likelihood value.  The ndims+nderived+2 column is the log-likelihood
value that the point was born at. They are is identical to the
[root]_phys_live.txt and [root]_dead.txt file, except for an additional column
including the birth contours


Visualization of PolyChord Output:

[root].txt file created by PolyChord is compatable with the format
required by getdist package which is part of the CosmoMC package.
Refer to the following website in order to download or get more
information about getdist:
http://getdist.readthedocs.org/en/latest/


Common Problems & FAQs:


Run time Issues
===============

1 Output files ([root].txt & [root]_equal_weights.dat) files have very few (of order tens) points. 

These files only become populated as the algorithm approaches the peak(s) of the posterior. Wait for the run to be closer to finishing.

2 MPI doesn't help

* Currently, the MPI parallelisation will only increase speed for 
  'slow' likelihoods, i.e. likelihoods where the slice sampling step
  is the dominant computational cost (compared to the organisation of
  live points and clustering steps). 
* Parallelisation is only effective up to ncores~O(nlive).


Compilation Issues
==================
Most issues are usually one associated with an out-of-date MPI library or
fortran compiler. Ideally you should be using:

* gfortran 4.8    or    ifort 14
* openMPI 1.6.5   or    Intel MPI 4.1
