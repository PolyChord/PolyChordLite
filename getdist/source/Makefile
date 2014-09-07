
OUTPUT_DIR ?= Release

BUILD ?= MPI
#set BUILD to MPI to force MPI, should be set in ../Makefile

ifortErr = $(shell which ifort >/dev/null; echo $$?)

#these settings for ifort 13 and higher
ifeq "$(ifortErr)" "0"
#ifort; Can remove -xHost if your cluster is not uniform, or specify specific processor optimizations -x...
F90C     = ifort
MPIF90C = mpif90
FFLAGS = -mkl -openmp -O3 -xHost -no-prec-div -fpp
DEBUGFLAGS = -mkl -openmp -g -check all -check noarg_temp_created -traceback -no-prec-div -fpp 
#add -fpe0 to check for floating point errors (think lowLike also throws these harmlessly)
MODOUT = -module $(OUTPUT_DIR)
LAPACKL =
else
#gfortran does not yet work yet, open GCC 4.9 compiler bugs are
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=60414
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=60334
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=60289
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=60255
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=51976
#and preferably
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=44672
F90C = gfortran
MPIF90C = mpif90
FFLAGS = -cpp -O3 -ffast-math
DEBUGFLAGS = -cpp -Og -g -fbounds-check -fbacktrace
MODOUT =  -J$(OUTPUT_DIR)
#LAPACKL = -Wl,-framework -Wl,accelerate
#commented above is (I think) for Mac; this is standard linux (sudo apt-get install liblapack-dev)
LAPACKL = -lblas -llapack
endif

#COSMOS: use "module load cosmolib latest"
#use "runCosmomc" (globally installed) to run, defining required memory usage
ifeq ($(COSMOHOST),cosmos)
F90C = ifort
FFLAGS = -openmp -fast -w -fpp2
LINKMPI = -lmpi
LAPACKL = -mkl=sequential -lmkl_lapack
cfitsio = $(CFITSIO)
WMAP = $(COSMOLIB)/WMAP9
GSLPATH = $(GSL_ROOT)
endif

ifeq ($(COSMOHOST),darwin)
WMAP = /scratch/aml1005/wmap_likelihood_v5
cfitsio = /usr/local/Cluster-Apps/cfitsio/intel/3.300
endif

ifeq ($(OUTPUT_DIR),Debug)
FFLAGS=$(DEBUGFLAGS)
endif

GETDISTFLAGS = $(FFLAGS)
F90FLAGS = $(FFLAGS) 
LINKFLAGS = $(LAPACKL) $(F90CRLINK)

UTILS = $(OUTPUT_DIR)/MiscUtils.o $(OUTPUT_DIR)/ArrayUtils.o $(OUTPUT_DIR)/StringUtils.o $(OUTPUT_DIR)/MpiUtils.o $(OUTPUT_DIR)/FileUtils.o \
    $(OUTPUT_DIR)/RandUtils.o $(OUTPUT_DIR)/ObjectLists.o $(OUTPUT_DIR)/Interpolation.o

DISTFILES = $(UTILS) $(OUTPUT_DIR)/IniObjects.o \
    $(OUTPUT_DIR)/ObjectParamNames.o $(OUTPUT_DIR)/Matrix_utils_new.o \
    $(OUTPUT_DIR)/settings.o $(OUTPUT_DIR)/samples.o $(OUTPUT_DIR)/IO.o $(OUTPUT_DIR)/GetDist.o



$(OUTPUT_DIR)/RandUtils.o: $(OUTPUT_DIR)/MpiUtils.o
$(OUTPUT_DIR)/FileUtils.o: $(OUTPUT_DIR)/MpiUtils.o $(OUTPUT_DIR)/MiscUtils.o $(OUTPUT_DIR)/StringUtils.o
$(OUTPUT_DIR)/ObjectLists.o: $(OUTPUT_DIR)/FileUtils.o
$(OUTPUT_DIR)/IniObjects.o: $(OUTPUT_DIR)/FileUtils.o
$(OUTPUT_DIR)/ObjectParamNames.o: $(UTILS)
$(OUTPUT_DIR)/settings.o:  $(OUTPUT_DIR)/IniObjects.o $(OUTPUT_DIR)/ObjectParamNames.o
$(OUTPUT_DIR)/propose.o: $(OUTPUT_DIR)/settings.o $(OUTPUT_DIR)/Matrix_utils_new.o
$(OUTPUT_DIR)/IO.o: $(OUTPUT_DIR)/settings.o
$(OUTPUT_DIR)/samples.o: $(OUTPUT_DIR)/ObjectLists.o $(OUTPUT_DIR)/settings.o $(OUTPUT_DIR)/Matrix_utils_new.o
$(OUTPUT_DIR)/GetDist.o: $(OUTPUT_DIR)/IO.o $(OUTPUT_DIR)/samples.o


export FFLAGS
export F90C
export OUTPUT_DIR

directories:
	mkdir -p $(OUTPUT_DIR)

$(OUTPUT_DIR)/%.o: %.c
	$(CC) $(GSLINC) -c $*.c -o $(OUTPUT_DIR)/$*.o

$(OUTPUT_DIR)/%.o: %.f90 Makefile
	$(F90C) $(F90FLAGS) -c $*.f90 -o $(OUTPUT_DIR)/$*.o

$(OUTPUT_DIR)/%.o: %.F90 Makefile
	$(F90C) $(F90FLAGS) -c $*.F90 -o $(OUTPUT_DIR)/$*.o

$(OUTPUT_DIR)/fpico_interface.o: $(PICO)/fpico_interface.f90 Makefile
	$(F90C) $(F90FLAGS) -c $(PICO)/fpico_interface.f90 -o $(OUTPUT_DIR)/fpico_interface.o

getdist: directories $(DISTFILES)
	$(F90C) -o ../getdist $(DISTFILES) $(LAPACKL) $(GETDISTFLAGS)

clean:
	rm -f $(OUTPUT_DIR)/*.o $(OUTPUT_DIR)/*.mod ./Debug/* ./Release/*  ./DebugMPI/* ./ReleaseMPI/* ../core
