# template for the Intel fortran compiler
# typical use with mkmf
# mkmf -t template.ifc -c"-Duse_libMPI -Duse_netCDF" path_names /usr/local/include
############
# commands #
############
FC = mpif90
CC = mpicc
CXX = CC
LD = mpif90
#########
# flags #
#########
DEBUG =
REPRO =
VERBOSE =
OPENMP =

##############################################
# Need to use at least GNU Make version 3.81 #
##############################################
need := 3.81
ok := $(filter $(need),$(firstword $(sort $(MAKE_VERSION) $(need))))
ifneq ($(need),$(ok))
$(error Need at least make version $(need).  Load module gmake/3.81)
endif

MAKEFLAGS += --jobs=8

NETCDF_ROOT = $(NETCDFDIR)
MPI_ROOT    = 
INCLUDE = -I$(NETCDFDIR)/include 

FPPFLAGS := -fpp -Wp,-w $(INCLUDE)

# -msse2 is added as a workaround for reproducibility on the c3 system.  We in the
# modeling systems group are looking for why this is needed to allow run-to-run
# reproducibility on the c3 system.
FFLAGS := -xHost -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -i4 -r8 -nowarn -sox -traceback $(INCLUDE)
FFLAGS_OPT = -O3 -debug minimal -fp-model source 
FFLAGS_DEBUG = -g -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -ftrapuv
FFLAGS_REPRO = -O2 -debug minimal -fp-model source
FFLAGS_OPENMP = -openmp
FFLAGS_VERBOSE = -v -V -what

CFLAGS := -D__IFC -xHost -sox -traceback
CFLAGS_OPT = -O3 -debug minimal
CFLAGS_REPRO = -O2 -debug minimal
CFLAGS_OPENMP = -openmp
CFLAGS_DEBUG = -O0 -g -ftrapuv

# Optional Testing compile flags.  Mutually exclusive from DEBUG, REPRO, and OPT
# *_TEST will match the production if no new option(s) is(are) to be tested.
FFLAGS_TEST = -O3 -debug minimal -fp-model source -override-limits
CFLAGS_TEST = -O2

LDFLAGS :=
LDFLAGS_OPENMP := -openmp
LDFLAGS_VERBOSE := -Wl,-V,--verbose,-cref,-M

# start with blank LIBS
LIBS :=

ifneq ($(REPRO),)
CFLAGS += $(CFLAGS_REPRO)
FFLAGS += $(FFLAGS_REPRO)
else ifneq ($(DEBUG),)
CFLAGS += $(CFLAGS_DEBUG)
FFLAGS += $(FFLAGS_DEBUG)
else ifneq ($(TEST),)
CFLAGS += $(CFLAGS_TEST)
FFLAGS += $(FFLAGS_TEST)
else
CFLAGS += $(CFLAGS_OPT)
FFLAGS += $(FFLAGS_OPT)
endif

ifneq ($(OPENMP),)
CFLAGS += $(CFLAGS_OPENMP)
FFLAGS += $(FFLAGS_OPENMP)
LDFLAGS += $(LDFLAGS_OPENMP)
# to correct a loader bug on gaea: envars below set by module load intel
#LIBS += -L$(INTEL_PATH)/$(INTEL_MAJOR_VERSION)/$(INTEL_MINOR_VERSION)/lib/intel64 -lifcoremt
endif

ifneq ($(VERBOSE),)
CFLAGS += $(CFLAGS_VERBOSE)
FFLAGS += $(FFLAGS_VERBOSE)
LDFLAGS += $(LDFLAGS_VERBOSE)
endif

ifeq ($(NETCDF),3)
  # add the use_LARGEFILE cppdef
  ifneq ($(findstring -Duse_netCDF,$(CPPDEFS)),)
    CPPDEFS += -Duse_LARGEFILE
  endif
endif

#ifneq ($(findstring netcdf-4.0.1,$(LOADEDMODULES)),)
#  LIBS += -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
#else
  LIBS += -L${NETCDFDIR}/lib64 -lnetcdf -lnetcdff -lhdf5
#endif

LIBS +=
LDFLAGS += $(LIBS)

#---------------------------------------------------------------------------
# you should never need to change any lines below.

# see the MIPSPro F90 manual for more details on some of the file extensions
# discussed here.
# this makefile template recognizes fortran sourcefiles with extensions
# .f, .f90, .F, .F90. Given a sourcefile <file>.<ext>, where <ext> is one of
# the above, this provides a number of default actions:

# make <file>.opt       create an optimization report
# make <file>.o         create an object file
# make <file>.s         create an assembly listing
# make <file>.x         create an executable file, assuming standalone
#                       source
# make <file>.i         create a preprocessed file (for .F)
# make <file>.i90       create a preprocessed file (for .F90)

# The macro TMPFILES is provided to slate files like the above for removal.

RM = rm -f
SHELL = /bin/csh -f
TMPFILES = .*.m *.B *.L *.i *.i90 *.l *.s *.mod *.opt

.SUFFIXES: .F .F90 .H .L .T .f .f90 .h .i .i90 .l .o .s .opt .x

