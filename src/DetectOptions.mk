######################################################
## DO NOT TOUCH THESE UNLESS YOU'VE ADDED A NEW FEATURE
## FOR PUBLIC RELEASE
## @author Charles Lena
## @date 2017-01-15
###
####################################
# use some form of FFTW3 bindings
####################################
ifeq ($(strip $(WITH_FFTW3)),1)
DEF += -DUSEFFTW3
#####################################
# FFTW3 has f2003+ bindings required
#####################################
ifeq ($(strip $(WITH_FFTW3_MODERN)),1)
DEF += -DCREATE_FFTW3_MODULE
endif
endif


####################################
# use advanced OMP4 stuff
# @bug CML:Intel's provided libomp5 has a bug as of
# version Intel Compiler Suite 17.0.1 regarding
# omp task and leaking memory. This means our kernel
# will eventually destroy itself using -DGOODOMP
# unless we update to a manually patched version.
# Intel expects this to be patched in "update 1"
####################################
ifeq ($(strip $(WITH_OMP4)),1)
DEF += 
endif

#############################################################
# \note Default MPI uses neighborhood collectives now
# \bug CML: It's been confirmed that or normal MPI no longer functions
#  as of the new kernel
##############################################################
ifeq ($(strip $(WITH_MPI)),1)
DEF += -DMPI
endif

#############################################################
# User wants ARPACK
##############################################################
ifeq ($(strip $(WITH_ARPACK)),1)
DEF += -DUSEARPACK
else
DEF +=
endif

#############################################################
# User wants ELPA
#############################################################
ifeq ($(strip $(WITH_ELPA)),1)
DEF += -DUSE_ELPA
# require blacs and scalapack 
WITH_SCALAPACK=1
ifeq ($(ELPA_LIBS),)
$(warn "elpa libs unspecified")
endif
ELPA_INCLUDES=$(MODULE_IN)$(ELPA_MOD_PATH) -I$(ELPA_INC_PATH)
ifneq ($(ELPA_LIB_PATH),)
ELPA_LINK = -Wl,-rpath=$(ELPA_LIB_PATH) -L$(ELPA_LIB_PATH)
endif
ifneq ($(ELPA_LIBS),)
ELPA_LINK += $(ELPA_LIBS) 
endif
endif

#############################################################
# User wants scalapack
##############################################################
ifeq ($(strip $(WITH_SCALAPACK)),1)
DEF += -DUSE_SCALAPACK
ifneq ($(strip $(WITH_MPI)),1)
$(error "SCALAPACK requires MPI : Please set both WITH_SCALAPACK and WITH_MPI to 1")
endif
ifneq ($(SCALAPACK_LINK),)
ELPA_LINK += $(SCALAPACK_LINK)
endif
ifneq ($(BLACS_LINK),)
ELPA_LINK += $(BLACS_LINK)
endif
endif

#############################################################
# User wants MKL_VSL : defaults to ON
##############################################################
ifneq ($(strip $(WITH_MKL_VSL)),0)
ifeq ($(strip $(MKLROOT)),)
$(error "WITH_MKL_VSL=1 requires defining MKLROOT : See config/make.$(MACH)")
endif
DEF += -DUSEMKL_VSL -I$(MKLROOT)/include
else
DEF +=
endif

### Allows for using C + glibc to check memory
### Requires cluster-style-linux (RHEL) etc
### or any distribution that uses /proc/meminfo
ifeq ($(strip $(WITH_CAUX)),1)
DEF += -DWITH_CAUX
else
DEF +=
endif

### VTUNE
ifeq ($(strip $(WITH_VTUNE)),1)
ifeq ($(strip $(WITH_ITAC)),1)
$(error "Mututally exclusive options (ITAC) - VTUNE ALREADY DEFINED")
endif
ifeq ($(strip $(VTUNE_AMPLIFIER_DIR)),)
$(error "VTUNE_AMPLIFIER_DIR must be defined")
endif
DEF += -DWITH_VTUNE -I$(VTUNE_AMPLIFIER_DIR)/include
LIBVTUNE?=$(VTUNE_AMPLIFIER_DIR)/lib64/libittnotify.a
else
DEF +=
endif

### ITAC
ifeq ($(strip $(WITH_ITAC)),1)
ifeq ($(strip $(WITH_VTUNE)),1)
$(error "Mututally exclusive options (VTUNE) - ITAC ALREADY DEFINED")
endif
DEF += -DITAC -I$(VT_ROOT)/include 
LIBMPI+= -L$(VT_LIB_DIR) -lVT  $(VT_ADD_LIBS)
else
DEF +=
endif

### Allows for using C + glibc to check memory
### Requires cluster-style-linux (RHEL) etc
### or any distribution that uses /proc/meminfo
ifeq ($(strip $(WITH_SDE)),1)
DEF += -DWITH_SDE
else
DEF +=
endif

############################################################
# User wants libxc
##############################################################
ifeq ($(strip $(WITH_LIBXC)),1)
DEF += -DUSE_LIBXC
else
DEF +=
endif


############################################################
# Use of HDF5
##############################################################
ifeq ($(strip $(WITH_HDF5)),1)
DEF += -DUSE_HDF5
else
endif
