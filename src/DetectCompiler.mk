
pathF90 := $(shell which $(firstword $(F90)))
rF90 := $(notdir $(abspath $(pathF90)))

$(info DEBUG: pathF90 is $(pathF90))
$(info DEBUG: rF90 is $(rF90))

# ================================
# Detect "mpif90" and decontstruct
# ================================
ifeq ($(findstring mpif90,$(notdir $(realpath $(pathF90)))),mpif90)
	tF90:=$(firstword $(shell $(pathF90) -show))
	TAG_COMPILER:=mpif90-
endif

# ================================
# Detect "mpiifx" and decontstruct - special intel family wrappers
# ================================
ifeq ($(findstring mpiifx,$(notdir $(realpath $(pathF90)))),mpiifx)
	tF90:=$(firstword $(shell $(pathF90) --version))
	TAG_COMPILER:=$(MACH)-
endif

# ================================
# Detect "mpiifort" and decontstruct - special intel family wrappers
# ================================
ifeq ($(findstring mpiifort,$(notdir $(realpath $(pathF90)))),mpiifort)
	tF90:=$(firstword $(shell $(pathF90) --version))
	TAG_COMPILER:=$(MACH)-
endif

# ================================
# Detect "mpifort" and decontstruct - special intel family wrappers
# ================================
ifeq ($(findstring mpifort,$(notdir $(realpath $(pathF90)))),mpifort)
	tF90:=$(firstword $(shell $(pathF90) --version))
	TAG_COMPILER:=$(MACH)-
endif

# ========================================================================================
# Detect "ftn" and deconstruct - special cray compiler wrapper - can be one of 3 compilers
# ========================================================================================
ifeq ($(findstring ftn,$(rF90)),ftn)
	tF90:=$(filter pgf90% ifort%, $(shell $(pathF90) --version))
	TAG_COMPILER:=$(MACH)-ftn-
endif


# ================================
# Detect "openmpi's mpif90"
# ================================
ifeq ($(findstring opal_wrapper,$(notdir $(realpath $(pathF90)))),opal_wrapper)
	tF90:=$(firstword $(shell $(pathF90) -show))
	TAG_COMPILER:=$(MACH)-
endif


# ===================================
# Neither worked - figure its serial
# ===================================
ifeq ($(TAG_COMPILER),)
	tF90:=$(firstword $(F90))
	TAG_COMPILER:=$(MACH)-
endif



$(info DEBUG: tF90 is $(tF90))



# =================
# Detect "icx"
# =================
#Check if $F90 filename has "icc" in it, but don't be fooled if path has "icc"
ifeq ($(shell echo $(tF90) | sed 's,^ifort.*$$,ifx,g'),ifx)
	IS_ICC := 1
	TAG_COMPILER_VERSION :=$(wordlist 3, 3, $(shell $(pathF90) --version))
	MODULE_IN :=-I
	MODULE_OUT :=-module  
	TAG_COMPILER :=$(TAG_COMPILER)$(tF90)-$(TAG_COMPILER_VERSION)
else
	IS_ICC := 0
endif

# =================
# Detect "icc"
# =================
#Check if $F90 filename has "icc" in it, but don't be fooled if path has "icc"
ifeq ($(shell echo $(tF90) | sed 's,^ifort.*$$,ifort,g'),ifort)
	IS_ICC := 1
	TAG_COMPILER_VERSION :=$(wordlist 3, 3, $(shell $(pathF90) --version))
	MODULE_IN :=-I
	MODULE_OUT :=-module  
	TAG_COMPILER :=$(TAG_COMPILER)$(tF90)-$(TAG_COMPILER_VERSION)
else
	IS_ICC := 0
endif

# =================
# Detect "gcc"
# =================
#Check if $F90 filename has "gcc" in it, but don't be fooled if path has "gcc"
ifeq ($(shell echo $(notdir $(tF90)) | sed 's,gfortran[-mp].*,gfortran,g'),gfortran)
	IS_GCC := 1
	TAG_COMPILER_VERSION := $(shell ./dumpgccversion.sh $(tF90))
	INCLUDE_VERB :=-I
	MODULE_IN :=-I
	MODULE_OUT :=-J
	TAG_COMPILER :=$(TAG_COMPILER)gfortran-$(value TAG_COMPILER_VERSION)
else
	IS_GCC := 0
endif



# =================
# Detect "GNU", do the same things as "gcc" option
# This part is needed for recognizing gnu compiler on Ubuntu
# =================
#Check if $F90 filename has "gcc" in it, but don't be fooled if path has "gcc"
ifeq ($(tF90),GNU)
	tF90:=gfortran
	IS_GCC := 1
	TAG_COMPILER_VERSION  := $(shell ./dumpgccversion.sh $(tF90))
	INCLUDE_VERB :=-I
	MODULE_IN :=-I
	MODULE_OUT :=-J
	TAG_COMPILER :=$(TAG_COMPILER)$(tF90)-$(value TAG_COMPILER_VERSION)
else
	IS_GCC := 0
endif


# =================
# Detect "pgf90"
# =================
#Check if $F90 filename /is/ "pgf90"
ifeq ($(tF90),pgf90)
	IS_PGF90 := 1
	TAG_COMPILER_VERSION :=$(wordlist 2, 2, $(shell pgf90 --version))
	MODULE_IN :=-I
	MODULE_OUT :=-module  
	TAG_COMPILER :=$(TAG_COMPILER)pgf90-$(TAG_COMPILER_VERSION)
else
	IS_PGF90 := 0
endif


# =================
# Detect "g95"
# =================
#Check if $F90 filename /is/ "g95", but don't be fooled if path has "g95"
# can't tell if this actually works, because g95 just chokes...
ifeq ($(tF90),g95)
	IS_G95 := 1
	TAG_COMPILER_VERSION :=$(wordlist 3, 3, $(shell g95 --version))
	MODULE_IN :=-I
	MODULE_OUT :=-fmod=
	TAG_COMPILER :=$(TAG_COMPILER)g95-$(TAG_COMPILER_VERSION)
else
	IS_G95 := 0
endif

F90:=$(pathF90) $(filter-out $(firstword $(F90)), $(F90))
