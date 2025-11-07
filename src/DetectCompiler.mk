# What`s the compiler called?
pathF90 := $(shell which $(firstword $(F90)))
# $(info DEBUG: pathF90 is $(pathF90))

# What`s the real path to it?
rF90 := $(notdir $(abspath $(pathF90)))
# $(info DEBUG: rF90 is $(rF90))


# ########################################
#
# Check which MPI wrapper we are using
#
# ########################################

TAG_COMPILER:=$(MACH)

# ================================
# Detect "mpif90" and decontstruct
# ================================
ifeq ($(findstring mpif90,$(notdir $(realpath $(pathF90)))),mpif90)
	tF90:=$(firstword $(shell $(pathF90) -show))
	TAG_COMPILER:=$(TAG_COMPILER)-mpif90
endif

# ================================
# Detect "mpiifx" and decontstruct
# ================================
ifeq ($(findstring mpiifx,$(notdir $(realpath $(pathF90)))),mpiifx)
	tF90:=$(firstword $(shell $(pathF90) --version))
	TAG_COMPILER:=$(MACH)-mpiifx
endif

# ================================
# Detect "mpiifort" and decontstruct
# ================================
ifeq ($(findstring mpiifort,$(notdir $(realpath $(pathF90)))),mpiifort)
	tF90:=$(firstword $(shell $(pathF90) --version))
	TAG_COMPILER:=$(MACH)-mpiifort
endif

# ================================
# Detect "mpifort" and decontstruct
# ================================
ifeq ($(findstring mpifort,$(notdir $(realpath $(pathF90)))),mpifort)
	tF90:=$(firstword $(shell $(pathF90) --version))
	TAG_COMPILER:=$(MACH)-mpifort
endif

# ========================================================================================
# Detect "ftn" and deconstruct - special cray compiler wrapper - can be one of 3 compilers
# ========================================================================================
ifeq ($(findstring ftn,$(rF90)),ftn)
	tF90:=$(filter pgf90% ifort% GNU%, $(shell $(pathF90) -craype-verbose -V 2>>/dev/null ))
	TAG_COMPILER:=$(MACH)-ftn
endif

# =====================================================================================================
# Detect "ftn" as "driver" and deconstruct - special cray compiler wrapper - can be one of 3 compilers
# =====================================================================================================
ifeq ($(findstring driver,$(rF90)),driver)
	tF90:=$(filter pgf90% ifort% GNU% ftn_driver.exe%, $(shell $(pathF90) -craype-verbose -V 2>>/dev/null ))
	TAG_COMPILER:=$(TAG_COMPILER)-ftn
endif

# ================================
# Detect "openmpi's mpif90"
# ================================
ifeq ($(findstring opal_wrapper,$(notdir $(realpath $(pathF90)))),opal_wrapper)
	tF90:=$(firstword $(shell $(pathF90) --version))
	TAG_COMPILER:=$(MACH)-$(notdir $(pathF90))
endif

# ===================================
# Neither worked - figure its serial
# ===================================
ifeq ($(TAG_COMPILER),$(MACH))
	tF90:=$(firstword $(F90))
endif

ifeq ($(tF90),GNU)
	tF90:=gfortran
endif

# (info DEBUG: tF90 is $(tF90))
# $(info DEBUG: TAG_COMPILER starts as $(TAG_COMPILER))


# ########################################
#
# Check which C compiler we are using
#
# ########################################

IS_GCC := 0

# =================
# Detect "ifx"
# =================
# Check if $F90 filename has "ifx" in it, but don't be fooled if path has "ifx"
ifeq ($(shell echo $(tF90) | sed 's,^ifx.*$$,ifx,g'),ifx)
	TAG_COMPILER_VERSION :=$(wordlist 3, 3, $(shell $(pathF90) --version))
	MODULE_IN :=-I
	MODULE_OUT :=-module 
	CFLAGS += -std=c99
	TAG_COMPILER :=$(TAG_COMPILER)-$(tF90)-$(TAG_COMPILER_VERSION)
endif

# =================
# Detect "ifort"
# =================
# Check if $F90 filename has "ifort" in it, but don't be fooled if path has "ifort"
ifeq ($(shell echo $(tF90) | sed 's,^ifort.*$$,ifort,g'),ifort)
	TAG_COMPILER_VERSION :=$(wordlist 3, 3, $(shell $(pathF90) --version))
	MODULE_IN :=-I
	MODULE_OUT :=-module 
	CFLAGS += -std=c99
	TAG_COMPILER :=$(TAG_COMPILER)-$(tF90)-$(TAG_COMPILER_VERSION)
endif

# =================
# Detect "gfortran"
# =================
# Check if $F90 filename has "gfortran" in it, but don't be fooled if path has "gfortran"
ifeq ($(shell echo $(notdir $(tF90)) | sed 's,gfortran[-mp].*,gfortran,g'),gfortran)
	TAG_COMPILER_VERSION := $(shell ./dumpgccversion.sh $(tF90))
	INCLUDE_VERB :=-I
	MODULE_IN :=-I
	MODULE_OUT :=-J
	CFLAGS += -std=c99
	TAG_COMPILER :=$(TAG_COMPILER)-gfortran-$(value TAG_COMPILER_VERSION)
	IS_GCC := 1
endif

# =================
# Detect "pgf90"
# =================
# Check if $F90 filename /is/ "pgf90"
ifeq ($(tF90),pgf90)
	TAG_COMPILER_VERSION :=$(wordlist 2, 2, $(shell pgf90 --version))
	MODULE_IN :=-I
	MODULE_OUT :=-module 
	TAG_COMPILER :=$(TAG_COMPILER)-pgf90-$(TAG_COMPILER_VERSION)
endif

# ======================
# Detect Cray Compiler
# ======================
ifeq ($(shell echo $(notdir $(tF90)) | sed 's,^ftn_driver.*,ftn_driver,g'),ftn_driver)
    TAG_COMPILER_VERSION :=$(wordlist 5, 5, $(shell $(pathF90) -V 2>&1))
    MODULE_IN :=-I
    MODULE_OUT :=-J
    TAG_COMPILER :=$(TAG_COMPILER)-cray-$(TAG_COMPILER_VERSION)
endif

# =================
# Detect "g95"
# =================
# Check if $F90 filename /is/ "g95"
# can't tell if this actually works, because g95 just chokes...
ifeq ($(tF90),g95)
	TAG_COMPILER_VERSION :=$(wordlist 3, 3, $(shell g95 --version))
	MODULE_IN :=-I
	MODULE_OUT :=-fmod=
	TAG_COMPILER :=$(TAG_COMPILER)-g95-$(TAG_COMPILER_VERSION)
endif

# =================
# Detect "nvfortran"
# =================
# Check if $F90 filename has "nvfortran" in it, but don't be fooled if path has "nvfortran"
ifeq ($(shell echo $(tF90) | sed 's,^nvfortran.*$$,nvfortran,g'),nvfortran)
	TAG_COMPILER_VERSION :=$(wordlist 2, 2, $(shell $(pathF90) --version))
	MODULE_IN :=-I
	MODULE_OUT :=-module 
	TAG_COMPILER :=$(TAG_COMPILER)-nvfortran-$(TAG_COMPILER_VERSION)
endif

F90:=$(pathF90) $(filter-out $(firstword $(F90)), $(F90))
