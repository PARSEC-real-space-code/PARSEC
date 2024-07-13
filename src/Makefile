# the included file contains the machine specific 
# compiler, flags, libraries, etc
# change to fit your environment!
# if you want to do a lot of different builds quickly
# get rid of this line and use:
# make  MACH=yourfilename-without-prefixhere -f Makefile
#
MACH ?= stampede3
include config/make.$(MACH)

# Detect Compiler - modified from original 762Studios under DWTFYWTPL2
# DON'T MODIFY THIS UNLESS YOU KNOW WHAT ITS DOING....
# You shouldn't need to change this unless one of the compilers
# goes through a major revamping... where they change and deprecate
# often used flags
include DetectCompiler.mk

#$(shell echo "$(TAG_COMPILER)")

#
# source files without preprocessing
#
SRCF	=  anderson.f90 cg.f90 corecd.f90 dipole.f90 domove.f90 \
 eigensave.f90 efvdw.f90 erfc.f90 ewald_sum.f90 exc_nspn.f90 exc_spn.f90 export_grid_data.f90 \
 flevel.f90 forceion.f90 forceion_wire.f90 force_loc_pb.f90 force_charged_sh.f90 forcept.f90 \
 forcloc.f90 forcloc_wire.f90 forcnloc.f90 fornberg.f90 forpbc.f90 \
 fourier_f.f90 genbro.f90 getsre.f90 grid_partition.f90 g_space.f90 \
 inipbc.f90 initchrg.f90 initchrg_f.f90 initial.f90 ionpbc.f90 ionpot.f90 ionpot_wire.f90 \
 irrbz.f90 isolat.f90 moldyn.f90 qinv.f90 \
 msecant1.f90 msecant2.f90 msecant3.f90 \
 pbc_grid_coefs.f90 pbe.f90 polar.f90 pot_local.f90 potshift.f90 \
 pseudo.f90 ptable.f90 spinmix.f90 spinunmix.f90 spline.f90 spn_pbe.f90 \
 subspace_sub.f90 symgen.f90 symmetries.f90 symm_forces.f90 \
 totnrg.f90 v_first.f90 zmatvec_so.f90 ewald_slab.f90 \
 ionpot_slab.f90 forcloc_slab.f90 calc_bands.f90 laplacian_test.f90
OBJF	= $(SRCF:.f90=.o) lbfgsb.o linpack.o timer.o 

#
# source files with preprocessing
#
PSRCF     = cfftw.F90 chkmemory.F90 comm_neigh.F90 eigval.F90 \
 exit_err.F90 hartset.F90 hartset_wire.F90 hpotcg.F90 hpot_g.F90 init_var.F90 \
 myflush.F90 newrho.F90 nonloc.F90 parsec.F90 pls.F90 restart.F90 \
 restart_run.F90 rho_hart.F90 setup.F90 soforc.F90 \
 topo_aware.F90 upot.F90 usrinputfile.F90 wfnsave.F90 \
 hartset_slab.F90 dos.F90 sort.F90 timing.F90 bxc.F90 spec_tet.F90 \
 calc_nscf.F90 gw_write.F90 h5_wfnsave.F90
POBJF = $(PSRCF:.F90=.o)

#
# source files with preprocessing and complex/real algebra
# they will be preprocessed twice
#
PSRCZ   =  arpk.f90z buffer_module.f90z matvecB.f90z calcden.f90z chebdav.f90z dgks.f90z diagla.f90z \
 fdo5scl.f90z lancz_bound.f90z \
 matvec.f90z orth_normal.f90z \
 peig_inter.f90z plusu.f90z psum.f90z subspace.f90z chebff.f90z
POBJZ   = $(PSRCZ:.f90z=.o)

#
# modules, with optimization
#
FOPTS?=$(FFLAGS)
MODULE	= const.o esdf_key_mod.o esdf_mod.o structures.o
BUILDDIR=obj/$(TAG_COMPILER)$(EXT)
TMPSRCDIR=tmp/$(TAG_COMPILER)$(EXT)
OBJ=$(addprefix $(BUILDDIR)/,$(POBJF) $(POBJZ) $(OBJF))

####
# Sets module path to find those pesky .mod files
####
MODULE_PATH=-I./ $(MODULE_IN)$(BUILDDIR) $(MODULE_OUT)$(BUILDDIR) $(TRLAN_MOD_PATH) \
	    $(BGW_MOD_PATH) $(HDF5_MOD_PATH)

#Makes temporary directories
$(shell mkdir -p $(BUILDDIR))
$(shell mkdir -p $(TMPSRCDIR))

Q = @

.PHONY: all vars

all: vars parsec-$(TAG_COMPILER)$(EXT)
	@echo $(TAG_COMPILER) finished

parsec-$(TAG_COMPILER)$(EXT): $(addprefix $(BUILDDIR)/,$(MODULE)) $(OBJ)
	$(FC) -o parsec-$(TAG_COMPILER)$(EXT) $(FOPTS) $(MEMFLAG) \
	$(OBJ) $(addprefix $(BUILDDIR)/,$(MODULE)) \
	$(LIBTRLAN) $(LIBARPACK2) $(LIBARPACK1) \
	$(LIBLAPA) $(LIBBLAS) $(LIBFFTW) $(LIBMPI) $(LIBBGW) $(LIBHDF5)

.SUFFIXES: .f90p .f90z $(SUFFIXES)


vars:
	@echo pathF90 $(pathF90)
	@echo rF90 $(rF90)
	@echo tF90 $(tF90) 
	@echo TAG_COMPILER_VERSION $(TAG_COMPILER_VERSION)
	@echo MODULE_IN $(MODULE_IN)
	@echo MODULE_OUT $(MODULE_OUT)
	@echo TAG_COMPILER $(TAG_COMPILER)
ifeq ($(TAG_COMPILER),)
	@echo "Compiler not autodetected - this is a problem. Fix it!"
	@exit 1
endif
	

$(OBJ) : $(addprefix $(BUILDDIR)/,$(MODULE))


$(BUILDDIR)/%.o: %.f
	$(Q)$(FC) -c $(FOPTS) $*.f $(MODULE_PATH) -o $(BUILDDIR)/$*.o
	$(Q)echo "$(FC) $<"

$(BUILDDIR)/%.o: %.f90
	$(Q)$(F90) -c $(FOPTS) $*.f90 $(MODULE_PATH) -o $(BUILDDIR)/$*.o
	$(Q)echo "$(F90) $<"

$(BUILDDIR)/usrinputfile.o: usrinputfile.F90
	$(Q)echo '#define GIT_COMMIT "$(shell git log -1 --pretty=format:"%h")"' > $(TMPSRCDIR)/tmp_$<
	$(Q)echo '#define MACH "$(MACH)"' >> $(TMPSRCDIR)/tmp_$<
	$(Q)cat $< >> $(TMPSRCDIR)/tmp_$<
	$(Q)$(F90) $(CPPOPT) -c $(FOPTS) $(MODULE_PATH) $(TMPSRCDIR)/tmp_$< -o $*.o
	$(Q)echo "$(F90) $<"

$(BUILDDIR)/%.o: %.F90
	$(Q)$(F90) $(CPPOPT) -c $(FOPTS) $(MODULE_PATH) $< -o $(BUILDDIR)/$*.o
	$(Q)echo "$(F90) $<"

$(BUILDDIR)/%.o: %.f90p
	$(CPP) $(CPPOPT) $< | sed '/^#pragma/d' > $(TMPSRCDIR)/$*.tmp.f90
	$(F90) -c $(FOPTS) $(MODULE_PATH) $(TMPSRCDIR)/$*.tmp.f90 -o $(BUILDDIR)/$*.o

$(BUILDDIR)/%.o: %.f90z
	$(Q)$(CPP) $(CPPOPT) $< | sed '/^#pragma/d' > $(TMPSRCDIR)/$*.tmp.f90
	$(Q)echo "$(CPP) $<"
	$(Q)$(CPP) $(CPPOPT) -DCPLX $< | sed '/^#pragma/d' >> $(TMPSRCDIR)/$*.tmp.f90
	$(Q)echo "$(CPP) -DCPLX $<"
	$(Q)$(F90) -c $(FOPTS) $(MODULE_PATH) $(TMPSRCDIR)/$*.tmp.f90 -o $(BUILDDIR)/$*.o
	$(Q)echo "$(F90) $<"

cleanall: clean
	/bin/rm -f parsec-*.ser* parsec-*.mpi*
	/bin/rm -f *~
	/bin/rm -rf obj/ tmp/

allclean: cleanall

cleanlib: clean

clean: cleantmp
	/bin/rm -f *.o
	/bin/rm -f *.mod
	/bin/rm -f parsec-$(TAG_COMPILER)$(EXT) parsec.ser* parsec.mpi*

cleantmp: 
	/bin/rm -f *.tmp.f *.tmp.f90 *__genmod.f90
	/bin/rm -rf $(TMPSRCDIR) *.tmp.f *.tmp.f90 *__genmod.f90
	/bin/rm -rf $(BUILDDIR)

tags:
	ctags --if0=yes --language-force=Fortran -R *.f90z *.F90 *.f90 *.f -h ".h,.inc" 
############################################
# Dependencies on "include" or "use" files #
############################################
$(BUILDDIR)/structures.o : $(BUILDDIR)/const.o
$(BUILDDIR)/esdf_mod.o: $(BUILDDIR)/esdf_key_mod.o
$(BUILDDIR)/parsec.o : def.h
$(BUILDDIR)/arpk.o : arpkdbg.h
# something wasn't quite right (shrugs)
$(BUILDDIR)/arpk.o $(BUILDDIR)/chebdav.o $(BUILDDIR)/chebff.o $(BUILDDIR)/diagla.o $(BUILDDIR)/lancz_bound.o $(BUILDDIR)/subspace.o $(BUILDDIR)/zmatvec_so.o: $(BUILDDIR)/matvecB.o
$(BUILDDIR)/arpk.o $(BUILDDIR)/chebdav.o $(BUILDDIR)/chebff.o $(BUILDDIR)/diagla.o $(BUILDDIR)/lancz_bound.o $(BUILDDIR)/subspace.o $(BUILDDIR)/zmatvec_so.o: $(BUILDDIR)/matvec.o 
$(BUILDDIR)/matvecB.o :$(BUILDDIR)/buffer_module.o
#$(addprefix $(BUILDDIR),trlanc_diag.o arpk.o chebdav.o diagla.o lancz_bound.o subspace.o zmatvec_so.o):$(BUILDDIR)/matvec.o

############################################################
# Special build rules for certain files compiler dependent
# ############################################################
ifeq ($(IS_GCC),1)
$(BUILDDIR)/myflush.o : myflush.F90
	$(F90) $(CPPOPT) -c $(FOPTS) -fall-intrinsics $(MODULE_PATH) $< -o $*.o
endif
