
! ==============================================================
!
!     Copyright (C) 2005 Finite Difference Research Group
!     This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! ==============================================================
!                  Intel Trace Analyzer setup 
! ==============================================================
     common /vtclasses/ calc_class, eig_class, chebdav_class, chebff_class, subspace_class, mvec_class,setup_class, post_class
     integer ::         calc_class, eig_class, chebdav_class, chebff_class, subspace_class, mvec_class,setup_class, post_class
     common /vtstates/ &
        vt_hart_xc, vt_lapmvs, vt_nonloc, vt_scf_state,&
        vt_post_scf_state, vt_arpack, vt_so,&
        vt_eigensort_state, vt_eigencomm_state,&
        vt_orthnormal, vt_dgks, vt_gemv,&
        vt_lancz, vt_matvec_ps, vt_matvec,&
        vt_bdxc_leftover, vt_bdxc, vt_chebff,&
        vt_chebff_update, vt_chebff_decomp,&
        vt_chebff_filter, vt_chebdav,&
        vt_chebdav_MM, vt_chebdav_decomp,&
        vt_chebdav_Iupdate, vt_chebdav_IIupdate,&
        vt_chebdav_reorder, vt_subspace_main,&
        vt_subspace_filter, vt_subspace_update,&
        vt_forces_state, vt_move_state,&
        vt_pbc_ion, vt_post_processing_state,&
        vt_grid_setup_state, vt_pre_scf_state, &
        vt_naupd,vt_seupd,vt_buffer
        integer :: &
        vt_hart_xc, vt_lapmvs, vt_nonloc,vt_scf_state,&
        vt_post_scf_state, vt_arpack, vt_so,&
        vt_eigensort_state, vt_eigencomm_state,&
        vt_orthnormal, vt_dgks, vt_gemv,&
        vt_lancz, vt_matvec_ps, vt_matvec,&
        vt_bdxc_leftover, vt_bdxc, vt_chebff,&
        vt_chebff_update, vt_chebff_decomp,&
        vt_chebff_filter, vt_chebdav,&
        vt_chebdav_MM, vt_chebdav_decomp,&
        vt_chebdav_Iupdate, vt_chebdav_IIupdate,&
        vt_chebdav_reorder, vt_subspace_main,&
        vt_subspace_filter, vt_subspace_update,&
        vt_forces_state, vt_move_state,&
        vt_pbc_ion, vt_post_processing_state,&
        vt_grid_setup_state, vt_pre_scf_state,&
        vt_naupd,vt_seupd,vt_buffer   
