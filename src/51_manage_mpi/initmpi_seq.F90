!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_seq
!! NAME
!!  initmpi_seq
!!
!! FUNCTION
!!  Initializes the MPI information for a sequential use of other routines.
!!
!! COPYRIGHT
!!  Copyright (C) 2004-2012 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!  mpi_enreg=informations about MPI parallelization
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! PARENTS
!!      bethe_salpeter,calc_sigc_me,calc_sigx_me,calc_vhxc_me,crho,debug_tools
!!      denfgr,fftprof,ks_ddiago,kss2wfk,linear_optics_paw,m_cprj_bspline
!!      m_fft_prof,m_gsphere,m_hamiltonian,m_io_kss,m_paw_pwij,m_ppmodel
!!      m_screening,m_wfs,mlwfovlp_qp,mpi_enreg_tools,mrggkk,mrgscr,optic
!!      paw_qpscgw,pawmkaewf,phfrq3,scfcv,screening,setshells,setup_bse
!!      setup_screening,setup_sigma,sigma,suscep_stat,susk,wfk_read_ene
!!      xc_kernel,xc_kernel_ADA
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initmpi_seq(mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_seq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(MPI_type),intent(out) :: mpi_enreg


!Local variables-------------------------------
 integer :: comm_self

! ***********************************************************************

 DBG_ENTER("COLL")

!$call nullify_mpi_enreg(MPI_enreg)

 mpi_enreg%paral_compil_kpt=0
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%paral_compil_respfn=0
 mpi_enreg%paral_spin=0
 mpi_enreg%flag_ind_kg_mpi_to_seq = 0
 mpi_enreg%paral_compil_mpio=0
 mpi_enreg%paral_level=0
 mpi_enreg%paralbd=0
 mpi_enreg%me=0
 mpi_enreg%nproc=0
 mpi_enreg%nproc_spin=1
 mpi_enreg%me_group=0
 mpi_enreg%nproc_group=0
 mpi_enreg%me_kpt=0
 mpi_enreg%me_fft=0
 mpi_enreg%me_spin=0
 mpi_enreg%me_cart_4d=0
 mpi_enreg%me_cart_3d=0
!mpi_enreg%nproc_fft=0
 mpi_enreg%nproc_fft=1    !changed to 1 by MG, see kpgsph
 mpi_enreg%paral_fft=0
 mpi_enreg%fft_option_lob=0
 mpi_enreg%me_g0=0
 mpi_enreg%num_group_fft=0
 mpi_enreg%num_group=0
 mpi_enreg%nproc_per_kpt=0
 mpi_enreg%world_group=0
 mpi_enreg%has_band_comm=0
 mpi_enreg%nproc_atom=1

!Set communicators to MPI_COMM_SELF to avoid problems in MPI calls
!MG: Not all the communicators defined in MPI_type are set to comm_self
!There are additional arrays that however should not be used in seq mode.

 comm_self = xmpi_self

 mpi_enreg%comm_fft            = comm_self
 mpi_enreg%comm_band           = comm_self
 mpi_enreg%comm_kpt            = comm_self
 mpi_enreg%commcart            = comm_self
 mpi_enreg%commcart_4d         = comm_self
 mpi_enreg%commcart_3d         = comm_self
 mpi_enreg%fft_master_comm     = comm_self
 mpi_enreg%respfn_master_comm  = comm_self
 mpi_enreg%comm_atom           = comm_self
 mpi_enreg%comm_bandspin       = comm_self
 mpi_enreg%comm_spin           = comm_self
 mpi_enreg%comm_spinfft        = comm_self

 mpi_enreg%mode_para=" " ! value for seq execution (anything but "b...")
 nullify(mpi_enreg%nscatterarr)
 nullify(mpi_enreg%ngatherarr )
 nullify(mpi_enreg%bandfft_kpt)
 nullify(mpi_enreg%tab_kpt_distrib)

 mpi_enreg%paral_img=0
 mpi_enreg%nimage=1
 mpi_enreg%me_one_img=mpi_enreg%me
 mpi_enreg%nproc_one_img=mpi_enreg%nproc
 mpi_enreg%comm_one_img=xmpi_world
 mpi_enreg%me_img=0
 mpi_enreg%nproc_img=1
 mpi_enreg%comm_img=comm_self

 DBG_EXIT("COLL")

end subroutine initmpi_seq
!!***
