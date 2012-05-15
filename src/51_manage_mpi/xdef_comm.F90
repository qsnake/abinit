!{\src2tex{textfont=tt}}
!!****f* ABINIT/xdef_comm
!! NAME
!!  xdef_comm
!!
!! FUNCTION
!!  Defines communicator and tools for MPI.
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2012 ABINIT group (MB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! NOTES
!!  Should become a module.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


!!***

!!****f* ABINIT/xcomm_world
!! NAME
!!  xcomm_world
!!
!! FUNCTION
!!  Defines the communicator with the world communicator
!!
!! INPUTS
!!  mpi_enreg= information about MPI parallelization
!!
!! OUTPUT
!!  spaceComm= MPI communicator
!!  === optional arguments ===
!!   myrank,mysize= rank and size for current proc
!!
!! SIDE EFFECTS
!!  None
!!
!! NOTES
!!  Used to distinguish serial execution when performing reduction operations.
!!
!! PARENTS
!!      afterscfloop,calc_cs,clnup1,datafordmft,dyfnl3,eig2tot,eltfrkin3
!!      eltfrnl3,fxphas,getshell,gstate,initberry,initmv,inwffil3,ioarr,mag_out
!!      mean_fftr,mklocl_wavelets,mkrho3,mover,mv_3dte,nstdy3,nstpaw3
!!      optics_paw,optics_paw_core,optics_vloc,outqmc,outscfcv,outwf,prtrhomxmn
!!      psolver_hartree,psolver_kernel,psolver_rhohxc,rdm,resp3dte,respfn
!!      rhofermi3,scfcv,scfcv3,scphon,sqnorm_v,subdiago,tddft,vtorho3
!!      wvl_memory,wvl_mkrho,wvl_newvtr,wvl_nl_gradient,wvl_rwwf,wvl_setngfft
!!      wvl_tail_corrections,wvl_vtorho,wvl_wfsinp_disk,wvl_wfsinp_reformat
!!      wvl_wfsinp_scratch,xdefineoff
!!
!! CHILDREN
!!
!! SOURCE

subroutine xcomm_world(mpi_enreg,spaceComm,myrank,mysize)

 use m_profiling

 use defs_abitypes
 use m_xmpi

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcomm_world'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 type(MPI_type),intent(in) :: mpi_enreg
 integer,intent(out) :: spaceComm
 integer,optional,intent(out) :: myrank,mysize

!Local variables-------------------
 integer :: ierr

! *********************************************************************

#if defined HAVE_MPI
 if (mpi_enreg%paral_compil_respfn == 1) then
   spaceComm = mpi_enreg%comm_respfn
 else if (mpi_enreg%paral_img == 1) then
   spaceComm = mpi_enreg%comm_one_img
 else
   spaceComm = xmpi_world
 end if
 if (present(myrank)) then
   call MPI_COMM_RANK(spaceComm,myrank,ierr)
 end if
 if (present(mysize)) then
   call MPI_COMM_SIZE(spaceComm,mysize,ierr)
 end if
#else
 spaceComm = abinit_comm_serial
 if (present(myrank)) myrank=0
 if (present(mysize)) mysize=1
 ierr=0
#endif
end subroutine xcomm_world
!!***


!!****f* ABINIT/xcomm_self
!! NAME
!!  xcomm_self
!!
!! FUNCTION
!!  Defines the communicator with the self communicator
!!
!! INPUTS
!!  None
!!
!! OUTPUT
!!  spaceComm= MPI communicator
!!
!! SIDE EFFECTS
!!  None
!!
!! NOTES
!!  Used to distinguish serial execution when performing reduction operations.
!!
!! PARENTS
!!      initmpi_img,inwffil,ioarr,outwf,outxfhist,wffopen
!!
!! CHILDREN
!!
!! SOURCE

subroutine xcomm_self(spaceComm)

 use m_profiling

 use defs_abitypes

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcomm_self'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(out) :: spaceComm

!Local variables-------------------

! *********************************************************************

#if defined HAVE_MPI
 spaceComm = MPI_COMM_SELF
#else
 spaceComm = abinit_comm_serial
#endif
end subroutine xcomm_self
!!***


!!****f* ABINIT/xcomm_init
!! NAME
!!  xcomm_init
!!
!! FUNCTION
!!  Initializes the communicator.
!!
!! INPUTS
!!  mpi_enreg= information about MPI parallelization
!!  spaceComm_bandfft= MPI communicator which has to be used
!!    if band-fft parallelism is activated (optional argument)
!!
!! OUTPUT
!!  spaceComm= MPI communicator
!!
!! SIDE EFFECTS
!!  None
!!
!! PARENTS
!!      atm2fft,atomden,back,back_wf,berryphase_new,bestwfs,calcdensph,cgwf
!!      chkexi,crho,ctocprj,datafordmft,dmft_solve,dotprod_g,dotprod_v
!!      dotprod_vn,eltfrhar3,eltfrloc3,eltfrxc3,energy,extrapwf,fftwfn
!!      forstrnps,forw,forw_wf,fourwf,fxphas,indirect_parallel_Fourier,inwffil
!!      ioarr,kpgio,ladielmt,lavnl,lobpcgccIIwf,lobpcgwf,loper3,m_green
!!      m_lobpcg,m_self,magcart,matrixelmt_g,meanvalue_g,mklocl_realspace
!!      mklocl_recipspace,mkrho,mlwfovlp,mlwfovlp_proj,mlwfovlp_pw,mlwfovlp_qp
!!      newrho,newvtr,newvtr3,nonlop_pl,nselt3,opernla_ylm,optics_paw
!!      optics_paw_core,optics_vloc,orthonormalize,outkss,outwf,pareigocc
!!      partial_dos_fractions,partial_dos_fractions_paw,pawdij,pawfrnhat
!!      pawgrnl,pawmkaewf,pawmkrhoij,poslifetime,prctfvw1,prctfvw2,pre_gather
!!      pre_scatter,precon,precon2,prep_fourwf,prep_getghc,prep_nonlop,projbd
!!      prtrhomxmn,pw_orthon,rhohxc,smatrix_pawinit,sqnorm_g,strhar,suscep_dyn
!!      suscep_kxc_dyn,suscep_stat,vtorho,vtorhotf,vtowfk,vtowfk3,wfconv
!!      wfkfermi3,wfsinp,zorthonormalize,zprecon3
!!
!! CHILDREN
!!
!! SOURCE

subroutine xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft)

 use m_profiling
 use defs_basis
 use defs_abitypes
 use m_xmpi

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcomm_init'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 type(MPI_type),intent(in) :: mpi_enreg
 integer,intent(out) :: spaceComm
 integer,intent(in),optional :: spaceComm_bandfft

!Local variables-------------------

! *********************************************************************

#if defined HAVE_MPI
 if (mpi_enreg%mode_para=='b'.and.present(spaceComm_bandfft)) then
   spaceComm=spaceComm_bandfft
 else if (mpi_enreg%paral_compil_respfn == 1) then
   if (mpi_enreg%paral_level == 2) then
     spaceComm = mpi_enreg%comm_respfn
   else
     if (mpi_enreg%num_group_fft /= 0) then
       spaceComm =  mpi_enreg%fft_comm(mpi_enreg%num_group_fft)
     else
       spaceComm = MPI_COMM_SELF
     end if
   end if
 elseif (mpi_enreg%paral_level > 1) then
   if (mpi_enreg%paral_level == 2) then
     if (mpi_enreg%paral_img==1) then
       spaceComm = mpi_enreg%comm_one_img
     else
       spaceComm = xmpi_world
     end if
   else
     if (mpi_enreg%num_group_fft /= 0) then
       spaceComm =  mpi_enreg%fft_comm(mpi_enreg%num_group_fft)
     else
       spaceComm = MPI_COMM_SELF
     end if
   end if
 else
   if (mpi_enreg%paral_img==1) then
     spaceComm = mpi_enreg%comm_one_img
   else
     spaceComm = xmpi_world
   end if
 end if
#else
 spaceComm = abinit_comm_serial
#endif

end subroutine xcomm_init
!!***


!!****f* ABINIT/xmaster_init
!! NAME
!!  xmaster_init
!!
!! FUNCTION
!!  Defines the master MPI node.
!!
!! INPUTS
!!  mpi_enreg= information about MPI parallelization
!!
!! OUTPUT
!!  master= master MPI node
!!
!! SIDE EFFECTS
!!  None
!!
!! PARENTS
!!      atomden,crho,energy,fftwfn,forstrnps,inwffil,ladielmt,lavnl,m_self
!!      mkrho,mkrho3,mlwfovlp,mlwfovlp_proj,mlwfovlp_qp,nselt3,outgkk,outwf
!!      pawmkaewf,prctfvw1,prctfvw2,rhofermi3,smatrix_pawinit,suscep_dyn
!!      suscep_kxc_dyn,suscep_stat,wfsinp
!!
!! CHILDREN
!!
!! SOURCE

subroutine xmaster_init(mpi_enreg,master)

 use m_profiling
 use defs_abitypes

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmaster_init'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 type(MPI_type),intent(in) :: mpi_enreg
 integer,intent(out) :: master

! *********************************************************************

 master = 0
 if(.false.)write(std_out,*)mpi_enreg%nproc

end subroutine xmaster_init
!!***

!!****f* ABINIT/xmaster_init_fft
!! NAME
!!  xmaster_init_fft
!!
!! FUNCTION
!!  Defines the master MPI-FFT node.
!!
!! INPUTS
!!  mpi_enreg= information about MPI parallelization
!!
!! OUTPUT
!!  master= master MPI node
!!
!! SIDE EFFECTS
!!  None
!!
!! PARENTS
!!      fxphas,prtrhomxmn
!!
!! CHILDREN
!!
!! SOURCE

subroutine xmaster_init_fft(mpi_enreg,master)

 use m_profiling
 use defs_abitypes

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmaster_init_fft'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 type(MPI_type),intent(in) :: mpi_enreg
 integer,intent(out) :: master

!Local variables-------------------

! *********************************************************************

#if defined HAVE_MPI
 if (mpi_enreg%paral_fft == 0) then
   if (mpi_enreg%paral_img==1) then
     master = mpi_enreg%me_one_img
   else
     master = mpi_enreg%me
   end if
 else
   master = mpi_enreg%master_fft
 end if
#else
 master = 0
#endif
end subroutine xmaster_init_fft
!!***

!!****f* ABINIT/xme_init
!! NAME
!!  xme_init
!!
!! FUNCTION
!!  Defines who is the current node ("me") in the world (or k-points) communicator
!!
!! INPUTS
!!  mpi_enreg= information about MPI parallelization
!!  option_comm= (optional argument) 1 = k-points communicator (default)
!!                                   2 = world communicator
!!
!! OUTPUT
!!  me= index of the current MPI node
!!
!! SIDE EFFECTS
!!  None
!!
!! PARENTS
!!      atomden,berryphase_new,cgwf,cprj_utils_mpi,crho,ctocprj,dmft_solve
!!      dyfnl3,eig2tot,eltfrkin3,eltfrnl3,energy,extrapwf,fftwfn,forstrnps
!!      getgsc,initylmg,inwffil,inwffil3,kpgio,ladielmt,lavnl,loper3,m_green
!!      m_qparticles,m_self,mklocl_realspace,mkrho,mkrho3,mksubham,mlwfovlp
!!      mlwfovlp_proj,mlwfovlp_pw,mlwfovlp_qp,newkpt,nselt3,nstdy3,nstpaw3
!!      nstwf3,nstwf4,optics_paw,optics_paw_core,optics_vloc,outgkk,outwf
!!      pareigocc,partial_dos_fractions_paw,pawmkaewf,pawmkrhoij,prctfvw1
!!      prctfvw2,rhofermi3,scfcv3,smatrix_pawinit,suscep,suscep_dyn
!!      suscep_kxc_dyn,suscep_stat,vtorho,vtorho3,vtowfk,vtowfk3,wfkfermi3
!!      wfsinp
!!
!! CHILDREN
!!
!! SOURCE

subroutine xme_init(mpi_enreg,me,option_comm)

 use m_profiling
 use defs_abitypes

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xme_init'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 type(MPI_type),intent(in) :: mpi_enreg
 integer,intent(out) :: me
 integer,intent(in),optional :: option_comm

!Local variables-------------------

! *********************************************************************

#if defined HAVE_MPI
 if(mpi_enreg%paral_compil_respfn == 1.and.mpi_enreg%me_respfn/=-1) then
   me = mpi_enreg%me_respfn
 else if ((mpi_enreg%paral_compil_kpt==1).and.(mpi_enreg%paral_compil_fft==1)) then
   me=mpi_enreg%me_kpt
   if (present(option_comm)) then
     if (option_comm==2) then
       if (mpi_enreg%paral_img==1) then
         me = mpi_enreg%me_one_img
       else
         me=mpi_enreg%me
       end if
     end if
   end if
 else
   if (mpi_enreg%paral_img==1) then
     me = mpi_enreg%me_one_img
   else
     me=mpi_enreg%me
   end if
 end if
#else
 me = 0
#endif
end subroutine xme_init
!!***

!!****f* ABINIT/xproc_init
!! NAME
!!  xproc_init
!!
!! FUNCTION
!!  Defines the total number of procesoors.
!!
!! INPUTS
!!  mpi_enreg= information about MPI parallelization
!!
!! OUTPUT
!!  nproc_max= total number of processors
!!
!! SIDE EFFECTS
!!  None
!!
!! PARENTS
!!      m_green,m_self,mlwfovlp,mlwfovlp_proj,mlwfovlp_pw,mlwfovlp_qp,newkpt
!!      smatrix_pawinit,wfsinp
!!
!! CHILDREN
!!
!! SOURCE

subroutine xproc_init(mpi_enreg,nproc_max)

 use m_profiling
 use defs_abitypes

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xproc_init'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(out) :: nproc_max
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------

! *********************************************************************

#if defined HAVE_MPI
 if (mpi_enreg%paral_compil_respfn == 1) then
   nproc_max=mpi_enreg%nproc_respfn
 else
   if (mpi_enreg%paral_img==1) then
     nproc_max = mpi_enreg%nproc_one_img
   else
     nproc_max=mpi_enreg%nproc
   end if
 end if

#else
 nproc_max = 1
#endif
end subroutine xproc_init
!!***
