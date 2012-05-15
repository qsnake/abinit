#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
!!***

!!****f* ABINIT/init_mpi_enreg
!! NAME
!! init_mpi_enreg
!!
!! FUNCTION
!!  Initialise a mpi_enreg structure with dataset independent values.
!!  Other values of mpi_enreg are dataset dependent, and should NOT be initialized
!!  inside abinit.F90 .
!!  XG 071118 : At present several other values are
!!  initialized temporarily inside invars1.F90, FROM THE DTSET
!!  VALUES. In order to releave the present constraint of having mpi_enreg
!!  equal for all datasets, they should be reinitialized from the dtset values
!!  inside invars2m.F90 (where there is a loop over datasets, and finally,
!!  reinitialized from the dataset values inside each big routine called by driver,
!!  according to the kind of parallelisation that is needed there.
!!  One should have one init_mpi_dtset routine (or another name) per big routine (well, there is also
!!  the problem of TDDFT ...). Also, one should have a clean_mpi_dtset called at the end
!!  of each big routine, as well as invars1.F90 or invars2m.F90 .
!!
!! INPUTS
!!  init_mpi= --optional, default=.true.-- if true, MPI is initialized before MPI_enreg
!!
!! SIDE EFFECTS
!!  MPI_enreg<MPI_type>=All pointer set to null().
!!
!! PARENTS
!!      abinit,m_ab6_invars_f90,m_results_out
!!
!! CHILDREN
!!
!! SOURCE
subroutine init_mpi_enreg(mpi_enreg,init_mpi)

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
#define ABI_FUNC 'init_mpi_enreg'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi, except_this_one => init_mpi_enreg
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 logical,optional,intent(in) :: init_mpi
 type(MPI_type),intent(inout) :: MPI_enreg

!Local variables-------------------------------
 logical :: init_mpi_
#if defined HAVE_MPI
 integer :: ierr
#endif
#if defined HAVE_MPI_IO
 character(len=500) :: message
#endif

! *********************************************************************

 init_mpi_=.true.;if (present(init_mpi)) init_mpi_=init_mpi
 if (init_mpi_) then
   call xmpi_init()
 end if

!Default for sequential use
 call nullify_mpi_enreg(mpi_enreg)

 mpi_enreg%world_comm=0
 mpi_enreg%world_group=0
 mpi_enreg%me=0
 mpi_enreg%nproc=1
 mpi_enreg%num_group_fft = 0 ! in some cases not initialized but referenced in xdef_comm.F90
 mpi_enreg%paral_compil=0
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%paral_compil_mpio=0
 mpi_enreg%mode_para="n"
 mpi_enreg%flag_ind_kg_mpi_to_seq = 0
 mpi_enreg%paral_spin=0
 mpi_enreg%paral_img=0

!MG080916 If we want to avoid MPI preprocessing options, %proc_distr should be always allocated and
!set to mpi_enreg%me. In such a way we can safely test its value inside loops parallelized over k-points
!For the time being, do not remove this line since it is needed in outkss.F90.
!nullify(mpi_enreg%proc_distrb)
!nullify(mpi_enreg%bandfft_kpt,mpi_enreg%tab_kpt_distrib)

!Initialize MPI
#if defined HAVE_MPI
 mpi_enreg%world_comm=xmpi_world
!mpi_enreg%world_group=MPI_GROUP_NULL
 call MPI_COMM_GROUP(mpi_enreg%world_comm,mpi_enreg%world_group,ierr)
 call MPI_COMM_RANK(xmpi_world,mpi_enreg%me,ierr)
 call MPI_COMM_SIZE(xmpi_world,mpi_enreg%nproc,ierr)
 mpi_enreg%paral_compil=1
#endif

!Signal MPI I/O compilation has been activated
#if defined HAVE_MPI_IO
 mpi_enreg%paral_compil_mpio=1
 if(mpi_enreg%paral_compil==0)then
   write(message,'(6a)') ch10,&
&   ' abinit : ERROR -',ch10,&
&   '  In order to use MPI_IO, you must compile with the MPI flag ',ch10,&
&   '  Action : recompile your code with different CPP flags.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
!Test the opening, writing, closing and deleting of a test file
!call MPI_FILE_OPEN(xmpi_world,"Testfile",MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)
!call print_ierr(ierr,"abinit","MPI_FILE_OPEN")
!Write 1 integer per proc
!data=mpi_enreg%me;offset=4*mpi_enreg%me
!call MPI_FILE_WRITE_AT(fh,offset,data,1,MPI_INTEGER,mpi_status,ierr)
!call print_ierr(ierr,"abinit","MPI_FILE_WRITE_AT")
!Close file
!call MPI_FILE_CLOSE(fh,ierr)
!call print_ierr(ierr,"abinit","MPI_FILE_CLOSE")
!Delete file
!call MPI_FILE_DELETE(fh,MPI_INFO_NULL,ierr)
!call print_ierr(ierr,"abinit","MPI_FILE_DELETE")
#endif

!Initialize comm_respfn, used in respfn
 mpi_enreg%comm_respfn=mpi_enreg%world_comm

!Initialize paral_compil_kpt, actually always equal to paral_compil
!(paral_compil_kpt should be suppressed after big cleaning)
 mpi_enreg%paral_compil_kpt=0
 if(mpi_enreg%paral_compil==1) mpi_enreg%paral_compil_kpt=1
end subroutine init_mpi_enreg
!!***

!!****f* ABINIT/nullify_mpi_enreg
!! NAME
!! nullify_mpi_enreg
!!
!! FUNCTION
!!  nullify a mpi_enreg datastructure
!!
!! SIDE EFFECTS
!!  MPI_enreg<MPI_type>=All pointer set to null().
!!
!! PARENTS
!!      anaddb,cut3d,fftprof,kss2wfk,lwf,m_fft_prof,m_wfs,mpi_enreg_tools
!!      mrgddb,mrggkk,mrgscr,newsp,optic,ujdet
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_mpi_enreg(MPI_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_mpi_enreg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: MPI_enreg

! *********************************************************************

 nullify(mpi_enreg%band_comm)
 nullify(mpi_enreg%fft_comm)
 nullify(mpi_enreg%proc_distrb)
 nullify(mpi_enreg%kpt_comm)
 nullify(mpi_enreg%kptdstrb)
 nullify(mpi_enreg%kptdstrbi)
 nullify(mpi_enreg%kpt_loc2fbz_sp)
 nullify(mpi_enreg%kpt_loc2ibz_sp)
 nullify(mpi_enreg%fmkmem)
 nullify(mpi_enreg%mkmem)
 nullify(mpi_enreg%nplanes_fft)
 nullify(mpi_enreg%ind_fft_planes)

 nullify(mpi_enreg%bandfft_kpt)
 nullify(mpi_enreg%tab_kpt_distrib)
 nullify(mpi_enreg%respfn_group)
 nullify(mpi_enreg%respfn_comm)
 nullify(mpi_enreg%nscatterarr)
 nullify(mpi_enreg%ngatherarr)
 nullify(mpi_enreg%sizecart)
 nullify(mpi_enreg%coords)
 nullify(mpi_enreg%atom_indx)
 nullify(mpi_enreg%keywp)
 nullify(mpi_enreg%trialproc)

 nullify(mpi_enreg%distrb_img)
 nullify(mpi_enreg%index_img)

end subroutine nullify_mpi_enreg
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/destroy_mpi_enreg
!! NAME
!! destroy_mpi_enreg
!!
!! FUNCTION
!!  Destroy a mpi_enreg datastructure
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2012 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SIDE EFFECTS
!!  MPI_enreg<MPI_type>=Datatype gathering information on the parallelism.
!!
!! PARENTS
!!      abinit,anaddb,cut3d,fftprof,inwffil,kss2wfk,lwf,m_ab6_invars_f90
!!      m_fft_prof,m_wfs,mrgddb,mrggkk,mrgscr,newsp,ujdet
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_mpi_enreg(MPI_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_mpi_enreg'
 use interfaces_51_manage_mpi, except_this_one => destroy_mpi_enreg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: MPI_enreg

!Local variables-------------------------------

! *********************************************************************

 if (associated(mpi_enreg%bandfft_kpt).or.associated(mpi_enreg%tab_kpt_distrib)) then
   call clnmpi_bandfft(mpi_enreg)
 end if

 if (associated(mpi_enreg%band_comm)) then
   ABI_DEALLOCATE(mpi_enreg%band_comm)
   nullify(mpi_enreg%band_comm)
 end if
 if (associated(mpi_enreg%fft_comm)) then
   ABI_DEALLOCATE(mpi_enreg%fft_comm)
   nullify(mpi_enreg%fft_comm)
 end if
 if (associated(mpi_enreg%proc_distrb)) then
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
   nullify(mpi_enreg%proc_distrb)
 end if
 if (associated(mpi_enreg%kpt_comm)) then
   ABI_DEALLOCATE(mpi_enreg%kpt_comm)
   nullify(mpi_enreg%kpt_comm)
 end if
 if (associated(mpi_enreg%kptdstrb)) then
   ABI_DEALLOCATE(mpi_enreg%kptdstrb)
   nullify(mpi_enreg%kptdstrb)
 end if
 if (associated(mpi_enreg%kptdstrbi)) then
   ABI_DEALLOCATE(mpi_enreg%kptdstrbi)
   nullify(mpi_enreg%kptdstrbi)
 end if
 if (associated(mpi_enreg%kpt_loc2fbz_sp)) then
   ABI_DEALLOCATE(mpi_enreg%kpt_loc2fbz_sp)
   nullify(mpi_enreg%kpt_loc2fbz_sp)
 end if
 if (associated(mpi_enreg%kpt_loc2ibz_sp)) then
   ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
   nullify(mpi_enreg%kpt_loc2ibz_sp)
 end if
 if (associated(mpi_enreg%fmkmem)) then
   ABI_DEALLOCATE(mpi_enreg%fmkmem)
   nullify(mpi_enreg%fmkmem)
 end if
 if (associated(mpi_enreg%mkmem)) then
   ABI_DEALLOCATE(mpi_enreg%mkmem)
   nullify(mpi_enreg%mkmem)
 end if
 if (associated(mpi_enreg%nplanes_fft)) then
   ABI_DEALLOCATE(mpi_enreg%nplanes_fft)
   nullify(mpi_enreg%nplanes_fft)
 end if
 if (associated(mpi_enreg%ind_fft_planes)) then
   ABI_DEALLOCATE(mpi_enreg%ind_fft_planes)
   nullify(mpi_enreg%ind_fft_planes)
 end if
 if (associated(mpi_enreg%bandfft_kpt)) then
!  FIXME here there is a memory leak.
   ABI_DEALLOCATE(mpi_enreg%bandfft_kpt)
   nullify(mpi_enreg%bandfft_kpt)
 end if
 if (associated(mpi_enreg%tab_kpt_distrib)) then
   ABI_DEALLOCATE(mpi_enreg%tab_kpt_distrib)
   nullify(mpi_enreg%tab_kpt_distrib)
 end if
 if (associated(mpi_enreg%respfn_group)) then
   ABI_DEALLOCATE(mpi_enreg%respfn_group)
   nullify(mpi_enreg%respfn_group)
 end if
 if (associated(mpi_enreg%respfn_comm)) then
   ABI_DEALLOCATE(mpi_enreg%respfn_comm)
   nullify(mpi_enreg%respfn_comm)
 end if
 if (associated(mpi_enreg%nscatterarr)) then
   ABI_DEALLOCATE(mpi_enreg%nscatterarr)
   nullify(mpi_enreg%nscatterarr)
 end if
 if (associated(mpi_enreg%ngatherarr)) then
   ABI_DEALLOCATE(mpi_enreg%ngatherarr)
   nullify(mpi_enreg%ngatherarr)
 end if
 if (associated(mpi_enreg%sizecart)) then
   ABI_DEALLOCATE(mpi_enreg%sizecart)
   nullify(mpi_enreg%sizecart)
 end if
 if (associated(mpi_enreg%coords)) then
   ABI_DEALLOCATE(mpi_enreg%coords)
   nullify(mpi_enreg%coords)
 end if
 if (associated(mpi_enreg%atom_indx)) then
   ABI_DEALLOCATE(mpi_enreg%atom_indx)
   nullify(mpi_enreg%atom_indx)
 end if
 if (associated(mpi_enreg%keywp)) then
   ABI_DEALLOCATE(mpi_enreg%keywp)
   nullify(mpi_enreg%keywp)
 end if
 if (associated(mpi_enreg%trialproc)) then
   ABI_DEALLOCATE(mpi_enreg%trialproc)
   nullify(mpi_enreg%trialproc)
 end if
 if (associated(mpi_enreg%distrb_img)) then
   ABI_DEALLOCATE(mpi_enreg%distrb_img)
   nullify(mpi_enreg%distrb_img)
 end if
 if (associated(mpi_enreg%index_img)) then
   ABI_DEALLOCATE(mpi_enreg%index_img)
   nullify(mpi_enreg%index_img)
 end if

 call initmpi_seq(mpi_enreg)

end subroutine destroy_mpi_enreg
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/copy_mpi_enreg
!! NAME
!! copy_mpi_enreg
!!
!! FUNCTION
!!  Copy a mpi_enreg datastructure into another
!!
!! INPUTS
!!  opt_bandfft=if 1, the mpi_enreg_bandfft_kpt field has to be copied
!!              else, it is ignored
!!  MPI_enreg1<MPI_type>=input mpi_enreg datastructure
!!
!! OUTPUT
!!  MPI_enreg2<MPI_type>=output mpi_enreg datastructure
!!
!! PARENTS
!!      inwffil,m_fft_prof,m_wfs
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_mpi_enreg(MPI_enreg1,MPI_enreg2,opt_bandfft)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_mpi_enreg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: opt_bandfft
 type(MPI_type),intent(inout) :: mpi_enreg1,MPI_enreg2

!Local variables-------------------------------
!scalars
 integer :: ikpt,isppol,jkpt
!arrays

! *********************************************************************

!scalars
 mpi_enreg2%world_comm=mpi_enreg1%world_comm
 mpi_enreg2%world_group=mpi_enreg1%world_group
 mpi_enreg2%me=mpi_enreg1%me
 mpi_enreg2%nproc=mpi_enreg1%nproc
 mpi_enreg2%paral_compil=mpi_enreg1%paral_compil
 mpi_enreg2%paral_compil_mpio=mpi_enreg1%paral_compil_mpio
 mpi_enreg2%paral_compil_kpt=mpi_enreg1%paral_compil_kpt
 mpi_enreg2%paral_compil_fft=mpi_enreg1%paral_compil_fft
 mpi_enreg2%paral_spin=mpi_enreg1%paral_spin
 mpi_enreg2%paral_level=mpi_enreg1%paral_level
 mpi_enreg2%paralbd=mpi_enreg1%paralbd
 mpi_enreg2%me_group=mpi_enreg1%me_group
 mpi_enreg2%nproc_group=mpi_enreg1%nproc_group
 mpi_enreg2%me_fft=mpi_enreg1%me_fft
 mpi_enreg2%me_band=mpi_enreg1%me_band
 mpi_enreg2%nproc_fft=mpi_enreg1%nproc_fft
 mpi_enreg2%master_fft=mpi_enreg1%master_fft
 mpi_enreg2%paral_fft=mpi_enreg1%paral_fft
 mpi_enreg2%me_g0=mpi_enreg1%me_g0
 mpi_enreg2%num_group_fft=mpi_enreg1%num_group_fft
 mpi_enreg2%num_group=mpi_enreg1%num_group
 mpi_enreg2%nproc_per_kpt=mpi_enreg1%nproc_per_kpt
 mpi_enreg2%fft_master_comm=mpi_enreg1%fft_master_comm
 mpi_enreg2%fft_option_lob=mpi_enreg1%fft_option_lob
 mpi_enreg2%has_band_comm=mpi_enreg1%has_band_comm
 mpi_enreg2%flag_ind_kg_mpi_to_seq=mpi_enreg1%flag_ind_kg_mpi_to_seq
 mpi_enreg2%paral_compil_respfn=mpi_enreg1%paral_compil_respfn
 mpi_enreg2%me_respfn=mpi_enreg1%me_respfn
 mpi_enreg2%nproc_respfn=mpi_enreg1%nproc_respfn
 mpi_enreg2%my_respfn_group=mpi_enreg1%my_respfn_group
 mpi_enreg2%my_respfn_comm=mpi_enreg1%my_respfn_comm
 mpi_enreg2%respfn_master_group=mpi_enreg1%respfn_master_group
 mpi_enreg2%respfn_master_comm=mpi_enreg1%respfn_master_comm
 mpi_enreg2%ngroup_respfn=mpi_enreg1%ngroup_respfn
 mpi_enreg2%comm_respfn=mpi_enreg1%comm_respfn
 mpi_enreg2%ngfft3_ionic=mpi_enreg1%ngfft3_ionic
 mpi_enreg2%mode_para=mpi_enreg1%mode_para
 mpi_enreg2%commcart=mpi_enreg1%commcart
 mpi_enreg2%comm_band=mpi_enreg1%comm_band
 mpi_enreg2%comm_fft=mpi_enreg1%comm_fft
 mpi_enreg2%dimcart=mpi_enreg1%dimcart
 mpi_enreg2%nproc_band=mpi_enreg1%nproc_band
 mpi_enreg2%commcart_3d=mpi_enreg1%commcart_3d
 mpi_enreg2%comm_kpt=mpi_enreg1%comm_kpt
 mpi_enreg2%me_kpt=mpi_enreg1%me_kpt
 mpi_enreg2%nproc_kpt=mpi_enreg1%nproc_kpt
 mpi_enreg2%me_cart_2d=mpi_enreg1%me_cart_2d
 mpi_enreg2%nproc_atom=mpi_enreg1%nproc_atom
 mpi_enreg2%natom=mpi_enreg1%natom
 mpi_enreg2%comm_atom=mpi_enreg1%comm_atom
 mpi_enreg2%bandpp=mpi_enreg1%bandpp
 mpi_enreg2%paral_img=mpi_enreg1%paral_img
 mpi_enreg2%comm_img=mpi_enreg1%comm_img
 mpi_enreg2%me_img=mpi_enreg1%me_img
 mpi_enreg2%nproc_img=mpi_enreg1%nproc_img
 mpi_enreg2%group_one_img=mpi_enreg1%group_one_img
 mpi_enreg2%comm_one_img=mpi_enreg1%comm_one_img
 mpi_enreg2%me_one_img=mpi_enreg1%me_one_img
 mpi_enreg2%nproc_one_img=mpi_enreg1%nproc_one_img
 mpi_enreg2%nproc_spin=mpi_enreg1%nproc_spin
 mpi_enreg2%me_spin=mpi_enreg1%me_spin
 mpi_enreg2%commcart_4d=mpi_enreg1%commcart_4d
 mpi_enreg2%comm_spinfft=mpi_enreg1%comm_spinfft
 mpi_enreg2%me_cart_4d=mpi_enreg1%me_cart_4d
 mpi_enreg2%me_cart_3d=mpi_enreg1%me_cart_3d

!pointers
 if (associated(mpi_enreg1%band_comm)) then
   ABI_ALLOCATE(mpi_enreg2%band_comm,(size(mpi_enreg1%band_comm)))
   mpi_enreg2%band_comm=mpi_enreg1%band_comm
 else
   nullify(mpi_enreg2%band_comm)
 end if
 if (associated(mpi_enreg1%fft_comm)) then
   ABI_ALLOCATE(mpi_enreg2%fft_comm,(size(mpi_enreg1%fft_comm)))
   mpi_enreg2%fft_comm=mpi_enreg1%fft_comm
 else
   nullify(mpi_enreg2%fft_comm)
 end if
 if (associated(mpi_enreg1%proc_distrb)) then
   sz1=size(mpi_enreg1%proc_distrb,1)
   sz2=size(mpi_enreg1%proc_distrb,2)
   sz3=size(mpi_enreg1%proc_distrb,3)
   ABI_ALLOCATE(mpi_enreg2%proc_distrb,(sz1,sz2,sz3))
   mpi_enreg2%proc_distrb=mpi_enreg1%proc_distrb
 else
   nullify(mpi_enreg2%proc_distrb)
 end if
 if (associated(mpi_enreg1%kpt_comm)) then
   ABI_ALLOCATE(mpi_enreg2%kpt_comm,(size(mpi_enreg1%kpt_comm)))
   mpi_enreg2%kpt_comm=mpi_enreg1%kpt_comm
 else
   nullify(mpi_enreg2%kpt_comm)
 end if
 if (associated(mpi_enreg1%kptdstrb)) then
   sz1=size(mpi_enreg1%kptdstrb,1)
   sz2=size(mpi_enreg1%kptdstrb,2)
   sz3=size(mpi_enreg1%kptdstrb,3)
   ABI_ALLOCATE(mpi_enreg2%kptdstrb,(sz1,sz2,sz3))
   mpi_enreg2%kptdstrb=mpi_enreg1%kptdstrb
 else
   nullify(mpi_enreg2%kptdstrb)
 end if
 if (associated(mpi_enreg1%kptdstrbi)) then
   sz1=size(mpi_enreg1%kptdstrbi,1)
   sz2=size(mpi_enreg1%kptdstrbi,2)
   sz3=size(mpi_enreg1%kptdstrbi,3)
   ABI_ALLOCATE(mpi_enreg2%kptdstrbi,(sz1,sz2,sz3))
   mpi_enreg2%kptdstrbi=mpi_enreg1%kptdstrbi
 else
   nullify(mpi_enreg2%kptdstrbi)
 end if
 if (associated(mpi_enreg1%kpt_loc2fbz_sp)) then
   sz1=size(mpi_enreg1%kpt_loc2fbz_sp,1)-1
   sz2=size(mpi_enreg1%kpt_loc2fbz_sp,2)
   sz3=size(mpi_enreg1%kpt_loc2fbz_sp,3)
   ABI_ALLOCATE(mpi_enreg2%kpt_loc2fbz_sp,(0:sz1,1:sz2,1:sz3))
   mpi_enreg2%kpt_loc2fbz_sp=mpi_enreg1%kpt_loc2fbz_sp
 else
   nullify(mpi_enreg2%kpt_loc2fbz_sp)
 end if
 if (associated(mpi_enreg1%kpt_loc2ibz_sp)) then
   sz1=size(mpi_enreg1%kpt_loc2ibz_sp,1)-1
   sz2=size(mpi_enreg1%kpt_loc2ibz_sp,2)
   sz3=size(mpi_enreg1%kpt_loc2ibz_sp,3)
   ABI_ALLOCATE(mpi_enreg2%kpt_loc2ibz_sp,(0:sz1,1:sz2,1:sz3))
   mpi_enreg2%kpt_loc2ibz_sp=mpi_enreg1%kpt_loc2ibz_sp
 else
   nullify(mpi_enreg2%kpt_loc2ibz_sp)
 end if
 if (associated(mpi_enreg1%fmkmem)) then
   ABI_ALLOCATE(mpi_enreg2%fmkmem,(0:size(mpi_enreg1%fmkmem,1)-1))
   mpi_enreg2%fmkmem=mpi_enreg1%fmkmem
 else
   nullify(mpi_enreg2%fmkmem)
 end if
 if (associated(mpi_enreg1%mkmem)) then
   ABI_ALLOCATE(mpi_enreg2%mkmem,(0:size(mpi_enreg1%mkmem,1)-1))
   mpi_enreg2%mkmem=mpi_enreg1%mkmem
 else
   nullify(mpi_enreg2%mkmem)
 end if
 if (associated(mpi_enreg1%nplanes_fft)) then
   ABI_ALLOCATE(mpi_enreg2%nplanes_fft,(size(mpi_enreg1%nplanes_fft)))
   mpi_enreg2%nplanes_fft=mpi_enreg1%nplanes_fft
 else
   nullify(mpi_enreg2%nplanes_fft)
 end if
 if (associated(mpi_enreg1%ind_fft_planes)) then
   sz1=size(mpi_enreg1%ind_fft_planes,1)
   sz2=size(mpi_enreg1%ind_fft_planes,2)
   ABI_ALLOCATE(mpi_enreg2%ind_fft_planes,(sz1,sz2))
   mpi_enreg2%ind_fft_planes=mpi_enreg1%ind_fft_planes
 else
   nullify(mpi_enreg2%ind_fft_planes)
 end if
 if (associated(mpi_enreg1%respfn_group)) then
   ABI_ALLOCATE(mpi_enreg2%respfn_group,(size(mpi_enreg1%respfn_group)))
   mpi_enreg2%respfn_group=mpi_enreg1%respfn_group
 else
   nullify(mpi_enreg2%respfn_group)
 end if
 if (associated(mpi_enreg1%respfn_comm)) then
   ABI_ALLOCATE(mpi_enreg2%respfn_comm,(size(mpi_enreg1%respfn_comm)))
   mpi_enreg2%respfn_comm=mpi_enreg1%respfn_comm
 else
   nullify(mpi_enreg2%respfn_comm)
 end if
 if (associated(mpi_enreg1%nscatterarr)) then
   ABI_ALLOCATE(mpi_enreg2%nscatterarr,(0:size(mpi_enreg1%nscatterarr,1)-1,size(mpi_enreg1%nscatterarr,2)))
   mpi_enreg2%nscatterarr=mpi_enreg1%nscatterarr
 else
   nullify(mpi_enreg2%nscatterarr)
 end if
 if (associated(mpi_enreg1%ngatherarr)) then
   ABI_ALLOCATE(mpi_enreg2%ngatherarr,(0:size(mpi_enreg1%ngatherarr,1)-1,size(mpi_enreg1%ngatherarr,2)))
   mpi_enreg2%ngatherarr=mpi_enreg1%ngatherarr
 else
   nullify(mpi_enreg2%ngatherarr)
 end if
 if (associated(mpi_enreg1%sizecart)) then
   ABI_ALLOCATE(mpi_enreg2%sizecart,(size(mpi_enreg1%sizecart)))
   mpi_enreg2%sizecart=mpi_enreg1%sizecart
 else
   nullify(mpi_enreg2%sizecart)
 end if
 if (associated(mpi_enreg1%coords)) then
   ABI_ALLOCATE(mpi_enreg2%coords,(size(mpi_enreg1%coords)))
   mpi_enreg2%coords=mpi_enreg1%coords
 else
   nullify(mpi_enreg2%coords)
 end if
 if (associated(mpi_enreg1%atom_indx)) then
   ABI_ALLOCATE(mpi_enreg2%atom_indx,(size(mpi_enreg1%atom_indx)))
   mpi_enreg2%atom_indx=mpi_enreg1%atom_indx
 else
   nullify(mpi_enreg2%atom_indx)
 end if
 if (associated(mpi_enreg1%keywp)) then
   ABI_ALLOCATE(mpi_enreg2%keywp,(size(mpi_enreg1%keywp,1),size(mpi_enreg1%keywp,2)))
   mpi_enreg2%keywp=mpi_enreg1%keywp
 else
   nullify(mpi_enreg2%keywp)
 end if
 if (associated(mpi_enreg1%trialproc)) then
   ABI_ALLOCATE(mpi_enreg2%trialproc,(size(mpi_enreg1%trialproc)))
   mpi_enreg2%trialproc=mpi_enreg1%trialproc
 else
   nullify(mpi_enreg2%trialproc)
 end if
 if (associated(mpi_enreg1%distrb_img)) then
   ABI_ALLOCATE(mpi_enreg2%distrb_img,(size(mpi_enreg1%distrb_img)))
   mpi_enreg2%distrb_img=mpi_enreg1%distrb_img
 else
   nullify(mpi_enreg2%distrb_img)
 end if
 if (associated(mpi_enreg1%index_img)) then
   ABI_ALLOCATE(mpi_enreg2%index_img,(size(mpi_enreg1%index_img)))
   mpi_enreg2%index_img=mpi_enreg1%index_img
 else
   nullify(mpi_enreg2%index_img)
 end if

!Optional pointers
 if (opt_bandfft==0) then
   nullify(mpi_enreg2%tab_kpt_distrib)
   nullify(mpi_enreg2%bandfft_kpt)
 else if (opt_bandfft==1) then
   if (associated(mpi_enreg1%tab_kpt_distrib)) then
     ABI_ALLOCATE(mpi_enreg2%tab_kpt_distrib,(size(mpi_enreg1%tab_kpt_distrib)))
     mpi_enreg2%tab_kpt_distrib=mpi_enreg1%tab_kpt_distrib
   end if
   if (associated(mpi_enreg1%bandfft_kpt)) then
     ABI_ALLOCATE(mpi_enreg2%bandfft_kpt,(size(mpi_enreg1%bandfft_kpt)))
     do isppol=1,size(mpi_enreg1%proc_distrb,3)
       do ikpt=1,size(mpi_enreg1%proc_distrb,1)
         if(minval(abs(mpi_enreg1%proc_distrb(ikpt,:,isppol)-mpi_enreg1%me_kpt))/=0) then
           cycle
         end if
         jkpt=mpi_enreg1%tab_kpt_distrib(ikpt)
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%ind_kg_mpi_to_seq)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%ind_kg_mpi_to_seq)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%ind_kg_mpi_to_seq,(sz1))
           mpi_enreg2%bandfft_kpt(jkpt)%ind_kg_mpi_to_seq= &
&           mpi_enreg1%bandfft_kpt(jkpt)%ind_kg_mpi_to_seq
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%ind_kg_mpi_to_seq)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%kg_k_gather)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%kg_k_gather,1)
           sz2=size(mpi_enreg1%bandfft_kpt(jkpt)%kg_k_gather,2)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%kg_k_gather,(sz1,sz2))
           mpi_enreg2%bandfft_kpt(jkpt)%kg_k_gather= &
&           mpi_enreg1%bandfft_kpt(jkpt)%kg_k_gather
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%kg_k_gather)
         end if
         mpi_enreg2%bandfft_kpt(jkpt)%flag1_is_allocated=mpi_enreg1%bandfft_kpt(jkpt)%flag1_is_allocated
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%gbound)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%gbound,1)
           sz2=size(mpi_enreg1%bandfft_kpt(jkpt)%gbound,2)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%gbound,(sz1,sz2))
           mpi_enreg2%bandfft_kpt(jkpt)%gbound= &
&           mpi_enreg1%bandfft_kpt(jkpt)%gbound
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%gbound)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%recvcounts)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%recvcounts)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%recvcounts,(sz1))
           mpi_enreg2%bandfft_kpt(jkpt)%recvcounts= &
&           mpi_enreg1%bandfft_kpt(jkpt)%recvcounts
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%recvcounts)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%sendcounts)) then
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%sendcounts,(size(mpi_enreg1%bandfft_kpt(jkpt)%sendcounts)))
           mpi_enreg2%bandfft_kpt(jkpt)%sendcounts= &
&           mpi_enreg1%bandfft_kpt(jkpt)%sendcounts
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%sendcounts)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%rdispls)) then
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%rdispls,(size(mpi_enreg1%bandfft_kpt(jkpt)%rdispls)))
           mpi_enreg2%bandfft_kpt(jkpt)%rdispls= &
&           mpi_enreg1%bandfft_kpt(jkpt)%rdispls
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%rdispls)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%sdispls)) then
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%sdispls,(size(mpi_enreg1%bandfft_kpt(jkpt)%sdispls)))
           mpi_enreg2%bandfft_kpt(jkpt)%sdispls= &
&           mpi_enreg1%bandfft_kpt(jkpt)%sdispls
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%sdispls)
         end if
         mpi_enreg2%bandfft_kpt(jkpt)%flag2_is_allocated=mpi_enreg1%bandfft_kpt(jkpt)%flag2_is_allocated
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%ffnl_gather)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%ffnl_gather,1)
           sz2=size(mpi_enreg1%bandfft_kpt(jkpt)%ffnl_gather,2)
           sz3=size(mpi_enreg1%bandfft_kpt(jkpt)%ffnl_gather,3)
           sz4=size(mpi_enreg1%bandfft_kpt(jkpt)%ffnl_gather,4)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%ffnl_gather,(sz1,sz2,sz3,sz4))
           mpi_enreg2%bandfft_kpt(jkpt)%ffnl_gather= &
&           mpi_enreg1%bandfft_kpt(jkpt)%ffnl_gather
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%ffnl_gather)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%kinpw_gather)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%kinpw_gather)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%kinpw_gather,(sz1))
           mpi_enreg2%bandfft_kpt(jkpt)%kinpw_gather= &
&           mpi_enreg1%bandfft_kpt(jkpt)%kinpw_gather
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%kinpw_gather)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%ph3d_gather)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%ph3d_gather,1)
           sz2=size(mpi_enreg1%bandfft_kpt(jkpt)%ph3d_gather,2)
           sz3=size(mpi_enreg1%bandfft_kpt(jkpt)%ph3d_gather,3)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%ph3d_gather,(sz1,sz2,sz3))
           mpi_enreg2%bandfft_kpt(jkpt)%ph3d_gather= &
&           mpi_enreg1%bandfft_kpt(jkpt)%ph3d_gather
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%ph3d_gather)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%kpg_k_gather)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%kpg_k_gather,1)
           sz2=size(mpi_enreg1%bandfft_kpt(jkpt)%kpg_k_gather,2)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%kpg_k_gather,(sz1,sz2))
           mpi_enreg2%bandfft_kpt(jkpt)%kpg_k_gather= &
&           mpi_enreg1%bandfft_kpt(jkpt)%kpg_k_gather
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%kpg_k_gather)
         end if
         mpi_enreg2%bandfft_kpt(jkpt)%flag3_is_allocated=mpi_enreg1%bandfft_kpt(jkpt)%flag3_is_allocated
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%kg_k_gather_sym)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%kg_k_gather_sym,1)
           sz2=size(mpi_enreg1%bandfft_kpt(jkpt)%kg_k_gather_sym,2)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%kg_k_gather_sym,(sz1,sz2))
           mpi_enreg2%bandfft_kpt(jkpt)%kg_k_gather_sym= &
&           mpi_enreg1%bandfft_kpt(jkpt)%kg_k_gather_sym
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%kg_k_gather_sym)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%rdispls_sym)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%rdispls_sym)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%rdispls_sym,(sz1))
           mpi_enreg2%bandfft_kpt(jkpt)%rdispls_sym= &
&           mpi_enreg1%bandfft_kpt(jkpt)%rdispls_sym
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%rdispls_sym)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%recvcounts_sym)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%recvcounts_sym)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%recvcounts_sym,(sz1))
           mpi_enreg2%bandfft_kpt(jkpt)%recvcounts_sym= &
&           mpi_enreg1%bandfft_kpt(jkpt)%recvcounts_sym
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%recvcounts_sym)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%recvcounts_sym_tot)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%recvcounts_sym_tot)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%recvcounts_sym_tot,(sz1))
           mpi_enreg2%bandfft_kpt(jkpt)%recvcounts_sym_tot= &
&           mpi_enreg1%bandfft_kpt(jkpt)%recvcounts_sym_tot
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%recvcounts_sym_tot)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%sdispls_sym)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%sdispls_sym)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%sdispls_sym,(sz1))
           mpi_enreg2%bandfft_kpt(jkpt)%sdispls_sym= &
&           mpi_enreg1%bandfft_kpt(jkpt)%sdispls_sym
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%sdispls_sym)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%sendcounts_sym)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%sendcounts_sym)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%sendcounts_sym,(sz1))
           mpi_enreg2%bandfft_kpt(jkpt)%sendcounts_sym= &
&           mpi_enreg1%bandfft_kpt(jkpt)%sendcounts_sym
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%sendcounts_sym)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%sendcounts_sym_all)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%sendcounts_sym_all)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%sendcounts_sym_all,(sz1))
           mpi_enreg2%bandfft_kpt(jkpt)%sendcounts_sym_all= &
&           mpi_enreg1%bandfft_kpt(jkpt)%sendcounts_sym_all
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%sendcounts_sym_all)
         end if
         if (associated(mpi_enreg1%bandfft_kpt(jkpt)%tab_proc)) then
           sz1=size(mpi_enreg1%bandfft_kpt(jkpt)%tab_proc)
           ABI_ALLOCATE(mpi_enreg2%bandfft_kpt(jkpt)%tab_proc,(sz1))
           mpi_enreg2%bandfft_kpt(jkpt)%tab_proc= &
           mpi_enreg1%bandfft_kpt(jkpt)%tab_proc
         else
           nullify(mpi_enreg2%bandfft_kpt(jkpt)%tab_proc)
         end if
       end do
     end do
   end if
 end if

end subroutine copy_mpi_enreg
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/my_indeces
!! NAME
!! my_indeces
!!
!! FUNCTION
!!  Helper function returning useful local indeces from the global indeces (ikpt,isppol).
!!  It works only in the case of k-point parallelism or sequential run.
!!
!! INPUTS
!!  MPI_enreg<MPI_type>=Datatype gathering information on the parallelism.
!!  ikpt=The global index of the k-point.
!!  isppol=The global index for the spin.
!!  nkpt=The total number of k-points (global)
!!  nsppol=The total number of spins (global)
!!  nspinor=The number of spinorial components (on current proc)
!!  npwarr(nkpt)=Global array storing the number of planewaves at each k-point.
!!  nband(nkpt*nsppol)=Global array stoting the number of bands treated at each k-point and spin.
!!
!! OUTPUT
!!  kindex= kindex+1 is the local index in the cg array defining the beginning of the block
!!  of the wavefunctions with (ikpt,isppol) indeces
!!  bdtot_index
!!  ibg=Local index
!!  ikg= ikg+1 is the local index in the array kg_k
!!  ierr=Status error.
!!    == 0  => This processor has this (k,s)
!!    /= 0  => This node is not treating this (k,s). (kindex,ibg,ikg) are set to HUGE(0)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine my_indeces(MPI_enreg,ikpt,isppol,nkpt,nsppol,nspinor,npwarr,nband,kindex,bdtot_index,ibg,ikg,ierr)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'my_indeces'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikpt,isppol,nkpt,nsppol,nspinor
 integer,intent(out) :: kindex,ibg,ikg,bdtot_index
 type(MPI_type),intent(in) :: MPI_enreg
!arrays
 integer,intent(in) :: nband(nkpt*nsppol),npwarr(nkpt)
 integer,intent(out) :: ierr


!Local variables-------------------------------
!scalars
 integer :: iktot,isp,my_rank,nband_k,nprocs,npw_k
!arrays

! *********************************************************************
 if (MPI_enreg%paral_compil_kpt/=1) then
   MSG_ERROR(" %paral_compil_kpt/=1")
 end if

 if (MPI_enreg%paral_img==1) then
   nprocs  = MPI_enreg%nproc_one_img
   my_rank = MPI_enreg%me_one_img
 else
   nprocs  = MPI_enreg%nproc
   my_rank = MPI_enreg%me
 end if

 ierr=1
 if (nprocs==1) then ! Calculate local indeces, cannot use %proc_distrb as it is not allocated.

   kindex=0; bdtot_index=0; ibg=0
   isp_loop1: do isp=1,nsppol
     ikg=0
     do iktot=1,nkpt
       if (iktot==ikpt.and.isp==isppol) then
         EXIT isp_loop1
         ierr=0
       end if
       nband_k = nband(iktot+(isp-1)*nkpt)
       npw_k   = npwarr(iktot)

       kindex  = kindex + npw_k*nspinor*nband_k
       bdtot_index=bdtot_index+nband_k
       ibg     = ibg + nspinor*nband_k
       ikg     = ikg+npw_k
     end do
   end do isp_loop1

 else ! parallel case: calculate local indeces.

   kindex=0; bdtot_index=0; ibg=0
   isp_loop2: do isp=1,nsppol
     ikg=0
     do iktot=1,nkpt
       if (MINVAL(ABS(MPI_enreg%proc_distrb(iktot,:,isp)-my_rank))==0) then  ! FIXME this cannot be tested in seq.
         if (iktot==ikpt.and.isp==isppol) then
           ierr=0
           EXIT isp_loop2
         end if
         nband_k = nband(iktot+(isp-1)*nkpt)
         npw_k   = npwarr(iktot)

         kindex  = kindex + npw_k*nspinor*nband_k
         bdtot_index=bdtot_index+nband_k
         ibg     = ibg + nspinor*nband_k
         ikg     = ikg+npw_k
       end if
     end do
   end do isp_loop2
 end if

 if (ierr/=0) then ! This will lead to a SIGFAULT if the indeces are used in the caller.
   kindex = HUGE(0)
   bdtot_index = HUGE(0)
   ibg    = HUGE(0)
   ikg    = HUGE(0)
 end if

end subroutine my_indeces
!!***
