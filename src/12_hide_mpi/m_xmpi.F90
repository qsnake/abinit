!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xmpi
!! NAME
!!  m_xmpi
!!
!! FUNCTION
!!  This module provides MPI named constants, tools for inquiring the MPI environment
!!  and a set of generic interfaces wrapping the most commonly used MPI primitives.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG, MB, XG, YP, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_xmpi


 use defs_basis,   only : i4b, i8b, sp, dp, spc, dpc, std_out, std_out_default, tmp_unit, fnlen

#if defined HAVE_MPI2
 use mpi
#endif

 implicit none

 private
!!***

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

#ifdef HAVE_MPI
 ! MPI constants used in abinit. Make sure that a corresponding fake value is provided for the sequential version.
 integer,public :: xmpi_world          = MPI_COMM_WORLD
 integer,public,parameter :: xmpi_self           = MPI_COMM_SELF
 integer,public,parameter :: xmpi_undefined      = MPI_UNDEFINED
 integer,public,parameter :: xmpi_undefined_rank = MPI_UNDEFINED  ! MPI_UNDEFINED_RANK is not portable.
 integer,public,parameter :: xmpi_comm_null      = MPI_COMM_NULL
 integer,public,parameter :: xmpi_group_null     = MPI_GROUP_NULL
 !integer,public,parameter :: xmpi_group_null     = MPI_GROUP_NULL ! invalid handle.
 !integer,public,parameter :: xmpi_group_empty    = MPI_GROUP_EMPTY ! valid handle.
 integer,public,parameter :: xmpi_any_source     = MPI_ANY_SOURCE
 integer,public,parameter :: xmpi_msg_len= MPI_MAX_ERROR_STRING ! Length of fortran string used to store MPI error strings.
#else
 ! Fake replacements for the sequential version.
 integer,public,parameter :: xmpi_world          = 0
 integer,public,parameter :: xmpi_self           = 0
 integer,public,parameter :: xmpi_undefined      =-32765
 integer,public,parameter :: xmpi_undefined_rank =-32766
 integer,public,parameter :: xmpi_comm_null      = 0
 integer,public,parameter :: xmpi_group_null     = 0
 !integer,public,parameter :: xmpi_group_null    =
 !integer,public,parameter :: xmpi_group_empty   =
 integer,public,parameter :: xmpi_any_source     = 0
 integer,public,parameter :: xmpi_msg_len=1000
#endif

 integer,save,private  :: xmpi_tag_ub=32767
 ! The tag upper bound value must be at least 32767. An MPI implementation is free to make
 ! the value of MPI_TAG_UB larger than this hence xmpi_tag_ub is redefined when MPI is init in xmpi_init.

 ! Size in bytes of the entries used in MPI datatypes.
 integer,save,public :: xmpi_bsize_ch =0
 integer,save,public :: xmpi_bsize_int=0
 integer,save,public :: xmpi_bsize_sp =0
 integer,save,public :: xmpi_bsize_dp =0
 integer,save,public :: xmpi_bsize_spc=0
 integer,save,public :: xmpi_bsize_dpc=0

 ! kind of the offset used for MPI-IO.
#ifdef HAVE_MPI_IO
 integer,public,parameter :: xmpi_offset_kind =MPI_OFFSET_KIND
 integer,public,parameter :: xmpi_address_kind=MPI_ADDRESS_KIND
#else
 integer,public,parameter :: xmpi_offset_kind =i4b
 integer,public,parameter :: xmpi_address_kind=i4b
#endif

 ! The byte size and the MPI type of the Fortran record marker.
 ! These quantities are compiler-dependent and are initalized in xmpi_init (only if MPI-IO is on).
 integer,save,public :: xmpio_bsize_frm   =0
 integer,save,public :: xmpio_mpi_type_frm=0

 ! Options used for the MPI-IO wrappers used in abinit.
 integer,public,parameter :: xmpio_at    =1
 integer,public,parameter :: xmpio_at_all=2

!----------------------------------------------------------------------
!!***

!!****t* m_xmpi/comm_t
!! NAME
!!  comm_t
!!
!! FUNCTION
!!  Datatype used to store data associated to an MPI communicator.
!!
!! SOURCE

 type,public :: comm_t

  integer :: id = xmpi_undefined
  ! The MPI communicator identifier.

  integer :: my_rank = xmpi_undefined_rank
  ! The rank of the node inside comm.

  integer :: master = 0
  ! The rank of master node.

  integer :: nprocs = xmpi_undefined
  ! The number of processors in the communicator.

  !integer,pointer :: ranks_in_world(:)   SET2NULL
  ! The MPI ranks in MPI_COMM_WORLD of the nodes beloging to the communicator.
 end type comm_t
!!***

! Public procedures.
 public :: xmpi_init      ! Initialize the MPI environment.
 public :: xmpi_end       ! Terminate the MPI environment.
 public :: xmpi_abort     ! Abort all tasks in a group.
 public :: xmpi_show_info ! Printout of the basic variables stored in this module (useful for debugging).
 public :: xcomm_rank     ! The rank of the node inside the communicator.
 public :: xcomm_size     ! The number of processors inside the communicator.
 public :: xcomm_free     ! Marks the communicator for deallocation.
 public :: xbarrier_mpi   ! Hides MPI_BARRIER from MPI library.
 public :: xmpi_name      ! Return the name of this processor (usually the hostname).
 public :: xerror_string  ! Return a string describing the error from ierr.
 public :: xmpi_split_work
 public :: xmpi_distab

 interface xcomm_free
   module procedure xcomm_free_0D
   module procedure xcomm_free_1D
   module procedure xcomm_free_2D
   module procedure xcomm_free_3D
 end interface xcomm_free

 interface xmpi_split_work
   module procedure xmpi_split_work_i4b
 end interface xmpi_split_work

 public :: xmpi_split_work2_i4b
 public :: xmpi_split_work2_i8b
 !public :: xmpi_split_work2
 !
 ! g95@green v0.93 is not able to resolve the interface.
 ! For the time being, the use of this generic interface has been disabled.
 !interface xmpi_split_work2
 !  module procedure xmpi_split_work2_i4b
 !  module procedure xmpi_split_work2_i8b
 !end interface xmpi_split_work2

 interface xmpi_distab
   module procedure xmpi_distab_4D
 end interface xmpi_distab

 ! MPI generic interfaces.
 public :: xallgather_mpi
 public :: xallgatherv_mpi
 public :: xalltoall_mpi
 public :: xalltoallv_mpi
 public :: xcast_mpi
 public :: xexch_mpi
 public :: xgather_mpi
 public :: xgatherv_mpi
 public :: xmax_mpi
 public :: xmin_mpi
 public :: xrecv_mpi
 public :: xscatterv_mpi
 public :: xsend_mpi
 public :: xsum_master
 public :: xsum_mpi

#ifdef HAVE_MPI_IO
 public :: xmpio_type_struct
 public :: xmpio_get_info_frm
 public :: xmpio_read_frm
 public :: xmpio_write_frm
 public :: xmpio_check_frmarkers
 public :: xmpio_max_address
 public :: xmpio_write_frmarkers
 public :: xmpio_create_fsubarray_2D
 public :: xmpio_create_fsubarray_3D
 public :: xmpio_create_fsubarray_4D
 public :: xmpio_create_fherm_packed
 public :: xmpio_create_coldistr_from_fpacked
 public :: xmpio_create_coldistr_from_fp3blocks
#endif

!----------------------------------------------------------------------

interface xallgather_mpi
  module procedure xallgather_mpi_int
  module procedure xallgather_mpi_char
  module procedure xallgather_mpi_int1d
  module procedure xallgather_mpi_dp1d
  module procedure xallgather_mpi_dp2d
  module procedure xallgather_mpi_dp3d
  module procedure xallgather_mpi_dp4d
end interface xallgather_mpi

!----------------------------------------------------------------------

interface xallgatherv_mpi
  module procedure xallgatherv_mpi_int2d
  module procedure xallgatherv_mpi_int
  module procedure xallgatherv_mpi_dp
  module procedure xallgatherv_mpi_dp2d
  module procedure xallgatherv_mpi_dp3d
  module procedure xallgatherv_mpi_dp4d
end interface xallgatherv_mpi

!----------------------------------------------------------------------

interface xalltoall_mpi
  module procedure xalltoall_mpi_dp2d
end interface xalltoall_mpi

!----------------------------------------------------------------------

interface xalltoallv_mpi
  module procedure xalltoallv_mpi_dp2d
  module procedure xalltoallv_mpi_int2d
  module procedure xalltoallv_mpi_dp1d
  module procedure xalltoallv_mpi_dp1d2
end interface xalltoallv_mpi

!----------------------------------------------------------------------

interface xcast_mpi
  module procedure xcast_mpi_intv
  module procedure xcast_mpi_int1d
  module procedure xcast_mpi_int2d
  module procedure xcast_mpi_int3d
  module procedure xcast_mpi_dpv
  module procedure xcast_mpi_dp1d
  module procedure xcast_mpi_dp2d
  module procedure xcast_mpi_dp3d
  module procedure xcast_mpi_dp4d
  module procedure xcast_mpi_spv
  module procedure xcast_mpi_sp1d
  module procedure xcast_mpi_sp2d
  module procedure xcast_mpi_sp3d
  module procedure xcast_mpi_sp4d
  module procedure xcast_mpi_cplxv
  module procedure xcast_mpi_cplx1d
  module procedure xcast_mpi_cplx2d
  module procedure xcast_mpi_cplx3d
  module procedure xcast_mpi_cplx4d
  module procedure xcast_mpi_dcv
  module procedure xcast_mpi_dc1d
  module procedure xcast_mpi_dc2d
  module procedure xcast_mpi_dc3d
  module procedure xcast_mpi_dc4d
  module procedure xcast_mpi_ch0d
  module procedure xcast_mpi_ch1d
end interface xcast_mpi

!----------------------------------------------------------------------

interface xexch_mpi
  module procedure xexch_mpi_intn
  module procedure xexch_mpi_int2d
  module procedure xexch_mpi_dpn
  module procedure xexch_mpi_dp2d
  module procedure xexch_mpi_dp3d
  module procedure xexch_mpi_dp4d_tag
  module procedure xexch_mpi_dp5d_tag
  module procedure xexch_mpi_spc_1d
  module procedure xexch_mpi_dpc_1d
  module procedure xexch_mpi_dpc_2d
end interface xexch_mpi

!----------------------------------------------------------------------

interface xgather_mpi
  module procedure xgather_mpi_int
  module procedure xgather_mpi_int2d
  module procedure xgather_mpi_dp
  module procedure xgather_mpi_dp2d
  module procedure xgather_mpi_dp3d
  module procedure xgather_mpi_dp4d
end interface xgather_mpi

!----------------------------------------------------------------------

interface xgatherv_mpi
  module procedure xgatherv_mpi_int
  module procedure xgatherv_mpi_int2d
  module procedure xgatherv_mpi_dp
  module procedure xgatherv_mpi_dp2d
  module procedure xgatherv_mpi_dp3d
  module procedure xgatherv_mpi_dp4d
end interface xgatherv_mpi

!----------------------------------------------------------------------

interface xmax_mpi
  module procedure xmax_mpi_intv
  module procedure xmax_mpi_dpv
end interface xmax_mpi

!----------------------------------------------------------------------

interface xmin_mpi
  module procedure xmin_mpi_intv
  module procedure xmin_mpi_dpv
end interface xmin_mpi

!----------------------------------------------------------------------

interface xrecv_mpi
  module procedure xrecv_mpi_intv
  module procedure xrecv_mpi_dp2d
  module procedure xrecv_mpi_dp3d
end interface xrecv_mpi

!----------------------------------------------------------------------

interface xscatterv_mpi
  module procedure xscatterv_mpi_int
  module procedure xscatterv_mpi_int2d
  module procedure xscatterv_mpi_dp
  module procedure xscatterv_mpi_dp2d
  module procedure xscatterv_mpi_dp3d
  module procedure xscatterv_mpi_dp4d
end interface xscatterv_mpi

!----------------------------------------------------------------------

interface xsend_mpi
  module procedure xsend_mpi_intv
  module procedure xsend_mpi_dp2d
  module procedure xsend_mpi_dp3d
end interface xsend_mpi

!----------------------------------------------------------------------

interface xsum_master
  module procedure xsum_master_int
  module procedure xsum_master_dp1d
  module procedure xsum_master_dp2d
  module procedure xsum_master_dp3d
  module procedure xsum_master_dp4d
  module procedure xsum_master_dp5d
  module procedure xsum_master_dp6d
  module procedure xsum_master_dp7d
  module procedure xsum_master_int4d
  module procedure xsum_master_c1cplx
  module procedure xsum_master_c2cplx
  module procedure xsum_master_c3cplx
  module procedure xsum_master_c4cplx
  module procedure xsum_master_c5cplx
  module procedure xsum_master_c1dpc
  module procedure xsum_master_c2dpc
  module procedure xsum_master_c3dpc
  module procedure xsum_master_c4dpc
  module procedure xsum_master_c5dpc
end interface xsum_master

!----------------------------------------------------------------------

interface xsum_mpi
  module procedure xsum_mpi_int
  module procedure xsum_mpi_intv
  module procedure xsum_mpi_intv2
  module procedure xsum_mpi_intn
  module procedure xsum_mpi_int2t
  module procedure xsum_mpi_int2d
  module procedure xsum_mpi_int3d
  module procedure xsum_mpi_int4d
  module procedure xsum_mpi_dp
  module procedure xsum_mpi_dpvt
  module procedure xsum_mpi_dpv
  module procedure xsum_mpi_dpn
  module procedure xsum_mpi_dp2d
  module procedure xsum_mpi_dp3d
  module procedure xsum_mpi_dp4d
  module procedure xsum_mpi_dp5d
  module procedure xsum_mpi_dp6d
  module procedure xsum_mpi_dp7d
  module procedure xsum_mpi_dp2t
  module procedure xsum_mpi_dp3d2t
  module procedure xsum_mpi_dp4d2t
  module procedure xsum_mpi_c0dc
  module procedure xsum_mpi_c1dc
  module procedure xsum_mpi_c2dc
  module procedure xsum_mpi_c3dc
  module procedure xsum_mpi_c4dc
  module procedure xsum_mpi_c5dc
  module procedure xsum_mpi_c6dc
  module procedure xsum_mpi_c1cplx
  module procedure xsum_mpi_c2cplx
  module procedure xsum_mpi_c3cplx
  module procedure xsum_mpi_c4cplx
  module procedure xsum_mpi_c5cplx
  module procedure xsum_mpi_c6cplx
  module procedure xsum_mpi_log1d
  module procedure xsum_mpi_log2d
  module procedure xsum_mpi_log3d
end interface xsum_mpi

!----------------------------------------------------------------------

!!****t* m_xmpi/mpi_ffh_t
!! NAME
!!  mpi_ffh_t
!!
!! FUNCTION
!!  Datatype used to (read|write) Fortran files that are record-based using
!!  MPI-IO routines that use C-streams.
!!
!! SOURCE

 type,public :: mpi_ffh_t

  integer :: comm
  ! MPI communicator.

  integer :: fh
   ! MPI file handler used to access the file with MPI-IO.

  integer(XMPI_OFFSET_KIND)  :: off_rw
   ! offset used to (read|write) data.

  integer(XMPI_OFFSET_KIND)  :: off_1frm
   ! offset of the first Fortran record marker.

  integer(XMPI_OFFSET_KIND) :: bsize_frec
   ! byte size of the the record.

  integer(XMPI_OFFSET_KIND) :: bsize_frm
   ! byte size of Fortran file record markers.

  integer(XMPI_OFFSET_KIND) :: bsize_EOF
   ! byte size of Fortran end-of-file.

  integer :: mpi_type_frm
   ! MPI Datatype for Fortran record markers....

  !integer :: offset_mpi_type
   ! MPI Datatype for INTEGER(kind=MPI_OFFSET_KIND).

  integer :: advancing
  ! 1 if IO is done in advance mode.

  logical :: is_eof=.FALSE.
  ! .TRUE. if the end-of-file has been reached.

  character(len=fnlen) :: fname
  ! Name of the file.

 end type mpi_ffh_t
!!***
 integer :: ABI_ALLOC_STAT, ABI_ALLOC_SIZE

!----------------------------------------------------------------------

CONTAINS  !===========================================================
!!***

!!****f* m_xmpi/xmpi_init
!! NAME
!!  xmpi_init
!!
!! FUNCTION
!!  Hides MPI_INIT from MPI library. Perform the initialization of some basic variables
!!  used by the MPI routines employed in abinit.
!!
!! INPUTS
!!  None
!!
!! PARENTS
!!      aim,anaddb,conducti,cut3d,fftprof,kss2wfk,lwf,mpi_enreg_tools,mrgddb
!!      mrggkk,mrgscr,newsp,optic,ujdet
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xmpi_init()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpi_init'
!End of the abilint section

 implicit none

!Local variables-------------------
 integer :: mpi_err
#if defined HAVE_MPI
 integer :: attribute_val
 logical :: lflag
#endif

! *************************************************************************

 mpi_err=0
#if defined HAVE_MPI
 call MPI_INIT(mpi_err)

 ! Deprecated in MPI2 but not all MPI2 implementations provide MPI_Comm_get_attr !
 call MPI_ATTR_GET(xmpi_world, MPI_TAG_UB, attribute_val, lflag, mpi_err)
 !call MPI_Comm_get_attr(xmpi_world, MPI_TAG_UB, attribute_val, lflag, mpi_err)

 if (lflag) xmpi_tag_ub = attribute_val

!  Define type values.
 call MPI_TYPE_SIZE(MPI_CHARACTER,xmpi_bsize_ch,mpi_err)

 call MPI_TYPE_SIZE(MPI_INTEGER,xmpi_bsize_int,mpi_err)

 call MPI_TYPE_SIZE(MPI_REAL,xmpi_bsize_sp,mpi_err)

 call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,xmpi_bsize_dp,mpi_err)

 call MPI_TYPE_SIZE(MPI_COMPLEX,xmpi_bsize_spc,mpi_err)

 call MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX,xmpi_bsize_dpc,mpi_err)
 !
 ! Find the byte size of Fortran record marker used in MPI-IO routines.
 ! FIXME This call makes the code crash on several buildbot slaves,
 !  likely due to some configuration problem in the MPI libraries.
 call xmpio_get_info_frm(xmpio_bsize_frm, xmpio_mpi_type_frm, xmpi_world)
#endif

end subroutine xmpi_init
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_end
!! NAME
!!  xmpi_end
!!
!! FUNCTION
!!  Hides MPI_FINALIZE from MPI library.
!!
!! INPUTS
!!  None
!!
!! PARENTS
!!      abinit,aim,anaddb,conducti,cut3d,fftprof,kss2wfk,leave_new,lwf,mrgddb
!!      mrggkk,mrgscr,newsp,optic,ujdet
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xmpi_end()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpi_end'
!End of the abilint section

 implicit none

!Local variables-------------------
 integer :: mpi_err

! *************************************************************************

 mpi_err=0
#if defined HAVE_MPI
 call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)  !  Needed by some HPC architectures (MT, 20110315)
 call MPI_FINALIZE(mpi_err)
#endif

end subroutine xmpi_end
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_abort
!! NAME
!!  xmpi_abort
!!
!! FUNCTION
!!  Hides MPI_ABORT from MPI library.
!!
!! INPUTS
!!  comm=communicator of tasks to abort.
!!  mpi_err=Error code to return to invoking environment.
!!
!! PARENTS
!!      initmpi_grid,leave_new,leave_test,m_initcuda,m_xmpi
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xmpi_abort(comm,mpi_err)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpi_abort'
!End of the abilint section

 implicit none

!Arguments-------------------------
 integer,optional,intent(in) :: comm,mpi_err

!Local variables-------------------
 integer :: ierr,my_comm,my_errorcode,ilen,ierr2
 logical :: testopen
 character(len=xmpi_msg_len) :: mpi_msg_error

! *************************************************************************

 ierr=0
 my_comm = xmpi_world; if (PRESENT(comm)) my_comm = comm

#if defined HAVE_MPI
 my_errorcode=MPI_ERR_UNKNOWN; if (PRESENT(mpi_err)) my_errorcode=mpi_err

 call MPI_ERROR_STRING(my_errorcode, mpi_msg_error, ilen, ierr2)

 if (ilen>xmpi_msg_len) write(std_out,*)" WARNING: MPI message has been truncated!"
 if (ierr2/=MPI_SUCCESS) then
   write(std_out,*)" WARNING: MPI_ERROR_STRING returned ierr2= ",ierr2
 else
   write(std_out,*)" MPI_ERROR_STRING: ",TRIM(mpi_msg_error)
 end if

 inquire(std_out,OPENED=testopen)
 if (testopen) close(std_out)

 call MPI_ABORT(my_comm,my_errorcode,ierr)
#endif

end subroutine xmpi_abort
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_show_info
!! NAME
!!  xmpi_show_info
!!
!! FUNCTION
!!  Printout of the most important variables stored in this module (useful for debugging).
!!
!! INPUTS
!!  unt=Unit number for formatted output.
!!
!! PARENTS
!!      leave_new
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xmpi_show_info(unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpi_show_info'
!End of the abilint section

 implicit none

!Arguments-------------------------
 integer,optional,intent(in) :: unit

!Local variables-------------------
 integer :: my_unt

! *************************************************************************

 !@m_xmpi
 my_unt = std_out; if (PRESENT(unit)) my_unt=unit

#if defined HAVE_MPI1
  write(my_unt,*)" ==== Using MPI-1 specifications ==== "
#endif
#if defined HAVE_MPI2
  write(my_unt,*)" ==== Using MPI-2 specifications ==== "
#endif

#if defined HAVE_MPI_IO
  write(my_unt,*)" MPI-IO support is ON"
#else
  write(my_unt,*)" MPI-IO support is OFF"
#endif

#if defined HAVE_MPI
 write(my_unt,*)" xmpi_tag_ub ................ ",xmpi_tag_ub
 write(my_unt,*)" xmpi_bsize_ch .............. ",xmpi_bsize_ch
 write(my_unt,*)" xmpi_bsize_int ............. ",xmpi_bsize_int
 write(my_unt,*)" xmpi_bsize_sp .............. ",xmpi_bsize_sp
 write(my_unt,*)" xmpi_bsize_dp .............. ",xmpi_bsize_dp
 write(my_unt,*)" xmpi_bsize_spc ............. ",xmpi_bsize_spc
 write(my_unt,*)" xmpi_bsize_dpc ............. ",xmpi_bsize_dpc
 write(my_unt,*)" xmpio_bsize_frm ............ ",xmpio_bsize_frm
 write(my_unt,*)" xmpi_address_kind .......... ",xmpi_address_kind
 write(my_unt,*)" xmpi_offset_kind ........... ",xmpi_offset_kind
#endif

end subroutine xmpi_show_info
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xcomm_rank
!! NAME
!!  xcomm_rank
!!
!! FUNCTION
!!  Hides MPI_COMM_RANK from MPI library.
!!
!! INPUTS
!!  spaceComm=MPI communicator.
!!
!! OUTPUT
!!  xcomm_rank=The rank of the node inside spaceComm
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function xcomm_rank(spaceComm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcomm_rank'
!End of the abilint section

 implicit none

!Arguments-------------------------
 integer,intent(in) :: spaceComm
 integer :: xcomm_rank

!Local variables-------------------
 integer :: mpi_err

! *************************************************************************

 mpi_err=0; xcomm_rank=0
#if defined HAVE_MPI
 if ( spaceComm/=MPI_COMM_NULL ) then
   call MPI_COMM_RANK(spaceComm,xcomm_rank,mpi_err)
 end if
#endif

end function xcomm_rank
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xcomm_size
!! NAME
!!  xcomm_size
!!
!! FUNCTION
!!  Hides MPI_COMM_SIZE from MPI library.
!!
!! INPUTS
!!  spaceComm=MPI communicator.
!!
!! OUTPUT
!!  xcomm_size=The number of processors inside spaceComm.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function xcomm_size(spaceComm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcomm_size'
!End of the abilint section

 implicit none

!Arguments-------------------------
 integer,intent(in) :: spaceComm
 integer :: xcomm_size

!Local variables-------------------------------
!scalars
 integer :: mpi_err

! *************************************************************************

 mpi_err=0; xcomm_size=1
#if defined HAVE_MPI
 if ( spaceComm/=MPI_COMM_NULL ) then
   call MPI_COMM_SIZE(spaceComm,xcomm_size,mpi_err)
 end if
#endif

end function xcomm_size
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xcomm_free_0D
!! NAME
!!  xcomm_free_0D
!!
!! FUNCTION
!!  Hides MPI_COMM_FREE from MPI library.
!!
!! INPUTS
!!  spaceComm=MPI communicator.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_xmpi
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xcomm_free_0D(spaceComm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcomm_free_0D'
!End of the abilint section

 implicit none

!Arguments-------------------------
 integer,intent(inout) :: spaceComm

!Local variables-------------------------------
!scalars
 integer :: mpi_err

! *************************************************************************

 mpi_err=0
#if defined HAVE_MPI
 if ( spaceComm/=MPI_COMM_NULL ) then
   call MPI_COMM_FREE(spaceComm,mpi_err)
 end if
#endif

end subroutine xcomm_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xcomm_free_1D
!! NAME
!!  xcomm_free_1D
!!
!! FUNCTION
!!  Hides MPI_COMM_FREE from MPI library. Target 1D arrays
!!
!! INPUTS
!!  comms(:)=MPI communicators
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xcomm_free_1D(comms)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcomm_free_1D'
!End of the abilint section

 implicit none

!Arguments-------------------------
 integer,intent(inout) :: comms(:)

!Local variables-------------------------------
!scalars
 integer :: ii

! *************************************************************************

 do ii=LBOUND(comms,DIM=1),UBOUND(comms,DIM=1)
   call xcomm_free_0D(comms(ii))
 end do

end subroutine xcomm_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xcomm_free_2D
!! NAME
!!  xcomm_free_2D
!!
!! FUNCTION
!!  Hides MPI_COMM_FREE from MPI library. Target 2D arrays
!!
!! INPUTS
!!  comms=MPI communicator.
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xcomm_free_2D(comms)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcomm_free_2D'
!End of the abilint section

 implicit none

!Arguments-------------------------
 integer,intent(inout) :: comms(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,jj

! *************************************************************************

 do jj=LBOUND(comms,DIM=2),UBOUND(comms,DIM=2)
   do ii=LBOUND(comms,DIM=1),UBOUND(comms,DIM=1)
     call xcomm_free_0D(comms(ii,jj))
   end do
 end do

end subroutine xcomm_free_2D
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xcomm_free_3D
!! NAME
!!  xcomm_free_3D
!!
!! FUNCTION
!!  Hides MPI_COMM_FREE from MPI library. Target 3D arrays
!!
!! INPUTS
!!  comms=MPI communicator.
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xcomm_free_3D(comms)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcomm_free_3D'
!End of the abilint section

 implicit none

!Arguments-------------------------
 integer,intent(inout) :: comms(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk

! *************************************************************************

 do kk=LBOUND(comms,DIM=3),UBOUND(comms,DIM=3)
   do jj=LBOUND(comms,DIM=2),UBOUND(comms,DIM=2)
     do ii=LBOUND(comms,DIM=1),UBOUND(comms,DIM=1)
       call xcomm_free_0D(comms(ii,jj,kk))
     end do
   end do
 end do

end subroutine xcomm_free_3D
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xbarrier_mpi
!! NAME
!!  xbarrier_mpi
!!
!! FUNCTION
!!  Hides MPI_BARRIER from MPI library.
!!
!! INPUTS
!!
!! PARENTS
!!      atomden,calc_optical_mels,calc_sigx_me,cohsex_me,datafordmft
!!      defs_scalapack,denfgr,dmft_solve,exc_build_block,exc_diago
!!      exc_diago_driver,exc_iterative_diago,exc_spectra,fermisolverec
!!      filnam_comm,gstateimg,haydock,haydock_psherm,ks_ddiago,m_bse_io,m_green
!!      m_io_kss,m_melemts,m_screen,m_screening,m_self,m_shexc,m_vcoul,m_wfs
!!      outkss,outwf,pawmkaewf,respfn,setup_bse,sigma,tddft,vtorho,vtorhorec
!!      wfk_read_ene
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xbarrier_mpi(spaceComm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xbarrier_mpi'
!End of the abilint section

 implicit none

!Arguments-------------------------
 integer,intent(in) :: spaceComm

!Local variables-------------------
 integer   :: ier
#if defined HAVE_MPI
 integer :: nprocs
#endif

! *************************************************************************

 ier = 0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nprocs,ier)
   if(nprocs>1)then
     call MPI_BARRIER(spaceComm,ier)
   end if
 end if
#endif

end subroutine xbarrier_mpi
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_name
!! NAME
!!  xmpi_name
!!
!! FUNCTION
!!  Hides MPI_GET_PROCESSOR_NAME from MPI library.
!!
!! OUTPUT
!!  name= the host name transformed to integer variable.
!!  mpi_err=Status error.
!!
!! PARENTS
!!      m_gpu_detect
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xmpi_name(name_ch, mpi_err)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpi_name'
!End of the abilint section

 implicit none

!Arguments-------------------------
 integer,intent(out) ::  mpi_err
 character(20),intent(out) :: name_ch

!Local variables-------------------
 integer :: name,len
! character(len=MPI_MAX_PROCESSOR_NAME) :: name_ch

! *************************************************************************
!Get the name of this processor (usually the hostname)

 name    = 0
 mpi_err = 0

#if defined HAVE_MPI
 call MPI_GET_PROCESSOR_NAME(name_ch, len, mpi_err)
 name_ch = trim(name_ch)

#else
 name_ch ='0'
#endif

end subroutine xmpi_name
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xerror_string
!! NAME
!!  xerror_string
!!
!! FUNCTION
!!  Hides MPI_ERROR_STRING from MPI library.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      print_ierr
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xerror_string(mpi_err,err_string,ilen,ierror)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xerror_string'
!End of the abilint section

 implicit none

!Arguments-------------------------
 integer,intent(in) :: mpi_err
 integer,intent(out) :: ilen,ierror
 character(len=*),intent(out) :: err_string

! *************************************************************************

 ilen=0
#if defined HAVE_MPI
 call MPI_Error_string(mpi_err,err_string,ilen,ierror)
#else
 ierror=1
 err_string="Sorry, no MPI_Error_string routine is available to interpret the error message"
#endif

end subroutine xerror_string
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_split_work_i4b
!! NAME
!!  split_work_i4b
!!
!! FUNCTION
!!  Splits the number of tasks, ntasks, among nprocs processors. Used for the MPI parallelization of simple loops.
!!
!! INPUTS
!!  ntasks=number of tasks
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  my_start,my_stop= indices defining the initial and final task for this processor
!!  warn_msg=String containing a possible warning message if the distribution is not optima.
!!  ierr=Error status
!!    +1 if ntasks is not divisible by nprocs.
!!    +2 if ntasks>nprocs.
!!
!! NOTES
!!  If nprocs>ntasks then :
!!    my_start=ntasks+1
!!    my_stop=ntask
!!
!!  In this particular case, loops of the form
!!
!!  do ii=my_start,my_stop
!!   ...
!!  end do
!!
!!  are not executed. Moreover allocation such as foo(my_start:my_stop) will generate a zero-sized array.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xmpi_split_work_i4b(ntasks,comm,my_start,my_stop,warn_msg,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpi_split_work_i4b'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)  :: ntasks,comm
 integer,intent(out) :: my_start,my_stop,ierr
 character(len=500) :: warn_msg

!Local variables-------------------------------
 integer :: res,nprocs,my_rank,block_p1,block

! *************************************************************************

 nprocs  = xcomm_size(comm)
 my_rank = xcomm_rank(comm)

 block   = ntasks/nprocs
 res     = MOD(ntasks,nprocs)
 block_p1= block+1

 warn_msg = ""; ierr=0
 if (res/=0) then
   write(warn_msg,'(4a,i5,a,i4)')ch10,&
&   ' xmpi_split_work : ',ch10,&
&   '  The number of tasks= ',ntasks,' is not divisible by nprocs= ',nprocs
   ierr=1
 end if
 if (block==0) then
   write(warn_msg,'(4a,i4,a,i5,3a)')ch10,&
&   ' xmpi_split_work : ',ch10,&
&   ' The number of processors= ',nprocs,' larger than number of tasks= ',ntasks,ch10,&
&   ' This is a waste ',ch10
    ierr=2
 end if

 if (my_rank<res) then
   my_start =  my_rank   *block_p1+1
   my_stop  = (my_rank+1)*block_p1
 else
   my_start = res*block_p1 + (my_rank-res  )*block + 1
   my_stop  = res*block_p1 + (my_rank-res+1)*block
 end if

end subroutine xmpi_split_work_i4b
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_split_work2_i4b
!! NAME
!!  xmpi_split_work2_i4b
!!
!! FUNCTION
!!  Splits a number of tasks, ntasks, among nprocs processors.
!!  The output arrays istart(1:nprocs) and istop(1:nprocs)
!!  report the starting and final task index for each CPU.
!!  Namely CPU with rank ii has to perform all the tasks between
!!  istart(ii+1) and istop(ii+1). Note the Fortran convention of using
!!  1 as first index of the array.
!!  Note, moreover, that if a proc has rank>ntasks then :
!!   istart(rank+1)=ntasks+1
!!   istop(rank+1)=ntask
!!
!!  In this particular case, loops of the form
!!
!!  do ii=istart(rank),istop(rank)
!!   ...
!!  end do
!!
!!  are not executed. Moreover allocation such as foo(istart(rank):istop(rank))
!!  will generate a zero-sized array
!!
!! INPUTS
!!  ntasks= number of tasks
!!  nprocs=Number of processors.
!!
!! OUTPUT
!!  istart(nprocs),istop(nprocs)= indices defining the initial and final task for each processor
!!  ierr=Error status.
!!  warn_msg=String containing the warning message.
!!    +1 if ntasks is not divisible by nprocs.
!!    +2 if ntasks>nprocs.
!!
!! PARENTS
!!      exc_build_block,m_screening,screening
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xmpi_split_work2_i4b(ntasks,nprocs,istart,istop,warn_msg,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpi_split_work2_i4b'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)  :: ntasks,nprocs
 integer,intent(out) :: ierr
 integer,intent(inout) :: istart(nprocs),istop(nprocs)
 character(len=500),intent(out) :: warn_msg

!Local variables-------------------------------
 integer :: res,irank,block,block_tmp

! *************************************************************************

 block_tmp = ntasks/nprocs
 res       = MOD(ntasks,nprocs)
 block     = block_tmp+1

 warn_msg = ""; ierr=0
 if (res/=0) then
   write(warn_msg,'(a,i0,a,i0,3a)')&
&   ' The number of tasks = ',ntasks,' is not divisible by nprocs = ',nprocs,ch10,&
&   ' parallelism is not efficient ',ch10
   ierr=+1
 end if
 !
 if (block_tmp==0) then
   write(warn_msg,'(a,i0,a,i0,3a)')&
&   ' The number of processors = ',nprocs,' is larger than number of tasks =',ntasks,ch10,&
&   ' This is a waste ',ch10
   ierr=+2
 end if

 do irank=0,nprocs-1
   if (irank<res) then
     istart(irank+1)= irank   *block+1
     istop (irank+1)=(irank+1)*block
   else
     istart(irank+1)=res*block+(irank-res  )*block_tmp+1
     istop (irank+1)=res*block+(irank-res+1)*block_tmp
   end if
 end do

end subroutine xmpi_split_work2_i4b
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_split_work2_i8b
!! NAME
!!  xmpi_split_work2_i8b
!!
!! FUNCTION
!!  Same as xmpi_split_work2_i8b but accepts 8 bytes integer.
!!
!! INPUTS
!!  ntasks= number of tasks
!!  nprocs=Number of processors.
!!
!! OUTPUT
!!  istart(nprocs),istop(nprocs)= indices defining the initial and final task for each processor
!!  ierr=Error status.
!!  warn_msg=String containing the warning message.
!!    +1 if ntasks is not divisible by nprocs.
!!    +2 if ntasks>nprocs.
!!
!! PARENTS
!!      exc_build_block
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xmpi_split_work2_i8b(ntasks,nprocs,istart,istop,warn_msg,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpi_split_work2_i8b'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)  :: nprocs
 integer(i8b),intent(in)  :: ntasks
 integer,intent(out) :: ierr
 integer(i8b),intent(inout) :: istart(nprocs),istop(nprocs)
 character(len=500),intent(out) :: warn_msg

!Local variables-------------------------------
 integer(i8b) :: res,irank,block,block_tmp

! *************************************************************************

 block_tmp = ntasks/nprocs
 res       = MOD(ntasks,INT(nprocs,KIND=i8b))
 block     = block_tmp+1

 warn_msg = ""; ierr=0
 if (res/=0) then
   write(warn_msg,'(a,i0,a,i0,3a)')&
&   ' The number of tasks = ',ntasks,' is not divisible by nprocs = ',nprocs,ch10,&
&   ' parallelism is not efficient ',ch10
   ierr=+1
 end if
 !
 if (block_tmp==0) then
   write(warn_msg,'(a,i0,a,i0,3a)')&
&   ' The number of processors = ',nprocs,' is larger than number of tasks =',ntasks,ch10,&
&   ' This is a waste ',ch10
   ierr=+2
 end if

 do irank=0,nprocs-1
   if (irank<res) then
     istart(irank+1)= irank   *block+1
     istop (irank+1)=(irank+1)*block
   else
     istart(irank+1)=res*block+(irank-res  )*block_tmp+1
     istop (irank+1)=res*block+(irank-res+1)*block_tmp
   end if
 end do

end subroutine xmpi_split_work2_i8b
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_distab_4D
!! NAME
!!  xmpi_distab_4D
!!
!! FUNCTION
!!  Fill table defining the distribution of the tasks according to the number of processors involved in the
!!  calculation. For each set of indeces, the table contains the rank of the node in the MPI communicator.
!!
!! INPUTS
!!  nprocs=The number of processors performing the calculation in parallel.
!!
!! OUTPUT
!!  task_distrib(:,:,:,:) = Contains the rank of the node that is taking care of this particular set of loop indeces.
!!  Tasks are distributed across the nodes in column-major order.
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

subroutine xmpi_distab_4D(nprocs,task_distrib)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpi_distab_4D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nprocs
!arrays
 integer,intent(inout) :: task_distrib(:,:,:,:)

!Local variables ------------------------------
!scalars
 integer :: ii,jj,n1,n2,n3,n4,ntasks,irank,remainder,ntpblock
 integer,allocatable :: list(:)

!************************************************************************

 n1= SIZE(task_distrib,DIM=1)
 n2= SIZE(task_distrib,DIM=2)
 n3= SIZE(task_distrib,DIM=3)
 n4= SIZE(task_distrib,DIM=4)
 ntasks = n1*n2*n3*n4

 ABI_ALLOCATE(list,(ntasks))
 list=-999

 ntpblock  = ntasks/nprocs
 remainder = MOD(ntasks,nprocs)

 if (ntpblock==0) then ! nprocs > ntasks
   do ii=1,ntasks
     list(ii) = ii-1
   end do
 else
   ii=1
   do irank=nprocs-1,0,-1 ! If remainder/=0, master will get less tasks.
     jj = ii+ntpblock-1
     if (remainder>0) then
       jj=jj+1
       remainder = remainder-1
     end if
     list(ii:jj)=irank
     ii=jj+1
   end do
 end if

 task_distrib = RESHAPE(list,(/n1,n2,n3,n4/))

 if ( ANY(task_distrib==-999) ) then
   write(std_out,*)" task_distrib == -999"
   call xmpi_abort()
 end if

 ABI_DEALLOCATE(list)

end subroutine xmpi_distab_4D
!!***

!----------------------------------------------------------------------

! Include files providing wrappers for some of the most commonly used MPI primitives.

#include "xallgather_mpi.F90"

#include "xallgatherv_mpi.F90"

#include "xalltoall_mpi.F90"

#include "xalltoallv_mpi.F90"

#include "xcast_mpi.F90"

#include "xexch_mpi.F90"

#include "xgather_mpi.F90"

#include "xgatherv_mpi.F90"

#include "xmax_mpi.F90"

#include "xmin_mpi.F90"

#include "xrecv_mpi.F90"

#include "xscatterv_mpi.F90"

#include "xsend_mpi.F90"

#include "xsum_master.F90"

#include "xsum_mpi.F90"

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_type_struct
!! NAME
!!  xmpio_type_struct
!!
!! FUNCTION
!!  Some highly non-standard MPI implementations support MPI-IO without
!!  implementing the full set of MPI-2 extensions.
!!  This wrapper will call the obsolete MPI_TYPE_STRUCT if MPI_TYPE_CREATE_STRUCT
!!  is not supported. Note that MPI_TYPE_STRUCT requires the displacement arrays
!!  to be an array of default integers whereas the argument block_displ is an array of kind XMPI_ADDRESS_KIND.
!!  The routine will abort if the displacement cannot be represented with a default integer.
!!
!! INPUTS
!! ncount= number of blocks (integer) --- also number of entries in arrays array_of_types, array_of_displacements and array_of_blocklengths
!! array_of_blocklength(ncount)=number of elements in each block (array of integer)
!! array_of_displacements(ncount)=byte displacement of each block (array of integer)
!! array_of_types(ncount)=type of elements in each block (array of handles to datatype objects)
!!
!! OUTPUT
!! new_type=new datatype (handle)
!! mpi_err=MPI status error
!!
!! PARENTS
!!      defs_scalapack,m_xmpi
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_type_struct(ncount,block_length,block_displ,block_type,new_type,mpi_err)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_type_struct'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncount
 integer,intent(out) :: new_type,mpi_err
!arrays
 integer,intent(in) :: block_length(ncount),block_type(ncount)
 integer(XMPI_ADDRESS_KIND),intent(in) :: block_displ(ncount)

!Local variables-------------------
#ifndef HAVE_MPI_TYPE_CREATE_STRUCT
 integer,allocatable :: tmp_displ(:)
#endif

!************************************************************************

#ifdef HAVE_MPI_TYPE_CREATE_STRUCT
 call MPI_TYPE_CREATE_STRUCT(ncount,block_length,block_displ,block_type,new_type,mpi_err)
#else

 ABI_ALLOCATE(tmp_displ,(ncount))
 tmp_displ = block_displ
 if (ANY(block_displ > HUGE(tmp_displ(1)) ))then
   write(std_out,*)" byte displacement cannot be represented with a default integer"
   call xmpi_abort()
 end if

 call MPI_TYPE_STRUCT(ncount,block_length,block_displ,block_type,new_type,mpi_err)
 ABI_DEALLOCATE(tmp_displ)
#endif

end subroutine xmpio_type_struct
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpio_get_info_frm
!! NAME
!!  xmpio_marker_info
!!
!! FUNCTION
!!  Return the byte size of the Fortran record and its corresponding MPI_type (compiler-dependent).
!!  These two values are needed to access sequential binary Fortran files with MPI/IO routines where
!!  C-streams are used.
!!
!! INPUTS
!! comm=MPI communicator. Only master will find the values for the record marker. The results
!! are then broadcast to all the other nodes in comm.
!!
!! OUTPUT
!!  bsize_frm=Byte size of the Fortran record marker.
!!  mpi_type_frm=MPI type of the marker.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xmpio_get_info_frm(bsize_frm,mpi_type_frm,comm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_get_info_frm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 integer,intent(out) :: mpi_type_frm,bsize_frm

!Local variables-------------------------------
 integer :: my_rank
#if defined HAVE_MPI_IO
!scalars
 integer,parameter :: master=0
 integer :: spt,ept,ii
 integer :: f90_unt,ierr,iimax,mpio_fh,bsize_int,mpi_err
 integer(XMPI_OFFSET_KIND) :: offset,rml
 character(len=fnlen) :: fname
!character(len=500) :: msg
 logical :: file_exists
!arrays
 integer :: xvals(2),ivals(100),read_5ivals(5),ref_5ivals(5)
 integer :: rm_lengths(4)=(/4,8,2,16/)
 integer :: statux(MPI_STATUS_SIZE)
 real(dp) :: xrand(fnlen)
#endif

!************************************************************************

 bsize_frm=0; mpi_type_frm=0

 my_rank = xcomm_rank(comm)

#if defined HAVE_MPI_IO
 if ( my_rank == master ) then
   ! Fortran scratch files cannot have a name so have to generate a random one.
   ! cannot use pick_aname since it is higher level.
   fname = "__MPI_IO_FRM__"
   spt=LEN(fname); ept=spt

   inquire(file=fname,exist=file_exists)

   do while (file_exists)
     call RANDOM_NUMBER(xrand(spt:ept))
     xrand = xrand*127
     do ii=spt,ept
      fname(ii:ii) = ACHAR(NINT(xrand(ii)))
     end do
     ept = MIN(ept+1,fnlen)
     inquire(file=fname,exist=file_exists)
   end do
   !
   ! Write five integers on the binary file open in Fortran mode, then try
   ! to reread the values with MPI-IO using different offsets for the record marker.
   !
   f90_unt = tmp_unit
   open(unit=f90_unt,file=fname,status="new",form="unformatted",iostat=ierr)

   ref_5ivals = (/(ii, ii=5,9)/)
   ivals = HUGE(1); ivals(5:9)=ref_5ivals
   write(f90_unt) ivals
   close(f90_unt)

   call MPI_FILE_OPEN(xmpi_self, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh,mpi_err)
   !ABI_CHECK_MPI(ierr,"Open")

   iimax=3 ! Define number of INTEGER types to be tested
#if defined HAVE_FC_INT_QUAD
   iimax=4
#endif
   !
   ! Try to read ivals(5:9) from file.
   ii=0; bsize_frm=-1
   call MPI_TYPE_SIZE(MPI_INTEGER,bsize_int,mpi_err)

   do while (bsize_frm<=0 .and. ii<iimax)
     ii=ii+1
     rml = rm_lengths(ii)
     offset = rml + 4 * bsize_int
     call MPI_FILE_READ_AT(mpio_fh,offset,read_5ivals,5,MPI_INTEGER,statux,mpi_err)
     !write(std_out,*)read_5ivals
     if ( mpi_err==MPI_SUCCESS .and. ALL(read_5ivals==ref_5ivals) ) bsize_frm=rml
   end do

   if (ii==iimax.and.bsize_frm<=0) then
     if (iimax>=4) then
       write(std_out,'(3a)') &
&        ' Your architecture is not able to handle 16, 8, 4 or 2-bytes FORTRAN file record markers! ',ch10,&
&        ' You cannot use ABINIT and MPI/IO.'
     else
       write(std_out,'(3a)') &
&        ' Your architecture is not able to handle 8, 4 or 2-bytes FORTRAN file record markers! ',ch10,&
&        ' You cannot use ABINIT and MPI/IO.'
     end if
     !MSG_ERROR(msg)
     call xmpi_abort()
   else
     !write(std_out,'(a,i0)')' Detected FORTRAN record mark length: ',bsize_frm
   end if

   call MPI_FILE_CLOSE(mpio_fh, mpi_err)
   !
   ! Select MPI datatype corresponding to the Fortran marker.
   SELECT CASE (bsize_frm)
   CASE (4)
     mpi_type_frm=MPI_INTEGER4
   CASE (8)
     mpi_type_frm=MPI_INTEGER8
#if defined HAVE_FC_INT_QUAD
   CASE (16)
     mpi_type_frm=MPI_INTEGER16
#endif
   CASE (2)
     mpi_type_frm=MPI_INTEGER2
   CASE DEFAULT
     write(std_out,'(a,i0)')" Wrong bsize_frm: ",bsize_frm
     !MSG_ERROR(msg)
     call xmpi_abort()
   END SELECT

   open(unit=f90_unt,file=fname)
   close(f90_unt,status="delete")
 end if
 !
 ! Broadcast data.
 xvals = (/bsize_frm,mpi_type_frm/)
 call xcast_mpi(xvals,master,comm,mpi_err)

 bsize_frm    = xvals(1)
 mpi_type_frm = xvals(2)
#endif

end subroutine xmpio_get_info_frm
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/xmpio_read_frm
!! NAME
!!  xmpio_read_frm
!!
!! FUNCTION
!!  Read the content of a single record marker in a FORTRAN file at a given offset using MPI-IO.
!!  the file pointer is modified according to the value of advance.
!!
!! INPUTS
!!  mpi_fh=MPI-IO file handler.
!!  at_option=
!!         xmpio_at     ==> for reading by current proc.
!!         xmpio_at_all ==> for collective reading.
!!  offset=MPI/IO file pointer
!!  [advance]=By default the routine will move the file pointer to the next record.
!!    advance=.FALSE. can be used so that the next read will continue picking information
!!    off of the currect record.
!!
!! OUTPUT
!!  fmarker=Content of the Fortran record marker.
!!  mpi_err= MPI error code
!!
!! SIDE EFFECTS
!!  offset=
!!     input: file pointer used to access the Fortran marker.
!!     output: new offset updated after the reading, depending on advance.
!!
!! PARENTS
!!      exc_diago,exc_iterative_diago,m_bse_io,m_header,m_io_screening
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

#if defined HAVE_MPI_IO

subroutine xmpio_read_frm(mpi_fh,offset,at_option,fmarker,mpi_err,advance)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_read_frm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpi_fh,at_option
 integer(XMPI_OFFSET_KIND),intent(inout) :: offset
 integer(XMPI_OFFSET_KIND),intent(out) :: fmarker
 integer,intent(out) :: mpi_err
 logical,optional,intent(in) :: advance

!Local variables-------------------------------
!scalars
 integer :: bsize_frm,mpi_type_frm
 integer*2 :: delim_record2
 integer*4 :: delim_record4
 integer*8 :: delim_record8
#if defined HAVE_FC_INT_QUAD
 integer*16 :: delim_record16
#endif
 character(len=500) :: msg
!arrays
 integer :: statux(MPI_STATUS_SIZE)

!************************************************************************

 bsize_frm    = xmpio_bsize_frm    ! Byte size of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

 SELECT CASE (at_option)

 CASE (xmpio_at)

   if (bsize_frm==4) then
     call MPI_FILE_READ_AT(mpi_fh,offset,delim_record4,1,mpi_type_frm,statux,mpi_err)
     fmarker = delim_record4
   else if (bsize_frm==8) then
     call MPI_FILE_READ_AT(mpi_fh,offset,delim_record8,1,mpi_type_frm,statux,mpi_err)
     fmarker = delim_record8
#if defined HAVE_FC_INT_QUAD
   else if (bsize_frm==16) then
     call MPI_FILE_READ_AT(mpi_fh,offset,delim_record16,1,mpi_type_frm,statux,mpi_err)
     fmarker = delim_record16
#endif
   else if (bsize_frm==2) then
     call MPI_FILE_READ_AT(mpi_fh,offset,delim_record2 ,1,mpi_type_frm,statux,mpi_err)
     fmarker = delim_record2
   else
     !MSG_BUG('Wrong record marker length!')
     call xmpi_abort()
   end if

 CASE (xmpio_at_all)

   if (bsize_frm==4) then
     call MPI_FILE_READ_AT_ALL(mpi_fh,offset,delim_record4 ,1,mpi_type_frm,statux,mpi_err)
     fmarker = delim_record4
   else if (bsize_frm==8) then
     call MPI_FILE_READ_AT_ALL(mpi_fh,offset,delim_record8 ,1,mpi_type_frm,statux,mpi_err)
     fmarker = delim_record8
#if defined HAVE_FC_INT_QUAD
   else if (bsize_frm==16) then
     call MPI_FILE_READ_AT_ALL(mpi_fh,offset,delim_record16,1,mpi_type_frm,statux,mpi_err)
     fmarker = delim_record16
#endif
   else if (bsize_frm==2) then
     call MPI_FILE_READ_AT_ALL(mpi_fh,offset,delim_record2 ,1,mpi_type_frm,statux,mpi_err)
     fmarker = delim_record2
   else
     !MSG_BUG('Wrong record marker length !')
     call xmpi_abort()
   end if

 CASE DEFAULT
   write(msg,"(a,i0)")" Wrong value for at_option: ",at_option
   !MSG_BUG(msg)
   call xmpi_abort()
 END SELECT

 if (PRESENT(advance)) then
   if (advance) then
     offset = offset + fmarker + 2*bsize_frm ! Move the file pointer to the next record.
   else
     offset = offset + bsize_frm  ! Move the pointer after the marker.
   end if
 else
   offset = offset + fmarker + 2*bsize_frm
 end if

end subroutine xmpio_read_frm
!!***

#endif

!------------------------------------------------------------------------------------

!!****f* m_wffile/xmpio_write_frm
!! NAME
!!  xmpio_write_frm
!!
!! FUNCTION
!!  Write a single record marker in a FORTRAN file at a given offset using MPI-IO.
!!  The file pointer is modified according to the value of advance.
!!
!! INPUTS
!!  mpi_fh=MPI-IO file handler.
!!  at_option=
!!         xmpio_at     ==> for reading by current proc.
!!         xmpio_at_all ==> for collective reading.
!!  offset=MPI/IO file pointer
!!  [advance]=By default the routine will move the file pointer to the next record.
!!    advance=.FALSE. can be used so that the next read will continue writing data
!!    on the currect record.
!!
!! OUTPUT
!!  mpi_err= error code
!!
!! SIDE EFFECTS
!!  offset=
!!     input: file pointer used to access the Fortran marker.
!!     output: new offset updated after the writing, depending on advance.
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

#if defined HAVE_MPI_IO

subroutine xmpio_write_frm(mpi_fh,offset,at_option,fmarker,mpi_err,advance)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_write_frm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) ::mpi_fh
 integer,intent(in) :: at_option
 integer(XMPI_OFFSET_KIND),intent(in) :: fmarker
 integer(XMPI_OFFSET_KIND),intent(inout) :: offset
 integer,intent(out) :: mpi_err
 logical,optional,intent(in) :: advance

!Local variables-------------------------------
!scalars
 integer :: bsize_frm,mpi_type_frm
 integer*2  :: delim_record2
 integer*4  :: delim_record4
 integer*8  :: delim_record8
#if defined HAVE_FC_INT_QUAD
 integer*16 :: delim_record16
#endif
 character(len=500) :: msg
!arrays
 integer :: statux(MPI_STATUS_SIZE)

!************************************************************************

 bsize_frm    = xmpio_bsize_frm    ! Byte size of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

 SELECT CASE (at_option)

 CASE (xmpio_at)
   if (bsize_frm==4) then
     delim_record4 = fmarker
     call MPI_FILE_WRITE_AT(mpi_fh,offset,delim_record4 ,1,mpi_type_frm,statux,mpi_err)
   else if (bsize_frm==8) then
     delim_record8 = fmarker
     call MPI_FILE_WRITE_AT(mpi_fh,offset,delim_record8 ,1,mpi_type_frm,statux,mpi_err)
#if defined HAVE_FC_INT_QUAD
   else if (bsize_frm==16) then
     delim_record16 = fmarker
     call MPI_FILE_WRITE_AT(mpi_fh,offset,delim_record16,1,mpi_type_frm,statux,mpi_err)
#endif
   else if (bsize_frm==2) then
     delim_record2 = fmarker
     call MPI_FILE_WRITE_AT(mpi_fh,offset,delim_record2, 1,mpi_type_frm,statux,mpi_err)
   else
     !MSG_BUG('Wrong record marker length!')
     call xmpi_abort()
   end if

 CASE (xmpio_at_all)
   if (bsize_frm==4) then
     delim_record4 = fmarker
     call MPI_FILE_WRITE_AT_ALL(mpi_fh,offset,delim_record4 ,1,mpi_type_frm,statux,mpi_err)
   else if (bsize_frm==8) then
     delim_record8 = fmarker
     call MPI_FILE_WRITE_AT_ALL(mpi_fh,offset,delim_record8 ,1,mpi_type_frm,statux,mpi_err)
#if defined HAVE_FC_INT_QUAD
   else if (bsize_frm==16) then
     delim_record16 = fmarker
     call MPI_FILE_WRITE_AT_ALL(mpi_fh,offset,delim_record16,1,mpi_type_frm,statux,mpi_err)
#endif
   else if (bsize_frm==2) then
     delim_record2 = fmarker
     call MPI_FILE_WRITE_AT_ALL(mpi_fh,offset,delim_record2 ,1,mpi_type_frm,statux,mpi_err)
   else
     !MSG_BUG('Wrong record marker length!')
     call xmpi_abort()
   end if

 CASE DEFAULT
   write(msg,"(a,i0)")" Wrong value for at_option: ",at_option
   !MSG_BUG(msg)
   call xmpi_abort()
 END SELECT

 if (PRESENT(advance)) then
   if (advance) then
     offset = offset + fmarker + 2*bsize_frm  ! Move the file pointer to the next record.
   else
     offset = offset + bsize_frm              ! Move the pointer after the marker.
   end if
 else
   offset = offset + fmarker + 2*bsize_frm
 end if

end subroutine xmpio_write_frm
!!***

#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_fsubarray_2D
!! NAME
!!  xmpio_create_fsubarray_2D
!!
!! FUNCTION
!!
!! INPUTS
!!  array_of_sizes(2)=number of elements of type old_type in each dimension of the full array (array of positive integers)
!!  array_of_subsizes(2)=number of elements of type old_type in each dimension of the subarray (array of positive integers)
!!  array_of_starts(2)=starting coordinates of the subarray in each dimension (array of nonnegative integers >=1, <=array_of_sizes)
!!  old_type=Old MPI type.
!!
!! OUTPUT
!!  my_offpad=Offset to be added to the file pointer giving the position of the first Fortran record
!!    marker individuating the beginning of the matrix. (lets call it "base").
!!    Each node should (read|write) using my_offset = base + my_offpad.
!!    my_offpad is used so that one can safely change the way the fileview is generated (for example
!!    to make it more efficient) without having to change the client code.
!!  new_type=New MPI type.
!!  mpi_err= MPI error code
!!
!! PARENTS
!!      exc_build_block,exc_iterative_diago
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

#if defined HAVE_MPI_IO

subroutine xmpio_create_fsubarray_2D(array_of_sizes,array_of_subsizes,array_of_starts,old_type,new_type,my_offpad,mpi_err)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_create_fsubarray_2D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: old_type
 integer(XMPI_OFFSET_KIND),intent(out) :: my_offpad
 integer,intent(out) :: mpi_err,new_type
!arrays
 integer,intent(in) :: array_of_sizes(2),array_of_subsizes(2),array_of_starts(2)
!Local variables-------------------------------
!scalars
 integer :: bsize_frm,bsize_old,nx,ny,stride_x
 integer :: column_type,ldx,st_x,st_y
 !character(len=500) :: msg

!************************************************************************

 ! Byte size of the Fortran record marker.
 bsize_frm = xmpio_bsize_frm

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpi_err)
 !
 ! Number of columns and rows of the submatrix.
 nx   = array_of_subsizes(1)
 ny   = array_of_subsizes(2)

 ldx = array_of_sizes(1)
 st_x = array_of_starts(1)
 st_y = array_of_starts(2)

 ! The view starts at the first element of the submatrix.
 my_offpad = (st_x-1)*bsize_old + (st_y-1)*(ldx*bsize_old+2*xmpio_bsize_frm) + xmpio_bsize_frm

 ! Byte size of the Fortran record + the two markers.
 stride_x = ldx*bsize_old + 2*bsize_frm

 call MPI_Type_contiguous(nx,old_type,column_type,mpi_err)

 call MPI_Type_hvector(ny,1,stride_x,column_type,new_type,mpi_err)

 call MPI_TYPE_COMMIT(new_type,mpi_err)

end subroutine xmpio_create_fsubarray_2D
!!***

#endif


!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_fsubarray_3D
!! NAME
!!  xmpio_create_fsubarray_3D
!!
!! FUNCTION
!!
!! INPUTS
!!  array_of_sizes(3)=number of elements of type old_type in each dimension of the full array (array of positive integers)
!!  array_of_subsizes(3)=number of elements of type old_type in each dimension of the subarray (array of positive integers)
!!  old_type=Old MPI type.
!!
!! OUTPUT
!!  mpi_err= MPI error code
!!  new_type=New MPI type.
!!
!! PARENTS
!!      m_io_screening
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

#if defined HAVE_MPI_IO

subroutine xmpio_create_fsubarray_3D(array_of_sizes,array_of_subsizes,old_type,new_type,mpi_err)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_create_fsubarray_3D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: old_type
 integer,intent(out) :: mpi_err,new_type
!arrays
 integer,intent(in) :: array_of_sizes(3),array_of_subsizes(3)
!Local variables-------------------------------
!scalars
 integer :: bsize_frm,bsize_old,nx,ny,nz,stride_x
 integer :: column_type,plane_type,ldx,ldy,ldz
 !character(len=500) :: msg
!arrays

!************************************************************************

 bsize_frm    = xmpio_bsize_frm    ! Byte size of the Fortran record marker.

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpi_err)
 !
 ! Number of columns and rows of the submatrix.
 nx = array_of_subsizes(1)
 ny = array_of_subsizes(2)
 nz = array_of_subsizes(3)

 ldx = array_of_sizes(1)
 ldy = array_of_sizes(2)
 ldz = array_of_sizes(3)

 ! Byte size of the Fortran record + the two markers.
 stride_x = ldx*bsize_old + 2*bsize_frm

 call MPI_Type_contiguous(nx,old_type,column_type,mpi_err)

 call MPI_Type_hvector(ny,1,stride_x,column_type,plane_type,mpi_err)

 call MPI_Type_hvector(nz,1,ldy*stride_x,plane_type,new_type,mpi_err)

 call MPI_TYPE_COMMIT(new_type,mpi_err)

end subroutine xmpio_create_fsubarray_3D
!!***

#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_fsubarray_4D
!! NAME
!!  xmpio_create_fsubarray_4D
!!
!! FUNCTION
!!
!! INPUTS
!!  array_of_sizes(4)=number of elements of type old_type in each dimension of the full array (array of positive integers)
!!  array_of_subsizes(4)=number of elements of type old_type in each dimension of the subarray (array of positive integers)
!!  old_type=Old MPI type.
!!
!! OUTPUT
!!  mpi_err= MPI error code
!!  new_type=New MPI type.
!!
!! PARENTS
!!      m_io_screening
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE
#if defined HAVE_MPI_IO

subroutine xmpio_create_fsubarray_4D(array_of_sizes,array_of_subsizes,old_type,new_type,mpi_err)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_create_fsubarray_4D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: old_type
 integer,intent(out) :: mpi_err,new_type
!arrays
 integer,intent(in) :: array_of_sizes(4),array_of_subsizes(4)

!Local variables-------------------------------
!scalars
 integer :: bsize_frm
 integer :: bsize_old,nx,ny,nz,na,stride_x
 integer :: column_type,plane_type,ldx,ldy,ldz,lda,vol_type

!************************************************************************

 bsize_frm    = xmpio_bsize_frm    ! Byte size of the Fortran record marker.

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpi_err)
 !
 ! Number of columns and rows of the submatrix.
 nx = array_of_subsizes(1)
 ny = array_of_subsizes(2)
 nz = array_of_subsizes(3)
 na = array_of_subsizes(4)

 ldx = array_of_sizes(1)
 ldy = array_of_sizes(2)
 ldz = array_of_sizes(3)
 lda = array_of_sizes(4)

 ! Byte size of the Fortran record + the two markers.
 stride_x = ldx*bsize_old + 2*bsize_frm

 call MPI_Type_contiguous(nx,old_type,column_type,mpi_err)

 call MPI_Type_hvector(ny,1,stride_x,column_type,plane_type,mpi_err)

 call MPI_Type_hvector(nz,1,ldy*stride_x,plane_type,vol_type,mpi_err)

 call MPI_Type_hvector(na,1,ldz*ldy*stride_x,vol_type,new_type,mpi_err)

 call MPI_TYPE_COMMIT(new_type,mpi_err)

end subroutine xmpio_create_fsubarray_4D
!!***

#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_check_frmarkers
!! NAME
!!  xmpio_check_frmarkers
!!
!! FUNCTION
!!  Check a set of Fortran record markers starting at a given offset using MPI-IO.
!!
!! INPUTS
!!  mpi_fh=MPI-IO file handler.
!!  offset=MPI-IO file pointer
!!  at_option=Option for individual or collective reading.
!!  nfrec=Number of Fortran records to be checked.
!!  bsize_frecord(nfrec)=Byte size of the Fortran records (markers are NOT included)
!!    These values will be compared with the markers reported in the file.
!!
!! OUTPUT
!!  ierr=A non-zero error code signals failure.
!!
!! PARENTS
!!      defs_scalapack,exc_iterative_diago,m_bse_io
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

#if defined HAVE_MPI_IO

subroutine xmpio_check_frmarkers(mpi_fh,offset,at_option,nfrec,bsize_frecord,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_check_frmarkers'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfrec,at_option
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(inout) :: mpi_fh
 integer,intent(out) :: ierr
!arrays
 integer(XMPI_OFFSET_KIND),intent(in) :: bsize_frecord(nfrec)

!Local variables-------------------------------
!scalars
 integer :: nb,irec,frmarkers_type,jj,bsize_frm,mpi_type_frm,mpi_err
 integer(XMPI_OFFSET_KIND) :: displ
!arrays
 integer*2,allocatable :: bufdelim2(:)
 integer*4,allocatable :: bufdelim4(:)
 integer*8,allocatable :: bufdelim8(:)
#ifdef HAVE_FC_INT_QUAD
 integer*16,allocatable :: bufdelim16(:)
#endif
!integer :: statux(MPI_STATUS_SIZE)
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)
 integer(XMPI_OFFSET_KIND),allocatable :: delim_record(:)

!************************************************************************

 ierr=0

 bsize_frm    = xmpio_bsize_frm     ! Byte size of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm  ! MPI type of the record marker.
 !
 ! Define the view for the file.
 nb=2*nfrec
 ABI_ALLOCATE(block_length,(nb+2))
 ABI_ALLOCATE(block_displ,(nb+2))
 ABI_ALLOCATE(block_type,(nb+2))
 block_length(1)=1
 block_displ (1)=0
 block_type  (1)=MPI_LB

 jj=2; displ=0
 do irec=1,nfrec
   block_type (jj:jj+1) =mpi_type_frm
   block_length(jj:jj+1)=1
   block_displ(jj  )     = displ
   block_displ(jj+1)     = bsize_frm + displ + bsize_frecord(irec)
   jj=jj+2
   displ = displ + bsize_frecord(irec) + 2*bsize_frm ! Move to the beginning of the next column.
   if (xmpio_max_address(displ)) ierr=-1  ! Check for wraparound.
 end do

 block_length(nb+2)=1
 block_displ (nb+2)=displ
 block_type  (nb+2)=MPI_UB

 call xmpio_type_struct(nb+2,block_length,block_displ,block_type,frmarkers_type,mpi_err)
 ABI_DEALLOCATE(block_length)
 ABI_DEALLOCATE(block_displ)
 ABI_DEALLOCATE(block_type)

 call MPI_TYPE_COMMIT(frmarkers_type,mpi_err)

 call MPI_FILE_SET_VIEW(mpi_fh,offset,MPI_BYTE,frmarkers_type,"native",MPI_INFO_NULL,mpi_err)
 !
 jj=1
 ABI_ALLOCATE(delim_record,(nb))
 do irec=1,nfrec
   delim_record(jj:jj+1)=bsize_frecord(irec)
   jj=jj+2
 end do
 !
 ! Read markers according to the MPI type of the Fortran marker.
 !
 SELECT CASE (bsize_frm)

 CASE (4)
   ABI_ALLOCATE(bufdelim4,(nb))
   if (at_option==xmpio_at) then
     call MPI_FILE_READ    (mpi_fh,bufdelim4,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else if (at_option==xmpio_at_all) then
     call MPI_FILE_READ_ALL(mpi_fh,bufdelim4,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else
     ierr=2
   end if
   if (ANY(bufdelim4/=delim_record)) ierr=1
   !if (ierr==1) then
   !  do irec=1,2*nfrec
   !    write(std_out,*)"irec, bufdelim4, delim_record: ",irec,bufdelim4(irec),delim_record(irec)
   !  end do
   !end if
   ABI_DEALLOCATE(bufdelim4)

 CASE (8)
   ABI_ALLOCATE(bufdelim8,(nb))
   if (at_option==xmpio_at) then
     call MPI_FILE_READ    (mpi_fh,bufdelim8,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else if (at_option==xmpio_at_all) then
     call MPI_FILE_READ_ALL(mpi_fh,bufdelim8,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else
     ierr=2
   end if
   if (ANY(bufdelim8/=delim_record)) ierr=1
   ABI_DEALLOCATE(bufdelim8)

#if defined HAVE_FC_INT_QUAD
 CASE (16)
   ABI_ALLOCATE(bufdelim16,(nb))
   if (at_option==xmpio_at) then
     call MPI_FILE_READ    (mpi_fh,bufdelim16,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else if (at_option==xmpio_at_all) then
     call MPI_FILE_READ_ALL(mpi_fh,bufdelim16,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else
     ierr=2
   end if
   if (ANY(bufdelim16/=delim_record)) ierr=1
   ABI_DEALLOCATE(bufdelim16)
#endif

 CASE (2)
   ABI_ALLOCATE(bufdelim2,(nb))
   if (at_option==xmpio_at) then
     call MPI_FILE_READ    (mpi_fh,bufdelim2,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else if (at_option==xmpio_at_all) then
     call MPI_FILE_READ_ALL(mpi_fh,bufdelim2,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else
     ierr=2
   end if
   if (ANY(bufdelim2/=delim_record)) ierr=1
   ABI_DEALLOCATE(bufdelim2)

 CASE DEFAULT
   ierr=-2
 END SELECT
 !
 ! Free memory
 call MPI_TYPE_FREE(frmarkers_type,mpi_err)
 ABI_DEALLOCATE(delim_record)

end subroutine xmpio_check_frmarkers
!!***

#endif
!!***

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_max_address
!! NAME
!!  xmpio_max_address
!!
!! FUNCTION
!!  Returns .TRUE. if offset cannot be stored in a Fortran integer of kind XMPI_ADDRESS_KIND.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_MPI_IO

function xmpio_max_address(offset)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_max_address'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical :: xmpio_max_address
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
!arrays

!Local variables-------------------------------
!scalars
 integer(XMPI_ADDRESS_KIND) :: address
 integer(XMPI_OFFSET_KIND),parameter :: max_address=HUGE(address)-100

!************************************************************************

 xmpio_max_address = (offset >= max_address)

end function xmpio_max_address
!!***

#endif
!!***

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_write_frmarkers
!! NAME
!!  xmpio_write_frmarkers
!!
!! FUNCTION
!!  Write a set of Fortran record markers starting at a given offset using MPI-IO.
!!
!! INPUTS
!!  mpi_fh=MPI-IO file handler.
!!  offset=MPI-IO file pointer
!!  at_option=Option for individual or collective reading.
!!  nfrec=Number of Fortran records to be read.
!!  bsize_frecord(nfrec)=Byte size of the Fortran records to be written (markers are NOT included)
!!
!! OUTPUT
!!  ierr=A non-zero error code signals failure.
!!
!! PARENTS
!!      defs_scalapack,exc_build_block,exc_iterative_diago
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

#if defined HAVE_MPI_IO

subroutine xmpio_write_frmarkers(mpi_fh,offset,at_option,nfrec,bsize_frecord,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_write_frmarkers'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfrec,at_option
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(inout) :: mpi_fh
 integer,intent(out) :: ierr
!arrays
 integer(XMPI_OFFSET_KIND),intent(in) :: bsize_frecord(nfrec)

!Local variables-------------------------------
!scalars
 integer :: nb,irec,frmarkers_type,jj,bsize_frm,mpi_type_frm,mpi_err
 integer(XMPI_OFFSET_KIND) :: displ
!character(len=500) :: msg
!arrays
 integer*2,allocatable :: bufdelim2(:)
 integer*4,allocatable :: bufdelim4(:)
 integer*8,allocatable :: bufdelim8(:)
#if defined HAVE_FC_INT_QUAD
 integer*16,allocatable :: bufdelim16(:)
#endif
!integer :: statux(MPI_STATUS_SIZE)
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)
 integer(XMPI_OFFSET_KIND),allocatable :: delim_record(:)

!************************************************************************

 ierr=0

 bsize_frm    = xmpio_bsize_frm     ! Byte size of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm  ! MPI type of the record marker.
 !
 ! Define the view for the file
 nb=2*nfrec
 ABI_ALLOCATE(block_length,(nb+2))
 ABI_ALLOCATE(block_displ,(nb+2))
 ABI_ALLOCATE(block_type,(nb+2))
 block_length(1)=1
 block_displ (1)=0
 block_type  (1)=MPI_LB

 jj=2; displ=0
 do irec=1,nfrec
   block_type (jj:jj+1) =mpi_type_frm
   block_length(jj:jj+1)=1
   block_displ(jj  )     = displ
   block_displ(jj+1)     = bsize_frm + displ + bsize_frecord(irec)
   jj=jj+2
   displ = displ + bsize_frecord(irec) + 2*bsize_frm ! Move to the beginning of the next column.
   if (xmpio_max_address(displ)) ierr=-1  ! Check for wraparound.
 end do

 block_length(nb+2)=1
 block_displ (nb+2)=displ
 block_type  (nb+2)=MPI_UB

 call xmpio_type_struct(nb+2,block_length,block_displ,block_type,frmarkers_type,mpi_err)
 ABI_DEALLOCATE(block_length)
 ABI_DEALLOCATE(block_displ)
 ABI_DEALLOCATE(block_type)

 call MPI_TYPE_COMMIT(frmarkers_type,mpi_err)

 call MPI_FILE_SET_VIEW(mpi_fh,offset,MPI_BYTE,frmarkers_type,"native",MPI_INFO_NULL,mpi_err)

 jj=1
 ABI_ALLOCATE(delim_record,(nb))
 do irec=1,nfrec
   delim_record(jj:jj+1)=bsize_frecord(irec)
   jj=jj+2
 end do
 !
 ! Write all markers according to the MPI type of the Fortran marker.
 !
 SELECT CASE (bsize_frm)

 CASE (4)
   ABI_ALLOCATE(bufdelim4,(nb))
   bufdelim4=delim_record
   if (at_option==xmpio_at) then
     call MPI_FILE_WRITE    (mpi_fh,bufdelim4,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else if (at_option==xmpio_at_all) then
     call MPI_FILE_WRITE_ALL(mpi_fh,bufdelim4,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else
     ierr=2
   end if
   ABI_DEALLOCATE(bufdelim4)

 CASE (8)
   ABI_ALLOCATE(bufdelim8,(nb))
   bufdelim8=delim_record
   if (at_option==xmpio_at) then
     call MPI_FILE_WRITE    (mpi_fh,bufdelim8,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else if (at_option==xmpio_at_all) then
     call MPI_FILE_WRITE_ALL(mpi_fh,bufdelim8,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else
     ierr=2
   end if
   ABI_DEALLOCATE(bufdelim8)

#if defined HAVE_FC_INT_QUAD
 CASE (16)
   ABI_ALLOCATE(bufdelim16,(nb))
   bufdelim16=delim_record
   if (at_option==xmpio_at) then
     call MPI_FILE_WRITE    (mpi_fh,bufdelim16,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else if (at_option==xmpio_at_all) then
     call MPI_FILE_WRITE_ALL(mpi_fh,bufdelim16,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else
     ierr=2
   end if
   ABI_DEALLOCATE(bufdelim16)
#endif

 CASE (2)
   ABI_ALLOCATE(bufdelim2,(nb))
   bufdelim2=delim_record
   if (at_option==xmpio_at) then
     call MPI_FILE_WRITE    (mpi_fh,bufdelim2,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else if (at_option==xmpio_at_all) then
     call MPI_FILE_WRITE_ALL(mpi_fh,bufdelim2,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpi_err)
   else
     ierr=2
   end if
   ABI_DEALLOCATE(bufdelim2)

 CASE DEFAULT
   ierr=-2
 END SELECT
 !
 ! Free memory
 call MPI_TYPE_FREE(frmarkers_type,mpi_err)
 ABI_DEALLOCATE(delim_record)

end subroutine xmpio_write_frmarkers
!!***

#endif
!!***

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_fherm_packed
!! NAME
!!  xmpio_create_fherm_packed
!!
!! FUNCTION
!!  Returns an MPI datatype that can be used to (read|write) with MPI-IO the columns of an
!!  Hermitian matrix whose upper triangle is written on a Fortran binary file.
!!  Note that the view assumes that the file pointer used to create the MPI-IO view
!!  points to the first element of the first column. In other words,the first Fortran record marker
!!  (if any) is not taken into account in the calculation of the displacements.
!!
!! INPUTS
!!  array_of_starts(2)=starting coordinates in the global Hermitian matrix
!!     (array of positive integers with jj>=ii, Fortran convention)
!!  array_of_ends(2)=final coordinates in the global Hermitian matrix
!!     (array of positive integers, jj>=ii, Fortran convention)
!!  is_fortran_file=.FALSE. is C stream is used. .TRUE. for writing Fortran binary files.
!!  old_type=MPI datatype of the elements of the matrix.
!!
!! OUTPUT
!!  my_offset=Offset relative to the beginning of the matrix in the file.
!!  hmat_type=New MPI type.
!!  offset_err= error code
!!
!! NOTES
!!  The matrix on file is written in the following FORTRAN format (let us assume a 3x3 matrix for simplicity)
!!
!!    m (1,1)             m
!!    m (1,2) (2,2)       m
!!    m (1,3) (2,3) (3,3) m
!!
!!  each Fortran record stores a column of the packed Hermitian matrix, "m" denotes the Fortran
!!  record marker that introduces holes in the MPI-IO file view.
!!  To read the columns from (1,2) up to (2,2) one should use array_of_starts=(1,2) and array_of_ends=(2,2).
!!  The MPI-IO file view should be created by moving the file pointer so that it points to the elements (1,2).
!!
!! NOTES
!!  File views for C-streams is not optimal since one can use a single slice of contigous data.
!!
!! PARENTS
!!      exc_build_block
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

#if defined HAVE_MPI_IO

subroutine xmpio_create_fherm_packed(array_of_starts,array_of_ends,is_fortran_file,my_offset,old_type,hmat_type,offset_err)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_create_fherm_packed'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: old_type
 integer,intent(out) :: offset_err,hmat_type
 integer(XMPI_OFFSET_KIND),intent(out) :: my_offset
 logical,intent(in) :: is_fortran_file
!arrays
 integer,intent(in) :: array_of_starts(2),array_of_ends(2)

!Local variables-------------------------------
!scalars
 integer :: nrow,my_ncol,ii,bsize_old,col,jj_glob,bsize_frm,prev_col,mpi_err
 integer(XMPI_OFFSET_KIND) :: col_displ
 !character(len=500) :: msg
!arrays
 integer,allocatable :: col_type(:),block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)

!************************************************************************

 offset_err=0

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpi_err)

 bsize_frm=0; if (is_fortran_file) bsize_frm = xmpio_bsize_frm

 my_ncol = array_of_ends(2) - array_of_starts(2) + 1
 !
 ! Calculate my offset relative to the beginning of the matrix in the file.
 prev_col = array_of_starts(2)-1
 my_offset = (prev_col*(prev_col+1)/2)*bsize_old + (array_of_starts(1)-1)*bsize_old + 2*prev_col*bsize_frm + bsize_frm
 !
 ! col_type(col) describes the col-th column of the packed matrix.
 ! block_displ(col+1) stores its displacement taking into account the Fortran marker.
 ABI_ALLOCATE(col_type,(my_ncol))
 ABI_ALLOCATE(block_displ,(my_ncol+2))

 if (my_ncol>1) then
   col_displ=0
   do col=1,my_ncol
    jj_glob = (col-1) + array_of_starts(2)
    nrow = jj_glob
    if (jj_glob==array_of_starts(2)) nrow = jj_glob - array_of_starts(1) + 1 ! First column treated by me.
    if (jj_glob==array_of_ends(2))   nrow = array_of_ends(1)                 ! Last column treated by me.
    call MPI_Type_contiguous(nrow,old_type,col_type(col),mpi_err)
    !
    if (xmpio_max_address(col_displ)) offset_err=1  ! Test for wraparounds
    block_displ(col+1) = col_displ
    col_displ = col_displ + nrow * bsize_old + 2 * bsize_frm  ! Move to the next column.
   end do

 else if (my_ncol==1) then  ! The case of a single column is treated separately.
    block_displ(2) = 0
    nrow = array_of_ends(1) - array_of_starts(1) + 1
    call MPI_Type_contiguous(nrow,old_type,col_type(2),mpi_err)
    col_displ= nrow*bsize_old
    if (xmpio_max_address(col_displ)) offset_err=1  ! Test for wraparounds
 else
   write(std_out,*)" my_ncol cannot be negative!"
   call xmpi_abort()
 end if

 ABI_ALLOCATE(block_length,(my_ncol+2))
 ABI_ALLOCATE(block_type,(my_ncol+2))

 block_length(1)=1
 block_displ (1)=0
 block_type  (1)=MPI_LB

 do ii=2,my_ncol+1
   block_length(ii)=1
   block_type(ii)  =col_type(ii-1)
   !write(std_out,*)" ii-1, depl, length, type: ",ii-1,block_displ(ii),block_length(ii),block_type(ii)
 end do

 block_length(my_ncol+2)= 1
 block_displ (my_ncol+2)= col_displ
 block_type  (my_ncol+2)= MPI_UB

 call xmpio_type_struct(my_ncol+2,block_length,block_displ,block_type,hmat_type,mpi_err)

 call MPI_TYPE_COMMIT(hmat_type,mpi_err)

 ABI_DEALLOCATE(block_length)
 ABI_DEALLOCATE(block_displ)
 ABI_DEALLOCATE(block_type)

 do col=1,my_ncol
   call MPI_TYPE_FREE(col_type(col),mpi_err)
 end do
 ABI_DEALLOCATE(col_type)

end subroutine xmpio_create_fherm_packed
!!***

#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_coldistr_from_fpacked
!! NAME
!!  xmpio_create_coldistr_from_fpacked
!!
!! FUNCTION
!!  Returns an MPI datatype that can be used to MPI-IO (read|write) the columns of an
!!  (Hermitian|Symmetric) matrix whose upper triangle is written on a Fortran binary file.
!!  Note that the view assumes that the file pointer used to instanciate the MPI-IO view
!!  points to the first element of the first column. In other words,the first Fortran record marker
!!  (if any) is not taken into account in the calculation of the displacements.
!!
!! INPUTS
!!  array_of_sizes(2)=Number of elements of type old_type in each dimension of the full array (array of positive integers)
!!  my_cols(2)=initial and final column to (read|write). Array of positive integers, Fortran convention.
!!  old_type=MPI datatype of the elements of the matrix.
!!
!! OUTPUT
!!  new_type=New MPI type that can be used to instanciate the MPI-IO view for the Fortran file.
!!  my_offpad=Offset to be added to the file pointer giving the position of the first Fortran record
!!    marker (lets call it "base"). Each node should (read|write) using my_offset = base + my_offpad.
!!    my_offpad is used so that one can safely change the way the fileview is generated (for example
!!    to make it more efficient) without having to change the client code.
!!  offset_err=Error code. A non-zero returned value signals that the global matrix is tool large
!!    for a single MPI-IO access (see notes below).
!!
!! NOTES
!!  1) The matrix on file is written in the following FORTRAN format (let us assume a 3x3 matrix for simplicity)
!!
!!      m (1,1)             m
!!      m (1,2) (2,2)       m
!!      m (1,3) (2,3) (3,3) m
!!
!!     each Fortran record stores a column of the packed matrix, "m" denotes the Fortran
!!     record marker that introduces holes in the file view.
!!
!!  2) With (signed) Fortran integers, the maximum size of the file that
!!     that can be read in one-shot is around 2Gb when etype is set to byte.
!!     Using a larger etype might create portability problems (real data on machines using
!!     integer*16 for the marker) since etype must be a multiple of the Fortran record marker
!!     Due to the above reason, block_displ is given in bytes but it has to be defined as Fortran
!!     integer. If the displacement cannot be stored in a Fortran integer, the routine returns
!!     offset_err=1 so that the caller will know that several MPI-IO reads are nedded to
!!     read the file.
!!
!! PARENTS
!!      m_bse_io
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

#if defined HAVE_MPI_IO

subroutine xmpio_create_coldistr_from_fpacked(array_of_sizes,my_cols,old_type,new_type,my_offpad,offset_err)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_create_coldistr_from_fpacked'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: old_type
 integer,intent(out) :: new_type,offset_err
 integer(XMPI_OFFSET_KIND),intent(out) :: my_offpad
!arrays
 integer,intent(in) :: array_of_sizes(2),my_cols(2)

!Local variables-------------------------------
!scalars
 integer :: my_ncol,bsize_old,my_col
 integer :: my_nels,my_el,row_glob,ii_hpk,jj_hpk,col_glob,bsize_frm,mpi_err
 integer(XMPI_OFFSET_KIND) :: my_offset,ijp_glob
 !character(len=500) :: msg
!arrays
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)

!************************************************************************

 ! Byte size of the Fortran record marker.
 bsize_frm = xmpio_bsize_frm

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpi_err)

 ! my number of columns and total numer of elements to be read.
 my_ncol = my_cols(2) - my_cols(1) + 1
 my_nels = my_ncol*array_of_sizes(1)
 !
 ! block_displ(el+1) stores the displacement of the local element el taking into account the Fortran marker.
 ABI_ALLOCATE(block_displ,(my_nels+2))
 ABI_ALLOCATE(block_length,(my_nels+2))
 ABI_ALLOCATE(block_type,(my_nels+2))

 block_length(1)=1
 block_displ (1)=0
 block_type  (1)=MPI_LB
 !
 ! * the view assumes that the file pointer used to instanciate the MPI-IO view
 !   points to the first element of the first column. In other words,the first Fortran record marker
 !   is not taken into account in the calculation of the displacements.
 my_offpad=xmpio_bsize_frm

 ! * Some matrix elements are read twice. This part has to be tested.
 offset_err=0; my_el=0
 do my_col=1,my_ncol
   col_glob = (my_col-1) + my_cols(1)
   do row_glob=1,array_of_sizes(1)
     if (col_glob>=row_glob) then
       ii_hpk = row_glob
       jj_hpk = col_glob
       ijp_glob = row_glob + col_glob*(col_glob-1)/2  ! Index for packed form
     else ! Exchange the indeces as (jj,ii) will be read.
       ii_hpk = col_glob
       jj_hpk = row_glob
       ijp_glob = col_glob + row_glob*(row_glob-1)/2  ! Index for packed form
     end if
     my_el = my_el+1
     my_offset = (ijp_glob-1)* bsize_old + (jj_hpk-1)*2*bsize_frm
     if (xmpio_max_address(my_offset)) offset_err=1   ! Check for wraparounds.
     block_displ (my_el+1)=my_offset
     block_length(my_el+1)=1
     block_type  (my_el+1)=old_type
     !write(std_out,*)" my_el, displ: ",my_el,block_displ(my_el+1)
   end do
 end do

 block_length(my_nels+2)=1
 block_displ (my_nels+2)=my_offset
 block_type  (my_nels+2)=MPI_UB

 call xmpio_type_struct(my_nels+2,block_length,block_displ,block_type,new_type,mpi_err)

 call MPI_TYPE_COMMIT(new_type,mpi_err)

 ABI_DEALLOCATE(block_length)
 ABI_DEALLOCATE(block_displ)
 ABI_DEALLOCATE(block_type)

end subroutine xmpio_create_coldistr_from_fpacked
!!***

#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_coldistr_from_fp3blocks
!! NAME
!!  xmpio_create_coldistr_from_fp3blocks
!!
!! FUNCTION
!!  Returns an MPI datatype that can be used to MPI-IO (read|write) the columns of a
!!  matrix of the form  M = (S1    F3)
!!                          (F3^H  S2)
!!  where S1 and S2 are square (symmetric|Hermitian) matrices whose upper triangle is stored on file
!!  while F3 is a generic matrix (not necessarily square) stored in full mode.
!!  The Fortran file contains the blocks in the following order.
!!      upper(S1)
!!      upper(S2)
!!      F3
!! INPUTS
!!  array_of_sizes(2)=Number of elements of type old_type in each dimension of the full array M (array of positive integers)
!!  my_cols(2)=initial and final column to (read|write). Array of positive integers, Fortran convention.
!!  block_sizes(2,3)=The sizes of S1, S2, F.
!!  old_type=MPI datatype of the elements of the matrix.
!!
!! OUTPUT
!!  new_type=New MPI type that can be used to instanciate the MPI-IO view for the Fortran file.
!!  my_offpad=Offset to be added to the file pointer giving the position of the first Fortran record
!!    marker (lets call it "base"). Each node should (read|write) using my_offset = base + my_offpad.
!!    my_offpad is used so that one can safely change the way the fileview is generated (for example
!!    to make it more efficient) without having to change the client code.
!!  offset_err=Error code. A non-zero returned value signals that the global matrix is tool large
!!    for a single MPI-IO access (see notes below).
!!
!! NOTES
!!  1) block_displ is given in bytes due to the presence of the marker.
!!     If the displacement of an element is too large, the routine returns
!!     offset_err=1 so that the caller knows that several MPI-IO reads are required to (read| write) the file.
!!
!! PARENTS
!!      m_bse_io
!!
!! CHILDREN
!!      mpi_type_commit,mpi_type_size,xmpi_abort,xmpio_type_struct
!!
!! SOURCE

#if defined HAVE_MPI_IO

subroutine xmpio_create_coldistr_from_fp3blocks(array_of_sizes,block_sizes,my_cols,old_type,new_type,my_offpad,offset_err)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmpio_create_coldistr_from_fp3blocks'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: old_type
 integer,intent(out) :: new_type,offset_err
 integer(XMPI_OFFSET_KIND),intent(out) :: my_offpad
!arrays
 integer,intent(in) :: array_of_sizes(2),my_cols(2),block_sizes(2,3)

!Local variables-------------------------------
!scalars
 integer :: my_ncol,bsize_old,my_col,which_block,uplo,swap
 integer :: my_nels,my_el,row_glob,ii_hpk,jj_hpk,ii,jj
 integer :: col_glob,bsize_frm,mpi_err,row_shift,col_shift,n1,n2
 integer(XMPI_OFFSET_KIND) :: my_offset,ijp,bsize_tot,max_displ,min_displ
 integer(XMPI_ADDRESS_KIND) :: address
!arrays
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)
 integer(XMPI_OFFSET_KIND) :: bsize_mat(2)

!************************************************************************

 if ( array_of_sizes(1) /= SUM(block_sizes(1,1:2)) .or. &
&     array_of_sizes(2) /= SUM(block_sizes(2,1:2)) ) then
   write(std_out,*)" xmpio_create_coldistr_from_fp3blocks: Inconsistency between block_sizes ans array_of_sizes "
   call xmpi_abort()
 end if

 if ( block_sizes(1,1)/=block_sizes(2,1) .or.&
&     block_sizes(1,2)/=block_sizes(2,2) ) then
   write(std_out,*)" xmpio_create_coldistr_from_fp3blocks: first two blocks must be square"
   call xmpi_abort()
 end if

 if ( block_sizes(2,3)/=block_sizes(2,2) .or.&
&     block_sizes(1,3)/=block_sizes(1,1) ) then
   write(std_out,*)" xmpio_create_coldistr_from_fp3blocks: Full matrix must be square"
   call xmpi_abort()
 end if

 write(std_out,*)" xmpio_create_coldistr_from_fp3blocks is still under testing"
 !call xmpi_abort()

 ! Byte size of the Fortran record marker.
 bsize_frm = xmpio_bsize_frm

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpi_err)

 ! my number of columns and total numer of elements to be read.
 my_ncol = my_cols(2) - my_cols(1) + 1
 my_nels = array_of_sizes(1)*my_ncol
 !
 ! block_displ(el+1) stores the displacement of the local element el taking into account the Fortran marker.
 ABI_ALLOCATE(block_displ,(my_nels+2))
 ABI_ALLOCATE(block_length,(my_nels+2))
 ABI_ALLOCATE(block_type,(my_nels+2))
 !
 ! * the view assumes that the file pointer used to instanciate the MPI-IO view
 !   points to the first element of the first column. In other words,the first Fortran record marker
 !   is not taken into account in the calculation of the displacements.
 my_offpad=xmpio_bsize_frm
 !
 ! Byte size of the first two blocks including the markers.
 n1=block_sizes(1,1)
 bsize_mat(1) = (n1*(n1+1)/2)*bsize_old + 2*n1*bsize_frm

 n2=block_sizes(1,2)
 bsize_mat(2) = (n2*(n2+1)/2)*bsize_old + 2*n2*bsize_frm

 bsize_tot=SUM(bsize_mat) +  PRODUCT(block_sizes(:,3))*bsize_old + block_sizes(2,3)*2*bsize_frm - bsize_frm
 write(std_out,*)"bsize_mat",bsize_mat,"bsize_tot",bsize_tot
 !
 ! * Some matrix elements are read twice. This part has to be tested.
 offset_err=0; my_el=0; max_displ=0; min_displ=HUGE(address)
 do my_col=1,my_ncol
   col_glob = (my_col-1) + my_cols(1)
   do row_glob=1,array_of_sizes(1)
     !
     which_block=3
     if (row_glob<=block_sizes(1,1).and.col_glob<=block_sizes(2,1)) which_block=1
     if (row_glob >block_sizes(1,1).and.col_glob >block_sizes(2,1)) which_block=2

     if ( ANY(which_block == (/1,2/)) ) then ! S1 or S2
       !
       row_shift=(which_block-1)*block_sizes(1,1)
       col_shift=(which_block-1)*block_sizes(2,1)

       ii_hpk = row_glob - row_shift
       jj_hpk = col_glob - col_shift
       if (jj_hpk<ii_hpk) then ! Exchange the indeces so that the symmetric is read.
         swap   = jj_hpk
         jj_hpk = ii_hpk
         ii_hpk = swap
       end if
       ijp = ii_hpk + jj_hpk*(jj_hpk-1)/2  ! Index for packed form
       my_offset = (ijp-1)*bsize_old + (jj_hpk-1)*2*bsize_frm
       if (which_block==2) my_offset=my_offset+bsize_mat(1)    ! Shift the offset to account for S1.
       !my_offset=4
       !
     else
       ! The element belongs either to F3 of F3^H.
       ! Now find whether it is the upper or the lower block since only F3 is stored on file.
       uplo=1; if (row_glob>block_sizes(1,1)) uplo=2

       if (uplo==1) then
         row_shift=0
         col_shift=block_sizes(2,1)
       else
         row_shift=block_sizes(1,1)
         col_shift=0
       end if
       ii = row_glob - row_shift
       jj = col_glob - col_shift

       if (uplo==2) then ! Exchange the indeces since the symmetric element will be read.
         swap=jj
         jj  =ii
         ii  =swap
       end if

       my_offset = (ii-1)*bsize_old + (jj-1)*block_sizes(1,3)*bsize_old + (jj-1)*2*bsize_frm
       my_offset = my_offset + SUM(bsize_mat)
       !if (uplo==1) my_offset=my_offset + bsize_mat(1)
       !my_offset=0
       !if (ii==1.and.jj==1) write(std_out,*)" (1,1) offset = ",my_offset
       !if (ii==block_sizes(1,3).and.jj==block_sizes(2,3)) write(std_out,*)" (n,n) offset =", my_offset
       if (my_offset>=bsize_tot-1*bsize_old) then
         write(std_out,*)"WARNING (my_offset>bsize_tot-bsize_old),",ii,jj,my_offset,bsize_tot
       end if
     end if

     if (xmpio_max_address(my_offset)) offset_err=1   ! Check for wraparounds.
     my_el = my_el+1
     block_displ (my_el+1)=my_offset
     block_length(my_el+1)=1
     block_type  (my_el+1)=old_type
     max_displ = MAX(max_displ,my_offset)
     min_displ = MIN(min_displ,my_offset)
     !if (which_block==3) write(std_out,*)" my_el, which, displ: ",my_el,which_block,block_displ(my_el+1)
   end do
 end do

 write(std_out,*)" MAX displ = ",max_displ," my_nels = ",my_nels
 write(std_out,*)" MIN displ = ",MINVAL(block_displ(2:my_nels+1))

 !block_displ (1)=max_displ ! Do not change this value.
 !if (min_displ>0) block_displ (1)=min_displ ! Do not change this value.

 block_displ (1)=min_displ
 block_displ (1)=0
 block_length(1)=0
 block_type  (1)=MPI_LB

 block_length(my_nels+2)=0
 !block_displ (my_nels+2)=bsize_tot
 block_displ (my_nels+2)=max_displ
 block_type  (my_nels+2)=MPI_UB

 call xmpio_type_struct(my_nels+2,block_length,block_displ,block_type,new_type,mpi_err)
 !call xmpio_type_struct(my_nels,block_length(2:),block_displ(2:),block_type(2:),new_type,mpi_err)

 !call MPI_TYPE_CREATE_INDEXED_BLOCK(my_nels, block_length(2:), block_displ(2:), old_type, new_type, mpi_err)

 call MPI_TYPE_COMMIT(new_type,mpi_err)

 ABI_DEALLOCATE(block_length)
 ABI_DEALLOCATE(block_displ)
 ABI_DEALLOCATE(block_type)

end subroutine xmpio_create_coldistr_from_fp3blocks
!!***

#endif

!------------------------------------------------------------------------------------

END MODULE m_xmpi
!!***
