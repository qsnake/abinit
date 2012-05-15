!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffOpen
!! NAME
!! WffOpen
!!
!! FUNCTION
!! This subroutine opens a Wf file. It might be accessed
!! by different mechanisms (usual F90 IO routines,
!!  MPI I/O, or, in the future, NetCDF). The routine
!! provides a file handler, wff (a data structure containing
!! all needed information).
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MB,MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! accesswff=access mode (0 means all procs access using usual F90
!!  routines ; -1 means only the master proc access, using usual
!!  F90 routines ; 1 means MPI I/O; 2 means netcdf I/O)
!! filename=name of the file
!! master=the number of the master proc (only needed in parallel)
!! me=my number (only needed in parallel)
!! spaceComm= the space communicator handler (only needed in MPI parallel I/O)
!! spaceWorld= the space communicator for the whole set of procs
!! unwff=the file unit number
!!
!! OUTPUT
!! ier=error code
!! wff= structured info about the wavefunction file
!!
!! PARENTS
!!      conducti_nc,conducti_paw,conducti_paw_core,cut3d,emispec_paw,gstate
!!      inwffil,inwffil3,ioarr,kss2wfk,linear_optics_paw,loper3,m_wfs,nstdy3
!!      nstpaw3,optic,optics_paw,optics_paw_core,optics_vloc,outwf,uderiv
!!      wfk_read_ene
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,etsf_io_low_open_modify,handle_ncerr,leave_new
!!      mpi_comm_rank,mpi_comm_size,mpi_file_open,mpi_type_size,wrtout
!!      xcomm_self
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine WffOpen(accesswff,spaceComm,filename,ier,wff,master,me,unwff,&
&                  spaceComm_mpiio) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_wffile
 use m_errors      
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffOpen'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi, except_this_one => WffOpen
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer, intent(in)  :: accesswff,spaceComm,master,me,unwff
 integer, intent(in),optional  :: spaceComm_mpiio
 integer, intent(out) :: ier
 character(len=fnlen), intent(in) :: filename
 type(wffile_type), intent(out) :: wff

!Local variables-------------------------------
 character(len=500) :: message
 character(len=fnlen) :: fildata
!no_abirules
#if defined HAVE_MPI_IO
 integer :: isize
#endif
#if defined HAVE_TRIO_NETCDF
 integer :: ncerr, ndim,nvar,natt,uid
#endif
#if defined HAVE_TRIO_ETSF_IO
 type(etsf_io_low_error) :: error
 logical                 :: lstat
 character(len = etsf_io_low_error_len)   :: errmess
#endif

! *************************************************************************

!Initialize the mandatory data of the wff datastructure
 wff%unwff    =unwff
 wff%accesswff=accesswff
 wff%fname    =filename

!Initialize info useful for parallel use
 wff%nproc    =1
 wff%master   =master
 wff%me       =me
 wff%me_mpiio =0
 wff%spaceComm=spaceComm
 call xcomm_self(wff%spaceComm_mpiio)
 if (spaceComm/=abinit_comm_serial) then
#if defined HAVE_MPI
!  This case occurs when wff is connected to a DENSITY file
!  abinit_comm_output is generally equal to MPI_COMM_WORLD (except if paral. over images)
   if (spaceComm==MPI_COMM_SELF) wff%spaceComm=abinit_comm_output 
!  if (spaceComm==MPI_COMM_SELF) wff%spaceComm=MPI_COMM_WORLD
   call MPI_COMM_SIZE(wff%spaceComm,wff%nproc,ier)
!  Redefine the default MPIIO communicator if MPI, although MPIIO features should not be used unless 
!  present(spaceComm_mpiio).and.accesswff==1
   wff%spaceComm_mpiio=wff%spaceComm
   wff%me_mpiio=wff%me
#endif
   if (present(spaceComm_mpiio).and.accesswff==IO_MODE_MPI) wff%spaceComm_mpiio=spaceComm_mpiio
#if defined HAVE_MPI
   call MPI_COMM_RANK(wff%spaceComm_mpiio,wff%me_mpiio,ier)
#endif
 end if

 ier=0
 if (accesswff==IO_MODE_FORTRAN) then !  All processors see a local file
   open (unit=unwff,file=filename,form='unformatted')
   rewind(unwff)

 else if(accesswff==IO_MODE_FORTRAN_MASTER) then !  Only the master processor see a local file
   if(master==me)then
     open (unit=unwff,file=filename,form='unformatted')
     rewind(unwff)
   end if

#if defined HAVE_MPI_IO
 else if(accesswff==IO_MODE_MPI)then ! In the parallel case, only the master open filename file
   if(master==me)then
     open(unit=unwff,file=filename,form='unformatted')
     rewind(unwff)
   end if
   call MPI_FILE_OPEN(wff%spaceComm,filename,&
&   MPI_MODE_CREATE + MPI_MODE_RDWR,&
&   MPI_INFO_NULL,wff%fhwff,ier)

   ABI_CHECK_MPI(ier,"WffOpen")

!  Define all type values
   call MPI_Type_size(MPI_INTEGER,isize,ier)
   wff%nbOct_int=isize
   call MPI_Type_size(MPI_DOUBLE_PRECISION,isize,ier)
   wff%nbOct_dp=isize
   call MPI_Type_size(MPI_CHARACTER,isize,ier)
   wff%nbOct_ch=isize
   wff%nbOct_recMarker=-1;wff%kgwff=-1;wff%formwff=-1
   wff%offwff=0;wff%off_recs=0;wff%lght_recs=0
   wff%marker_mpi_type=MPI_INTEGER ! Default value
   if (MPI_OFFSET_KIND==4) then
     wff%offset_mpi_type=MPI_INTEGER4
   else  if (MPI_OFFSET_KIND==8) then
     wff%offset_mpi_type=MPI_INTEGER8
#if defined HAVE_FC_INT_QUAD
   else  if (MPI_OFFSET_KIND==16) then
     wff%offset_mpi_type=MPI_INTEGER16
#endif
   else  if (MPI_OFFSET_KIND==2) then
     wff%offset_mpi_type=MPI_INTEGER2
   end if
#endif
#if defined HAVE_TRIO_NETCDF
 else if (accesswff==IO_MODE_NETCDF)then
!  here we modify wff%unwff to save the netCDF identifier in it
   ncerr = nf90_open(path=filename, mode=NF90_WRITE, ncid=wff%unwff)
   call handle_ncerr(ncerr," open netcdf wavefunction file")
   write(std_out,*) ' WffOpen : open a netCDF file ', trim(filename), wff%unwff
!  DEBUG
   ncerr = nf90_Inquire(ncid=wff%unwff,nDimensions=ndim,nVariables=nvar,nAttributes=natt,unlimitedDimId=uid)
   call handle_ncerr(ncerr, " general Inquire ")
   write(std_out,*) 'WffOpen : found ndim,nvar,natt,uid = ', ndim,nvar,natt,uid
!  ENDDEBUG
#endif
#if defined HAVE_TRIO_ETSF_IO
 else if (accesswff==IO_MODE_ETSF)then
   write(fildata, "(A,A)") trim(filename), "-etsf.nc"
   call etsf_io_low_open_modify(wff%unwff, fildata, lstat, error_data = error)
   if (.not. lstat) then
     write(fildata, "(A,A,I0,A)") trim(filename), "_", me, "-etsf.nc"
     call etsf_io_low_open_modify(wff%unwff, fildata, lstat, error_data = error)
   end if
   if (.not. lstat) then
     call etsf_io_low_error_to_str(errmess, error)
     write(message, "(A,A,A,A)") ch10, " WffOpen: ERROR -", ch10, &
&     errmess(1:min(475, len(errmess)))
     call wrtout(std_out, message, 'COLL')
     call leave_new('COLL')
   end if
   write(message, '(3A,I0)' ) ' WffOpen: opening ', trim(wff%fname), &
&   "-etsf.nc on unit ", wff%unwff
   call wrtout(std_out, message, 'COLL')
#endif
 else

   write(message, '(12a,i6,3a)' ) ch10, &
&   ' WffOpen : ERROR -',ch10,&
&   '  For the time being the input variable accesswff is restricted ',ch10,&
&   '  to 0 (all cases), 1 (in case MPI is enabled),',ch10,&
&   '  2 (only sequential, and if the NetCDF library has been enabled),',ch10,&
&   '  or 3 (only sequential, and if the NetCDF and ETSF_IO libraries have been enabled).',ch10,&
&   '  Its value is accesswff=',accesswff,'.',ch10,&
&   '  Action : change accesswff or use ABINIT in parallel or enable NetCDF and/or ETSF_IO.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')

 end if

end subroutine WffOpen
!!***
