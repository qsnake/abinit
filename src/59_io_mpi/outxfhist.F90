!{\src2tex{textfont=tt}}
!!****f* ABINIT/outxfhist
!! NAME
!! outxfhist
!!
!! FUNCTION
!!  read/write xfhist
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  option =
!!   1: write
!!   2: read only nxfh
!!   3: read xfhist
!!  response =
!!   0: GS wavefunctions
!!   1: RF wavefunctions
!!  natom = number of atoms in unit cell
!!  mxfh = last dimension of the xfhist array
!!
!! OUTPUT
!!  ios = error code returned by read operations
!!
!! SIDE EFFECTS
!!  nxfh = actual number of (x,f) history pairs, see xfhist array
!!  wff2 = structured info for wavefunctions
!!  xfhist(3,natom+4,2,ab_xfh%mxfh) = (x,f) history array, also including
!!   rprim and stress
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      handle_ncerr,leave_new,wrtout,xcomm_self,xderiveread,xderiverrecend
!!      xderiverrecinit,xderivewrecend,xderivewrecinit,xderivewrite
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outxfhist(ab_xfh,natom,option,wff2,ios)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_mover
 use m_xmpi
 use m_wffile
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outxfhist'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi, except_this_one => outxfhist
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer          ,intent(in)    :: natom,option
 integer          ,intent(out)   :: ios
 type(wffile_type),intent(inout)    :: wff2
 type(ab_xfh_type),intent(inout) :: ab_xfh

!Local variables-------------------------------
 integer :: ierr,ixfh,ncid_hdr,spaceComm,xfdim2
 real(dp),allocatable :: xfhist_tmp(:)
 character(len=500) :: message
!no_abirules
#if defined HAVE_TRIO_NETCDF
 integer :: ncerr
 integer :: nxfh_id, mxfh_id, xfdim2_id, dim2inout_id, dimr3_id,xfhist_id
 integer :: nxfh_tmp,mxfh_tmp,xfdim2_tmp,dim2inout_tmp
#endif


! *************************************************************************

!DEBUG
!write(std_out,*)'outxfhist  : enter, option = ', option
!ENDDEBUG
 ncid_hdr = wff2%unwff
 xfdim2 = natom+4

 ios = 0

!### (Option=1) Write out content of all iterations
!#####################################################################
 if ( option == 1 ) then

!  Write the (x,f) history
   if (wff2%accesswff == IO_MODE_FORTRAN) then
     write(unit=wff2%unwff)ab_xfh%nxfh
     do ixfh=1,ab_xfh%nxfh
       write(unit=wff2%unwff)ab_xfh%xfhist(:,:,:,ixfh)
     end do

   else if (wff2%accesswff == IO_MODE_FORTRAN_MASTER) then
!    FIXME: should copy the xfhist to other processors, and check that we are on the master to read in this case
!    if node is master
     write(message, "(A,A,A,A)") ch10, " outxfhist: ERROR -", ch10, &
&     'accesswff == -1 (localrdwf ) has not been coded yet for xfhist rereading.'
     call wrtout(std_out, message, 'COLL')
     call leave_new('COLL')
     write(unit=wff2%unwff)ab_xfh%nxfh
     do ixfh=1,ab_xfh%nxfh
       write(unit=wff2%unwff)ab_xfh%xfhist(:,:,:,ixfh)
     end do

!    insert mpi broadcast here

   else if(wff2%accesswff==IO_MODE_MPI)then
     ABI_ALLOCATE(xfhist_tmp,(3*(natom+4)*2))
     call xcomm_self(spaceComm)
     call xderiveWRecInit(wff2,ierr)
     call xderiveWrite(wff2,ab_xfh%nxfh,ierr)
     call xderiveWRecEnd(wff2,ierr)
     do ixfh=1,ab_xfh%nxfh
       xfhist_tmp(:)=reshape(ab_xfh%xfhist(:,:,:,ixfh),(/3*(natom+4)*2/))
       call xderiveWRecInit(wff2,ierr)
       call xderiveWrite(wff2,xfhist_tmp,3*(natom+4)*2,spaceComm,ierr)
       call xderiveWRecEnd(wff2,ierr)
     end do
     ABI_DEALLOCATE(xfhist_tmp)

#if defined HAVE_TRIO_NETCDF
   else if (wff2%accesswff == IO_MODE_NETCDF) then
!    check if nxfh and xfhist are defined
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)

     if (ncerr /= NF90_NOERR) then
!      need to define everything
       ncerr = nf90_redef (ncid=ncid_hdr)
       call handle_ncerr(ncerr," outxfhist : going to define mode ")

       ncerr = nf90_def_dim(ncid=ncid_hdr,name="dim2inout",len=2,dimid=dim2inout_id)
       call handle_ncerr(ncerr," outxfhist : define dim2inout")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="mxfh",len=ab_xfh%mxfh,dimid=mxfh_id)
       call handle_ncerr(ncerr," outxfhist : define mxfh")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="nxfh",len=ab_xfh%nxfh,dimid=nxfh_id)
       call handle_ncerr(ncerr," outxfhist : define nxfh")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="xfdim2",len=xfdim2,dimid=xfdim2_id)
       call handle_ncerr(ncerr," outxfhist : define xfdim2")

       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dimr3",dimid=dimr3_id)
       call handle_ncerr(ncerr," outxfhist : inquire dimr3")

!      ab_xfh%xfhist(3,natom+4,2,ab_xfh%mxfh)
       ncerr = nf90_def_var(ncid=ncid_hdr,name="xfhist",xtype=NF90_DOUBLE,&
&       dimids=(/dimr3_id,xfdim2_id,dim2inout_id,mxfh_id/),varid=xfhist_id)
       call handle_ncerr(ncerr," outxfhist : define xfhist")

!      End define mode and go to data mode
       ncerr = nf90_enddef(ncid=ncid_hdr)
       call handle_ncerr(ncerr," outxfhist : enddef call ")
     else
!      check that the dimensions are correct
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
       call handle_ncerr(ncerr," outxfhist : inquire nxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
&       len=nxfh_tmp)
       call handle_ncerr(ncerr,"  outxfhist : get nxfh")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="xfdim2",dimid=xfdim2_id)
       call handle_ncerr(ncerr," outxfhist : inquire xfdim2")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=xfdim2_id,&
&       len=xfdim2_tmp)
       call handle_ncerr(ncerr,"  outxfhist : get xfdim2")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="mxfh",dimid=mxfh_id)
       call handle_ncerr(ncerr," outxfhist : inquire mxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=mxfh_id,&
&       len=mxfh_tmp)
       call handle_ncerr(ncerr,"  outxfhist : get mxfh")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dim2inout",dimid=dim2inout_id)
       call handle_ncerr(ncerr," outxfhist : inquire dim2inout")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=dim2inout_id,&
&       len=dim2inout_tmp)
       call handle_ncerr(ncerr,"  outxfhist : get dim2inout")

       ncerr = nf90_inq_varid(ncid=ncid_hdr,name="xfhist",varid=xfhist_id)
       call handle_ncerr(ncerr," outxfhist : inquire xfhist")

       if (mxfh_tmp /= ab_xfh%mxfh .or. dim2inout_tmp /= 2 .or. xfdim2_tmp /= xfdim2) then
         write (message,"(A)") 'outxfhist : ERROR xfhist has bad dimensions in NetCDF file. Can not re-write it.'
         call wrtout(std_out, message, 'COLL')
         call leave_new('COLL')
       end if

     end if

!    Now fill the data
     ncerr = nf90_put_var(ncid=ncid_hdr,varid=xfhist_id,values=ab_xfh%xfhist)
     call handle_ncerr(ncerr," outxfhist : fill xfhist")

!    end NETCDF definition ifdef
#endif
   end if  ! end accesswff if

!  ### (Option=2) Read in number of iterations
!  #####################################################################
 else if ( option == 2 ) then

   if (wff2%accesswff == IO_MODE_FORTRAN) then
     read(unit=wff2%unwff,iostat=ios)ab_xfh%nxfh

   else if (wff2%accesswff == IO_MODE_FORTRAN_MASTER) then
!    FIXME: should copy the xfhist to other processors, and check that we are on the master to read in this case
!    if node is master
     write(message, "(A,A,A,A)") ch10, " outxfhist: ERROR -", ch10, &
&     'accesswff == -1 (localrdwf ) has not been coded yet for xfhist rereading.'
     call wrtout(std_out, message, 'COLL')
     call leave_new('COLL')
     read(unit=wff2%unwff,iostat=ios)ab_xfh%nxfh

   else if (wff2%accesswff == IO_MODE_MPI) then
     call xderiveRRecInit(wff2,ierr)
     call xderiveRead(wff2,ab_xfh%nxfh,ierr)
     call xderiveRRecEnd(wff2,ierr)

#if defined HAVE_TRIO_NETCDF
   else if (wff2%accesswff == IO_MODE_NETCDF) then
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
     call handle_ncerr(ncerr," outxfhist : inquire nxfh")
     ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
&     len=ab_xfh%nxfh)
     call handle_ncerr(ncerr,"  outxfhist : get nxfh")
#endif
   end if

!  ### (Option=3) Read in iteration content
!  #####################################################################
 else if ( option == 3 ) then
   if (wff2%accesswff == IO_MODE_FORTRAN) then
     do ixfh=1,ab_xfh%nxfhr
       read(unit=wff2%unwff,iostat=ios)ab_xfh%xfhist(:,:,:,ixfh)
     end do
   else if (wff2%accesswff == IO_MODE_FORTRAN_MASTER) then
!    FIXME: should copy the xfhist to other processors, and check that we are on the master to read in this case
!    if node is master
     write(message, "(A,A,A,A)") ch10, " outxfhist: ERROR -", ch10, &
&     'accesswff == -1 (localrdwf ) has not been coded yet for xfhist rereading.'
     call wrtout(std_out, message, 'COLL')
     call leave_new('COLL')
     do ixfh=1,ab_xfh%nxfhr
       read(unit=wff2%unwff,iostat=ios)ab_xfh%xfhist(:,:,:,ixfh)
     end do

   else if (wff2%accesswff == IO_MODE_MPI) then
     ABI_ALLOCATE(xfhist_tmp,(3*(natom+4)*2))
     call xcomm_self(spaceComm)
     do ixfh=1,ab_xfh%nxfhr
       call xderiveRRecInit(wff2,ierr)
       call xderiveRead(wff2,xfhist_tmp,3*(natom+4)*2,spaceComm,ierr)
       call xderiveRRecEnd(wff2,ierr)
       xfhist_tmp(:)=xfhist_tmp(:)
     end do
     ABI_DEALLOCATE(xfhist_tmp)
   end if

!  FIXME: should this be inside the if not mpi as above for options 1 and 2?
!  it is placed here because the netcdf read is a single operation
#if defined HAVE_TRIO_NETCDF
   if (wff2%accesswff == IO_MODE_NETCDF) then
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
     call handle_ncerr(ncerr," outxfhist : inquire nxfh")
     ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
&     len=ab_xfh%nxfhr)
     call handle_ncerr(ncerr,"  outxfhist : get nxfh")

     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=xfhist_id,name="xfhist")
     call handle_ncerr(ncerr," outxfhist : inquire xfhist")
     ncerr = nf90_get_var(ncid=ncid_hdr,varid=xfhist_id,values=ab_xfh%xfhist,&
&     start=(/1,1,1,1/),count=(/3,natom+4,2,ab_xfh%nxfhr/))
     call handle_ncerr(ncerr," outxfhist : read xfhist")
   end if
#endif

 else
!  write(std_out,*)' outxfhist : option ', option , ' not available '
   write(message, "(A,A,A,A,I3,A)") ch10, "outxfhist: ERROR -", ch10, &
&   "option ", option, " not available."
   call wrtout(std_out, message, "COLL")
   call leave_new("COLL")

 end if

!DEBUG
!write(std_out,*)' outxfhist : exit'
!ENDDEBUG

end subroutine outxfhist
!!***
