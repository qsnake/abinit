!{\src2tex{textfont=tt}}
!!****f* ABINIT/ioarr
!!
!! NAME
!! ioarr
!!
!! FUNCTION
!! Read or write rho(r) or v(r), either ground-state or response-functions.
!! If ground-state, these arrays are real, if response-functions, these arrays are complex.
!! (in general, an array stored in unformatted form on a real space fft grid).
!! rdwr=1 to read, 2 to write
!!
!! This subroutine should be called only by one processor
!! in the writing mode
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! (some may be output)
!! accessfil=1 if netcdf file is to be in/out-put
!! dtset <type(dataset_type)>=all input variables for this dataset
!! fform=integer specification for data type:
!!   2 for wf; 52 for density; 102 for potential
!!   old format (prior to ABINITv2.0): 1, 51 and 101.
!! fildata=file name
!! hdr <type(hdr_type)>=the header of wf, den and pot files
!!  if rdwr=1 , used to compare with the hdr of the read disk file
!!  if rdwr=2 , used as the header of the written disk file
!! mpi_enreg=informations about MPI parallelization
!! rdwr=choice parameter, see above
!! rdwrpaw=1 only if rhoij PAW quantities have to be read (if rdwr=1)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! arr(ncplxfft,dtset%nspden)=array on real space grid, returned for rdwr=1, input for rdwr=2
!! etotal=total energy (Ha), returned for rdwr=1
!! === if rdwrpaw/=0 ===
!!  pawrhoij(hdr%natom*hdr%usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!
!! PARENTS
!!      gstate,gw_tools,loop3dte,loper3,nonlinear,outscfcv,respfn,scfcv,scfcv3
!!      setup_positron,sigma,suscep
!!
!! CHILDREN
!!      etsf_io_low_close,etsf_io_low_error_to_str,etsf_io_low_open_modify
!!      etsf_io_low_open_read,etsf_io_main_get,etsf_io_main_put,handle_ncerr
!!      hdr_check,hdr_clean,hdr_io,hdr_io_etsf,hdr_io_netcdf,leave_new
!!      rhoij_copy,wffclose,wffopen,wrtout,xcomm_init,xcomm_self,xcomm_world
!!      xderiveread,xderiverrecend,xderiverrecinit,xderivewrecend
!!      xderivewrecinit,xderivewrite
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ioarr(accessfil,arr,dtset,etotal,fform,fildata,hdr,mpi_enreg, &
&                ncplxfft,pawrhoij,rdwr,rdwrpaw,wvl)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_wffile
 use m_errors
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 use m_header,  only : hdr_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ioarr'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accessfil,ncplxfft,rdwr,rdwrpaw
 integer,intent(inout) :: fform
 real(dp),intent(inout) :: etotal
 character(len=fnlen),intent(in) :: fildata
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 real(dp),intent(inout),target :: arr(ncplxfft,dtset%nspden)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*dtset%usepaw*rdwrpaw)

!Local variables-------------------------------
#if defined HAVE_TRIO_ETSF_IO
 character(len = fnlen) :: file_etsf
 type(etsf_main), target :: main_folder
 logical :: lstat
 type(etsf_io_low_error) :: error
 character(len = etsf_io_low_error_len) :: errmess
 integer :: ncid
#endif
!scalars
 integer :: accesswff,arr_id,fform_dum,i,i1,i2,i3,ia,iarr,ierr,ind,ispden,me,me_fft,ncerr
 integer :: ncid_hdr,ncplxfft_id,nspden_id,restart,restartpaw,spaceComm,spaceComm_io
 integer :: zindex,zstart,zstop
 character(len=500) :: message
 type(hdr_type) :: hdr0
 type(wffile_type) :: wff
!arrays
 real(dp),pointer :: my_density(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 restartpaw=0

!Check validity of arguments--only rho(r) (51,52) and V(r) (101,102) are presently supported

 if ( (fform-1)/2 /=25 .and. (fform-1)/2 /=50 ) then
   write(message,'(a,i10,a)')' Input fform=',fform,' not allowed.'
   MSG_BUG(message)
 end if

!Print input fform
 if ( (fform-1)/2==25 .and. rdwr==1) then
   write(message, '(a)' ) ' ioarr: reading density data '
 else if ( (fform-1)/2==25 .and. rdwr==2) then
   write(message, '(a)' ) ' ioarr: writing density data'
 else if ( (fform-1)/2==50 .and. rdwr==1) then
   write(message, '(a)' ) ' ioarr: reading potential data'
 else if ( (fform-1)/2==50 .and. rdwr==2) then
   write(message, '(a)' ) ' ioarr: writing potential data'
 end if
 call wrtout(std_out,message,'COLL')

 write(message, '(a,a)' ) ' ioarr: file name is ',trim(fildata)
 call wrtout(std_out,message,'COLL')

#if defined HAVE_TRIO_ETSF_IO
 if (accessfil == 3) then ! Initialize filename in case of ETSF file.
   file_etsf = trim(fildata) // '-etsf.nc'
   write(message, '(a,a)' ) ' ioarr: created file name for ETSF access ', trim(file_etsf)
   call wrtout(std_out, message, 'COLL')
 end if
#endif

!Some definitions for MPI-IO access
 if (accessfil == 4) then
   accesswff=IO_MODE_MPI
   if (rdwr==1) then
     call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_4d)
   else
     spaceComm=mpi_enreg%comm_fft
   end if
   me=xcomm_rank(spaceComm)
   if (mpi_enreg%nproc_fft>1) then
     me_fft=mpi_enreg%me_fft
     spaceComm_io=mpi_enreg%comm_fft
   else
     me_fft=0
     call xcomm_self(spaceComm_io)
   end if
 end if
 if (dtset%usewvl==1) then
   call xcomm_world(mpi_enreg,spaceComm,myrank=me)
 end if

!DEBUG
!write(std_out,*)' ioarr: here '
!ENDDEBUG

!=======================================
!Handle input from disk file
!=======================================

 if (rdwr==1) then
   if (accessfil == 0 .or. accessfil == 4) then
     if(accessfil == 4) then
       call WffOpen(accesswff,spaceComm,fildata,ierr,wff,0,me,tmp_unit,spaceComm_io)
       call hdr_io(fform_dum,hdr0,rdwr,wff)
!      Compare the internal header and the header from the file
       call hdr_check(fform,fform_dum,hdr,hdr0,'COLL',restart,restartpaw)

     else
       open (unit=tmp_unit,file=fildata,form='unformatted',status='old')
!      Initialize hdr0, thanks to reading of unwff1
       call hdr_io(fform_dum,hdr0,rdwr,tmp_unit)
!      Compare the internal header and the header from the file
       call hdr_check(fform,fform_dum,hdr,hdr0,'COLL',restart,restartpaw)
     end if
     etotal=hdr0%etot

!    NOTE : should check that restart is possible !!

!    Read data
     do ispden=1,dtset%nspden
       if(accessfil == 4) then
         call xderiveRRecInit(wff,ierr)
         call xderiveRead(wff,arr(1:ncplxfft,ispden),ncplxfft,spaceComm_io,ierr)
         call xderiveRRecEnd(wff,ierr)
       else
         read (tmp_unit) (arr(iarr,ispden),iarr=1,ncplxfft)
       end if
     end do

     if(accessfil == 4) then
       call wffclose(wff,ierr)
     else
       close (unit=tmp_unit)
     end if

#if defined HAVE_TRIO_NETCDF
!    CASE netcdf file input
   else if (accessfil == 1) then

     ncerr = nf90_open(path=fildata, mode=NF90_NOWRITE, ncid=ncid_hdr)
     call handle_ncerr(ncerr," open netcdf in ioarr ")

!    Initialize hdr0, thanks to reading of unwff1
     call hdr_io_netcdf(fform_dum,hdr0,rdwr,ncid_hdr)

!    DEBUG
!    write(std_out,*) 'hdr0%nband = ',hdr0%nband
!    ENDDEBUG

!    Compare the internal header and the header from the file
     call hdr_check(fform,fform_dum,hdr,hdr0,'COLL',restart,restartpaw)

     etotal=hdr0%etot

!    NOTE : should check that restart is possible !!

!    Read data: all dimensions at once. Can change through the
!    way things are written (lower down)
     ncerr = nf90_inq_varid(ncid=ncid_hdr,name="arr",varid=arr_id)
     call handle_ncerr(ncerr," inquire arr ")
     ncerr = nf90_get_var(ncid=ncid_hdr,varid=arr_id,values=arr)
     call handle_ncerr(ncerr," get arr ")

!    close netCDF file
     ncerr = nf90_close(ncid_hdr)
     call handle_ncerr(ncerr," close netcdf file")
#endif

#if defined HAVE_TRIO_ETSF_IO
   else if ( accessfil == 3 ) then

!    We open the file
     call etsf_io_low_open_read(ncid, trim(file_etsf), lstat, error_data = error)

     if(lstat)then

!      Read the header
       call hdr_io_etsf(fform_dum, hdr0, rdwr, ncid)

!      Compare the internal header and the header from the file
       call hdr_check(fform, fform_dum, hdr, hdr0, 'COLL', restart, restartpaw)

!      Read the array
       if (fform==52) then !for density
         main_folder%density%data2D => arr
       else if (fform==102) then ! for all potential forms!!!!
         main_folder%exchange_correlation_potential%data2D => arr
       end if
       call etsf_io_main_get(ncid, main_folder, lstat, error)

     end if

     if(lstat)then ! close the file
       call etsf_io_low_close(ncid, lstat, error_data = error)
     end if

     if (.not. lstat) then ! handle the error
       call etsf_io_low_error_to_str(errmess, error)
       MSG_ERROR(errmess)
     end if
#endif

   else
     write(message,'(a,i0,a)')'Bad value for accessfil', accessfil, ' on read '
     MSG_BUG(message)
   end if
!  end accessfil if

   write(message, '(a,a)' ) ' ioarr: data read from disk file ',trim(fildata)
   call wrtout(std_out,message,'COLL')

   etotal=hdr0%etot
!  Eventually copy (or distribute) PAW data
   if (rdwrpaw==1.and.restartpaw/=0) then
     call rhoij_copy(hdr0%pawrhoij,pawrhoij,mpi_enreg=mpi_enreg)
   end if

   if (accessfil == 0 .or. accessfil == 1 .or. accessfil == 3) then
     call hdr_clean(hdr0)
   end if


!  =======================================
!  Set up for writing data
!  =======================================
 else if (rdwr==2) then

!  In the wavelet case (isolated boundary counditions), the
!  arr array has a buffer that we need to remove.
   if (dtset%usewvl == 1) then
#if defined HAVE_DFT_BIGDFT
     zindex = mpi_enreg%nscatterarr(me, 3)
     zstart = max(15 - zindex, 0)
     zstop  = mpi_enreg%nscatterarr(me, 2) + &
&     mpi_enreg%nscatterarr(me, 4) - &
&     max(zindex + mpi_enreg%nscatterarr(me, 2) &
&     - 2 * wvl%Glr%d%n3 - 15, 0)
     if (zstop - zstart + 1 > 0) then
!      Our slab contains (zstop - zstart + 1) elements
       ABI_ALLOCATE(my_density,((wvl%Glr%d%n1*2)*(wvl%Glr%d%n2*2)*(zstop-zstart),dtset%nspden))
!      We copy the data except the buffer to my_density
       ind = 0

       do i3 = zstart, zstop - 1, 1
         ia = (i3 - 1) * wvl%Glr%d%n1i * wvl%Glr%d%n2i
         do i2 = 0, 2 * wvl%Glr%d%n2 - 1, 1
           i = ia + (i2 + 14) * wvl%Glr%d%n1i + 14
           do i1 = 0, 2 * wvl%Glr%d%n1 - 1, 1
             i   = i + 1
             ind = ind + 1
             my_density(ind, :) = arr(i, :)
           end do
         end do
       end do
     else
       nullify(my_density)
     end if
#else
     write(message, '(a,a,a,a)' ) ch10,&
&     ' ioarr : BigDFT library is not compiled.', ch10, &
&     '   Action, used the flag --enable-bigdft when configuring.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
#endif
   end if

   if (accessfil == 0 .or. accessfil == 4) then
     if(accessfil == 4) then
       call WffOpen(accesswff,spaceComm,fildata,ierr,wff,0,me,tmp_unit)
       call hdr_io(fform,hdr,rdwr,wff)
     else
       open(unit=tmp_unit,file=fildata,form='unformatted',status='unknown')
!      Write header
       call hdr_io(fform,hdr,rdwr,tmp_unit)
     end if

!    Write actual data
     do ispden=1,dtset%nspden
       if(accessfil == 4) then
         call xderiveWRecInit(wff,ierr,me_fft)
         call xderiveWrite(wff,arr(1:ncplxfft,ispden),ncplxfft,spaceComm_io,ierr)
         call xderiveWRecEnd(wff,ierr,me_fft)
       else
         if (dtset%usewvl == 0) then
           write(tmp_unit) (arr(iarr,ispden),iarr=1,ncplxfft)
         else
           write(tmp_unit) (my_density(iarr,ispden),iarr=1,size(my_density, 1))
         end if
       end if
     end do

     if(accessfil == 4) then
       call WffClose(wff,ierr)
     else
       close (tmp_unit)
     end if

#if defined HAVE_TRIO_NETCDF
   else if (accessfil == 1) then

!    Create empty netCDF file
     ncerr = nf90_create(path=fildata, cmode=NF90_CLOBBER, ncid=ncid_hdr)
     call handle_ncerr(ncerr," create netcdf file")
     ncerr = nf90_close(ncid_hdr)
     call handle_ncerr(ncerr," close netcdf file")

!    open the file. By default in data mode.
     ncerr = nf90_open(path=fildata, mode=NF90_WRITE, ncid=ncid_hdr)
     call handle_ncerr(ncerr," open netcdf in ioarr ")

     call hdr_io_netcdf(fform,hdr,rdwr,ncid_hdr)

!    get id for nspden dimension
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nspden",dimid=nspden_id)
     call handle_ncerr(ncerr," inquire nspden")
!    get id for other dimension
!    
!    WARNING !! This is a very ugly way to write netcdf, and is used as a test
!    and for compatibility with the standard dumb way of writing !!!!
!    Cleaner would be to pass the dimensions you want the arr to have, and
!    a descriptive string "density...", and to use dimids created in hdr_io,
!    but this would change the interface to ioarr a lot.
!    
!    go to definition mode
     ncerr = nf90_redef(ncid=ncid_hdr)
     call handle_ncerr(ncerr," redef mode ")
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="ncplxfft",len=ncplxfft,dimid=ncplxfft_id)
     call handle_ncerr(ncerr," define ncplxfft")


     ncerr = nf90_def_var(ncid=ncid_hdr,name="arr",xtype=NF90_DOUBLE,&
&     dimids=(/ncplxfft_id,nspden_id/),varid=arr_id)
     call handle_ncerr(ncerr," define arr")

!    end netCDF definition mode
     ncerr = nf90_enddef(ncid=ncid_hdr)
     call handle_ncerr(ncerr," enddef")

!    Write data: all dimensions at once.
!    dimensions should have been created in hdr_io_netcdf
     ncerr = nf90_put_var(ncid=ncid_hdr,varid=arr_id,values=arr)
     call handle_ncerr(ncerr," fill arr")


!    close netCDF file
     ncerr = nf90_close(ncid_hdr)
     call handle_ncerr(ncerr," close netcdf file")
#endif

#if defined HAVE_TRIO_ETSF_IO
   else if ( accessfil == 3 ) then
!    We open the file
     call etsf_io_low_open_modify(ncid, trim(file_etsf), lstat, error_data = error)
     if (lstat) then

!      Write the header
       call hdr_io_etsf(fform, hdr, rdwr, ncid)

!      Write the array
       if (fform==52) then !for density
         if (dtset%usewvl == 0) then
           main_folder%density%data2D => arr
         else
           main_folder%density%data2D => my_density
         end if
       else if (fform==102) then ! for all potential forms!!!!
         main_folder%exchange_correlation_potential%data2D => arr
       end if
       call etsf_io_main_put(ncid, main_folder, lstat, error)

     end if

     if (lstat) then ! Close the file
       call etsf_io_low_close(ncid, lstat, error_data = error)
     end if

     if (.not. lstat) then ! Handle the error
       call etsf_io_low_error_to_str(errmess, error)
       MSG_ERROR(errmess)
     end if
#endif

   else
     write(message,'(a,i0,a)')'Bad value for accessfil', accessfil, ' on write '
     MSG_ERROR(message)
   end if

   if (dtset%usewvl == 1) then
     if (associated(my_density))  then
       ABI_DEALLOCATE(my_density)
     end if
   end if

   write(message, '(a,a)' ) ' ioarr: data written to disk file ',trim(fildata)
   call wrtout(std_out,message,'COLL')

 else
   write(message,'(a,i12,a)')'Called with rdwr=',rdwr,' :  not allowed.'
   MSG_BUG(message)
 end if

 DBG_EXIT("COLL")

end subroutine ioarr
!!***
