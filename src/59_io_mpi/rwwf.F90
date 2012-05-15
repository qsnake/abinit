!{\src2tex{textfont=tt}}
!!****f* ABINIT/rwwf
!! NAME
!! rwwf
!!
!! FUNCTION
!!  This subroutine reads (different options) or write (option=2) the block of records
!!  related to one k point, and one spin-polarization, that
!!  contains the wavefunctions (as well as the eigenvalues and occupations).
!!  If called with option -1, the records will be skipped.
!!  If called with option -2, only the wavefunctions are read.
!!  The disk file unitwf should have been prepared
!!  outside of this routine, in order to read or write the correct records.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA,XG,GMR,MVer,MB,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  formeig=format of the eigenvalues
!!     0 => vector of eigenvalues
!!     1 => hermitian matrix of eigenvalues
!!  headform=format of the header of the wf file, also governing the k block format
!!    in case headform=0, use the default (current) format and headform
!!  icg=shift to be given to the location of the cg array
!!  ikpt=index of current k point (only needed for error message)
!!  isppol=spin polarization currently treated (only needed for error message)
!!  mband=maximum number of bands (dimension of cg, eigen and occ)
!!  mcg=dimention of cg
!!  nband=number of bands actually in cg, eigen and occ
!!   (if writing mode : must be larger or equal to nband_disk, only nband_disk bands are written ;
!!    if reading mode : can be equal, larger or smaller than nband_disk, but
!!     cg, eigen and occ will not be completely filled if nband>nband_disk)
!!  nband_disk=number of bands on the disk file
!!  npw=number of plane waves
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  option= 2 for writing cg, eigen and occ,
!!          1 for reading cg and eigen,
!!         -1 for reading/skipping,
!!         -2 for reading cg only
!!          3 for reading the eigenvalues only
!!          4 for writing a file containing only eigenvalues and occupations
!!         -4 for reading a file written with 4
!!                 (different from 3 which reads a normal option 2 file)
!!  optkg= if 1 , read or write kg_k ; if 0, do not care about kg_k
!!  tim_rwwf=timing code of the calling routine (set to 0 if not attributed)
!!  wff=struct info for wavefunction
!!   | unitwf=unit number for wavefunction
!!
!! SIDE EFFECTS
!!  cg(2,npw*nspinor*mband)=planewave coefficients of wavefunctions,
!!    input if option=2; output if option=1 or -2
!!  eigen((2*mband)**formeig *mband)=array for holding eigenvalues (hartree)
!!    input if option=2 or 4; output if option=1
!!  kg_k(3,optkg*npw)=k+g data  (only if optkg==1)
!!    input if option=2; output if option=1 or -2
!!  nband_disk=number of bands on disk
!!    input if option=2 or 4; output in the other cases
!!  occ(mband)=array for holding eigenvalues (hartree)
!!    input if option=2 or 4; output if option=1
!!    no meaning if frmeig/=0
!!
!! NOTES
!!  WARNING : occ is not read in the present status of this routine
!!  WARNING : skipping k-blocks is also done in the randac subroutine
!!  WARNING : reading the two first records is also done in the rdnpw routine
!!  WARNING : writing the two first records is also done in the vtowfk3 routine
!!
!! TODO
!!  Some arguments are contained in the wff datastructure, and should be eliminated.
!!  option 3 should be called -3 (reading -> negative option) and others (-1,1) re-shuffled.
!!
!! PARENTS
!!      WffReadEigK,WffReadSkipK,berryphase,compare_interpol,ctocprj,dyfnl3
!!      eltfrkin3,eltfrnl3,energy,forstrnps,initwf,kss2wfk,ladielmt,lavnl,m_wfs
!!      mkrho,mkrho3,mrggkk,newkpt,optics_paw,optics_vloc,outkss,outwant,outwf
!!      overlap_wf,partial_dos_fractions,pawmkaewf,prctfvw1,prctfvw2,rhofermi3
!!      suscep_dyn,suscep_kxc_dyn,suscep_stat,tddft,uderiv,vtorho,vtorho3
!!      wffile,wfk_read_ene,wfread,wfsinp
!!
!! CHILDREN
!!      etsf_io_basisdata_put,etsf_io_electrons_put,etsf_io_main_put
!!      handle_ncerr,mpi_bcast,wffreadwrite_mpio,wffwritenpwrec,xderivewrecend
!!      xderivewrecinit,xderivewrite,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rwwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&
&               nband,nband_disk,npw,nspinor,occ,option,optkg,tim_rwwf,wff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_wffile
#if defined HAVE_MPI2 && !defined FC_G95
 use mpi
#endif
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rwwf'
 use interfaces_18_timing
 use interfaces_59_io_mpi, except_this_one => rwwf
!End of the abilint section

 implicit none
#if defined HAVE_MPI1 || (defined HAVE_MPI && defined FC_G95)
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: formeig,headform,icg,ikpt,isppol,mband,mcg,nband,npw
 integer,intent(inout) :: nband_disk
 integer,intent(in) :: nspinor,option,optkg,tim_rwwf
 integer,intent(inout),target :: kg_k(3,optkg*npw)
 real(dp),intent(inout),target :: cg(2,mcg),eigen((2*mband)**formeig*mband),occ(mband)
 type(wffile_type),intent(inout) :: wff
 type(MPI_type), intent(inout) :: mpi_enreg

!Local variables-------------------------------
 character(len=500) :: msg
 real(dp) :: tsec(2)

! *************************************************************************

 call timab(270+tim_rwwf,1,tsec)

!Might check that icg+npw*nband*nspinor is smaller than mcg

!Check that nband is smaller than mband, if one will not skip the records.
 if(nband>mband .and. option/=-1)then
   write(msg,'(a,i0,a,i0,a)')' One should have nband<=mband. However, nband=',nband,', and mband=',mband,'.'
   MSG_PERS_BUG(msg)
 end if

!Check that formeig is 0 or 1.
 if ( ALL(formeig/=(/0,1/)) ) then
   write(msg,'(a,i0,a)')' The argument formeig should be 0 or 1. However, formeig=',formeig,'.'
   MSG_PERS_BUG(msg)
 end if

!Check the value of option
 if ( ALL( option /= (/1,2,3,4,-1,-2,-4/) )) then
   write(msg,'(a,i0,a)')' The argument option should be 1, -1, 2, -2, 3, 4 or -4. However, option=',option,'.'
   MSG_PERS_BUG(msg)
 end if

 if (option/=2.and.option/=4) then ! read
   call readwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&
&   nband,nband_disk,npw,nspinor,occ,option,optkg,wff)
 else                                ! write
   call writewf(cg,eigen,formeig,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&
&   nband,nband_disk,npw,nspinor,occ,option,optkg,wff)
 end if

 call timab(270+tim_rwwf,2,tsec)

end subroutine rwwf
!!***

! -------------------------------------------------------------------------------------------------

!!****f* ABINIT/readwf
!! NAME
!! readwf
!!
!! FUNCTION
!!  This subroutine reads the block of records related to one k point, and one spin-polarization, that
!!  contains the wavefunctions (as well as the eigenvalues).
!!  The disk file unitwf should have been prepared outside of this routine, in order to read the correct records.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA,XG,GMR,MVer,MB,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  formeig=format of the eigenvalues
!!    0 => vector of eigenvalues
!!    1 => hermitian matrix of eigenvalues
!!  icg=shift to be given to the location of the cg array
!!  ikpt=index of current k point (only needed for error message)
!!  isppol=spin polarization currently treated (only needed for error message)
!!  kg_k(3,optkg*npw)=k+g data  (only if optkg==1)
!!  mband=maximum number of bands (dimension of cg, eigen and occ)
!!  mcg=dimension of cg
!!  nband=number of bands actually in cg, eigen and occ
!!  nband_disk=number of bands on the disk file
!!  npw=number of plane waves (on current proc)
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  option= 1 for reading cg, eigen and occ,
!!         -1 for reading/skipping,
!!         -2 for reading cg only
!!          3 for reading the eigenvalues only
!!         -4 for reading a file written with 4
!!  optkg= if 1 , read or write kg_k ; if 0, do not care about kg_k
!!  wff=struct info for wavefunction
!!
!! SIDE EFFECTS
!!  Current kpt and spin updated
!!  cg(2,npw*nspinor*mband)=planewave coefficients of wavefunctions,
!!  eigen((2*mband)**formeig *mband)=array for holding eigenvalues (hartree)
!!  occ(mband)=array for holding electronic occupations
!!
!! NOTES
!!  WARNING : occ is not read in the present status of this routine
!!  WARNING : skipping k-blocks is also done in the randac subroutine
!!  WARNING : reading the two first records is also done in the rdnpw routine
!!
!! PARENTS
!!      rwwf
!!
!! CHILDREN
!!      etsf_io_basisdata_put,etsf_io_electrons_put,etsf_io_main_put
!!      handle_ncerr,mpi_bcast,wffreadwrite_mpio,wffwritenpwrec,xderivewrecend
!!      xderivewrecinit,xderivewrite,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine readwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&
&                 nband,nband_disk,npw,nspinor,occ,option,optkg,wff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wffile
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
#define ABI_FUNC 'readwf'
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi, except_this_one => readwf
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: formeig,headform,icg,ikpt,isppol,mband,mcg,nband,npw
 integer,intent(in) :: nspinor,option,optkg
 integer,intent(inout) :: nband_disk
 integer,intent(inout),target :: kg_k(3,optkg*npw)
 real(dp),intent(inout),target :: cg(2,mcg),eigen((2*mband)**formeig*mband),occ(mband)
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
 integer :: iband,indxx,ios,ipw,ispinor_index,nband1,ncid_hdr
 integer :: npw1,npwso,npwso1,npwsotot,npwtot,nrec,nspinor1,nspinortot,unitwf,use_f90
 character(len=500) :: msg
 character(len=fnlen) :: fname
#if defined HAVE_MPI_IO
 integer :: ikpt_this_proc,ispinor
 integer,allocatable :: ind_cg_mpi_to_seq(:)
#endif
#if defined HAVE_TRIO_NETCDF
 integer :: cg_id,eigen_id,kg_id,ncerr,occ_id
#endif
#if defined HAVE_TRIO_ETSF_IO
 logical :: lstat
 type(etsf_basisdata),target :: wave_folder
 type(etsf_electrons),target :: electrons_folder
 type(etsf_io_low_error) :: error
 type(etsf_main) :: main_folder
#endif

! *********************************************************************

!Check the options
 if ( ALL(option/=(/1,3,-1,-2,-4/) )) then
   write(msg,'(a,i0,a)') '  The argument option should be -4, -2, -1, 1 or 3.  However, option=',option,'.'
   MSG_PERS_BUG(msg)
 end if

 npwtot=npw; npwso=npw*nspinor
 npwsotot=npwso
 nspinortot=min(2,(1+mpi_enreg%paral_spin)*nspinor)

 unitwf=wff%unwff;ncid_hdr=unitwf
 use_f90=0
 if (wff%accesswff==IO_MODE_FORTRAN.or.(wff%accesswff==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) use_f90=1
!$use_f90 = wff_usef90(wff)

#if defined HAVE_MPI_IO

 if (option==1.or.option==-2) then

   if (wff%accesswff==IO_MODE_MPI.and.nband>0) then
     call xsum_mpi(npwsotot,wff%spaceComm_mpiio,ios)
     npwtot=npwsotot/nspinortot
     ABI_ALLOCATE(ind_cg_mpi_to_seq,(npwso))
     if (mpi_enreg%flag_ind_kg_mpi_to_seq==1) then
       ikpt_this_proc=mpi_enreg%tab_kpt_distrib(ikpt)
       if ( ikpt_this_proc <= 0  ) then
         msg='rwwf Bug '
         MSG_PERS_BUG(msg)
       end if
       do ispinor=1,nspinor
         ispinor_index=ispinor
         if (mpi_enreg%nproc_spin>1) ispinor_index=mpi_enreg%me_spin + 1
         ind_cg_mpi_to_seq(1+npw*(ispinor-1):npw*ispinor)=npwtot*(ispinor_index-1) &
&         + mpi_enreg%bandfft_kpt(ikpt_this_proc)%ind_kg_mpi_to_seq(1:npw)
       end do
     else
       ind_cg_mpi_to_seq(1:npwso) = (/(ipw,ipw=1,npwso)/)
     end if
   end if
 end if
#endif

!---------------------------------------------------------------------------
!Read the first record: npw, nspinor, nband_disk
!---------------------------------------------------------------------------

 if (headform>=40.or.headform==0) then ! headform==0 refers to the current headform
   call WffReadNpwRec(ios,ikpt,isppol,nband_disk,npw1,nspinor1,wff)
   npwso1=npw1*nspinor1
   if(ios/=0)then
     inquire (unit=unitwf, NAME=fname)
     write(msg,'(3a,i4,2a,i4,5a)') &
&     '  Reading option of rwwf. Trying to read',ch10,&
&     '  the (npw,nspinor,nband) record of a wf file, unit=',unitwf,ch10,&
&     '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&     '  Action: check your input wf file:',ch10,&
&     trim(fname)
     MSG_PERS_ERROR(msg)
   end if
 else ! Old format
   if(use_f90==1)then
     read (unitwf,iostat=ios) npwso1,nband_disk
#if defined HAVE_MPI_IO
   else if(wff%accesswff==IO_MODE_MPI)then
     call xderiveRRecInit(wff,ios)
     call xderiveRead(wff,npwso1,ios)
     call xderiveRead(wff,nband_disk,ios)
     call xderiveRRecEnd(wff,ios)
#endif
   end if
   if(ios/=0)then
     inquire (unit=unitwf, NAME=fname)
     write(msg,'(3a,i4,2a,i4,5a)') &
&     ' Reading option of rwwf. Trying to read',ch10,&
&     ' the (npw,nband) record of a wf file with headform <40 , unit=',unitwf,ch10,&
&     ' gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&     ' Action: check your input wf file:',ch10,&
&     trim(fname)
     MSG_PERS_ERROR(msg)
   end if
 end if ! headform

 if (option==1.or.option==-2) then !  Will read the wavefunction and/or kg data, so check npw and nspinor
!  
   if (headform>=40.or.headform==0) then ! New format. headform==0 refers to the current headform
     if (npwtot/=npw1) then
       write(msg,'(3a,i0,a,i0,a)') &
&       ' Reading option of rwwf. One should have npw=npw1',ch10,&
&       ' However, npw=',npwtot,', and npw1=',npw1,'.'
       MSG_PERS_BUG(msg)
     end if
     if(nspinortot/=nspinor1)then
       write(msg,'(3a,i0,a,i0,a)') &
&       ' Reading option of rwwf. One should have nspinor=nspinor1',ch10,&
&       ' However, nspinor=',nspinortot,', and nspinor1=',nspinor1,'.'
       MSG_PERS_BUG(msg)
     end if
   else ! Treat the Old format.
     if(npwsotot/=npwso1)then
       write(msg,'(3a,i0,a,i0,a)') &
&       ' Reading option of rwwf. One should have npwso=npwso1',ch10,&
&       ' However, npwso=',npwsotot,', and npwso1=',npwso1,'.'
       MSG_PERS_BUG(msg)
     end if
   end if ! headform
!  
 end if ! option==1.or.option==2

!---------------------------------------------------------------------------
!Read the second record: (k+G) vectors
!---------------------------------------------------------------------------

 if (headform>=40.or.headform==0) then ! headform==0 refers to the current headform
   if ((option==1.or.option==-2.or.option==3).and.optkg/=0 )then
     if(use_f90==1)then
       read(unitwf,iostat=ios) kg_k(1:3,1:npw)
#if defined HAVE_MPI_IO
     else if(wff%accesswff==IO_MODE_MPI)then
       call xderiveRRecInit(wff,ios)
       if (mpi_enreg%flag_ind_kg_mpi_to_seq==1) then
         ikpt_this_proc=mpi_enreg%tab_kpt_distrib(ikpt)
         call xderiveRead(wff,kg_k(1:3,1:npw),3,npw,wff%spaceComm_mpiio,&
&         mpi_enreg%bandfft_kpt(ikpt_this_proc)%ind_kg_mpi_to_seq(1:npw),ios)
       else
!        MG The call below uses MPI_SCAN but here we want to read the full set of G.
!        call xderiveRead(wff,kg_k(1:3,1:npw),3,npw,wff%spaceComm_mpiio,ios)
         call xmpi_read_int2d(wff,kg_k(1:3,1:npw),wff%spaceComm_mpiio,xmpio_at_all,ios)
       end if
       call xderiveRRecEnd(wff,ios)
#endif
#if defined HAVE_TRIO_NETCDF
     else if (wff%accesswff==IO_MODE_NETCDF) then
!      In netcdf file saved as: kg(3,mpw,nkpt)
       ncerr = nf90_inq_varid(ncid=ncid_hdr,name="kg",varid=kg_id)
       call handle_ncerr(ncerr," inquire kg ")
!      only get part of kg for this kpoint ikpt and sppol isppol
       ncerr = nf90_get_var(ncid=ncid_hdr,varid=kg_id,values=kg_k,start=(/1,1,ikpt/),count=(/3,npw,1/))
       call handle_ncerr(ncerr," get kg_k ")
#endif
#if defined HAVE_TRIO_ETSF_IO
     else if (wff%accesswff == IO_MODE_ETSF) then ! Read reduced_coordinates_of_plane_waves for this k point.
       wave_folder%reduced_coordinates_of_plane_waves%data2D => kg_k
       wave_folder%red_coord_pw__kpoint_access = ikpt
       call etsf_io_basisdata_get(wff%unwff, wave_folder, lstat, error)
       ETSF_CHECK_MYERROR(lstat,error)
#endif
     end if
   else ! option
     call WffReadSkipRec(ios,1,wff) ! Skip the record
   end if
   if(ios/=0)then
     inquire (unit=unitwf, NAME=fname)
     write(msg,'(3a,i4,2a,i4,5a)')  &
&     '  Reading option of rwwf. Trying to read',ch10,&
&     '  the k+g record of a wf file, unit=',unitwf,ch10,&
&     '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&     '  Action: check your input wf file:',ch10,&
&     trim(fname)
     MSG_PERS_ERROR(msg)
   end if
 end if ! headform

!---------------------------------------------------------------------------
!Read the third record: eigenvalues
!---------------------------------------------------------------------------
!The reading of occ should be enabled, BUT taking into account
!of headform of the disk file : occ was NOT present in the disk files with headform=22

 nband1 = min(nband,nband_disk)

!===== Case formeig=0: read eigenvalues =====
 if (formeig==0) then

   if (option==1.or.option==3.or.option==-4) then
     if(use_f90==1)then
       read (unitwf,iostat=ios) eigen(1:nband1)
#if defined HAVE_MPI_IO
     else if(wff%accesswff==IO_MODE_MPI)then
       call xderiveRRecInit(wff,ios)
       call xderiveRead(wff,eigen,nband1,MPI_COMM_SELF,ios)
       call xderiveRRecEnd(wff,ios)
#endif
     end if
   else
     call WffReadSkipRec(ios,1,wff)
   end if ! option

#if defined HAVE_TRIO_NETCDF
   if(wff%accesswff==IO_MODE_NETCDF.and.(option== 1.or.option==3.or.option==-4)) then
!    In netcdf file saved as: eigen((2*mband)**formeig*mband,nkpt,nsppol)
     ncerr = nf90_inq_varid(ncid=ncid_hdr,name="eigen",varid=eigen_id)
     call handle_ncerr(ncerr," inquire eigen ")
!    Only get part of eigen for this kpoint ikpt and sppol isppol
     ncerr = nf90_get_var(ncid=ncid_hdr,varid=eigen_id,values=eigen,start=(/1,ikpt,isppol/),&
&     count=(/(2*nband1)**formeig*nband1,1,1/))
     call handle_ncerr(ncerr," get eigen ")
!    In netcdf file saved as: occ(mband,nkpt,nsppol)
     ncerr = nf90_inq_varid(ncid=ncid_hdr,name="occ",varid=occ_id)
     call handle_ncerr(ncerr," inquire occ ")
!    Only get part of occ for this kpoint ikpt and sppol isppol
     ncerr = nf90_get_var(ncid=ncid_hdr,varid=occ_id,values=occ,start=(/1,ikpt,isppol/),&
&     count=(/nband1,1,1/))
     call handle_ncerr(ncerr," get occ ")
#if defined HAVE_TRIO_ETSF_IO
   else if (wff%accesswff == IO_MODE_ETSF) then
!    We get eigenvalues and occupations
     electrons_folder%eigenvalues%data1D => eigen
     electrons_folder%eigenvalues__kpoint_access = ikpt
     electrons_folder%eigenvalues__spin_access = isppol
     electrons_folder%occupations%data1D => occ
     electrons_folder%occupations__kpoint_access = ikpt
     electrons_folder%occupations__spin_access = isppol
     call etsf_io_electrons_get(wff%unwff, electrons_folder, lstat, error)
     ETSF_CHECK_MYERROR(lstat,error)
#endif
   end if
   if(ios/=0)then
     inquire (unit=unitwf, NAME=fname)
     write(msg,'(3a,i4,2a,i4,5a)') &
&     '  Reading option of rwwf. Trying to read',ch10,&
&     '  an eigenvalue record of a wf file, unit=',unitwf,ch10,&
&     '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&     '  Action: check your input wf file.',ch10,&
&     trim(fname)
     MSG_PERS_ERROR(msg)
   end if
#endif

!  ===== Case formeig=1: read matrix of eigenvalues =====
!  Will be written later (together with wave-functions)
 else if(formeig==1)then
 end if ! formeig

!---------------------------------------------------------------------------
!Read the wave-function coefficients
!---------------------------------------------------------------------------

!Select bands
 nband1=min(nband,nband_disk)
 if(nband1>0.and.option/=-1)then

!  ===== Case formeig=0: read only wave-functions =====
   if (formeig==0) then

     if (option==1.or.option==-2) then
       if (use_f90==1) then
         do iband=1,nband1
           ipw=(iband-1)*npwso+icg
           read(unitwf,iostat=ios) cg(1:2,ipw+1:ipw+npwso)
           if (ios/=0) exit
         end do
#if defined HAVE_MPI_IO
       else if (wff%accesswff==IO_MODE_MPI) then
         call WffReadWrite_mpio(wff,1,cg,mcg,icg,nband1,npwso,npwsotot,ind_cg_mpi_to_seq,ios)
#endif
       end if
     else if (option/=-4) then
       do iband=1,nband1
         call WffReadSkipRec(ios,1,wff) ! Skip the record
       end do
     end if ! option
     if(ios/=0)then
       inquire (unit=unitwf, NAME=fname)
       write(msg,'(3a,i4,2a,i4,5a)') &
&       '  Reading option of rwwf. Trying to read',ch10,&
&       '  a RF wf record of a wf file, unit=',unitwf,ch10,&
&       '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&       '  Action: check your input wf file.',ch10,&
&       trim(fname)
       MSG_PERS_ERROR(msg)
     end if
#if defined HAVE_TRIO_NETCDF
!    Here read all wf_k as a block
!    For the moment use icg offset and save wfk as one big array, like in abinit.
!    Later could separate dimensions and make a sliced extraction of the wfk
     if(wff%accesswff==IO_MODE_NETCDF.and.(option==1.or.option==-2))then
!      In netcdf file saved as: cg(2,mpw,nspinor,mband,nkpt,nsppol)
       ncerr = nf90_inq_varid(ncid=ncid_hdr,name="cg",varid=cg_id)
       call handle_ncerr(ncerr," inquire cg ")
!      Only get part of cg for this kpoint ikpt and sppol isppol
       ncerr = nf90_get_var(ncid=ncid_hdr,varid=cg_id,values=cg(:,icg+1:icg+npw1*nspinor1*nband1),&
&       start=(/1,1,1,1,ikpt,isppol/),count=(/2,npw1,nspinor1,nband1,1,1/))
       call handle_ncerr(ncerr," get cg ")
#if defined HAVE_TRIO_ETSF_IO
     else if (wff%accesswff == IO_MODE_ETSF) then
!      We get the coefficients_of_wavefunctions
!      !!     main_folder%coefficients_of_wavefunctions%data2D => &
!      !!          & cg(:, icg + 1:icg + npw * nspinor * nband)
!      With g95, the association done above sometime leads to segfaults.
!      So we allocate a temporary array to store the wfs of our kpt.
       ABI_ALLOCATE(main_folder%coefficients_of_wavefunctions%data2D,(2,npw*nspinor*nband))
       main_folder%wfs_coeff__kpoint_access = ikpt
       main_folder%wfs_coeff__spin_access = isppol
       main_folder%wfs_coeff__number_of_states = nband
       main_folder%wfs_coeff__number_of_coefficients = npw * nspinor
       call etsf_io_main_get(wff%unwff, main_folder, lstat, error)
       ETSF_CHECK_MYERROR(lstat,error)
!      Now we copy our values and deallocate the temporary array.
!      cg(:,icg+1:icg+npw*nspinor*nband)=main_folder%coefficients_of_wavefunctions%data2D
!      this is better than the previous instruction
!      to optimize virtual memory.
       do iband=1,nband
         ipw=(iband-1)*npwso
         cg(:,icg+ipw+1:icg+ipw+npwso)=main_folder%coefficients_of_wavefunctions%data2D(:,ipw+1:ipw+npwso)
       end do
!      
       ABI_DEALLOCATE(main_folder%coefficients_of_wavefunctions%data2D)
#endif
     end if
#endif

!    ===== Case formeig=1: read eigenvalues, occupations and wave-functions =====
   else if(formeig==1)then
!    Not available for NETCDF and ETSF_IO
     indxx=0
     do iband=1,nband1
       if (option==1.or.option==3.or.option==-4) then
         if(use_f90==1)then
           read (unitwf,iostat=ios) eigen(1+indxx:2*nband1+indxx)
#if defined HAVE_MPI_IO
         else if(wff%accesswff==IO_MODE_MPI)then
!          Should use an access with a "view"
           call xderiveRRecInit(wff,ios)
           call xderiveRead(wff,eigen(1+indxx:2*nband1+indxx),2*nband1,MPI_COMM_SELF,ios)
           call xderiveRRecEnd(wff,ios)
#endif
         end if
         indxx=indxx+2*nband1
       else
         call WffReadSkipRec(ios,1,wff) ! Skip the record
       end if
       if(ios/=0)then
         inquire (unit=unitwf, NAME=fname)
         write(msg,'(3a,i4,2a,i4,5a)') &
&         '  Reading option of rwwf. Trying to read',ch10,&
&         '  a RF eigenvalue record of a wf file, unit=',unitwf,ch10,&
&         '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&         '  Action: check your input wf file.',ch10,&
&         trim(fname)
         MSG_PERS_ERROR(msg)
       end if
       if(option==1.or.option==-2)then
         ipw=(iband-1)*npwso+icg
         if(use_f90==1)then
           ipw=(iband-1)*npwso+icg
           read(unitwf,iostat=ios) cg(1:2,ipw+1:ipw+npwso)
#if defined HAVE_MPI_IO
         else if(wff%accesswff==IO_MODE_MPI)then
!          Should use an access with a "view"
           call xderiveRRecInit(wff,ios)
           call xderiveRead(wff,cg(1:2,ipw+1:ipw+npwso),2,npwso,wff%spaceComm_mpiio,ind_cg_mpi_to_seq,ios)
           call xderiveRRecEnd(wff,ios)
#endif
         end if
       else if (option/=-4) then
         call WffReadSkipRec(ios,1,wff) ! Skip the record
       end if ! option
       if(ios/=0)then
         inquire (unit=unitwf, NAME=fname)
         write(msg,'(3a,i4,2a,i4,5a)') &
&         '  Reading option of rwwf. Trying to read',ch10,&
&         '  a RF wf record of a wf file, unit=',unitwf,ch10,&
&         '  gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&         '  Action: check your input wf file.',ch10,&
&         trim(fname)
         MSG_PERS_ERROR(msg)
       end if
     end do ! iband
   end if ! formeig

 end if ! nband >0

!If fewer than all bands were read wind disk file forward to end of bands for this k point.
!Will have to fill the non-filled bands outside of this routine ...
 if (nband<nband_disk .or. option==-1) then
   nrec=(formeig+1)*(nband_disk-nband)
   if(option==-1)nrec=(formeig+1)*nband_disk
   call WffReadSkipRec(ios,nrec,wff)
 end if

!---------------------------------------------------------------------------
!Final statements
!---------------------------------------------------------------------------

#if defined HAVE_MPI_IO
 if(wff%accesswff==IO_MODE_MPI.and.nband>0.and.(option==1.or.option==-2))   then
   ABI_DEALLOCATE(ind_cg_mpi_to_seq)
 end if
#endif

 RETURN
 ABI_UNUSED(mpi_enreg%flag_ind_kg_mpi_to_seq) ! just to keep mpi_enreg as an argument when HAVE_MPI_IO is not true

end subroutine readwf
!!***

! -------------------------------------------------------------------------------------------------

!!****f* ABINIT/writewf
!! NAME
!! writewf
!!
!! FUNCTION
!!  This subroutine writes the block of records
!!  related to one k point, and one spin-polarization, that
!!  contains the wavefunctions (as well as the eigenvalues and occupations).
!!  The disk file unitwf should have been prepared
!!  outside of this routine, in order to write the correct records.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA,XG,GMR,MVer,MB,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,npw*nspinor*mband)=planewave coefficients of wavefunctions,
!!  eigen((2*mband)**formeig *mband)=array for holding eigenvalues (hartree)
!!  formeig=format of the eigenvalues
!!   0 => vector of eigenvalues
!!   1 => hermitian matrix of eigenvalues
!!  icg=shift to be given to the location of the cg array
!!  ikpt=index of current k point (only needed for error message)
!!  isppol=spin polarization currently treated (only needed for error message)
!!  kg_k(3,optkg*npw)=k+g data  (only if optkg==1)
!!  mband=maximum number of bands (dimension of cg, eigen and occ)
!!  mcg=dimension of cg
!!  nband=number of bands actually in cg, eigen and occ
!!   (must be larger or equal to nband_disk, only nband_disk bands are written)
!!  nband_disk=number of bands on the disk file
!!  npw=number of plane waves
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  occ(mband)=array for holding electronic occupations
!!  option= 2 for writing cg, eigen and occ,
!!          4 for writing a file containing only eigenvalues and occupations
!!  optkg= if 1 , read or write kg_k ; if 0, do not care about kg_k
!!  wff=struct info for wavefunction
!!
!! OUTPUT
!! (none, only writing)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  WARNING : skipping k-blocks is also done in the randac subroutine
!!  WARNING : writing the two first records is also done in the vtowfk3 routine
!!
!! PARENTS
!!      m_io_kss,rwwf
!!
!! CHILDREN
!!      etsf_io_basisdata_put,etsf_io_electrons_put,etsf_io_main_put
!!      handle_ncerr,mpi_bcast,wffreadwrite_mpio,wffwritenpwrec,xderivewrecend
!!      xderivewrecinit,xderivewrite,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine writewf(cg,eigen,formeig,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&
&                  nband,nband_disk,npw,nspinor,occ,option,optkg,wff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wffile
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
#define ABI_FUNC 'writewf'
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi, except_this_one => writewf
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: formeig,icg,ikpt,isppol,mband,mcg,nband,nband_disk,npw,nspinor,option,optkg
 integer,intent(in),target :: kg_k(3,optkg*npw)
 real(dp),intent(in),target :: cg(2,mcg),eigen((2*mband)**formeig*mband),occ(mband)
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
 integer :: iband,ii,ios,ipw,ispinor_index,jj,nband2,ncid_hdr,npwso,npwsotot,npwtot,nspinortot
 integer :: unitwf,use_f90
 character(len=500) :: msg
#if defined HAVE_MPI_IO
 integer :: ikpt_this_proc,ispinor
 integer,allocatable :: ind_cg_mpi_to_seq(:)
 real(dp),pointer :: cg_ptr(:,:)
#endif
#if defined HAVE_TRIO_NETCDF
 integer :: cg_id,eigen_id,kg_id,ncerr,occ_id
#endif
#if defined HAVE_TRIO_ETSF_IO
 logical :: lstat
 type(etsf_basisdata) :: wave_folder
 type(etsf_electrons) :: electrons_folder
 type(etsf_io_low_error) :: error
 type(etsf_main) :: main_folder
#endif

! *********************************************************************

!Check the options
 if ( ALL(option /= (/2,4/)) ) then
   write(msg,'(a,i0)')' The argument option should be 2 or 4. However, option=',option
   MSG_PERS_BUG(msg)
 end if

!Check that nband_disk is not larger than nband (only for writing)
 if (nband<nband_disk) then
   write(msg,'(3a,i5,a,i5,a)') &
&   ' Writing option of rwwf. One should have nband<=nband_disk',ch10,&
&   ' However, nband=',nband,', and nband_disk=',nband_disk,'.'
   MSG_PERS_BUG(msg)
 end if

 npwtot=npw; npwso=npw*nspinor
 unitwf=wff%unwff; ncid_hdr=unitwf
 npwsotot=npwso
 nspinortot=min(2,(1+mpi_enreg%paral_spin)*nspinor)

 use_f90=0
 if (wff%accesswff==IO_MODE_FORTRAN.or.(wff%accesswff ==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) use_f90=1
!use_f90 = wff_usef90(wff)

#if defined HAVE_MPI_IO
 if (wff%accesswff==IO_MODE_MPI) then

   call xsum_mpi(npwsotot,wff%spaceComm_mpiio,ios)
   npwtot=npwsotot/nspinortot

   if (option/=4) then
     ABI_ALLOCATE(ind_cg_mpi_to_seq,(npwso))
     if (mpi_enreg%flag_ind_kg_mpi_to_seq==1) then
       ikpt_this_proc=mpi_enreg%tab_kpt_distrib(ikpt)
       do ispinor=1,nspinor
         ispinor_index=ispinor
         if (mpi_enreg%nproc_spin > 1) ispinor_index = mpi_enreg%me_spin + 1
         ind_cg_mpi_to_seq(1+npw*(ispinor-1):npw*ispinor)=npwtot*(ispinor_index-1) &
&         + mpi_enreg%bandfft_kpt(ikpt_this_proc)%ind_kg_mpi_to_seq(1:npw)
       end do
     else
       ind_cg_mpi_to_seq(1:npwso) = (/(ipw,ipw=1,npwso)/)
     end if
   end if
 end if
#endif

!---------------------------------------------------------------------------
!Write the first record: npw, nspinor, nband_disk
!---------------------------------------------------------------------------
!Not modified for netCDF: no need to add writing of nband_disk,npw,nspinor

 call WffWriteNpwRec(ios,nband_disk,npwtot,nspinortot,wff,opt_paral=2)

!---------------------------------------------------------------------------
!Write the second record: (k+G) vectors
!---------------------------------------------------------------------------

 if (optkg/=0.and.option/=4) then
   if(use_f90==1)then
     write(unitwf) kg_k(1:3,1:optkg*npw)
#if defined HAVE_MPI_IO
   else if (wff%accesswff==IO_MODE_MPI) then

     if (mpi_enreg%flag_ind_kg_mpi_to_seq==1) then
       ikpt_this_proc=mpi_enreg%tab_kpt_distrib(ikpt)
       call xderiveWRecInit(wff,ios,mpi_enreg%me_cart_3d)
       if (mpi_enreg%me_spin==0) then
         call xderiveWrite(wff,kg_k,3,npw,mpi_enreg%commcart, &
&         mpi_enreg%bandfft_kpt(ikpt_this_proc)%ind_kg_mpi_to_seq(1:npw),ios)
       end if
       call xderiveWRecEnd(wff,ios,mpi_enreg%me_cart_3d)
     else
!      MG does it work if we are not using FFT distribution ?
       call xderiveWRecInit(wff,ios )
       if (mpi_enreg%me_spin==0) then
         call xderiveWrite(wff,kg_k,3,optkg*npw,Wff%spaceComm_mpiio,ios)
       end if
       call xderiveWRecEnd(wff,ios)
     end if
#endif
#if defined HAVE_TRIO_NETCDF
   else if (wff%accesswff==IO_MODE_NETCDF)then
!    In netcdf file saved as: kg(3,mpw,nkpt)
!    Write data: all dimensions at once.
!    Dimensions should have been created in hdr_io_netcdf
     ncerr = nf90_inq_varid(ncid=ncid_hdr,name="kg",varid=kg_id)
     call handle_ncerr(ncerr," inquire kg ")
     ncerr = nf90_put_var(ncid=ncid_hdr,varid=kg_id,values=kg_k,start=(/1,1,ikpt/),count=(/3,npw,1/))
     call handle_ncerr(ncerr," fill kg")
#endif
#if defined HAVE_TRIO_ETSF_IO
   else if (wff%accesswff == IO_MODE_ETSF) then
!    We write reduced_coordinates_of_plane_waves
     wave_folder%reduced_coordinates_of_plane_waves%data2D => kg_k
     wave_folder%red_coord_pw__kpoint_access = ikpt
     wave_folder%red_coord_pw__number_of_coefficients = npw
     call etsf_io_basisdata_put(wff%unwff, wave_folder, lstat, error)
     ETSF_CHECK_MYERROR(lstat,error)
#endif
   end if ! end if wff%accesswff
 else ! Still skip the record
   if (use_f90==1) then
     write(unitwf)
#if defined HAVE_MPI_IO
   else if (wff%accesswff==IO_MODE_MPI) then
     call xderiveWRecInit(wff,wff%spaceComm_mpiio,ios)
     call xderiveWRecEnd(wff,wff%spaceComm_mpiio,ios)
#endif
   end if
 end if

!---------------------------------------------------------------------------
!Write the third record: eigenvalues and occupations
!---------------------------------------------------------------------------

!===== Case formeig=0: write eigenvalues and occupations =====
 if (formeig==0) then
   if (use_f90==1) then
     write(unitwf) (eigen(iband),iband=1,nband_disk),(occ(iband),iband=1,nband_disk)
#if defined HAVE_MPI_IO
   else if(wff%accesswff==IO_MODE_MPI) then
     if (wff%me_mpiio==0) then
       call xderiveWRecInit(wff,ios)
       call xderiveWrite(wff,eigen,nband_disk,MPI_COMM_SELF,ios)
       call xderiveWrite(wff,occ,nband_disk,MPI_COMM_SELF,ios)
       call xderiveWRecEnd(wff,ios)
     end if
     call MPI_BCAST(wff%offwff,1,wff%offset_mpi_type,0,wff%spaceComm_mpiio,ios)
#endif
#if defined HAVE_TRIO_NETCDF
   else if (wff%accesswff==IO_MODE_NETCDF) then
!    In netcdf file saved as: eigen((2*mband)**formeig*mband,nkpt,nsppol)
     ncerr = nf90_inq_varid(ncid=ncid_hdr,name="eigen",varid=eigen_id)
     call handle_ncerr(ncerr," inquire eigen ")
!    Only get part of eigen for this kpoint ikpt and sppol isppol
     ncerr = nf90_put_var(ncid=ncid_hdr,varid=eigen_id,values=eigen,start=(/1,ikpt,isppol/),&
&     count=(/(2*nband_disk)**formeig*nband_disk,1,1/))
     call handle_ncerr(ncerr," put eigen ")
!    In netcdf file saved as: occ(mband,nkpt,nsppol)
     ncerr = nf90_inq_varid(ncid=ncid_hdr,name="occ",varid=occ_id)
     call handle_ncerr(ncerr," inquire occ ")
!    Only get part of occ for this kpoint ikpt and sppol isppol
     ncerr = nf90_put_var(ncid=ncid_hdr,varid=occ_id,values=occ,start=(/1,ikpt,isppol/),&
&     count=(/nband_disk,1,1/))
     call handle_ncerr(ncerr," put occ ")
#endif
#if defined HAVE_TRIO_ETSF_IO
   else if (wff%accesswff == IO_MODE_ETSF) then
     electrons_folder%eigenvalues%data1D => eigen
     electrons_folder%eigenvalues__kpoint_access = ikpt
     electrons_folder%eigenvalues__spin_access = isppol
     electrons_folder%eigenvalues__number_of_states = mband
     electrons_folder%occupations%data1D => occ
     electrons_folder%occupations__kpoint_access = ikpt
     electrons_folder%occupations__spin_access = isppol
     electrons_folder%occupations__number_of_states = mband
     call etsf_io_electrons_put(wff%unwff, electrons_folder, lstat, error)
     ETSF_CHECK_MYERROR(lstat,error)
#endif
   end if

!  ===== Case formeig=1: write matrix of eigenvalues =====
!  Will be written later (together with wave-functions)
 else if(formeig==1)then
 end if ! formeig

!---------------------------------------------------------------------------
!Write the wave-function coefficients
!---------------------------------------------------------------------------

!===== Case formeig=0: write only wave-functions =====
 if (formeig==0) then
!  If option=4, do not write wave functions
   if (option/=4) then
     if (use_f90==1) then
       do iband=1,nband_disk
         ipw=(iband-1)*npwso+icg
         write(unitwf) cg(1:2,ipw+1:ipw+npwso) ! VALGRIND complains some elements of cg are not initialized, but written
       end do
#if defined HAVE_MPI_IO
     else if(wff%accesswff==IO_MODE_MPI)then
       cg_ptr => cg ! Need pointer to bypass "inout" intent attribute
       call WffReadWrite_mpio(wff,2,cg_ptr,mcg,icg,nband_disk,npwso,npwsotot,ind_cg_mpi_to_seq,ios)
       nullify(cg_ptr)
#endif
#if defined HAVE_TRIO_NETCDF
     else if (wff%accesswff==IO_MODE_NETCDF) then
!      In netcdf file saved as: cg(2,mpw,nspinor,mband,nkpt,nsppol)
       ncerr = nf90_inq_varid(ncid=ncid_hdr,name="cg",varid=cg_id)
       call handle_ncerr(ncerr," inquire cg ")
!      Only get part of cg for this kpoint ikpt and sppol isppol
       ncerr = nf90_put_var(ncid=ncid_hdr,varid=cg_id,values=cg(:,icg+1:icg+npw*nspinor*mband),&
&       start=(/1,1,1,1,ikpt,isppol/),count=(/2,npw,nspinor,nband_disk,1,1/))
       call handle_ncerr(ncerr," put cg ")
#endif
#if defined HAVE_TRIO_ETSF_IO
     else if (wff%accesswff == IO_MODE_ETSF) then
!      !$main_folder%coefficients_of_wavefunctions%data2D => &
!      !! & cg(:, icg + 1:icg + npw * nspinor * nband)
!      See the read access, this direct association sometime leads to
!      segfaults with g95, so we use a temporary array.
       ABI_ALLOCATE(main_folder%coefficients_of_wavefunctions%data2D,(2,npw*nspinor*nband))
!      The following instruction leads to segmentation faults sometimes
!      I changed it for a longer expression
!      Tonatiuh Rangel 2009
!      main_folder%coefficients_of_wavefunctions%data2D = &
!      & cg(:, icg + 1:icg + npw * nspinor * nband)
       jj=0
       do ii=icg+1,icg+npw*nspinor*nband
         jj=jj+1
         main_folder%coefficients_of_wavefunctions%data2D(:,jj)=cg(:,ii)
       end do
       main_folder%wfs_coeff__kpoint_access = ikpt
       main_folder%wfs_coeff__spin_access = isppol
       main_folder%wfs_coeff__number_of_states = nband
       main_folder%wfs_coeff__number_of_coefficients = npw * nspinor
       call etsf_io_main_put(wff%unwff, main_folder, lstat, error)
       ETSF_CHECK_MYERROR(lstat,error)
       ABI_DEALLOCATE(main_folder%coefficients_of_wavefunctions%data2D)
#endif
     end if
   end if ! option/=4

!  ===== Case formeig=1: write eigenvalues and wave-functions =====
 else if(formeig==1)then
!  Not available for NETCDF and ETSF_IO
   nband2=2*nband_disk
   do iband=1,nband_disk
     ipw=(iband-1)*npwso+icg
     ii=(iband-1)*nband2
     if(use_f90==1)then
       write(unitwf) eigen(1+ii:nband2+ii)
       if (option/=4) then
         write(unitwf) cg(1:2,1+ipw:npwso+ipw)
       end if
#if defined HAVE_MPI_IO
     else if(wff%accesswff==IO_MODE_MPI)then
!      Should use an access with a "view"
       call xderiveWRecInit(wff,ios)
       call xderiveWrite(wff,eigen(ii:ii+nband2),nband2,wff%spaceComm_mpiio,ios)
       call xderiveWRecEnd(wff,ios)
       if (option/=4) then
         call xderiveWRecInit(wff,ios)
         call xderiveWrite(wff,cg(1:2,ipw+1:ipw+npwso),2,npwso,wff%spaceComm_mpiio,ios)
         call xderiveWRecEnd(wff,ios)
       end if
#endif
     end if
   end do
 end if ! formeig

!---------------------------------------------------------------------------
!Final statements
!---------------------------------------------------------------------------

#if defined HAVE_MPI_IO
 if (wff%accesswff==IO_MODE_MPI.and.option/=4)  then
   ABI_DEALLOCATE(ind_cg_mpi_to_seq)
 end if
#endif

 RETURN

!Silence compiler warning.
 ABI_UNUSED((/ii,mpi_enreg%me/))

end subroutine writewf
!!***
