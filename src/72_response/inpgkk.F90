!{\src2tex{textfont=tt}}
!!****f* ABINIT/inpgkk
!! NAME
!! inpgkk
!!
!! FUNCTION
!! read in gkk file and return eigenvalue matrix
!!   Only works for a single gkk matrix (1 perturbation and qpoint) in the file
!!   like the files produced by outgkk
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  
!!  filegkk= filename
!!  mpi_enreg=informations about MPI parallelization
!!
!! OUTPUT
!!  eigen1 = response function 1st order eigenvalue matrix
!!
!! PARENTS
!!      read_el_veloc
!!
!! CHILDREN
!!      hdr_clean,hdr_io,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine inpgkk(bantot1,eigen1,filegkk,hdr1,mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

 use m_header,          only : hdr_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inpgkk'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bantot1
 character(len=fnlen),intent(in) :: filegkk
 type(MPI_type),intent(in) :: mpi_enreg
 type(hdr_type), intent(out) :: hdr1
!arrays
 real(dp),intent(out) :: eigen1(bantot1)

!Local variables-------------------------------
!scalars
 integer :: isppol, ikpt, rdwr, mband, ikb, ios
 integer :: unitgkk, fform, ierr, n1wf, i1wf
 type(hdr_type) :: hdr0

 real(dp), allocatable :: eigen(:)
 character(len=500) :: message

! *************************************************************************

!Only for not having an unused arg
 if(.false.)write(std_out,*)mpi_enreg%me

 rdwr = 5

 unitgkk = 35
 open(unit=unitgkk,file=filegkk,form='unformatted',status='old',iostat=ios)
 if (ios/=0) then
   write(message,'(5a)')&
&   ' inpgkk: ERROR- ',ch10,&
&   ' opening file: ',trim(filegkk),' as old'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 rewind (unitgkk)

!read in header of GS file and eigenvalues
 call hdr_io(fform,hdr0,rdwr,unitgkk)
 
 mband = maxval(hdr0%nband(:))
 ABI_ALLOCATE(eigen,(mband))
 write(message,'(a)')'inpgkk : try to reread GS eigenvalues'
 call wrtout(std_out,message,'COLL')

 do isppol=1,hdr0%nsppol
   do ikpt=1,hdr0%nkpt
     read (unitgkk,IOSTAT=ierr) eigen(1:hdr0%nband(ikpt))
     if (ierr /= 0) write(std_out,*) 'error reading eigen from gkk file'
   end do
 end do
 read(unitgkk,IOSTAT=ierr) n1wf
 if (ierr /= 0) write(std_out,*) 'error reading eigen from gkk file'
 ABI_DEALLOCATE(eigen)
 call hdr_clean(hdr0)

 if (n1wf > 1) then
   write(message,'(6a)')ch10, &
&   ' inpgkk : ERROR- ',ch10, &
&   ' several 1wf records were found in the file,',ch10, &
&   ' which is not allowed for reading with this routine'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!read in header of 1WF file
 call hdr_io(fform,hdr1,rdwr,unitgkk)
 if (fform == 0) then
   write(message,'(4a,i4,a)')ch10,&
&   ' inpgkk : ERROR- ',ch10,&
&   ' 1WF header number ',i1wf,' was mis-read. fform == 0'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if (bantot1 < 2*hdr1%nsppol*hdr1%nkpt*mband**2) then
   write(message,'(4a,2I10)')ch10,&
&   ' inpgkk : ERROR- ',ch10,&
&   ' input size for eigenvalue matrix is not large enough ', bantot1, 2*hdr1%nsppol*hdr1%nkpt*mband**2
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!retrieve 1WF <psi_k+q | H | psi_k> from gkk file and echo to output
 ikb = 0
 do isppol=1,hdr1%nsppol
   do ikpt=1,hdr1%nkpt
     read (unitgkk,IOSTAT=ierr) eigen1(ikb+1:ikb+2*hdr1%nband(ikpt)**2)
     ikb = ikb + 2*hdr1%nband(ikpt)**2
     if (ierr /= 0) write(std_out,*) 'error reading eigen1 from gkk file',isppol,ikpt
   end do
 end do

end subroutine inpgkk
!!***
