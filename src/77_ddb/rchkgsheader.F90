!{\src2tex{textfont=tt}}
!!****f* ABINIT/rchkgsheader
!!
!! NAME 
!! rchkgsheader
!!
!! FUNCTION
!! This routine reads the GS header information in the GKK file and checks it
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  natom = number of atoms from DDB, for check
!!  kptirr_phon = coordinates of the irreducible kpoints close to the FS
!!
!! OUTPUT
!!  hdr = header information
!!  nband = number of bands for rest of calculation
!!          should be the same for all kpts
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      hdr_io
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rchkGSheader (hdr,natom,nband,unitgkk)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rchkGSheader'
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,unitgkk
 integer,intent(out) :: nband
 type(hdr_type),intent(out) :: hdr

!Local variables-------------------------------
!scalars
 integer :: fform,ikpt,rdwr
 character(len=500) :: message

! *************************************************************************
!
!read in general header of _GKK file
!this is where we get nkpt, ngkpt(:,:)... which are also read in
!rdddb9 and inprep8. Probably should do some checking to avoid
!using ddb files from other configurations
!
 rewind(unitgkk)
 rdwr = 5 !read in header of file without rewinding it
 call hdr_io(fform,hdr,rdwr,unitgkk)
 ABI_CHECK(fform/=0," GKK header mis-read. fform == 0")

 if (hdr%natom /= natom) then
   MSG_ERROR('natom in gkk file is different from anaddb input')
 end if

 do ikpt=1,hdr%nkpt
   if (hdr%nband(ikpt) /= hdr%nband(1)) then
     write (message,'(3a)')&
&     ' Use the same number of bands for all kpts : ',ch10,&
&     ' could have spurious effects if efermi is too close to the last band '
     MSG_ERROR(message)
   end if
 end do

 rdwr = 4 ! echo header to screen
 call hdr_io(fform,hdr,rdwr,6)

 nband=hdr%nband(1)

end subroutine rchkGSheader
!!***
