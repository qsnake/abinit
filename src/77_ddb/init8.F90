!{\src2tex{textfont=tt}}
!!****f* ABINIT/init8
!!
!! NAME
!! init8
!!
!! FUNCTION
!! Initialize the code : write heading and make the first i/os
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dscrpt=character string that describe the derivative database
!! filnam(3)=character strings giving file names
!! nddb=(=1 => will initialize the ddb, using an input GS file)
!!  (>1 => will merge the whole set of ddbs listed)
!!
!! OUTPUT
!!  (None)
!!
!! NOTES
!! 1. To be executed by one processor only.
!! 2. File names refer to following files, in order:
!!     (1) Output Derivative Database
!!    if nddb==1,
!!     (2) Formatted input file for the Corning ground-state code
!!    if nddb>1,
!!     (2 ... nddb+1) Derivative Databases to be added
!!
!! PARENTS
!!      mrgddb
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine init8(dscrpt,filnam,mddb,nddb)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init8'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mddb
 integer,intent(out) :: nddb
 character(len=fnlen),intent(out) :: dscrpt
!arrays
 character(len=fnlen),intent(out) :: filnam(mddb+1)

!Local variables -------------------------
!scalars
 integer :: iddb
 character(len=500) :: message

! *********************************************************************

!Read the name of the output ddb
 write(std_out,*)' Give name for output derivative database : '
 read(05, '(a)' ) filnam(1)
 write(std_out,'(a,a)' )' ',trim(filnam(1))

!Read the description of the derivative database
 write(std_out,*)' Give short description of the derivative database :'
 read(05, '(a)' )dscrpt
 write(std_out,'(a,a)' )' ',trim(dscrpt)

!Read the number of input ddbs, and check its value
 write(std_out,*)' Give number of input ddbs, or 1 if input GS file'
 read(05,*)nddb
 write(std_out,*)nddb
 if(nddb<=0.or.nddb>mddb)then
   write(message, '(a,a,a,a,i8,a,i8,a,a,a)' )&
&   ' init8 : ERROR -',ch10,&
&   '  nddb should be positive, >1 , and lower',&
&   '  than mddb =',mddb,' while the input nddb is ',nddb,'.',ch10,&
&   '  Action : correct the input nddb.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Read the file names
 if(nddb==1)then
   write(std_out,*)' Give name for ABINIT input file : '
   read(05, '(a)' ) filnam(2)
   write(std_out,'(a,a)' )' ',trim(filnam(2))
 else
   do iddb=1,nddb
     write(std_out,*)' Give name for derivative database number',&
&     iddb,' : '
     read(05, '(a)' ) filnam(iddb+1)
     write(std_out,'(a,a)' )' ',trim(filnam(iddb+1))
   end do
 end if

end subroutine init8
!!***
