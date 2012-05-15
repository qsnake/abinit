!{\src2tex{textfont=tt}}
!!****f* ABINIT/int2char
!! NAME
!! int2char
!!
!! FUNCTION
!! Convert a positive integer number (zero included) to ("2") a character(len=10),
!! with blanks to COMPLETE the string.
!! Exemple : 1234 will be mapped to "1234      "
!! Makes sure that the integer is between 0 and 9 999 999 999
!! Should be enough for integer*4
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  iint=integer to be converted
!!
!! OUTPUT
!!  string=character string
!!
!! TODO
!!  Should be included in m_fstrings
!! 
!! PARENTS
!!      gw_driver,m_bands_sym,m_dyson_solver,m_errors,m_qparticles,m_wfs
!!      prt_cif,wffile
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine int2char(iint,string)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'int2char'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iint
 character(len=10),intent(out) :: string

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *************************************************************************

!Note the use of floating numbers instead of large integers, for portability
 if(iint<0 .or. iint>=1.d10)then
   write(message, '(6a,i10)' ) ch10,&
&   ' int2char: ERROR -',ch10,&
&   '  The integer argument should be between 0 and 9999999999, while',ch10,&
&   '  it is ',iint
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(iint<10)then
   write(string,'(i1,9x)')iint
 else if(iint<100)then
   write(string,'(i2,8x)')iint
 else if(iint<1.0d3)then
   write(string,'(i3,7x)')iint
 else if(iint<1.0d4)then
   write(string,'(i4,6x)')iint
 else if(iint<1.0d5)then
   write(string,'(i5,5x)')iint
 else if(iint<1.0d6)then
   write(string,'(i6,4x)')iint
 else if(iint<1.0d7)then
   write(string,'(i7,3x)')iint
 else if(iint<1.0d8)then
   write(string,'(i8,2x)')iint
 else if(iint<1.0d9)then
   write(string,'(i9,1x)')iint
 else
   write(string,'(i10)')iint
 end if

end subroutine int2char
!!***
