!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkr8
!!
!! NAME
!! chkr8
!!
!! FUNCTION
!! This small subroutine check the identity of reali and realt,
!! who are integers, and eventually send a message and stop
!! if they are found unequal by more than tol
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! reali=first real number
!! intt=second  real number
!! character(len=6) name=name of the variable in the calling routine, to be echoed
!! tol=tolerance
!!
!! OUTPUT
!!  (only checking)
!!
!! PARENTS
!!      cmpar8
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chkr8(reali,realt,name,tol)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkr8'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 real(dp),intent(in) :: reali,realt,tol
 character(len=6),intent(in) :: name

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *********************************************************************

 if(abs(reali-realt)>tol) then
   write(message, '(a,a,a,a,a,a,a,es16.6,a,a,a,es16.6,a,a,a)' )&
&   ' chkr8 : ERROR -',ch10,&
&   '  Comparing reals for variable',name,'.',ch10,&
&   '  Value from input DDB is',reali,' and',ch10,&
&   '        from transfer DDB is',realt,'.',ch10,&
&   '  Action : check your DDBs.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

end subroutine chkr8
!!***
