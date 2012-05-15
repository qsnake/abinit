!{\src2tex{textfont=tt}}
!!****f* ABINIT/contract_int_le_val
!! NAME
!! contract_int_le_val
!!
!! FUNCTION
!!  "Design by contract" routine, ensuring that the
!!  target integer argument is lower or equal to some limiting value.
!!  Might be used to test composite quantities :
!!  integer_name might be "arg1-arg2*arg3", with corresponding input value
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  character(len=*) calling_routine = name of the routine of which the argument is checked
!!  character(len=*) integer_name = name of the integer variable that is checked
!!  integer_value = value that is checked
!!  limit = larger admitted limit for the integer
!!
!! OUTPUT
!!  (only checking, writing and stopping)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  "Design by contract" routines are simpler than the routines
!!  that check the input variables. Indeed, the error message
!!  can be much primitive.
!!
!! PARENTS
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine contract_int_le_val(calling_routine,integer_name,integer_value,limit)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'contract_int_le_val'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: integer_value,limit
 character(len=*),intent(in) :: calling_routine,integer_name

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*)' contract_int_le_val : enter '
!ENDDEBUG

 if(integer_value>limit)then
   write(message,'(10a,i6,5a,i6,a)') ch10,&
&   ' contract_int_le_val: BUG -',ch10,&
&   '  For the routine ',trim(calling_routine),',',ch10,&
&   '  the value of "',trim(integer_name),'" should be lower or equal to ',limit,'.',ch10,&
&   '  However, ',trim(integer_name),' = ',integer_value,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!DEBUG
!write(std_out,*)' contract_int_le_val : exit'
!stop
!ENDDEBUG

end subroutine contract_int_le_val
!!***
