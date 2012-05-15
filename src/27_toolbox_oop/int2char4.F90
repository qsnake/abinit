!{\src2tex{textfont=tt}}
!!****f* ABINIT/int2char4
!! NAME
!! int2char4
!!
!! FUNCTION
!! Convert an integer number to ("2") a character(len=4)
!! Makes sure that the integer is between 0 and 9999
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
!!      contract_int_list,driver,dtfil_init1,gaus_dos,gstateimg,iofn1,m_atprj
!!      m_green,m_phonon_supercell,m_self,mrgscr,optic,pawmkaewf,prtfatbands
!!      read_wfrspa,scfcv,tddft,tetrahedron
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine int2char4(iint,string)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'int2char4'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iint
 character(len=4),intent(out) :: string

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *************************************************************************

 if(iint<0 .or. iint>9999)then
   write(message, '(a,a,a,a,a,a,i10)' ) ch10,&
&   ' int2char4: ERROR -',ch10,&
&   '  The integer argument should be between 0 and 9999, while',ch10,&
&   '  it is ',iint
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(iint<10)then
   write(string,'("000",i1)')iint
 else if(iint<100)then
   write(string,'("00",i2)')iint
 else if(iint<1000)then
   write(string,'("0",i3)')iint
 else
   write(string,'(i4)')iint
 end if

end subroutine int2char4
!!***
