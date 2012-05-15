!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkilwf
!!
!! NAME
!! chkilwf
!!
!! FUNCTION
!! Open input file for the lwf code, then
!! reads or echoes the input information.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! grdsize(3)= size of the grid of q points = limit of the shells of the LWF
!!
!! OUTPUT
!! All the other arguments are outputs
!! and are read from the input file. NOT TRUE ! SOME ARGUMENTS ARE DUMMY AND ARE DECLARED AS INOUT.
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      lwf
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chkilwf(alpha,decflg,enwdmax,enwdmin,frozflg,grdsize,irwfl,&
& nqpt,subwdmax,subwdmin)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkilwf'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: decflg,frozflg,irwfl,nqpt
 real(dp),intent(inout) :: alpha,enwdmax,enwdmin,subwdmax,subwdmin
!arrays
 integer,intent(in) :: grdsize(3)

!Local variables ------------------------------
!scalars
 character(len=500) :: message
!arrays

!*********************************************************************

!DEBUG
!write(std_out,*) ' chkilwf : enter'
!ENDDEBUG

 if (alpha>one) then
   write(std_out,'(a,f12.7,a)') ' chkilwf : WARNING : alpha is ',alpha,' and it should be less or equal to 1'
   write(std_out,*) ' Assume experienced user !'
 end if
 if (alpha<zero) then
   write(message, '(a,a,a,a,f12.7,a,a,a,a)' ) ch10,&
&   ' chkilwf : ERROR -',ch10,&
&   '  alpha is negative: ',alpha,ch10,&
&   '  This is not allowed.  ',ch10,&
&   '  Action : modify the energy windows in the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

 if (enwdmax<enwdmin) then
   write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&   ' chkilwf : ERROR -',ch10,&
&   '  The global window has negative width.',ch10,&
&   '  This is not allowed.  ',ch10,&
&   '  Action : modify the energy windows in the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if
 if (frozflg==1) then
   if (subwdmax<subwdmin) then
     write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&     ' chkilwf  : ERROR -',ch10,&
&     '  The frozen states  window has negative width.',ch10,&
&     '  This is not allowed.  ',ch10,&
&     '  Action : modify the energy windows in the input file.'
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if
   if (subwdmax>enwdmax) then
     write(message, '(a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&     ' chkilwf  : ERROR -',ch10,&
&     '  The frozen-states energy window has an upper limit',ch10,&
&     '  that is greater than the global energy window.',ch10,&
&     '  This is not allowed.  ',ch10,&
&     '  Action : modify the energy windows in the input file.'
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if
   if (subwdmin<enwdmin) then
     write(message, '(a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&     ' chkilwf  : ERROR -',ch10,&
&     '  The frozen-states energy window has a lower limit',ch10,&
&     '  that is smaller than the global energy window.',ch10,&
&     '  This is not allowed.  ',ch10,&
&     '  Action : modify the energy windows in the input file.'
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if
 end if

 if (grdsize(1)*grdsize(2)*grdsize(3)>nqpt) then
   write(message, '(a,a,a,a,i7,a,a,i7,a,a)' ) ch10,&
&   ' chkilwf  : WARNING -',ch10,&
&   '  The extent of the lattice Wannier function',grdsize(1)*grdsize(2)*grdsize(3),ch10,&
&   '  is larger than the grid of Q points',nqpt,ch10,&
&   '  Assume experienced user !'
   call wrtout(std_out,  message,'COLL')
 end if

 if (.not.( (irwfl==1) .or. (irwfl==-1) )) then
   write(message, '(a,a,a,a,i5,a,a,a,a,a,a)' ) ch10,&
&   ' chkilwf  : ERROR -',ch10,&
&   '  The read/write integer flag is',irwfl,ch10,&
&   '  and it may have only +1/-1 values',ch10,&
&   '  This is not allowed.  ',ch10,&
&   '  Action : modify irwfl in the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

 if (.not.( (decflg==0) .or. (decflg==1))) then
   write(message, '(a,a,a,a,i5,a,a,a,a,a,a)' ) ch10,&
&   ' chkilwf  : ERROR -',ch10,&
&   '  The flag for the forced stop at the minimum of Omega_I ',decflg,ch10,&
&   '  and it may have only 0 or 1 values',ch10,&
&   '  This is not allowed.  ',ch10,&
&   '  Action : modify decflg in the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if


!DEBUG
!write(std_out,*) ' chkilwf: exit'
!ENDDEBUG

end subroutine chkilwf
!!***
