!{\src2tex{textfont=tt}}
!!****f* ABINIT/normev
!! NAME
!! normev
!!
!! FUNCTION
!! Normalize a set of num eigenvectors of complex length ndim
!! (real length 2*ndim) and set phases to make evec(i,i) real and positive.
!! Near convergence, evec(i,j) approaches delta(i,j).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  evec(2*ndim,num)=num unnormalized eigenvectors
!!  ndim=dimension of evec as shown
!!  num=number of eigenvectors and complex length thereof.
!!
!! OUTPUT
!!  evec(2*ndim,num)=num normalized eigenvectors
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      subdiago
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine normev(evec,ndim,num)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'normev'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndim,num
!arrays
 real(dp),intent(inout) :: evec(2*ndim,num)

!Local variables-------------------------------
!scalars
 integer :: ii,jj
 real(dp) :: den,evim,evre,phim,phre,xnorm
 character(len=500) :: message
!DEBUG
!integer :: kk
!real(dp) :: scalprod
!ENDDEBUG

! *************************************************************************
!
!Loop over vectors
 do ii=1,num
!  find norm
   xnorm=0.0d0
   do jj=1,2*ndim
     xnorm=xnorm+evec(jj,ii)**2
   end do
   if((xnorm-one)**2>tol6)then
     write(message,'(6a,i6,a,es16.6,3a)' )ch10,&
&     ' normev : BUG ',ch10,&
&     '   Starting xnorm should be close to one (tol is 0.001).',ch10,&
&     '   However, for state number',ii,', xnorm=',xnorm,ch10,&
&     '   It might be that your LAPACK library has not been correctly installed.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  DEBUG
!  if(ii<7)then
!  write(std_out,*)' normev : ii,xnorm=',ii,xnorm
!  do kk=1,6
!  scalprod=zero
!  do jj=1,2*ndim,2
!  scalprod=scalprod+evec(jj,kk)*evec(jj,ii)+evec(jj+1,kk)*evec(jj+1,ii) &
!  &                         -evec(jj+1,kk)*evec(jj,ii)+evec(jj,kk)*evec(jj+1,ii)
!  enddo
!  write(std_out,*)' normev : kk,ii,scalar product=',kk,ii,scalprod
!  enddo
!  endif
!  ENDDEBUG

   xnorm=1.0d0/sqrt(xnorm)
!  Set up phase
   phre=evec(2*ii-1,ii)
   phim=evec(2*ii,ii)
   if (phim/=0.0d0) then
     den=1.0d0/sqrt(phre**2+phim**2)
     phre=phre*xnorm*den
     phim=phim*xnorm*den
   else
!    give xnorm the same sign as phre (negate if negative)
     phre=sign(xnorm,phre)
     phim=0.0d0
   end if
!  normalize with phase change
   do jj=1,2*ndim,2
     evre=evec(jj,ii)
     evim=evec(jj+1,ii)
     evec(jj,ii)=phre*evre+phim*evim
     evec(jj+1,ii)=phre*evim-phim*evre
   end do
 end do

!DEBUG
!write(std_out,*)' normev : exit '
!stop
!ENDDEBUG

end subroutine normev
!!***
