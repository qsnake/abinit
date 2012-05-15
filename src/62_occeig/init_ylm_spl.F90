!{\src2tex{textfont=tt}}
!!****f* ABINIT/init_ylm_spl
!! NAME
!! init_ylm_spl
!!
!! FUNCTION
!! Pre-calculate the integral of y*y*j_v(y) for recip_ylm
!!     NOTE: spherical Bessel function small j!
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mbessint=  max number of points on grid for integral
!!  bessargmax=max point to which well integrate
!!  bessint_delta = space between integral arguments
!!  mlang=     max angular momentum
!!
!! OUTPUT
!!  spl_bessint(mbessint,mlang)= array of integrals
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      besjm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine init_ylm_spl(mbessint,bessargmax,bessint_delta,mlang,spl_bessint)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_ylm_spl'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mbessint,mlang
 real(dp),intent(in) :: bessargmax,bessint_delta
!arrays
 real(dp),intent(out) :: spl_bessint(mbessint,mlang)

!Local variables -------------------------
! function calls (NumRec)
!scalars
 integer :: ii,iintarg,iloop,iloop2,ll,nbess_int,nbesspoints,nfiner,nloops
 integer :: nrest,nstep
 real(dp) :: arg_bess,bess_int,coef1,coef2,coef3,cosdelta_bess,delta_bess
 real(dp) :: i_coef1,i_coef2,i_coef3,sindelta_bess,tol
!arrays
 real(dp),allocatable :: bessj_arr(:,:),cosbessx(:),sinbessx(:),x_bess(:)

! *********************************************************************

 tol = 1.0d-10

!DEBUG
!write(std_out,*) 'init_bess_ylm enter '
!ENDDEBUG

!-----------------------------------------------------------------
!Bessel function into array
!-----------------------------------------------------------------

!integration grid is nfiner times finer than the interpolation grid
 nfiner = 7
 delta_bess = bessint_delta / nfiner
 nbesspoints = mbessint * nfiner
 ABI_ALLOCATE(bessj_arr,(nbesspoints,mlang))
 ABI_ALLOCATE(x_bess,(nbesspoints))
 ABI_ALLOCATE(sinbessx,(nbesspoints))
 ABI_ALLOCATE(cosbessx,(nbesspoints))

 sindelta_bess = sin(delta_bess)
 cosdelta_bess = cos(delta_bess)
 nstep = 40
 nloops = int((nbesspoints)/nstep)
 nrest  = nbesspoints - nstep*nloops
 iintarg = 0
 do iloop=1,nloops
   iintarg = iintarg + 1
   x_bess(iintarg) = (iintarg-1)*delta_bess
   sinbessx(iintarg) = sin(x_bess(iintarg))
   cosbessx(iintarg) = cos(x_bess(iintarg))
   do iloop2=2,nstep
     iintarg = iintarg + 1
     x_bess(iintarg) = x_bess(iintarg-1)+delta_bess
!    get sin and cos of x_bess arguments
     sinbessx(iintarg) = sinbessx(iintarg-1)*cosdelta_bess &
&     + cosbessx(iintarg-1)*sindelta_bess
     cosbessx(iintarg) = cosbessx(iintarg-1)*cosdelta_bess &
&     - sinbessx(iintarg-1)*sindelta_bess
   end do
 end do
 do iloop2=1,nrest
   iintarg = iintarg + 1
   x_bess(iintarg) = x_bess(iintarg-1)+delta_bess
   sinbessx(iintarg) = sinbessx(iintarg-1)*cosdelta_bess &
&   + cosbessx(iintarg-1)*sindelta_bess
   cosbessx(iintarg) = cosbessx(iintarg-1)*cosdelta_bess &
&   - sinbessx(iintarg-1)*sindelta_bess
 end do

 do ll=1,mlang

   bessj_arr(1,ll) = zero
   call besjm(one,bessj_arr(2:nbesspoints,ll),cosbessx(2:nbesspoints),   &
&   ll,nbesspoints-1,sinbessx(2:nbesspoints),x_bess(2:nbesspoints))

   do iintarg=1,nbesspoints
     bessj_arr(iintarg,ll) = x_bess(iintarg)*x_bess(iintarg)*bessj_arr(iintarg,ll)
!    DEBUG
!    write(std_out,*) ' bess x2 ', ll, x_bess, bessj_arr(iintarg,ll)
!    ENDDEBUG
   end do
 end do

!-----------------------------------------------------------------
!Bessel function integral
!-----------------------------------------------------------------

 coef1 = 9.0_dp  / 24.0_dp
 coef2 = 28.0_dp / 24.0_dp
 coef3 = 23.0_dp / 24.0_dp
 i_coef1 = one - coef1
 i_coef2 = one - coef2
 i_coef3 = one - coef3
!first point is 0, second is calculated here
 spl_bessint(1, :) = zero
!DEBUG
!write(std_out,*) ' bess int  1  0.0 0.0'
!write(std_out,*) ' bess int  2  0.0 0.0'
!write(std_out,*) ' bess int  3  0.0 0.0'
!ENDDEBUG

 nbess_int = nfiner+1
 do ll=1,mlang
   bess_int =            coef1*bessj_arr(1,ll)
   bess_int = bess_int + coef2*bessj_arr(2,ll)
   bess_int = bess_int + coef3*bessj_arr(3,ll)

   bess_int = bess_int + coef3*bessj_arr(nbess_int-2,ll)
   bess_int = bess_int + coef2*bessj_arr(nbess_int-1,ll)
   bess_int = bess_int + coef1*bessj_arr(nbess_int  ,ll)

   do ii=4, nbess_int-3
     bess_int = bess_int + bessj_arr(ii,ll)
   end do

   spl_bessint(2, ll) = bess_int
!  DEBUG
!  write(std_out,*) ' bess int ', ll, arg_bess, bess_int * delta_bess
!  ENDDEBUG
 end do


!remaining points: add last terms to integral sum
 do iintarg=3,mbessint
   arg_bess = (iintarg-1)*bessargmax/mbessint
   nbess_int = (iintarg-1)*nfiner + 1
!  DEBUG
!  write(std_out,'(a,2F16.8)') ' arg_bess ', &
!  &                           arg_bess
!  ENDDEBUG
!  fill bessel function arrays
   do ll=1, mlang

!    do Simpson O(1/N^4)  from NumRec in C p 134
     bess_int = spl_bessint(iintarg-1, ll)
     bess_int = bess_int + i_coef3*bessj_arr(nbess_int-nfiner-2,ll)
     bess_int = bess_int + i_coef2*bessj_arr(nbess_int-nfiner-1,ll)
     bess_int = bess_int + i_coef1*bessj_arr(nbess_int-nfiner  ,ll)

     do ii=nbess_int-nfiner+1, nbess_int-3
       bess_int = bess_int + bessj_arr(ii,ll)
     end do

     bess_int = bess_int + coef3*bessj_arr(nbess_int-2,ll)
     bess_int = bess_int + coef2*bessj_arr(nbess_int-1,ll)
     bess_int = bess_int + coef1*bessj_arr(nbess_int  ,ll)
!    DEBUG
!    write(std_out,'(a,2I6,3F16.8)') ' ipw, ll, arg_bess, delta_bess, bess_int = ', &
!    &                               ipw, ll, arg_bess, delta_bess, bess_int
!    ENDDEBUG

     spl_bessint(iintarg, ll) = bess_int
!    DEBUG
!    write(std_out,*) ' bess int ', ll, arg_bess, bess_int * delta_bess
!    ENDDEBUG

   end do

 end do

 do iintarg=1,mbessint
   do ll=1,mlang
     spl_bessint(iintarg, ll) = spl_bessint(iintarg, ll) * delta_bess
   end do
 end do


end subroutine init_ylm_spl
!!***
