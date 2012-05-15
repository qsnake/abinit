!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxcsum
!! NAME
!! pawxcsum
!!
!! FUNCTION
!! PAW only:
!! Compute useful sums of moments of densities needed to compute on-site contributions
!! to XC energy and potential
!!  First order sums:
!!    Sum1(1)=Sum_L{Rho1_L(r)**2}
!!    Sum1(2)=Sum_L{Rho1_L(r)*Rho2_L(r)}
!!    Sum1(3)=Sum_L{Rho2_L(r)**2}
!!    With L>0
!!  Second order sums:
!!    Sum2(L,1)=Sum_L1_L2{Rho1_L1(r)*Rho1_L2(r)*Gaunt_(L,L1,L2)}
!!    Sum2(L,2)=Sum_L1_L2{Rho1_L1(r)*Rho2_L2(r)*Gaunt_(L,L1,L2)}
!!    Sum2(L,3)=Sum_L1_L2{Rho2_L1(r)*Rho2_L2(r)*Gaunt_(L,L1,L2)}
!!    With L1>0, L2>0
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!! This routine has been written from rhohxc
!!
!! INPUTS
!!  lmselect1(lm_size)=select the non-zero LM-moments of input density Rho1
!!  lmselect2(lm_size)=select the non-zero LM-moments of input density Rho2
!!                     ignored if ndens=1
!!  lm_size=number of moments of the density
!!  ndens=number of input densities (1 or 2)
!!  nrad=number of radial points
!!  option= 1: compute first order sums
!!          2: compute first and second order sums
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  rho1(nrad,lm_size)=moments of first density on each radial point
!!  rho2(nrad,lm_size)=moments of 2nd density on each radial point
!!                     ignored if ndens=1
!!
!! OUTPUT
!!  sum1(nrad,2*ndens-1)=first order sums
!!  === if option>=2
!!    sum2(nrad,lm_size,2*ndens-1)=second order sums
!!
!! PARENTS
!!      pawxcm,pawxcmpositron,poslifetime
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawxcsum(lmselect1,lmselect2,lm_size,ndens,nrad,option,pawang,rho1,rho2,sum1,sum2)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawxcsum'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lm_size,ndens,nrad,option
!arrays
 logical,intent(in) :: lmselect1(lm_size),lmselect2(lm_size)
 real(dp),intent(in) :: rho1(nrad,lm_size),rho2(nrad,lm_size)
 real(dp),intent(out) :: sum1(nrad,2*ndens-1),sum2(nrad,lm_size,(2*ndens-1)*(option/2))
 type(pawang_type),intent(in) :: pawang

!Local variables-------------------------------
!scalars
 integer :: ilm,ilm1,ilm2,isel
 real(dp) :: fact
 character(len=500) :: msg
!arrays

!************************************************************************

 DBG_ENTER("COLL")

 if(ndens/=1.and.ndens/=2) then
   msg='  Only 1 or 2 densities !'
 end if
 if(pawang%gnt_option==0) then
   msg='  pawang%gnt_option=0 !'
   MSG_BUG(msg)
 end if

 if (option>=1) then

!  SUM1(r)= Sum_L{Rho1_L(r)*Rho2_L(r)} (L>0)
!  --------------------------------------------------
   sum1=zero

   if (ndens==1) then

!    Case: one density
     do ilm=2,lm_size
       if (lmselect1(ilm)) then
         sum1(:,1)=sum1(:,1)+rho1(:,ilm)**2
       end if
     end do
   else

!    Case: two densities
     do ilm=2,lm_size
       if (lmselect1(ilm)) then
         sum1(:,1)=sum1(:,1)+rho1(:,ilm)**2
         if (lmselect2(ilm)) sum1(:,2)=sum1(:,2)+rho1(:,ilm)*rho2(:,ilm)
       end if
       if (lmselect2(ilm)) sum1(:,3)=sum1(:,3)+rho2(:,ilm)**2
     end do
   end if

 end if !option

 if (option>=2) then

!  SUM2(r,L)= Sum_L1_L2{Rho1_L1(r)*Rho2_L2(r)*Gaunt_(L,L1,L2)}  (L1>0, L2>0)
!  --------------------------------------------------
   sum2=zero

   do ilm=1,lm_size
     do ilm1=2,lm_size
       if (lmselect1(ilm1)) then
         do ilm2=2,ilm1
           if (lmselect1(ilm2)) then
             isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
             if (isel>0) then
               fact=pawang%realgnt(isel);if (ilm1/=ilm2) fact=two*fact
               sum2(:,ilm,1)=sum2(:,ilm,1)+fact*rho1(:,ilm1)*rho1(:,ilm2)
             end if
           end if
         end do
       end if
     end do
   end do

!  Case: two densities
   if (ndens==2) then
     do ilm=1,lm_size
       do ilm1=2,lm_size
         if (lmselect2(ilm1)) then
           do ilm2=2,ilm1
             if (lmselect2(ilm2)) then
               isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
               if (isel>0) then
                 fact=pawang%realgnt(isel);if (ilm1/=ilm2) fact=two*fact
                 sum2(:,ilm,3)=sum2(:,ilm,3)+fact*rho2(:,ilm1)*rho2(:,ilm2)
               end if
             end if
           end do
         end if
       end do
       do ilm1=2,lm_size
         if (lmselect1(ilm1)) then
           do ilm2=2,ilm1
             if (lmselect2(ilm2)) then
               isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
               if (isel>0) then
                 fact=pawang%realgnt(isel)
                 sum2(:,ilm,2)=sum2(:,ilm,2)+fact*rho1(:,ilm1)*rho2(:,ilm2)
               end if
             end if
           end do
           if (ilm1<lm_size) then
             do ilm2=ilm1+1,lm_size
               if (lmselect2(ilm2)) then
                 isel=pawang%gntselect(ilm,ilm1+ilm2*(ilm2-1)/2)
                 if (isel>0) then
                   fact=pawang%realgnt(isel)
                   sum2(:,ilm,2)=sum2(:,ilm,2)+fact*rho1(:,ilm1)*rho2(:,ilm2)
                 end if
               end if
             end do
           end if
         end if
       end do
     end do
   end if ! ndens==2

 end if !option

 DBG_EXIT("COLL")

 end subroutine pawxcsum
!!***
