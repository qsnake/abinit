!{\src2tex{textfont=tt}}
!!****f* ABINIT/axial9
!!
!! NAME
!! axial9
!!
!! FUNCTION
!! Generates the local coordinates system from the
!! knowledge of the first vector (longitudinal) and
!! the ifc matrix in cartesian coordinates
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! ifccar(3,3)= matrix of interatomic force constants in cartesian coordinates
!! vect1(3)= cartesian coordinates of the first local vector
!!
!! OUTPUT
!! vect2(3)= cartesian coordinates of the second local vector
!! vect3(3)= cartesian coordinates of the third local vector
!!
!! PARENTS
!!      rsiaf9
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine axial9(ifccar,vect1,vect2,vect3)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'axial9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!arrays
 real(dp),intent(in) :: ifccar(3,3),vect1(3)
 real(dp),intent(out) :: vect2(3),vect3(3)

!Local variables -------------------------
!scalars
 integer :: flag,ii,itrial,jj
 real(dp) :: innorm,scprod
!arrays
 real(dp) :: work(3)

! *********************************************************************

 do jj=1,3
   work(jj)=0.0_dp
   do ii=1,3
     work(jj)=work(jj)+ifccar(jj,ii)*vect1(ii)
   end do
 end do

 flag=0
 do itrial=1,4
   scprod=0.0_dp
   do ii=1,3
     scprod=scprod+work(ii)*vect1(ii)
   end do

   do ii=1,3
     work(ii)=work(ii)-vect1(ii)*scprod
   end do

   scprod=0.0_dp
   do ii=1,3
     scprod=scprod+work(ii)**2
   end do

   if(scprod<1.0d-10)then
     work(1:3)=0.0_dp
     if(itrial>1)work(itrial-1)=1.0_dp
   else
     flag=1
   end if

   if(flag==1)exit
 end do

 innorm=scprod**(-0.5_dp)
 do ii=1,3
   vect2(ii)=work(ii)*innorm
 end do

 vect3(1)=vect1(2)*vect2(3)-vect1(3)*vect2(2)
 vect3(2)=vect1(3)*vect2(1)-vect1(1)*vect2(3)
 vect3(3)=vect1(1)*vect2(2)-vect1(2)*vect2(1)

end subroutine axial9
!!***
