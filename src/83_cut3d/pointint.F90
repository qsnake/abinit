!{\src2tex{textfont=tt}}
!!****f* ABINIT/pointint
!! NAME
!! pointint
!!
!! FUNCTION
!! Computes the values at any point rr (this point is input from keyboard)
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (GMR,RC,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gridtt(nr1,nr2,nr3)=Total density
!! gridux(nr1,nr2,nr3)=spin-Up density, or magnetization density in X direction
!! griddy(nr1,nr2,nr3)=spin-Down density, or magnetization density in Y direction
!! gridmz(nr1,nr2,nr3)=spin-polarization density or magnetization density in Z direction
!! nr1=grid size along x
!! nr2=grid size along y
!! nr3=grid size along z
!! nspden=number of spin-density components
!! rprimd(3,3)=orientation of the unit cell in 3D
!!
!! OUTPUT
!!   only writing
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!      interpol3d,reduce
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine pointint(gridt,gridu,gridd,gridm,nr1,nr2,nr3,nspden,rprimd)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pointint'
 use interfaces_32_util
 use interfaces_83_cut3d, except_this_one => pointint
!End of the abilint section

 implicit none

!Arguments--------------------------------------------------------------
!scalars
 integer,intent(in) :: nr1,nr2,nr3,nspden
!arrays
 real(dp),intent(in) :: gridd(nr1,nr2,nr3),gridm(nr1,nr2,nr3)
 real(dp),intent(in) :: gridt(nr1,nr2,nr3),gridu(nr1,nr2,nr3),rprimd(3,3)

!Local variables--------------------------------------------------------
!scalars
 integer :: inpopt,mu,okinp
 real(dp) :: denvaldy,denvalmz,denvaltt,denvalux
!arrays
 real(dp) :: rcart(3),rr(3)

! *************************************************************************

 okinp=0
 do while (okinp==0)
   write(std_out,*) ' Select the coordinate system:'
   write(std_out,*) ' Type 1) for cartesian coordinates'
   write(std_out,*) '  or 2) for crystallographic coordinates'
   read(*,*) inpopt
   if (inpopt==1 .or. inpopt==2) okinp=1
 end do

 if (inpopt==1) then

   write(std_out,*) ' Input point Cartesian Coord:  X  Y  Z'
   read(*,*) rcart(1),rcart(2),rcart(3)
   call reduce(rr,rcart,rprimd)
   write(std_out,'(a,3es16.6)' ) ' Crystallographic coordinates: ',rr(1:3)

 else

   write(std_out,*) ' Input point Crystallographic Coord:  X  Y  Z'
   read(*,*) rr(1),rr(2),rr(3)

   do mu=1,3
     rcart(mu)=rprimd(mu,1)*rr(1)+rprimd(mu,2)*rr(2)+rprimd(mu,3)*rr(3)
   end do

   write(std_out,*) ' Cartesian coordinates : '
   write(std_out,'(3es16.6)' ) rcart(1),rcart(2),rcart(3)

 end if

!At this moment the code knows everything needed about the geometric input
!It will further proceed to calculate the interpolation at the demanded point

 rr(1)=mod(mod(rr(1),1._dp)+1._dp,1._dp)
 rr(2)=mod(mod(rr(2),1._dp)+1._dp,1._dp)
 rr(3)=mod(mod(rr(3),1._dp)+1._dp,1._dp)

 write(std_out,'(a,es16.6)' ) ' X coordinate, r1 is:',rr(1)
 write(std_out,'(a,es16.6)' ) ' Y coordinate, r2 is:',rr(2)
 write(std_out,'(a,es16.6)' ) ' Z coordinate, r3 is:',rr(3)

!devalt = total density value
!devalu = spin-up density value
!devald = spin-down density value
!devalm = magnetization density value
 call interpol3d(rr,nr1,nr2,nr3,denvaltt,gridt)
 if(nspden==2 .or. nspden==4)then
   call interpol3d(rr,nr1,nr2,nr3,denvalux,gridu)
   call interpol3d(rr,nr1,nr2,nr3,denvaldy,gridd)
   call interpol3d(rr,nr1,nr2,nr3,denvalmz,gridm)
 end if
 write(std_out,*)
 write(std_out,*)'---------------------------------------------'
 write(std_out,'(a,es16.6)') ' Non-spin-polarized value= ',denvaltt
 if(nspden==2)then
   write(std_out,'(a,es16.6)')' Spin-up value           = ',denvalux
   write(std_out,'(a,es16.6)')' Spin-down value         = ',denvaldy
   write(std_out,'(a,es16.6)')' Spin difference value   = ',denvalmz
 else if(nspden==4)then
   write(std_out,'(a,es16.6)')' x component             = ',denvalux
   write(std_out,'(a,es16.6)')' y component             = ',denvaldy
   write(std_out,'(a,es16.6)')' z component             = ',denvalmz
 end if
 write(std_out,*)'---------------------------------------------'

end subroutine pointint
!!***
