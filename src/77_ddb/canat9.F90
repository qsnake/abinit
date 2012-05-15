!{\src2tex{textfont=tt}}
!!****f* ABINIT/canat9
!!
!! NAME
!! canat9
!!
!! FUNCTION
!! Transforms an atom whose coordinates (xred*rprim) would not be
!! in the chosen unit cell used to generate the interatomic forces
!! to its correspondent (rcan) in canonical coordinates.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (JCC,XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! brav= Bravais Lattice (1=S.C.;2=F.C.C.;3=BCC;4=Hex.)
!! natom= Number of atoms in the unit cell
!! rprim(3,3)= Normalized coordinates  of primitive vectors
!!
!! OUTPUT
!! rcan(3,natom)  = Atomic position in canonical coordinates
!! trans(3,natom) = Atomic translations : xred = rcan + trans
!!
!! PARENTS
!!      mkifc9
!!
!! CHILDREN
!!      leave_new,wrap2_pmhalf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine canat9(brav,natom,rcan,rprim,trans,xred)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'canat9'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav,natom
!arrays
 real(dp),intent(in) :: rprim(3,3),xred(3,natom)
 real(dp),intent(out) :: rcan(3,natom),trans(3,natom)

!Local variables -------------------------
!scalars
 integer :: found,iatom,ii
 character(len=500) :: message
!arrays
 real(dp) :: dontno(3,4),rec(3),rok(3),shift(3),tt(3)

! *********************************************************************

!Normalization of the cartesian atomic coordinates
!If not normalized : rcan(i) <- rcan(i) * acell(i)

 do iatom=1,natom
   rcan(:,iatom)=xred(1,iatom)*rprim(:,1)+&
&   xred(2,iatom)*rprim(:,2)+&
&   xred(3,iatom)*rprim(:,3)
 end do

!Study of the different cases for the Bravais lattice :

!Simple Cubic Lattice
 if (brav==1) then

   do iatom=1,natom

!    Canon will produces these coordinate transformations
!    (Note : here we still use reduced coordinates )
     call wrap2_pmhalf(xred(1,iatom),rok(1),shift(1))
     call wrap2_pmhalf(xred(2,iatom),rok(2),shift(2))
     call wrap2_pmhalf(xred(3,iatom),rok(3),shift(3))

!    New coordinates : rcan
     rcan(:,iatom)=rok(1)*rprim(:,1)+&
&     rok(2)*rprim(:,2)+&
&     rok(3)*rprim(:,3)
!    Translations between New and Old coordinates
     tt(:)=xred(1,iatom)*rprim(:,1)+&
&     xred(2,iatom)*rprim(:,2)+&
&     xred(3,iatom)*rprim(:,3)
     trans(:,iatom)=tt(:)-rcan(:,iatom)
   end do

 else if (brav==2) then
!  Face Centered Lattice
!  Special possible translations in the F.C.C. case
   dontno(:,:)=0.0_dp
   dontno(2,2)=0.5_dp
   dontno(3,2)=0.5_dp
   dontno(1,3)=0.5_dp
   dontno(3,3)=0.5_dp
   dontno(1,4)=0.5_dp
   dontno(2,4)=0.5_dp
   do iatom=1,natom
     found=0
     do ii=1,4
       if (found==1) exit
!      Canon will produces these coordinate transformations
       call wrap2_pmhalf(rcan(1,iatom)+dontno(1,ii),rok(1),shift(1))
       call wrap2_pmhalf(rcan(2,iatom)+dontno(2,ii),rok(2),shift(2))
       call wrap2_pmhalf(rcan(3,iatom)+dontno(3,ii),rok(3),shift(3))
!      In the F.C.C., ABS[ Ri ] + ABS[ Rj ] < or = 1/2
!      The equal signs has been treated using a tolerance parameter
!      not to have twice the same point in the unit cell !
       rok(1)=rok(1)-1.0d-10
       rok(2)=rok(2)-2.0d-10
       rok(3)=rok(3)-5.0d-10
       if (abs(rok(1))+abs(rok(2))<=0.5_dp) then
         if (abs(rok(1))+abs(rok(3))<=0.5_dp) then
           if (abs(rok(2))+abs(rok(3))<=0.5_dp) then
             tt(:)=rcan(:,iatom)
!            New coordinates : rcan
             rcan(1,iatom)=rok(1)+1.0d-10
             rcan(2,iatom)=rok(2)+2.0d-10
             rcan(3,iatom)=rok(3)+5.0d-10
!            Translations between New and Old coordinates
             trans(:,iatom)=tt(:)-rcan(:,iatom)
             found=1
           end if
         end if
       end if
     end do
   end do

 else if (brav==3) then
!  Body Centered Cubic Lattice
!  Special possible translations in the B.C.C. case

   dontno(:,1)=0.0_dp
   dontno(:,2)=0.5_dp
   do iatom=1,natom
     found=0
     do ii=1,2
       if (found==1) exit
!      Canon will produce these coordinate transformations
       call wrap2_pmhalf(rcan(1,iatom)+dontno(1,ii),rok(1),shift(1))
       call wrap2_pmhalf(rcan(2,iatom)+dontno(2,ii),rok(2),shift(2))
       call wrap2_pmhalf(rcan(3,iatom)+dontno(3,ii),rok(3),shift(3))
!      In the F.C.C., ABS[ Ri ] < or = 1/2
!      and    ABS[ R1 ] + ABS[ R2 ] + ABS[ R3 ] < or = 3/4
!      The equal signs has been treated using a tolerance parameter
!      not to have twice the same point in the unit cell !
       rok(1)=rok(1)-1.0d-10
       rok(2)=rok(2)-2.0d-10
       rok(3)=rok(3)-5.0d-10
       if(abs(rok(1))+abs(rok(2))+abs(rok(3))<=0.75_dp)then
         if ( abs(rok(1))<=0.5_dp .and.&
&         abs(rok(2))<=0.5_dp .and.&
&         abs(rok(3))<=0.5_dp       ) then
           tt(:)=rcan(:,iatom)
!          New coordinates : rcan
           rcan(1,iatom)=rok(1)+1.0d-10
           rcan(2,iatom)=rok(2)+2.0d-10
           rcan(3,iatom)=rok(3)+5.0d-10
!          Translations between New and Old coordinates
           trans(:,iatom)=tt(:)-rcan(:,iatom)
           found=1
         end if
       end if
     end do
   end do

 else if (brav==4) then
!  Hexagonal Lattice
!  In this case, it is easier first to work in reduced coor-
!  dinates space !
   do iatom=1,natom
!    Passage from the reduced space to the "lozenge" cell
     rec(1)=xred(1,iatom)-0.5_dp
     rec(2)=xred(2,iatom)-0.5_dp
     rec(3)=xred(3,iatom)
!    Canon will produces these coordinate transformations
     call wrap2_pmhalf(rec(1),rok(1),shift(1))
     call wrap2_pmhalf(rec(2),rok(2),shift(2))
     call wrap2_pmhalf(rec(3),rok(3),shift(3))
     rec(1)=rok(1)+0.5_dp
     rec(2)=rok(2)+0.5_dp
     rec(3)=rok(3)
!    Passage in Cartesian Normalized Coordinates
     rcan(:,iatom)=rec(1)*rprim(:,1)+&
&     rec(2)*rprim(:,2)+&
&     rec(3)*rprim(:,3)
!    Use of a tolerance parameter not to have twice the same point
!    in the unit cell !
     rcan(1,iatom)=rcan(1,iatom)-1.0d-10
     rcan(2,iatom)=rcan(2,iatom)-2.0d-10
!    Passage to the honey-com hexagonal unit cell !
     if (rcan(1,iatom)>0.5_dp) then
       rcan(1,iatom)=rcan(1,iatom)-1.0_dp
     end if
     if (rcan(1,iatom)>0.0_dp.and.rcan(1,iatom)+sqrt(3.0_dp)*rcan&
&     (2,iatom)>1.0_dp) then
       rcan(1,iatom)=rcan(1,iatom)-0.5_dp
       rcan(2,iatom)=rcan(2,iatom)-sqrt(3.0_dp)*0.5_dp
     end if
     if (rcan(1,iatom)<=0.0_dp.and.sqrt(3.0_dp)*rcan(2,iatom)-rcan&
&     (1,iatom)>1.0_dp) then
       rcan(1,iatom)=rcan(1,iatom)+0.5_dp
       rcan(2,iatom)=rcan(2,iatom)-sqrt(3.0_dp)*0.5_dp
     end if
!    Translations between New and Old coordinates
     tt(:)=xred(1,iatom)*rprim(:,1)+xred(2,iatom)*rprim(:,2)&
&     +xred(3,iatom)*rprim(:,3)
     trans(:,iatom)=tt(:)-rcan(:,iatom)
   end do

!  End of the possible cases for brav : 1, 2, 4.
 else

   write(message, '(a,a,a,i4,a,a,a)' )&
&   ' canat9 : BUG -',ch10,&
&   '  The required value of brav=',brav,' is not available.',ch10,&
&   '  It should be 1,2 or 4 .'
   call wrtout(std_out,message,'(COLL')
   call leave_new('COLL')

 end if

 write(message, '(a)' )' Canonical Atomic Coordinates '
 call wrtout(std_out,message,'COLL')

 do iatom=1,natom
   write(message, '(a,i5,3es18.8)' )&
&   ' atom',iatom,rcan(1,iatom),rcan(2,iatom),rcan(3,iatom)
   call wrtout(std_out,message,'COLL')
 end do

end subroutine canat9
!!***
