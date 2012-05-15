!{\src2tex{textfont=tt}}
!!****f* ABINIT/dist9
!!
!! NAME
!! dist9
!!
!! FUNCTION
!! Compute the distance between atoms
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales by which rprim is to be multiplied
!! dist(natom,natom,nrpt)=distances between atoms
!! gprim(3,3)=dimensionless primitive translations in reciprocal space
!! natom=number of atoms in unit cell
!! nrpt= Number of R points in the Big Box
!! rcan(3,natom)=canonical coordinates of atoms
!! rprim(3,3)=dimensionless primitive translations in real space
!! rpt(3,nrpt)=cartesian coordinates of the points in the BigBox.
!!
!! OUTPUT
!! dist(natom,natom,nrpt)=distances between atoms
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


subroutine dist9(acell,dist,gprim,natom,nrpt,rcan,rprim,rpt)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dist9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
! real(dp),intent(in) :: ucvol
!scalars
 integer,intent(in) :: natom,nrpt
!arrays
 real(dp),intent(in) :: acell(3),gprim(3,3),rcan(3,natom),rprim(3,3)
 real(dp),intent(in) :: rpt(3,nrpt)
 real(dp),intent(out) :: dist(natom,natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,ii,irpt
!arrays
 real(dp) :: ra(3),rb(3),rdiff(3),red(3),rptcar(3),xred(3)

! *********************************************************************

!DEBUG
!write(std_out,*)' dist9 : enter '
!write(std_out,*)' natom,nrpt,ucvol=',natom,nrpt,ucvol
!write(std_out,*)'  rcan(1:3,1:natom)='
!do ia=1,natom
!write(std_out,'(i4,3es16.6)')ia,rcan(:,ia)
!end do
!write(std_out,*)' acell=',acell(:)
!write(std_out,*)' gprim=',gprim(:,:)
!write(std_out,*)' rprim=',rprim(:,:)
!ENDDEBUG

!BIG loop on all generic atoms
 do ia=1,natom

!  First transform canonical coordinates to reduced coordinates
   do ii=1,3
     xred(ii)=gprim(1,ii)*rcan(1,ia)+gprim(2,ii)*rcan(2,ia)&
&     +gprim(3,ii)*rcan(3,ia)
   end do

!  Then to cartesian coordinates
   ra(:)=xred(1)*acell(1)*rprim(:,1)+xred(2)*acell(2)*rprim(:,2)+&
&   xred(3)*acell(3)*rprim(:,3)

   do ib=1,natom
     do ii=1,3
       xred(ii)=gprim(1,ii)*rcan(1,ib)+gprim(2,ii)*rcan(2,ib)&
&       +gprim(3,ii)*rcan(3,ib)
     end do
     do ii=1,3
       rb(ii)=xred(1)*acell(1)*rprim(ii,1)+&
&       xred(2)*acell(2)*rprim(ii,2)+&
&       xred(3)*acell(3)*rprim(ii,3)
     end do

     do irpt=1,nrpt
!      First transform it to reduced coordinates
       do ii=1,3
         red(ii)=gprim(1,ii)*rpt(1,irpt)+gprim(2,ii)*rpt(2,irpt)&
&         +gprim(3,ii)*rpt(3,irpt)
       end do
!      Then to cartesian coordinates
       do ii=1,3
         rptcar(ii)=red(1)*acell(1)*rprim(ii,1)+&
&         red(2)*acell(2)*rprim(ii,2)+&
&         red(3)*acell(3)*rprim(ii,3)
       end do
       do ii=1,3
         rdiff(ii)=-rptcar(ii)+ra(ii)-rb(ii)
       end do
       dist(ia,ib,irpt)=(rdiff(1)**2+rdiff(2)**2+rdiff(3)**2)**0.5

     end do

   end do
 end do

!DEBUG
!write(std_out,*)' dist9 : exit '
!ENDDEBUG

end subroutine dist9
!!***
