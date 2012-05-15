!{\src2tex{textfont=tt}}
!!****f* ABINIT/canct9
!!
!! NAME
!! canct9
!!
!! FUNCTION
!! Convert from canonical coordinates to cartesian coordinates
!! a vector defined by its index=ib+natom*(irpt-1)
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
!! gprim(3,3)=dimensionless primitive translations in reciprocal space
!! index= index of the atom
!! natom=number of atoms in unit cell
!! nrpt= Number of R points in the Big Box
!! rcan(3,natom)=canonical coordinates of atoms
!! rprim(3,3)=dimensionless primitive translations in real space
!! rpt(3,nrpt)=canonical coordinates of the points in the BigBox.
!!
!! OUTPUT
!! ib=number of the atom in the unit cell
!! irpt= number of the unit cell to which belong the atom
!! rcart(3)=cartesian coordinate of the atom indexed by index.
!!
!! PARENTS
!!      hybrid9,rsiaf9
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine canct9(acell,gprim,ib,index,irpt,natom,nrpt,&
&       rcan,rcart,rprim,rpt)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'canct9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: index,natom,nrpt
 integer,intent(out) :: ib,irpt
!arrays
 real(dp),intent(in) :: acell(3),gprim(3,3),rcan(3,natom),rprim(3,3)
 real(dp),intent(in) :: rpt(3,nrpt)
 real(dp),intent(out) :: rcart(3)

!Local variables -------------------------
!scalars
 integer :: jj
!arrays
 real(dp) :: xred(3)

! *********************************************************************

 irpt=(index-1)/natom+1
 ib=index-natom*(irpt-1)

!Transform the canonical coordinates to reduced coord.
 do jj=1,3
   xred(jj)=gprim(1,jj)*(rpt(1,irpt)+rcan(1,ib))&
&   +gprim(2,jj)*(rpt(2,irpt)+rcan(2,ib))&
&   +gprim(3,jj)*(rpt(3,irpt)+rcan(3,ib))
 end do

!Then to cartesian coordinates (here the position of
!the atom b)

 do jj=1,3
   rcart(jj)=xred(1)*acell(1)*rprim(jj,1)+&
&   xred(2)*acell(2)*rprim(jj,2)+&
&   xred(3)*acell(3)*rprim(jj,3)
 end do

end subroutine canct9
!!***
