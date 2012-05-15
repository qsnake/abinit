!{\src2tex{textfont=tt}}
!!****f* ABINIT/dtech9
!!
!! NAME
!! dtech9
!!
!! FUNCTION
!! Reads the Dielectric Tensor and the Effective Charges in the
!! Gamma Block coming from the Derivative Data Base.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (JCC,XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! natom= number of atoms in unit cell
!! iblok= index of the Gamma block
!! mpert =maximum number of ipert
!! nblok= number of blocks in the DDB
!! blkval(2,3*mpert*3*mpert,nblok)=  dynamical matrices
!!  In our case, the nblok is restricted to iblok
!!
!! OUTPUT
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement. Note the following convention:
!!  zeff(electric field direction, atomic direction, atom number)
!! dielt(3,3)=dielectric tensor
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dtech9(blkval,dielt,iblok,mpert,natom,nblok,zeff)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtech9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 real(dp),intent(out) :: dielt(3,3),zeff(3,3,natom)

!Local variables -------------------------
!scalars
 integer :: depl,elec,elec1,elec2,iatom

! *********************************************************************

!Extration of effectives charges
 do iatom=1,natom
   do elec=1,3
     do depl=1,3
       zeff(elec,depl,iatom)=0.5*&
&       (blkval(1,depl,iatom,elec,natom+2,iblok)+&
&       blkval(1,elec,natom+2,depl,iatom,iblok))
     end do
   end do
 end do

!Extration of dielectric tensor
 do elec1=1,3
   do elec2=1,3
     dielt(elec1,elec2)=blkval(1,elec1,natom+2,&
&     elec2,natom+2,iblok)
   end do
 end do

 write(std_out,'(/,a,/,3es16.6,/,3es16.6,/,3es16.6)' )&
& ' Dielectric Tensor ',&
& dielt(1,1),dielt(1,2),dielt(1,3),&
& dielt(2,1),dielt(2,2),dielt(2,3),&
& dielt(3,1),dielt(3,2),dielt(3,3)

 write(std_out,'(/,a)' ) ' Effectives Charges '
 do iatom=1,natom
   write(std_out,'(a,i4,3es16.6,/,3es16.6,/,3es16.6)' )' atom ',iatom,&
&   zeff(1,1,iatom),zeff(1,2,iatom),zeff(1,3,iatom),&
&   zeff(2,1,iatom),zeff(2,2,iatom),zeff(2,3,iatom),&
&   zeff(3,1,iatom),zeff(3,2,iatom),zeff(3,3,iatom)
 end do

end subroutine dtech9
!!***
