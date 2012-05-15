!{\src2tex{textfont=tt}}
!!****f* ABINIT/dymfz9
!!
!! NAME
!! dymfz9
!!
!! FUNCTION
!! As the subroutine canatm has transformed the coordinates of the
!! atoms in normalized canonical coordinates, the corresponding
!! dynamical matrix should be multiplied by a phase shift corresponding
!! to the translation between New and Old coordinates of its two
!! corresponding atoms.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (JCC,XG,MM)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dynmat = non-phase shifted dynamical matrices
!! natom = number of atoms
!! nqpt = number of qpoints
!! gprim = reciprocal lattice vectors (cartesian but dimensionless)
!! option=1 : the matrices are transformed from the old (tn)
!!  coordinate system to the new (normalized canonical)
!!        2 : the matrices are restored from the normalized
!!  canonical coordinate system to the usual (tn) one...
!! rcan = canonical coordinates of atoms
!! spqpt = qpoint coordinates (reduced reciprocal)
!! trans = Atomic translations : xred = rcan + trans
!!
!! OUTPUT
!! dynmat = phase shifted dynamical matrices
!!
!! PARENTS
!!      gtdyn9,mkifc9
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dymfz9(dynmat,natom,nqpt,gprim,option,&
&    spqpt,trans)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dymfz9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nqpt,option
!arrays
 real(dp),intent(in) :: gprim(3,3),spqpt(3,nqpt),trans(3,natom)
 real(dp),intent(inout) :: dynmat(2,3,natom,3,natom,nqpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,iqpt,mu,nu
 real(dp) :: im,ktrans,re
!arrays
 real(dp) :: kk(3)

! *********************************************************************

 do iqpt=1,nqpt

!  Definition of q in normalized reciprocal space
   kk(1)=spqpt(1,iqpt)*gprim(1,1)+spqpt(2,iqpt)*gprim(1,2)+&
&   spqpt(3,iqpt)*gprim(1,3)
   kk(2)=spqpt(1,iqpt)*gprim(2,1)+spqpt(2,iqpt)*gprim(2,2)+&
&   spqpt(3,iqpt)*gprim(2,3)
   kk(3)=spqpt(1,iqpt)*gprim(3,1)+spqpt(2,iqpt)*gprim(3,2)+&
&   spqpt(3,iqpt)*gprim(3,3)
   if(option==1)then
     kk(1)=-kk(1)
     kk(2)=-kk(2)
     kk(3)=-kk(3)
   end if
   do ia=1,natom
     do ib=1,natom
!      Product of q with the differences between the two atomic translations
       ktrans=kk(1)*(trans(1,ia)-trans(1,ib))+kk(2)*(trans(2,ia)-&
&       trans(2,ib))+kk(3)*(trans(3,ia)-trans(3,ib))

!      DEBUG
!      if(ia==1 .and. (ib==2 .or. ib==3) )then
!      write(std_out,'(a,3i3,es16.8)' )'iqpt,ia,ib,ktrans',iqpt,ia,ib,ktrans
!      write(std_out,'(a,3es16.8)' )'           kk    ',kk(1:3)
!      write(std_out,'(a,3es16.8)' )'           trans ',trans(1:3,ia)-trans(1:3,ib)
!      write(std_out,'(a,3es16.8)' )'           transa',trans(1:3,ia)
!      write(std_out,'(a,3es16.8)' )'           transb',trans(1:3,ib)
!      end if
!      ENDDEBUG

!      "!OCL SCALAR" needed for VPP Fujitsu machine, inserted by MM 19990722
!      OCL SCALAR
       do mu=1,3
         do nu=1,3
           re=dynmat(1,mu,ia,nu,ib,iqpt)
           im=dynmat(2,mu,ia,nu,ib,iqpt)
!          Transformation of the Old dynamical matrices by New ones by multi-
!          plication by a phase shift
           dynmat(1,mu,ia,nu,ib,iqpt)=re*cos(two_pi*ktrans)-&
&           im*sin(two_pi*ktrans)
           dynmat(2,mu,ia,nu,ib,iqpt)=re*sin(two_pi*ktrans)+&
&           im*cos(two_pi*ktrans)

!          DEBUG
!          if((ia==2 .or. ia==3) .and. ib==1)then
!          write(std_out,'(5i3,2es16.8)' )&
!          &       mu,ia,nu,ib,iqpt,dynmat(1:2,mu,ia,nu,ib,iqpt)
!          end if
!          ENDDEBUG

         end do
       end do
     end do
   end do
 end do

end subroutine dymfz9
!!***
