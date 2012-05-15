!{\src2tex{textfont=tt}}
!!****f* ABINIT/vel2hist
!!
!! NAME
!! vel2hist
!!
!! FUNCTION
!! Set the values of the history "hist"
!! related with velocities, ie 
!! The array of velocities and Kinetic Energy
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! amass = mass of the atoms
!! natom = number of atoms
!! vel(3,natom)= Velocities of the atoms
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist<type ab_movehistory>=Historical record of positions, forces
!!      |                    acell, stresses, and energies,
!!      |                    contains:
!!      | mxhist:  Maximun number of records
!!      | ihist:   Index of present record
!!      | histA:   Historical record of acell(A) and rprimd(R)
!!      | histE:   Historical record of energy(E)
!!      | histEk:  Historical record of Ionic kinetic energy(Ek)
!!      | histT:   Historical record of time(T) (For MD or iteration for GO)
!!      | histR:   Historical record of rprimd(R)
!!      | histS:   Historical record of strten(S)
!!      | histV:   Historical record of velocity(V)
!!      | histXF:  Historical record of positions(X) and forces(F)
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine vel2hist(amass,hist,natom,vel)

 use m_profiling

! define dp,sixth,third,etc...
use defs_basis
! type(ab_movetype), type(ab_movehistory)
use defs_mover

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vel2hist'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom
type(ab_movehistory),intent(inout) :: hist
!arrays
real(dp),intent(in) :: vel(3,natom)
real(dp),intent(in) :: amass(natom)

!Local variables-------------------------------
!scalars
integer :: ii,jj
real(dp) :: ekin 

! *************************************************************

!Store the velocities
 hist%histV(:,:,hist%ihist)=vel(:,:)

!Compute the Ionic Kinetic energy
 ekin=0.0
 do ii=1,natom
   do jj=1,3
     ekin=ekin+0.5_dp*amass(ii)*vel(jj,ii)**2
   end do
 end do

!Store the Ionic Kinetic Energy
 hist%histEk(hist%ihist)=ekin

end subroutine vel2hist
!!***
