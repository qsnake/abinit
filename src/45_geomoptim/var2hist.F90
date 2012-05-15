!{\src2tex{textfont=tt}}
!!****f* ABINIT/var2hist
!!
!! NAME
!! var2hist
!!
!! FUNCTION
!! Set the values of the history "hist"
!! with the values of xcart, xred, acell and rprimd
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
!! natom = number of atoms
!! xred(3,natom) = reduced dimensionless atomic
!!                           coordinates
!! xcart(3,natom)= cartesian dimensional atomic
!!                           coordinates (bohr)
!! acell(3)    = length scales of primitive translations (bohr)
!! rprimd(3,3) = dimensionlal real space primitive translations
!!               (bohr)
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
!!      | histEk:  Historical record of kinetic energy(Ek)
!!      | histT:   Historical record of time(T)
!!      | histR:   Historical record of rprimd(R)
!!      | histS:   Historical record of strten(S)
!!      | histV:   Historical record of velocity(V)
!!      | histXF:  Historical record of positions(X) and forces(F)
!!
!! PARENTS
!!      mover,pred_bfgs,pred_delocint,pred_diisrelax,pred_isokinetic
!!      pred_isothermal,pred_langevin,pred_moldyn,pred_nose,pred_srkna14
!!      pred_steepdesc,pred_verlet
!!
!! CHILDREN
!!      chkrprimd
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine var2hist(acell,hist,natom,rprim,rprimd,xcart,xred,zDEBUG)

 use m_profiling

! define dp,sixth,third,etc...
use defs_basis
! type(ab_movetype), type(ab_movehistory)
use defs_mover

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'var2hist'
 use interfaces_42_geometry
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom
type(ab_movehistory),intent(inout) :: hist
logical,intent(in) :: zDEBUG
!arrays
real(dp),intent(in) :: acell(3)
real(dp),intent(in) :: rprimd(3,3)
real(dp),intent(in) :: rprim(3,3)
real(dp),intent(in) :: xcart(3,natom)
real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
!scalars
integer :: kk

! *************************************************************

 hist%histXF(:,:,1,hist%ihist)=xcart(:,:)
 hist%histXF(:,:,2,hist%ihist)=xred(:,:)
 hist%histA(:,hist%ihist)     =acell(:)
 hist%histR(:,:,hist%ihist)   =rprimd(:,:)

 if(zDEBUG)then
   write (std_out,*) 'Atom positions and cell parameters '
   write (std_out,*) 'ihist: ',hist%ihist
   write (std_out,*) 'xcart:'
   do kk=1,natom
     write (std_out,*) xcart(:,kk)
   end do
   write (std_out,*) 'xred:'
   do kk=1,natom
     write (std_out,*) xred(:,kk)
   end do
   write(std_out,*) 'rprim:'
   do kk=1,3
     write(std_out,*) rprim(:,kk)
   end do
   write(std_out,*) 'rprimd:'
   do kk=1,3
     write(std_out,*) rprimd(:,kk)
   end do
   write(std_out,*) 'acell:'
   write(std_out,*) acell(:)
   call chkrprimd(acell,rprim,rprimd,std_out)
 end if

end subroutine var2hist
!!***
