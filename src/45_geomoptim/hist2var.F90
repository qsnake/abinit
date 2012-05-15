!{\src2tex{textfont=tt}}
!!****f* ABINIT/hist2var
!!
!! NAME
!! hist2var
!!
!! FUNCTION
!! Set the values of xcart, xred, acell and rprimd
!! with the values that comes from the history "hist"
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
!! natom = Number of atoms
!! hist<type ab_movehistory>=Historical record of positions, forces
!!      |                    acell, stresses, and energies,
!!      |                    contains:
!!      | mxhist:  Maximun number of records
!!      | histA:   Historical record of acell(A) and rprimd(R)
!!      | histE:   Historical record of energy(E)
!!      | histEk:  Historical record of Ionic kinetic energy(Ek)
!!      | histT:   Historical record of time(T) (For MD or iteration for GO)
!!      | histR:   Historical record of rprimd(R)
!!      | histS:   Historical record of strten(S)
!!      | histV:   Historical record of velocity(V)
!!      | histXF:  Historical record of positions(X) and forces(F)
!! zDebug = If true some output will be printed
!!
!! OUTPUT
!!  xred(3,natom) = reduced dimensionless atomic
!!                           coordinates
!!  xcart(3,natom)= cartesian dimensional atomic
!!                           coordinates (bohr)
!!  acell(3)    = length scales of primitive translations (bohr)
!!  rprim(3,3) = dimensionless real space primitive translations
!!  rprimd(3,3) = dimensional real space primitive translations
!!                (bohr)
!!
!! SIDE EFFECTS
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


subroutine hist2var(acell,hist,natom,rprim,rprimd,xcart,xred,zDEBUG)

 use m_profiling

! define dp,sixth,third,etc...
use defs_basis
! type(ab_movetype), type(ab_movehistory)
use defs_mover

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hist2var'
 use interfaces_42_geometry
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom
type(ab_movehistory),intent(in) :: hist
logical,intent(in) :: zDEBUG
!arrays
real(dp),intent(out) :: acell(3)
real(dp),intent(out) :: rprimd(3,3)
real(dp),intent(out) :: rprim(3,3)
real(dp),intent(out) :: xcart(3,natom)
real(dp),intent(out) :: xred(3,natom)

!Local variables-------------------------------
!scalars
integer :: jj,kk

! *************************************************************

 xcart(:,:) =hist%histXF(:,:,1,hist%ihist)
 xred(:,:)  =hist%histXF(:,:,2,hist%ihist)
 acell(:)   =hist%histA(:,hist%ihist)
 rprimd(:,:)=hist%histR(:,:,hist%ihist)

!Compute rprim from rprimd and acell
 do kk=1,3
   do jj=1,3
     rprim(jj,kk)=rprimd(jj,kk)/acell(kk)
   end do
 end do

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

end subroutine hist2var
!!***
