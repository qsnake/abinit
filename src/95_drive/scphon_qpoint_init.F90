!!****f* ABINIT/scphon_qpoint_init
!! NAME
!! scphon_qpoint_init
!!
!! FUNCTION
!! Initialize the qpoint grid which should be read in for the frequencies
!! and eigenvectors of the equilibrium primitive unit cell, and used for FT
!! of supercell quantities.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nphononq=number of phonon q-vectors input from anaddb run at equilibrium
!!   geometry
!! supercell_multiplicity=number of times the primitive unit cell is repeated
!!
!! OUTPUT
!! phononq= phonon q vectors, should be the same as those used in anaddb run
!!   (reduced coordinates)
!! FIXME: add a check to make sure this is the case.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!      smpbz,wrap2_pmhalf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_qpoint_init (nphononq,phononq,supercell_multiplicity)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_qpoint_init'
 use interfaces_32_util
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nphononq
!arrays
 integer,intent(in) :: supercell_multiplicity(3)
 real(dp),intent(out) :: phononq(3,nphononq)

!Local variables-------------------------------
!scalars
 integer :: brav,iout,iqpt,mqpt,nqpt,nqshft,option
 real(dp) :: res
!character(len=500) :: message
!arrays
 integer :: qptrlatt(3,3)
 real(dp) :: kpt(3),qshift(3)

! *************************************************************************

!always impose no shift in q
 nqshft=1
 qshift=(/zero,zero,zero/)

 qptrlatt = 0
 qptrlatt(1,1) = supercell_multiplicity(1)
 qptrlatt(2,2) = supercell_multiplicity(2)
 qptrlatt(3,3) = supercell_multiplicity(3)

 mqpt= qptrlatt(1,1)*qptrlatt(2,2)*qptrlatt(3,3) &
& +qptrlatt(1,2)*qptrlatt(2,3)*qptrlatt(3,1) &
& +qptrlatt(1,3)*qptrlatt(2,1)*qptrlatt(3,2) &
& -qptrlatt(1,2)*qptrlatt(2,1)*qptrlatt(3,3) &
& -qptrlatt(1,3)*qptrlatt(2,2)*qptrlatt(3,1) &
& -qptrlatt(1,1)*qptrlatt(2,3)*qptrlatt(3,2)

 iout = 6
 option=1
 brav=1
 call smpbz(brav,iout,qptrlatt,mqpt,nqpt,&
& nqshft,option,qshift,phononq)

 if (mqpt /= nqpt .or. nqpt /= nphononq) then
   write(std_out,*) 'scphon_qpoint_init: Error in num of qpt: mqpt,nqpt,nphononq = ',&
&   mqpt,nqpt,nphononq
   stop
 end if

!FIXME: should check these qpt are in the same order as the frequencies and
!eigenvectors!!!!
!
!reduce spqpt to correct zone
 do iqpt=1,nphononq
   call wrap2_pmhalf(phononq(1,iqpt),kpt(1),res)
   call wrap2_pmhalf(phononq(2,iqpt),kpt(2),res)
   call wrap2_pmhalf(phononq(3,iqpt),kpt(3),res)
   phononq(:,iqpt) = kpt
 end do

end subroutine scphon_qpoint_init
!!***

