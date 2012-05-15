!!****f* ABINIT/scphon_freq_to_normmode
!! NAME
!! scphon_freq_to_normmode
!!
!! FUNCTION
!! From updated phonon frequencies, and temperature, calculate displacement
!! of normal modes of phonons, adding an arbitrary sign (+- displacement)
!!
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! natom_primitive_cell=number of atoms in primitive cell (not supercell used
!!   for SC phonon calculation)
!! nphononq=number of phonon q-vectors input from anaddb run at equilibrium
!!   geometry
!! phonon_eigval=phonon eigenfrequencies, updated inside SC phonon run
!! scphon_temp= phononic temperature in Ha
!!
!! OUTPUT
!! normal_mode_displacements= calculated displacements of canonical coordinates
!!   (normal modes) of phonons at desired temperature.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!
!! SOURCE
! initialize the first normal mode displacements

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_freq_to_normmode (minusq_map,natom_primitive_cell,normal_mode_displacements,&
&   nphononq,phonon_eigval,scphon_temp)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_freq_to_normmode'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom_primitive_cell,nphononq
 real(dp),intent(in) :: scphon_temp
!arrays
 integer,intent(in) :: minusq_map(nphononq)
 real(dp),intent(in) :: phonon_eigval(3*natom_primitive_cell,nphononq)
 real(dp),intent(out) :: normal_mode_displacements(3*natom_primitive_cell,nphononq)

!Local variables-------------------------------
!scalars
 integer :: iminusq,imode_primitive_cell,iq,iseed=-10
 real(dp) :: bose_factor,random_sign

! *************************************************************************

!NOTE: what do we do about negative modes?

!NOTE: the only subtle part is choosing the sign of the displacement randomly
!Also presumes the mode is symmetric, which is only true for purely harmonic
!modes. As we have real displacements here, it may be an issue.
!One could also use a continuous variable instead of just +-1 for the
!amplitude of the normal coordinates...

 normal_mode_displacements=zero
 do iq=1, nphononq
!  
!  skip -q if it has already been filled
!  
   if ( abs(normal_mode_displacements(1,iq)) > tol10) cycle

   do imode_primitive_cell=1,3*natom_primitive_cell
!    skip gamma point acoustic modes
     if (abs(phonon_eigval(imode_primitive_cell,iq)) < tol10) cycle

!    for negative modes, use absolute value, until they eventually stabilize
!    also used (but not stated) in PRL
     bose_factor=one/(exp(abs(phonon_eigval(imode_primitive_cell,iq))/scphon_temp)-one)

!    always reseed random function
     random_sign=-one
     if (uniformrandom(iseed) > half) random_sign = one
     normal_mode_displacements(imode_primitive_cell,iq) = random_sign *&
&     sqrt( (half+bose_factor) / abs(phonon_eigval(imode_primitive_cell,iq)) )

!    TODO: MJV 24 11 2010: introduce cutoff ~ 1 Angstrom for using displacement
!    if so, use interpolated freq as for acoustic mode: q/q_max * w_Debye
!    q_max = Brillouin zone edge along direction of q
!    
!    NOTE: this is probably wrong in the PRL eq 6, as an atom mass is associated with a true
!    phonon mode, which moves all atoms a priori. The mass should be inserted
!    elsewhere, when the eigenvector comes in, which does have components by idir,
!    iatom
!    the sqrt(M) should be associated to the components of the eigenvectors, when
!    the cartesian displacements are calculated.
!    
   end do
!  
!  add normal_mode_displacement for -q
!  
   iminusq = minusq_map(iq)
   if (iminusq /= 0) then
     normal_mode_displacements(:,iminusq) = normal_mode_displacements(:,iq)
   end if

 end do

end subroutine scphon_freq_to_normmode
!!***



