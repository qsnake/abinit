!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawvhnzc
!! NAME
!! pawvhnzc
!!
!! FUNCTION
!! PAW: compute Hartree potential for n_{Zc}
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (JWZ, MT, GJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  ncore(radmesh_core%mesh_size)=atomic core density
!!  radmesh_core <type(pawrad_type)>=radial mesh (and related data) for the core densities
!!  znucl= valence and total charge of the atomic species
!!
!! OUTPUT
!!  vhnzc(radmesh_core%mesh_size)=Hartree potential due to Z_nc
!!
!! PARENTS
!!      pawdij0,pawkij,psp7in
!!
!! CHILDREN
!!      deducer0,poisson
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawvhnzc(ncore,radmesh_core,vhnzc,znucl)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

 use m_radmesh,   only : poisson, deducer0

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawvhnzc'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 real(dp),intent(in) :: znucl
 type(pawrad_type),intent(in) :: radmesh_core
!arrays
 real(dp),intent(in) :: ncore(radmesh_core%mesh_size)
 real(dp), intent(out) :: vhnzc(radmesh_core%mesh_size)

!Local variables ---------------------------------------
  real(dp) :: intg
  real(dp),allocatable :: nwk(:)
! *********************************************************************

 DBG_ENTER("COLL")

 ABI_ALLOCATE(nwk,(radmesh_core%mesh_size))

 nwk(:)=ncore(:)*four_pi*radmesh_core%rad(:)**2
 call poisson(nwk,0,intg,radmesh_core,vhnzc)
 vhnzc(2:radmesh_core%mesh_size)=(vhnzc(2:radmesh_core%mesh_size)-znucl)/radmesh_core%rad(2:radmesh_core%mesh_size)
 call deducer0(vhnzc,radmesh_core%mesh_size,radmesh_core)

 ABI_DEALLOCATE(nwk)

 DBG_EXIT("COLL")

 end subroutine pawvhnzc
!!***
