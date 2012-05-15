!!****m* ABINIT/interfaces_64_atompaw
!! NAME
!! interfaces_64_atompaw
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/64_atompaw
!!
!! COPYRIGHT
!! Copyright (C) 2010-2011 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_64_atompaw

 implicit none

interface
 subroutine pawdij0(indlmn,kij,lmnmax,ncore,opt_init,pawtab,radmesh,radmesh_core,radmesh_vloc,vhtnzc,znucl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: opt_init
  type(pawtab_type),intent(inout) :: pawtab
  type(pawrad_type),intent(in) :: radmesh
  type(pawrad_type),intent(in) :: radmesh_core
  type(pawrad_type),intent(in) :: radmesh_vloc
  real(dp),intent(in) :: znucl
  integer,intent(in) :: indlmn(6,lmnmax)
  real(dp),intent(in) :: kij(pawtab%lmn2_size)
  real(dp),intent(in) :: ncore(radmesh_core%mesh_size)
  real(dp),intent(in) :: vhtnzc(radmesh_vloc%mesh_size)
 end subroutine pawdij0
end interface

interface
 subroutine pawkij(indlmn,kij,lmnmax,ncore,opt_init,opt_vhnzc,pawtab,radmesh,radmesh_core,radmesh_vloc,vhtnzc,znucl)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: opt_init
  integer,intent(in) :: opt_vhnzc
  type(pawtab_type),intent(in) :: pawtab
  type(pawrad_type),intent(in) :: radmesh
  type(pawrad_type),intent(in) :: radmesh_core
  type(pawrad_type),intent(in) :: radmesh_vloc
  real(dp),intent(in) :: znucl
  integer,intent(in) :: indlmn(6,lmnmax)
  real(dp),intent(out) :: kij(pawtab%lmn2_size)
  real(dp),intent(in) :: ncore(radmesh_core%mesh_size)
  real(dp),intent(in) :: vhtnzc(radmesh_vloc%mesh_size)
 end subroutine pawkij
end interface

interface
 subroutine pawshpfun(ll,mesh,norm,pawtab,shapefunc)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ll
  type(pawrad_type),intent(in) :: mesh
  real(dp),intent(out) :: norm
  type(pawtab_type),intent(in) :: pawtab
  real(dp),intent(inout) :: shapefunc(mesh%mesh_size)
 end subroutine pawshpfun
end interface

interface
 subroutine pawvhnzc(ncore,radmesh_core,vhnzc,znucl)
  use defs_basis
  use defs_datatypes
  implicit none
  type(pawrad_type),intent(in) :: radmesh_core
  real(dp),intent(in) :: znucl
  real(dp),intent(in) :: ncore(radmesh_core%mesh_size)
  real(dp), intent(out) :: vhnzc(radmesh_core%mesh_size)
 end subroutine pawvhnzc
end interface

interface
 subroutine shapebes(al,ql,ll,rc)
  use defs_basis
  implicit none
  integer :: ll
  real(dp) :: rc
  real(dp) :: al(2)
  real(dp) :: ql(2)
 end subroutine shapebes
end interface

end module interfaces_64_atompaw
!!***
