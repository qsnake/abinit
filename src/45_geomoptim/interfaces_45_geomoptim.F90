!!****m* ABINIT/interfaces_45_geomoptim
!! NAME
!! interfaces_45_geomoptim
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/45_geomoptim
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

module interfaces_45_geomoptim

 implicit none

interface
 subroutine brdene(etotal,etotal_prev,hessin,ndim,vin,vin_prev,vout,vout_prev)
  use defs_basis
  implicit none
  integer,intent(in) :: ndim
  real(dp),intent(in) :: etotal
  real(dp),intent(inout) :: etotal_prev
  real(dp),intent(in) :: hessin(ndim,ndim)
  real(dp),intent(inout) :: vin(ndim)
  real(dp),intent(inout) :: vin_prev(ndim)
  real(dp),intent(in) :: vout(ndim)
  real(dp),intent(inout) :: vout_prev(ndim)
 end subroutine brdene
end interface

interface
 subroutine calc_b_matrix(deloc,natom,rprimd,xcart,b_matrix)
  use defs_basis
  use m_delocint
  implicit none
  integer,intent(in) :: natom
  type(ab_delocint),intent(inout) :: deloc
  real(dp),intent(out) :: b_matrix(deloc%ninternal,3*natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine calc_b_matrix
end interface

interface
 subroutine dbond_length_d1(r1,r2,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
 end subroutine dbond_length_d1
end interface

interface
 subroutine dang_d1(r1,r2,r3,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
 end subroutine dang_d1
end interface

interface
 subroutine dang_d2(r1,r2,r3,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
 end subroutine dang_d2
end interface

interface
 subroutine ddihedral_d1(r1,r2,r3,r4,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
  real(dp),intent(in) :: r4(3)
 end subroutine ddihedral_d1
end interface

interface
 subroutine ddihedral_d2(r1,r2,r3,r4,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
  real(dp),intent(in) :: r4(3)
 end subroutine ddihedral_d2
end interface

interface
 subroutine calc_prim_int(deloc,natom,rprimd,xcart,prim_int)
  use defs_basis
  use m_delocint
  implicit none
  integer,intent(in) :: natom
  type(ab_delocint),intent(inout) :: deloc
  real(dp),intent(out) :: prim_int(deloc%ninternal)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine calc_prim_int
end interface

interface
 function bond_length(r1,r2)
  use defs_basis
  implicit none
  real(dp) :: bond_length
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
 end function bond_length
end interface

interface
 function angle_ang(r1,r2,r3)
  use defs_basis
  implicit none
  real(dp) :: angle_ang
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
 end function angle_ang
end interface

interface
 function angle_dihedral(r1,r2,r3,r4)
  use defs_basis
  implicit none
  real(dp) :: angle_dihedral
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
  real(dp),intent(in) :: r4(3)
 end function angle_dihedral
end interface

interface
 subroutine deloc2xcart(deloc,natom,rprimd,xcart,&  
  &  deloc_int,btinv,u_matrix)
  use defs_basis
  use m_delocint
  implicit none
  integer,intent(in) :: natom
  type(ab_delocint),intent(inout) :: deloc
  real(dp),intent(out) :: btinv(3*(natom-1),3*natom)
  real(dp),intent(in) :: deloc_int(3*(natom-1))
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: u_matrix(deloc%ninternal,3*(natom-1))
  real(dp),intent(inout) :: xcart(3,natom)
 end subroutine deloc2xcart
end interface

interface
 subroutine fappnd(filapp,filnam,iapp,&  
  &  suff) ! optional argument
  use defs_basis
  implicit none
  integer,intent(in) :: iapp
  character(len=fnlen),intent(out) :: filapp
  character(len=fnlen),intent(in) :: filnam
  character(len=3),optional,intent(in) :: suff
 end subroutine fappnd
end interface

interface
 subroutine fcart2fred(fcart,fred,rprimd,natom)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  real(dp),intent(in) :: fcart(3,natom)
  real(dp),intent(out) :: fred(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine fcart2fred
end interface

interface
 subroutine fred2fcart(favg,fcart,fred,gprimd,jellslab,natom)
  use defs_basis
  implicit none
  integer,intent(in) :: jellslab
  integer,intent(in) :: natom
  real(dp),intent(out) :: favg(3)
  real(dp),intent(out) :: fcart(3,natom)
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine fred2fcart
end interface

interface
 subroutine fred2fdeloc(btinv,deloc_force,fred,natom,gprimd)
  use defs_basis
  implicit none
  integer, intent(in) :: natom
  real(dp),intent(in) :: btinv(3*(natom-1),3*natom)
  real(dp),intent(out) :: deloc_force(3*(natom-1))
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine fred2fdeloc
end interface

interface
 subroutine hessinit(ab_mover, hessin, init_matrix, ndim, ucvol)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: ndim
  type(ab_movetype),intent(in) :: ab_mover
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: hessin(ndim,ndim)
  real(dp),intent(in) :: init_matrix(3,3)
 end subroutine hessinit
end interface

interface
 subroutine hessupdt(hessin,iatfix,natom,ndim,vin,vin_prev,vout,vout_prev)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndim
  real(dp),intent(inout) :: hessin(ndim,ndim)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: vin(ndim)
  real(dp),intent(in) :: vin_prev(ndim)
  real(dp),intent(in) :: vout(ndim)
  real(dp),intent(in) :: vout_prev(ndim)
 end subroutine hessupdt
end interface

interface
 subroutine hist2var(acell,hist,natom,rprim,rprimd,xcart,xred,zDEBUG)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  type(ab_movehistory),intent(in) :: hist
  logical,intent(in) :: zDEBUG
  real(dp),intent(out) :: acell(3)
  real(dp),intent(out) :: rprim(3,3)
  real(dp),intent(out) :: rprimd(3,3)
  real(dp),intent(out) :: xcart(3,natom)
  real(dp),intent(out) :: xred(3,natom)
 end subroutine hist2var
end interface

interface
 subroutine hist_compare(hist_in,hist_out,natom,similar,tolerance)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(out) :: similar
  type(ab_movehistory),intent(in) :: hist_in
  type(ab_movehistory),intent(inout) :: hist_out
  real(dp),intent(in) :: tolerance
 end subroutine hist_compare
end interface

interface
 subroutine isotemp(amass,dtion,ekin,iatfix,ktemp,mttk_vars,natom,nnos,qmass,vel)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nnos
  real(dp),intent(in) :: dtion
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  integer,intent(in) :: iatfix(:,:)
  real(dp),intent(in) :: amass(natom)
  real(dp),intent(in) :: qmass(:)
  real(dp),intent(inout) :: vel(3,natom)
 end subroutine isotemp
end interface

interface
 subroutine isopress(amass,bmass,dtion,ekin,iatfix,ktemp,mttk_vars,natom,nnos,qmass,&  
  &  strten,strtarget,ucvol,vel,vlogv)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nnos
  real(dp),intent(in) :: bmass
  real(dp),intent(in) :: dtion
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  real(dp),intent(in) :: ucvol
  real(dp),intent(inout) :: vlogv
  integer,intent(in) :: iatfix(:,:)
  real(dp),intent(in) :: amass(natom)
  real(dp),intent(in) :: qmass(:)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(inout) :: vel(3,natom)
 end subroutine isopress
end interface

interface
 subroutine isostress(amass,bmass,dtion,ekin,iatfix,ktemp,mttk_vars,natom,nnos,&  
  &  qmass,strten,strtarget,ucvol,vel)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nnos
  real(dp),intent(in) :: bmass
  real(dp),intent(in) :: dtion
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  real(dp),intent(in) :: ucvol
  integer, intent(in) :: iatfix(:,:)
  real(dp),intent(in) :: amass(natom)
  real(dp),intent(in) :: qmass(:)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(inout) :: vel(3,natom)
 end subroutine isostress
end interface

interface
 subroutine make_angles(deloc,icenter,natom)
  use m_delocint
  implicit none
  integer,intent(in) :: icenter
  integer,intent(in) :: natom
  type(ab_delocint),intent(inout) :: deloc
 end subroutine make_angles
end interface

interface
 subroutine make_angles_new(angles,bonds,natom,ntypat,rprimd,typat,xcart,znucl)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(go_angles),intent(inout) :: angles
  type(go_bonds),intent(in) :: bonds
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine make_angles_new
end interface

interface
 subroutine make_bonds(deloc,natom,ntypat,icenter,rprimd,typat,xcart,znucl)
  use defs_basis
  use m_delocint
  implicit none
  integer,intent(in) :: icenter
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(ab_delocint),intent(inout) :: deloc
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(in) :: znucl(:)
 end subroutine make_bonds
end interface

interface
 subroutine make_bonds_new(bonds,natom,ntypat,rprimd,typat,xcart,znucl)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(go_bonds),intent(inout) :: bonds
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine make_bonds_new
end interface

interface
 subroutine make_dihedrals(badangles,deloc,icenter)
  use m_delocint
  implicit none
  integer,intent(in) :: icenter
  type(ab_delocint),intent(inout) :: deloc
  integer,intent(in) :: badangles(deloc%nang)
 end subroutine make_dihedrals
end interface

interface
 subroutine make_prim_internals(deloc,icenter,natom,ntypat,rprimd,&  
  &  typat,xcart,znucl)
  use defs_basis
  use m_delocint
  implicit none
  integer,intent(in) :: icenter
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(ab_delocint),intent(inout) :: deloc
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(in) :: znucl(:)
 end subroutine make_prim_internals
end interface

interface
 subroutine pimd_langevin_npt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&  
  &  rprimd,rprimd_next,rprimd_prev,trotter,vel,volume,xred,xred_next,xred_prev)
  use m_pimd
  use defs_basis
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: prtvolimg
  integer,intent(in) :: trotter
  type(pimd_type),intent(in) :: pimd_param
  real(dp),intent(in) :: volume
  real(dp),intent(in) :: etotal(trotter)
  real(dp),intent(inout) :: forces(3,natom,trotter)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: rprimd_next(3,3)
  real(dp),intent(in) :: rprimd_prev(3,3)
  real(dp),intent(inout) :: vel(3,natom,trotter)
  real(dp),intent(in),target :: xred(3,natom,trotter)
  real(dp),intent(out) :: xred_next(3,natom,trotter)
  real(dp),intent(in),target :: xred_prev(3,natom,trotter)
 end subroutine pimd_langevin_npt
end interface

interface
 subroutine pimd_langevin_nvt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&  
  &  rprimd,trotter,vel,volume,xred,xred_next,xred_prev)
  use m_pimd
  use defs_basis
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: prtvolimg
  integer,intent(in) :: trotter
  type(pimd_type),intent(in) :: pimd_param
  real(dp),intent(in) :: volume
  real(dp),intent(in) :: etotal(trotter)
  real(dp),intent(inout) :: forces(3,natom,trotter)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: vel(3,natom,trotter)
  real(dp),intent(in),target :: xred(3,natom,trotter)
  real(dp),intent(out) :: xred_next(3,natom,trotter)
  real(dp),intent(in),target :: xred_prev(3,natom,trotter)
 end subroutine pimd_langevin_nvt
end interface

interface
 subroutine pimd_nosehoover_npt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&  
  &  rprimd,rprimd_next,rprimd_prev,trotter,vel,volume,xred,xred_next,xred_prev)
  use m_pimd
  use defs_basis
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: prtvolimg
  integer,intent(in) :: trotter
  type(pimd_type),intent(in) :: pimd_param
  real(dp),intent(in) :: volume
  real(dp),intent(in) :: etotal(trotter)
  real(dp),intent(inout) :: forces(3,natom,trotter)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: rprimd_next(3,3)
  real(dp),intent(in) :: rprimd_prev(3,3)
  real(dp),intent(inout) :: vel(3,natom,trotter)
  real(dp),intent(in),target :: xred(3,natom,trotter)
  real(dp),intent(out) :: xred_next(3,natom,trotter)
  real(dp),intent(in),target :: xred_prev(3,natom,trotter)
 end subroutine pimd_nosehoover_npt
end interface

interface
 subroutine pimd_nosehoover_nvt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&  
  &  rprimd,trotter,vel,volume,xred,xred_next,xred_prev)
  use m_pimd
  use defs_basis
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: prtvolimg
  integer,intent(in) :: trotter
  type(pimd_type),intent(in) :: pimd_param
  real(dp),intent(in) :: volume
  real(dp),intent(in) :: etotal(trotter)
  real(dp),intent(inout) :: forces(3,natom,trotter)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: vel(3,natom,trotter)
  real(dp),intent(in),target :: xred(3,natom,trotter)
  real(dp),intent(out) :: xred_next(3,natom,trotter)
  real(dp),intent(in),target :: xred_prev(3,natom,trotter)
 end subroutine pimd_nosehoover_nvt
end interface

interface
 subroutine prec_simple(ab_mover,forstr_precon,hist,icycle,itime,iexit)
  use defs_mover
  implicit none
  integer,intent(in) :: icycle
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_forstr) :: forstr_precon
  type(ab_movehistory),intent(inout) :: hist
 end subroutine prec_simple
end interface

interface
 subroutine pred_bfgs(ab_mover,ab_xfh,forstr_precon,hist,ionmov,itime,zDEBUG,iexit)
  use defs_mover
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: ionmov
  integer,intent(in) :: itime
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_xfh_type),intent(inout) :: ab_xfh
  type(ab_forstr),intent(in) :: forstr_precon
  type(ab_movehistory),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_bfgs
end interface

interface
 subroutine pred_delocint(ab_mover,ab_xfh,forstr_precon,hist,ionmov,itime,zDEBUG,iexit)
  use defs_mover
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: ionmov
  integer,intent(in) :: itime
  type(ab_movetype),intent(inout) :: ab_mover
  type(ab_xfh_type),intent(inout) :: ab_xfh
  type(ab_forstr),intent(in) :: forstr_precon
  type(ab_movehistory),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_delocint
end interface

interface
 subroutine pred_diisrelax(ab_mover,hist,itime,ntime,zDEBUG,iexit)
  use defs_mover
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(in) :: ntime
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_movehistory),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_diisrelax
end interface

interface
 subroutine pred_isokinetic(ab_mover,hist,itime,ntime,zDEBUG,iexit)
  use defs_mover
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(in) :: ntime
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_movehistory),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_isokinetic
end interface

interface
 subroutine pred_isothermal(ab_mover,hist,itime,mttk_vars,ntime,zDEBUG,iexit)
  use defs_mover
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(in) :: ntime
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_movehistory),intent(inout) :: hist
  type(mttk_type),intent(inout) :: mttk_vars
  logical,intent(in) :: zDEBUG
 end subroutine pred_isothermal
end interface

interface
 subroutine pred_langevin(ab_mover,hist,icycle,itime,ncycle,ntime,zDEBUG,iexit)
  use defs_mover
  implicit none
  integer,intent(inout) :: icycle
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(inout) :: ncycle
  integer,intent(in) :: ntime
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_movehistory),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_langevin
end interface

interface
 subroutine pred_moldyn(ab_mover,hist,icycle,itime,ncycle,ntime,zDEBUG,iexit)
  use defs_mover
  implicit none
  integer,intent(in) :: icycle
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(inout) :: ncycle
  integer,intent(in) :: ntime
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_movehistory),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_moldyn
end interface

interface
 function fdtion(ab_mover,hist,itime)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: itime
  type(ab_movetype),intent(in) :: ab_mover
  real(dp) :: fdtion
  type(ab_movehistory),intent(in) :: hist
 end function fdtion
end interface

interface
 subroutine pred_nose(ab_mover,hist,itime,ntime,zDEBUG,iexit)
  use defs_mover
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  integer,intent(in) :: ntime
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_movehistory),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_nose
end interface

interface
 subroutine pred_simple(ab_mover,hist,iexit)
  use defs_mover
  implicit none
  integer,intent(in) :: iexit
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_movehistory),intent(inout) :: hist
 end subroutine pred_simple
end interface

interface
 subroutine pred_srkna14(ab_mover,hist,icycle,ncycle,zDEBUG,iexit)
  use defs_mover
  implicit none
  integer,intent(inout) :: icycle
  integer,intent(in) :: iexit
  integer,intent(in) :: ncycle
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_movehistory),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_srkna14
end interface

interface
 subroutine pred_steepdesc(ab_mover,forstr_precon,hist,itime,zDEBUG,iexit)
  use defs_mover
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: itime
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_forstr),intent(in) :: forstr_precon
  type(ab_movehistory),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_steepdesc
end interface

interface
 subroutine pred_verlet(ab_mover,hist,ionmov,itime,ntime,zDEBUG,iexit)
  use defs_mover
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: ionmov
  integer,intent(in) :: itime
  integer,intent(in) :: ntime
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_movehistory),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
 end subroutine pred_verlet
end interface

interface
 subroutine predict_copy(itimimage,list_dynimage,ndynimage,nimage,ntimimage,results_img)
  use m_results_img
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: ndynimage
  integer,intent(in) :: nimage
  integer,intent(in) :: ntimimage
  integer,intent(in) :: list_dynimage(ndynimage)
  type(results_img_type) :: results_img(nimage,ntimimage)
 end subroutine predict_copy
end interface

interface
 subroutine predict_pimd(imgmov,itimimage,mpi_enreg,natom,nimage,nimage_tot,&  
  &  ntimimage,pimd_param,prtvolimg,results_img)
  use m_pimd
  use m_results_img
  use defs_abitypes
  implicit none
  integer,intent(in) :: imgmov
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: nimage
  integer,intent(in) :: nimage_tot
  integer,intent(in) :: ntimimage
  integer,intent(in) :: prtvolimg
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pimd_type),intent(in) :: pimd_param
  type(results_img_type) :: results_img(nimage,ntimimage)
 end subroutine predict_pimd
end interface

interface
 subroutine predict_steepest(fxcartfactor,itimimage,list_dynimage,natom,ndynimage,nimage,&  
  &  ntimimage,results_img)
  use defs_basis
  use m_results_img
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: ndynimage
  integer,intent(in) :: nimage
  integer,intent(in) :: ntimimage
  real(dp),intent(in) :: fxcartfactor
  integer,intent(in) :: list_dynimage(ndynimage)
  type(results_img_type) :: results_img(nimage,ntimimage)
 end subroutine predict_steepest
end interface

interface
 subroutine predict_string(fxcartfactor,iatfix,itimimage,list_dynimage,mpi_enreg,natom,&  
  &  ndynimage,nimage,nimage_tot,ntimimage,results_img)
  use m_results_img
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: ndynimage
  integer,intent(in) :: nimage
  integer,intent(in) :: nimage_tot
  integer,intent(in) :: ntimimage
  real(dp),intent(in) :: fxcartfactor
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: iatfix(3,natom)
  integer,intent(in) :: list_dynimage(ndynimage)
  type(results_img_type) :: results_img(nimage,ntimimage)
 end subroutine predict_string
end interface

interface
 subroutine predictimg(deltae,fxcartfactor,iatfix,imagealgo_str,imgmov,itimimage,list_dynimage,&  
  &  mpi_enreg,natom,ndynimage,nimage,nimage_tot,ntimimage,&  
  &  pimd_param,prtvolimg,results_img)
  use m_pimd
  use defs_basis
  use defs_abitypes
  use m_results_img
  implicit none
  integer,intent(in) :: imgmov
  integer,intent(in) :: itimimage
  integer,intent(in) :: natom
  integer,intent(in) :: ndynimage
  integer,intent(in) :: nimage
  integer,intent(in) :: nimage_tot
  integer,intent(in) :: ntimimage
  integer,intent(in) :: prtvolimg
  real(dp),intent(in) :: deltae
  real(dp),intent(in) :: fxcartfactor
  character(len=60),intent(in) :: imagealgo_str
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pimd_type),intent(in) :: pimd_param
  integer,intent(in) :: iatfix(3,natom)
  integer,intent(in) :: list_dynimage(ndynimage)
  type(results_img_type) :: results_img(nimage,ntimimage)
 end subroutine predictimg
end interface

interface
 subroutine print_bonds(amu,bonds,natom,ntypat,symbol,typat,znucl)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  type(go_bonds),intent(in) :: bonds
  real(dp) :: amu(ntypat)
  character(len=2) :: symbol(ntypat)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine print_bonds
end interface

interface
 subroutine print_phonfreq(istep,natom_primitive_cell,nphononq,phonon_eigval)
  use defs_basis
  implicit none
  integer,intent(in) :: istep
  integer,intent(in) :: natom_primitive_cell
  integer,intent(in) :: nphononq
  real(dp),intent(in) :: phonon_eigval(3*natom_primitive_cell,nphononq)
 end subroutine print_phonfreq
end interface

interface
 subroutine prtxfase(ab_mover,hist,iout,pos)
  use defs_mover
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: pos
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_movehistory),intent(in) :: hist
 end subroutine prtxfase
end interface

interface
 subroutine gettag(atlist,index,natom,prtallatoms,tag)
  implicit none
  integer,intent(in) :: index
  integer,intent(in) :: natom
  logical,intent(in) :: prtallatoms
  character(len=7),intent(out) :: tag
  logical,intent(in) :: atlist(natom)
 end subroutine gettag
end interface

interface
 subroutine prtnatom(atlist,iout,message,natom,prtallatoms,thearray)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  character(len=80*(max(natom,3)+1)) :: message
  logical,intent(in) :: prtallatoms
  logical,intent(in) :: atlist(natom)
  real(dp) :: thearray(3,natom)
 end subroutine prtnatom
end interface

interface
 subroutine strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rprimd0(3,3)
  real(dp),intent(out) :: rprimd_symm(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine strainsym
end interface

interface
 subroutine var2hist(acell,hist,natom,rprim,rprimd,xcart,xred,zDEBUG)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  type(ab_movehistory),intent(inout) :: hist
  logical,intent(in) :: zDEBUG
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine var2hist
end interface

interface
 subroutine vel2hist(amass,hist,natom,vel)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  type(ab_movehistory),intent(inout) :: hist
  real(dp),intent(in) :: amass(natom)
  real(dp),intent(in) :: vel(3,natom)
 end subroutine vel2hist
end interface

interface
 subroutine xcart2deloc(deloc,natom,rprimd,xcart,&  
  &  bt_inv_matrix,u_matrix,deloc_int,prim_int)
  use defs_basis
  use m_delocint
  implicit none
  integer,intent(in) :: natom
  type(ab_delocint),intent(inout) :: deloc
  real(dp),intent(out) :: bt_inv_matrix(3*(natom-1),3*natom)
  real(dp),intent(out) :: deloc_int(3*(natom-1))
  real(dp),intent(out) :: prim_int(deloc%ninternal)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: u_matrix(deloc%ninternal,3*(natom-1))
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine xcart2deloc
end interface

interface
 subroutine calc_btinv_matrix(b_matrix,natom,&  
  &  ninternal,bt_inv_matrix,u_matrix)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ninternal
  real(dp),intent(in) :: b_matrix(ninternal,3*natom)
  real(dp),intent(out) :: bt_inv_matrix(3*(natom-1),3*natom)
  real(dp),intent(inout) :: u_matrix(ninternal,3*(natom-1))
 end subroutine calc_btinv_matrix
end interface

interface
 subroutine align_u_matrices(natom,ninternal,u_matrix,u_matrix_old,s_matrix,f_eigs)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ninternal
  real(dp),intent(inout) :: f_eigs(3*natom)
  real(dp),intent(inout) :: s_matrix(3*natom,3*natom)
  real(dp),intent(inout) :: u_matrix(ninternal,3*(natom-1))
  real(dp),intent(in) :: u_matrix_old(ninternal,3*(natom-1))
 end subroutine align_u_matrices
end interface

interface
 subroutine xfh_recover_deloc(ab_xfh,ab_mover,acell,acell0,cycl_main,&  
  &  fred,hessin,ndim,rprim,rprimd0,strten,ucvol,ucvol0,vin,vin_prev,&  
  &  vout,vout_prev,xred,deloc,deloc_int,deloc_force,btinv,gprimd,prim_int,&  
  &  u_matrix)
  use defs_mover
  use defs_basis
  use m_delocint
  implicit none
  integer,intent(out) :: cycl_main
  integer,intent(in) :: ndim
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_xfh_type),intent(inout) :: ab_xfh
  type(ab_delocint),intent(inout) :: deloc
  real(dp),intent(inout) :: ucvol
  real(dp),intent(inout) :: ucvol0
  real(dp),intent(inout) :: acell(3)
  real(dp),intent(in) :: acell0(3)
  real(dp),intent(inout) :: btinv(3*(ab_mover%natom-1),3*ab_mover%natom)
  real(dp),intent(inout) :: deloc_force(3*(ab_mover%natom-1))
  real(dp),intent(inout) :: deloc_int(3*(ab_mover%natom-1))
  real(dp),intent(inout) :: fred(3,ab_mover%natom)
  real(dp),intent(inout) :: gprimd(3,3)
  real(dp),intent(inout) :: hessin(:,:)
  real(dp),intent(inout) :: prim_int(:)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(inout) :: rprimd0(3,3)
  real(dp),intent(inout) :: strten(6)
  real(dp),intent(inout) :: u_matrix(:,:)
  real(dp),intent(inout) :: vin(:)
  real(dp),intent(inout) :: vin_prev(:)
  real(dp),intent(inout) :: vout(:)
  real(dp),intent(inout) :: vout_prev(:)
  real(dp),intent(inout) :: xred(3,ab_mover%natom)
 end subroutine xfh_recover_deloc
end interface

interface
 subroutine xfh_recover_new(ab_xfh,ab_mover,acell,acell0,cycl_main,fred,&  
  &  hessin,ndim,rprim,rprimd0,strten,ucvol,ucvol0,vin,vin_prev,vout,&  
  &  vout_prev,xred)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(out) :: cycl_main
  integer,intent(in) :: ndim
  type(ab_movetype),intent(in) :: ab_mover
  type(ab_xfh_type),intent(inout) :: ab_xfh
  real(dp),intent(inout) :: ucvol
  real(dp),intent(inout) :: ucvol0
  real(dp),intent(inout) :: acell(3)
  real(dp),intent(in) :: acell0(3)
  real(dp),intent(inout) :: fred(3,ab_mover%natom)
  real(dp),intent(inout) :: hessin(:,:)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(inout) :: rprimd0(3,3)
  real(dp),intent(inout) :: strten(6)
  real(dp),intent(inout) :: vin(:)
  real(dp),intent(inout) :: vin_prev(:)
  real(dp),intent(inout) :: vout(:)
  real(dp),intent(inout) :: vout_prev(:)
  real(dp),intent(inout) :: xred(3,ab_mover%natom)
 end subroutine xfh_recover_new
end interface

interface
 subroutine xfh_update(ab_xfh,acell,fred_corrected,natom,rprim,strten,xred)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  type(ab_xfh_type),intent(inout) :: ab_xfh
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: fred_corrected(3,natom)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine xfh_update
end interface

interface
 subroutine xfpack_f2vout(fred,natom,ndim,optcell,&  
  &  strtarget,strten,ucvol,vout)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndim
  integer,intent(in) :: optcell
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(out) :: vout(ndim)
 end subroutine xfpack_f2vout
end interface

interface
 subroutine xfpack_vin2x(acell,acell0,natom,ndim,nsym,optcell,&  
  &  rprim,rprimd0,symrel,ucvol,ucvol0,vin,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndim
  integer,intent(in) :: nsym
  integer,intent(in) :: optcell
  real(dp),intent(out) :: ucvol
  real(dp),intent(in) :: ucvol0
  real(dp),intent(out) :: acell(3)
  real(dp),intent(in) :: acell0(3)
  real(dp),intent(out) :: rprim(3,3)
  real(dp),intent(in) :: rprimd0(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: vin(ndim)
  real(dp),intent(out) :: xred(3,natom)
 end subroutine xfpack_vin2x
end interface

interface
 subroutine xfpack_x2vin(acell,acell0,natom,ndim,nsym,optcell,&  
  &  rprim,rprimd0,symrel,ucvol,ucvol0,vin,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndim
  integer,intent(in) :: nsym
  integer,intent(in) :: optcell
  real(dp),intent(out) :: ucvol
  real(dp),intent(in) :: ucvol0
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: acell0(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rprimd0(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(out) :: vin(ndim)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine xfpack_x2vin
end interface

end module interfaces_45_geomoptim
!!***
