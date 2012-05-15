!!****m* ABINIT/interfaces_61_ionetcdf
!! NAME
!! interfaces_61_ionetcdf
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/61_ionetcdf
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

module interfaces_61_ionetcdf

 implicit none

interface
 subroutine ab_define_var(ncid, var_dim_id, var_id, var_type, var_name, var_mnemo, var_units)
  implicit none
  integer, intent(in) :: ncid
  integer, intent(out) :: var_id
  integer,intent(in) :: var_type
  character(len=*), intent(in) :: var_mnemo
  character(len=*), intent(in) :: var_name
  character(len=*), intent(in) :: var_units
  integer,intent(in) :: var_dim_id(:)
 end subroutine ab_define_var
end interface

interface
 subroutine abi_etsf_electrons_put(dtset, filapp)
  use defs_basis
  use defs_abitypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  character(len=fnlen),intent(in) :: filapp
 end subroutine abi_etsf_electrons_put
end interface

interface
 subroutine abi_etsf_geo_put(dtset, filapp, psps)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  character(len=fnlen),intent(in) :: filapp
  type(pseudopotential_type),intent(in) :: psps
 end subroutine abi_etsf_geo_put
end interface

interface
 subroutine abi_etsf_init(dtset, filapp, itype, kdep, lmn_size, psps, wfs)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: itype
  type(dataset_type),intent(in) :: dtset
  character(len=fnlen),intent(in) :: filapp
  logical,intent(in) :: kdep
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_wf_type),intent(in) :: wfs
  integer,intent(in) :: lmn_size(psps%npsp)
 end subroutine abi_etsf_init
end interface

interface
 subroutine ini_wf_etsf(dtset, lmn_size, npsp, ntypat, unwff)
  use defs_abitypes
  implicit none
  integer,intent(in) :: npsp
  integer,intent(in) :: ntypat
  integer,intent(in) :: unwff
  type(dataset_type),intent(in) :: dtset
  integer,intent(in) :: lmn_size(npsp)
 end subroutine ini_wf_etsf
end interface

interface
 subroutine ini_wf_netcdf(mpw,ncid_hdr,response)
  implicit none
  integer,intent(in) :: mpw
  integer,intent(in) :: ncid_hdr
  integer,intent(in) :: response
 end subroutine ini_wf_netcdf
end interface

interface
 subroutine read_md_hist(filename,hist)
  use defs_mover
  use defs_basis
  implicit none
  character(len=fnlen),intent(in) :: filename
  type(ab_movehistory),intent(out) :: hist
 end subroutine read_md_hist
end interface

interface
 subroutine write_eig(eigen,filename,kptns,mband,nband,nkpt,nsppol)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  character(len=fnlen),intent(in) :: filename
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine write_eig
end interface

interface
 subroutine write_md_hist(filename,hist,icycle,itime,natom)
  use defs_mover
  use defs_basis
  implicit none
  integer,intent(in) :: icycle
  integer,intent(in) :: itime
  integer,intent(in) :: natom
  character(len=fnlen),intent(in) :: filename
  type(ab_movehistory),intent(in) :: hist
 end subroutine write_md_hist
end interface

interface
 subroutine wrt_moldyn_netcdf(amass,dtset,itime,option,moldyn_file,mpi_enreg,&  
  &  results_gs,rprimd,unpos,vel,xcart,xred)
  use defs_basis
  use defs_abitypes
  use m_results_gs
  implicit none
  integer,intent(in) :: itime
  integer,intent(in) :: option
  integer,intent(in) :: unpos
  type(dataset_type),intent(in) :: dtset
  character(fnlen),intent(in) :: moldyn_file
  type(mpi_type),intent(in) :: mpi_enreg
  type(results_gs_type),intent(in) :: results_gs
  real(dp),intent(in) :: amass(dtset%natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in),target :: vel(3,dtset%natom)
  real(dp),intent(in) :: xcart(3,dtset%natom)
  real(dp),intent(in) :: xred(3,dtset%natom)
 end subroutine wrt_moldyn_netcdf
end interface

end module interfaces_61_ionetcdf
!!***
