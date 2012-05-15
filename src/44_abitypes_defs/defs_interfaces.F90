!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_interfaces
!! NAME
!! defs_interfaces
!!
!! FUNCTION
!! This module contains interfaces for using several procedures with the
!! same generic name, with optional arguments or with assumed-shape arrays
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (XG, MVer, TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module defs_interfaces

 use m_profiling

 implicit none


 interface
  subroutine deloc2xcart(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,nrshift,&
&  dtset,rprimd,rshift,xcart,&
&  deloc_int,btinv,u_matrix)
   use defs_basis
   use defs_abitypes
   integer,intent(in) :: nbond
   integer,intent(in) :: nang
   integer,intent(in) :: ndihed
   integer,intent(in) :: ncart
   integer,intent(in) :: ninternal
   integer,intent(in) :: nrshift
   type(dataset_type),intent(in) :: dtset
   integer,pointer :: angs(:,:,:)
   integer,pointer :: bonds(:,:,:)
   integer,pointer :: carts(:,:)
   integer,pointer :: dihedrals(:,:,:)
   real(dp),intent(in) :: rprimd(3,3)
   real(dp),intent(in) :: rshift(3,nrshift)
   real(dp),intent(inout) :: xcart(3,dtset%natom)
   real(dp),intent(in) :: deloc_int(3*(dtset%natom-1))
   real(dp),intent(out) :: btinv(3*(dtset%natom-1),3*dtset%natom)
   real(dp),intent(inout) :: u_matrix(ninternal,3*(dtset%natom-1))
  end subroutine deloc2xcart
 end interface


!Generic interface of the routines hdr_skip
 interface hdr_skip
  subroutine hdr_skip_int(unitfi,ierr)
   integer, intent(in) :: unitfi
   integer, intent(out) :: ierr
  end subroutine hdr_skip_int
 end interface
!End of the generic interface of hdr_skip

!Generic interface of the routines hdr_io
 interface hdr_io

  subroutine hdr_io_int(fform,hdr,rdwr,unitfi)
   use defs_abitypes
   integer,intent(inout) :: fform
   integer,intent(in) :: rdwr
   integer,intent(in) :: unitfi
   type(hdr_type),intent(inout) :: hdr
  end subroutine hdr_io_int
 end interface
!End of the generic interface of hdr_io

!Generic interface of the routines hdr_io_netcdf
 interface hdr_io_netcdf

  subroutine hdr_io_netcdf_int(fform,hdr,rdwr,unitfi)
   use defs_abitypes
   integer,intent(inout) :: fform
   integer,intent(in) :: rdwr
   integer,intent(in) :: unitfi
   type(hdr_type),intent(inout) :: hdr
  end subroutine hdr_io_netcdf_int

end interface
!End of the generic interface of hdr_io


!Optional arguments

 !Subroutine belongs to 03ionetcdf
 interface
  subroutine netcdf_file_create(ncid,filename,blocks,dims)
   use defs_datatypes
   integer,intent(inout) :: ncid
   integer,intent(in) :: blocks
   character(len=*),intent(in) :: filename
   type(vardims_type),intent(in) :: dims
  end subroutine netcdf_file_create
 end interface

 interface
  subroutine make_prim_internals(angs,bonds,carts,dihedrals,icenter,nbond,nang,ndihed,ncart,&
&  dtset,nrshift,rprimd,rshift,xcart)
   use defs_basis
   use defs_abitypes
   integer,intent(in) :: icenter
   integer,intent(out) :: nbond
   integer,intent(out) :: nang
   integer,intent(out) :: ndihed
   integer,intent(out) :: ncart
   integer,intent(in) :: nrshift
   type(dataset_type),intent(in) :: dtset
   integer,pointer :: angs(:,:,:)
   integer,pointer :: bonds(:,:,:)
   integer,pointer :: carts(:,:)
   integer,pointer :: dihedrals(:,:,:)
   real(dp),intent(in) :: rprimd(3,3)
   real(dp),intent(in) :: rshift(3,nrshift)
   real(dp),intent(in) :: xcart(3,dtset%natom)
  end subroutine make_prim_internals
 end interface

 interface
  subroutine xcart2deloc(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
&  dtset,nrshift,rprimd,rshift,xcart,&
&  bt_inv_matrix,u_matrix,deloc_int,prim_int)
   use defs_basis
   use defs_abitypes
   integer,intent(in) :: nbond
   integer,intent(in) :: nang
   integer,intent(in) :: ndihed
   integer,intent(in) :: ncart
   integer,intent(in) :: ninternal
   integer,intent(in) :: nrshift
   type(dataset_type),intent(in) :: dtset
   integer,pointer :: angs(:,:,:)
   integer,pointer :: bonds(:,:,:)
   integer,pointer :: carts(:,:)
   integer,pointer :: dihedrals(:,:,:)
   real(dp),intent(in) :: rprimd(3,3)
   real(dp),intent(in) :: rshift(3,nrshift)
   real(dp),intent(in) :: xcart(3,dtset%natom)
   real(dp),intent(out) :: bt_inv_matrix(3*(dtset%natom-1),3*dtset%natom)
   real(dp),intent(inout) :: u_matrix(ninternal,3*(dtset%natom-1))
   real(dp),intent(out) :: deloc_int(3*(dtset%natom-1))
   real(dp),intent(out) :: prim_int(ninternal)
  end subroutine xcart2deloc
 end interface

 interface
  subroutine print_ij(a_ij,adim,cplex,ndim,opt_io,opt_l,opt_l_index,opt_pack,opt_prtvol,&
&                     pack2ij,test_value,unt, &
&                     opt_sym,asym_ij)    !Optional arguments
   use defs_basis
   integer,intent(in) :: adim,cplex,ndim,opt_io,opt_l,opt_pack,opt_prtvol,unt
   integer,intent(in),optional :: opt_sym
   real(dp),intent(in) :: test_value
   integer,intent(in) :: opt_l_index(ndim*min(1+opt_l,1)),pack2ij(adim*opt_pack)
   real(dp),intent(in) :: a_ij(cplex*adim)
   real(dp),intent(in),optional :: asym_ij(cplex*adim)
  end subroutine print_ij
 end interface

!Assumed-shape arrays

 interface
  subroutine cprj_alloc(cprj,ncpgr,nlmn)
   use defs_datatypes
   integer,intent(in) :: ncpgr
   type(cprj_type),intent(inout) :: cprj(:,:)
   integer,intent(in) :: nlmn(:)
  end subroutine cprj_alloc
 end interface

 interface
  subroutine cprj_free(cprj)
   use defs_datatypes
   type(cprj_type),intent(inout) :: cprj(:,:)
  end subroutine cprj_free
 end interface

 interface
  subroutine cprj_nullify(cprj)
   use defs_datatypes
   type(cprj_type),intent(inout) :: cprj(:,:)
  end subroutine cprj_nullify
 end interface

 interface
  subroutine cprj_copy(cprj_in,cprj_out)
   use defs_basis
   use defs_datatypes
   type(cprj_type),intent(in) :: cprj_in(:,:)
   type(cprj_type),intent(inout) :: cprj_out(:,:)
  end subroutine cprj_copy
 end interface

 interface
  subroutine cprj_axpby(alpha,beta,cprjx,cprjy)
   use defs_basis
   use defs_datatypes
   real(dp),intent(in) :: alpha,beta
   type(cprj_type),intent(in) :: cprjx(:,:)
   type(cprj_type),intent(inout) :: cprjy(:,:)
  end subroutine cprj_axpby
 end interface

 interface
  subroutine rhoij_alloc(nlmn,nspden,nsppol,pawrhoij,typat,&      ! Mandatory arguments
&                        ngrhoij,nlmnmix,use_rhoij_,use_rhoijres) ! Optional arguments
   use defs_basis
   use defs_datatypes
   integer,intent(in) :: nspden,nsppol
   integer,intent(in),optional :: ngrhoij,nlmnmix,use_rhoij_,use_rhoijres
   integer,intent(in) :: nlmn(:),typat(:)
   type(pawrhoij_type),intent(inout) :: pawrhoij(:)
  end subroutine rhoij_alloc
 end interface

 interface
  subroutine rhoij_free(pawrhoij)
   use defs_basis
   use defs_datatypes
   type(pawrhoij_type),intent(inout) :: pawrhoij(:)
  end subroutine rhoij_free
 end interface

 interface
  subroutine rhoij_copy(pawrhoij_in,pawrhoij_out,&
&                       keep_cplex,keep_nspden)
   use defs_basis
   use defs_datatypes
   logical,intent(in),optional :: keep_cplex,keep_nspden
   type(pawrhoij_type),intent(in) :: pawrhoij_in(:)
   type(pawrhoij_type),intent(inout) :: pawrhoij_out(:)
  end subroutine rhoij_copy
 end interface

 interface
  subroutine hdr_update(bantot,etot,fermie,hdr,natom,residm,rprimd,occ,pawrhoij,usepaw,xred)
   use defs_basis
   use defs_datatypes
   use defs_abitypes
   integer,intent(in) :: bantot,natom,usepaw
   real(dp),intent(in) :: etot,fermie,residm
   type(hdr_type),intent(out) :: hdr
   real(dp),intent(in) :: occ(bantot),rprimd(3,3),xred(3,natom)
   type(pawrhoij_type),intent(in) :: pawrhoij(:)
  end subroutine hdr_update
 end interface

end module defs_interfaces
!!***
