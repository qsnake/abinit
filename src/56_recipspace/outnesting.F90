!!****f* ABINIT/outnesting
!! NAME
!! outnesting
!!
!! FUNCTION
!!  Write ou the nesting factors calculated in mknesting
!!  Data on file in the X-Y format (prtnest 1) or
!!  in the XCrysden format (XSF)   (prtnest 2)
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  base_name = prefix of the output file
!!  gmet = metric in reciprocal space
!!  gprimd(3,3) dimensional reciprocal lattice vectors
!!  kptrlatt(3,3) basis vectors for k-grid
!!  nestordered = nesting function on full grid, points ordered in x, then y, then z
!!  nkpt = number of k points
!!  nqpath = number of points requested along the trajectory
!!  prtnest = flags governing the format of the output file
!!  qpath_vertices = vertices of the reciprocal space trajectory
!!
!! OUTPUT
!!  only write to file
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      mknesting
!!
!! CHILDREN
!!      interpol3d,make_path,printxsf,wrap2_zero_one
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outnesting(base_name,gmet,gprimd,kptrlatt,nestordered,nkpt,nqpath,prtnest,qpath_vertices)

 use m_profiling

 use defs_basis
 use m_io_tools
 use m_bz_mesh
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outnesting'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nqpath,prtnest,nkpt
 character(len=fnlen),intent(in) :: base_name
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: qpath_vertices(3,nqpath)
 real(dp),intent(in) :: nestordered(nkpt)


!local
 integer :: unit_nest, nkx,nky,nkz
 integer :: iost, indx, ii, ipoint
 integer :: npt_tot
 integer :: realrecip

 character(len=fnlen) :: fname
 character(len=500) :: message
 integer :: ndiv(nqpath-1)
 real(dp), pointer :: finepath(:,:)
 real(dp) :: tmpkpt(3)
 real(dp) :: origin(3),qpt(3), res, kval
 
! dummy variables for call to printxsf
 integer :: natom, ntypat, typat(1)
 real(dp) :: xcart (3,1), znucl(1)

! *************************************************************************

!===================================================================
!Definition of the q path along which ph linwid will be interpolated
!===================================================================
 nullify(finepath)
 call make_path(nqpath,qpath_vertices,gmet,'G',20,ndiv,npt_tot,finepath)

 nkx=kptrlatt(1,1)
 nky=kptrlatt(2,2)
 nkz=kptrlatt(3,3)

 if (nkpt /= nkx*nky*nkz) then
   write(message,'(a,9(i0,1x),2x,i0)')' Wrong input value for kptrlatt  ',kptrlatt, nkpt
   MSG_BUG(message)
 end if

!open output file and write header
 unit_nest=get_unit()
 fname=trim(base_name)

 open (unit=unit_nest,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
   MSG_ERROR(" opening file "//trim(fname))
 end if

 write (unit_nest,'(a)')'#'
 write (unit_nest,'(a)')'# ABINIT package : Nesting factor file'
 write (unit_nest,'(a)')'#'
 write (unit_nest,'(a,i10,a)')'# Nesting factor calculated on ',npt_tot,' Q-points'
 write (unit_nest,'(a)')'# Description of the Q-path :'
 write (unit_nest,'(a,i10)')'# Number of line segments = ',nqpath-1
 write (unit_nest,'(a)')'# Vertices of the Q-path and corresponding index = '
 indx=1
 do ii=1,nqpath
   write (unit_nest,'(a,3(E16.6,1x),i8)')'#  ',qpath_vertices(:,ii),indx
   if(ii<nqpath) indx=indx+ndiv(ii)
 end do
 write (unit_nest,'(a)')'#'

!Get qpoint along the q-path from finepath and interpolate the nesting factor
 indx=1

 do ipoint=1, npt_tot
   qpt(:) = finepath(:,ipoint)
   call wrap2_zero_one(qpt(1),tmpkpt(1),res)
   call wrap2_zero_one(qpt(2),tmpkpt(2),res)
   call wrap2_zero_one(qpt(3),tmpkpt(3),res)

   call interpol3d(tmpkpt,nkx,nky,nkz,kval,nestordered)

   write(unit_nest,'(i5,18e16.5)')indx,kval
   indx = indx+1
 end do !end ipoint do
 close (unit_nest)
 ABI_DEALLOCATE(finepath)

 if (prtnest==2) then !write also the nest factor in the XSF format
   fname=trim(base_name) // '_NEST_XSF'
   open (unit=unit_nest,file=fname,status='unknown',iostat=iost)
   if (iost /= 0) then
     MSG_ERROR("opening file "//trim(fname))
   end if

   origin(:)=zero
   realrecip=1 !reciprocal space
   natom = 1
   ntypat = 1
   typat = (/1/)
   xcart = reshape ((/zero, zero, zero/), (/3,1/)) 
   znucl = (/one/)
   call printxsf(nkx,nky,nkz,nestordered,gprimd,origin,natom, ntypat, typat, xcart, znucl, unit_nest,realrecip)

   close (unit_nest)
 end if

end subroutine outnesting
!!***
