!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_kptrank
!! NAME
!! m_kptrank
!!
!! FUNCTION
!! This module deals with rank objects for hashing k-point vector lists
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_kptrank

 use m_profiling

 use defs_basis
 use m_errors

 implicit none

 public
!!***

!!****t* m_kptrank/kptrank_type
!! NAME
!! kptrank_type
!! 
!! FUNCTION
!!  structure to contain a rank/inverse rank pair of arrays, with dimensions
!! 
!! SOURCE

 type,public :: kptrank_type
   integer :: max_linear_density
   integer :: max_rank
   integer :: npoints
   integer, pointer :: invrank(:)
   integer, pointer :: rank(:)
   integer, pointer :: multipl(:)
 end type kptrank_type

contains
!!***

!!****f* m_kptrank/mkkptrank
!!
!! NAME
!! mkkptrank
!!
!! FUNCTION
!! This routine sets up the kpt ranks for comparing kpts
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npt = number of kpoints
!!  kpt = coordinates of kpoints
!!
!! OUTPUT
!!  kptrank_t = object containing ranking and inverse ranking
!!
!! NOTES
!!
!! PARENTS
!!      get_full_kgrid,m_gamma,m_tetrahedron,mkfskgrid,mknesting,mkqptequiv
!!      order_fs_kpts,outelph,read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkkptrank (kpt,nkpt,kptrank_t, nsym, symrec)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkkptrank'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 integer,intent(in), optional :: nsym
!arrays
 type(kptrank_type), intent(out) :: kptrank_t
 real(dp),intent(in) :: kpt(3,nkpt)
 integer,intent(in), optional :: symrec(3,3, *)

!Local variables -------------------------
!scalars
 integer :: ikpt,istat, isym, symkptrank, irank
 real(dp) :: smallestlen
 character(len=500) :: msg
!arrays
 real(dp) :: symkpt(3)

! *********************************************************************

! find smallest linear length
 smallestlen = one
 do ikpt=1, nkpt
   if (abs(kpt(1,ikpt)) > tol10) &
&     smallestlen = min(smallestlen, abs(kpt(1,ikpt)))
   if (abs(kpt(2,ikpt)) > tol10) &
&     smallestlen = min(smallestlen, abs(kpt(2,ikpt)))
   if (abs(kpt(3,ikpt)) > tol10) &
&     smallestlen = min(smallestlen, abs(kpt(3,ikpt)))
 end do

 kptrank_t%max_linear_density = int(one/smallestlen)+1
 kptrank_t%max_rank = 2*kptrank_t%max_linear_density**3
 kptrank_t%npoints = nkpt

 ABI_ALLOCATE(kptrank_t%rank,(nkpt))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory %rank")

 ABI_ALLOCATE(kptrank_t%invrank,(kptrank_t%max_rank))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory %invrank")
 kptrank_t%invrank(:) = -1

!Ensure kpt(i)+one is positive, and the smallest
!difference between kpts should be larger than 1/100
!ie ngkpt < 100.
 do ikpt=1,nkpt
   call get_rank_1kpt (kpt(:,ikpt), kptrank_t%rank(ikpt), kptrank_t)

   if (kptrank_t%rank(ikpt) > kptrank_t%max_rank .or. kptrank_t%rank(ikpt) < 1) then
     write(msg,'(a,2i0)')" max rank exceeded or < 1, ikpt, rank ", ikpt, kptrank_t%rank(ikpt)
     MSG_ERROR(msg)
   end if
   kptrank_t%invrank(kptrank_t%rank(ikpt)) = ikpt
 end do
 
! if symrec is provided, fill invrank with appropriate irred kpt indices
! for symmetry completion: kptrank_t%invrank points to the irred k-point
! equivalent to the k-point whose rank is provided
 if (present(symrec)) then
   ABI_CHECK(present(nsym), "need both symrec and nsym arguments together")
   do ikpt=1,nkpt
     do isym = 1, nsym
       symkpt = matmul(symrec(:,:,isym), kpt(:, ikpt))
       
       call get_rank_1kpt (symkpt(:), symkptrank, kptrank_t)

       kptrank_t%invrank(symkptrank) = ikpt
     end do
   end do
 end if 

 ABI_ALLOCATE(kptrank_t%multipl,(nkpt))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory %multipl")

 kptrank_t%multipl = 0
! find multiplicity of ikpt
 do ikpt = 1, nkpt
   do irank = 1, kptrank_t%max_rank
     if (kptrank_t%invrank(irank)==ikpt) kptrank_t%multipl(ikpt) = kptrank_t%multipl(ikpt) + 1
   end do
 end do

end subroutine mkkptrank
!!***

!!****f* m_kptrank/get_rank_1kpt
!!
!! NAME
!! get_rank_1kpt
!!
!! FUNCTION
!! This routine calculates the rank for one kpt
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kpt = coordinates of kpoints
!!  kptrank_t = rank object for the k-grid we are using
!!
!! OUTPUT
!!  rank = rank of the kpoint
!!
!! NOTES
!!
!! PARENTS
!!      bfactor,elphon,get_full_kgrid,integrate_gamma,integrate_gamma_alt
!!      k_neighbors,m_gamma,m_kptrank,m_tetrahedron,mkfskgrid,mkqptequiv
!!      read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_rank_1kpt (kpt,rank,kptrank_t)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_rank_1kpt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: rank
 type(kptrank_type), intent(in) :: kptrank_t
!arrays
 real(dp),intent(in) :: kpt(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg
!arrays    
 real(dp) :: redkpt(3)

! *************************************************************************

! wrap to [0, 1[ -> replaced call to wrap2_zeroone inline, to encapsulate this module
 if (kpt(1)>zero) then
   redkpt(1)=mod((kpt(1)+tol12),one)-tol12
 else
   redkpt(1)=-mod(-(kpt(1)-one+tol12),one)+one-tol12
 end if
 if(abs(redkpt(1))<tol12)redkpt(1)=0.0_dp

 if (kpt(2)>zero) then
   redkpt(2)=mod((kpt(2)+tol12),one)-tol12
 else
   redkpt(2)=-mod(-(kpt(2)-one+tol12),one)+one-tol12
 end if
 if(abs(redkpt(2))<tol12)redkpt(2)=0.0_dp

 if (kpt(3)>zero) then
   redkpt(3)=mod((kpt(3)+tol12),one)-tol12
 else
   redkpt(3)=-mod(-(kpt(3)-one+tol12),one)+one-tol12
 end if
 if(abs(redkpt(3))<tol12)redkpt(3)=0.0_dp



 rank = int(real(kptrank_t%max_linear_density)*(redkpt(3)+half+tol8 +&
&           real(kptrank_t%max_linear_density)*(redkpt(2)+half+tol8 +&
&           real(kptrank_t%max_linear_density)*(redkpt(1)+half+tol8))))

 if (rank > kptrank_t%max_rank) then
   write(msg,'(a,i0)') ' rank should be inferior to ', kptrank_t%max_rank
   MSG_ERROR(msg)
 end if

end subroutine get_rank_1kpt
!!***

!----------------------------------------------------------------------

!!****f* m_kptrank/copy_kptrank
!!
!! NAME
!! copy_kptrank
!!
!! FUNCTION
!! This routine deallocates the arrays in a kptrank_type structure
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  kptrank_t = object containing ranking and inverse ranking, to be deallocated
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_kptrank (kptrank_t_in, kptrank_t_out)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_kptrank'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(kptrank_type), intent(in) :: kptrank_t_in
 type(kptrank_type), intent(out) :: kptrank_t_out

!Local variables -------------------------

! *********************************************************************
 kptrank_t_out%max_linear_density = kptrank_t_in%max_linear_density
 kptrank_t_out%max_rank = kptrank_t_in%max_rank
 kptrank_t_out%npoints = kptrank_t_in%npoints
 
 ABI_ALLOCATE(kptrank_t_out%rank,(kptrank_t_out%npoints))
 kptrank_t_out%rank = kptrank_t_in%rank
 
 ABI_ALLOCATE(kptrank_t_out%invrank,(kptrank_t_out%max_rank))
 kptrank_t_out%invrank = kptrank_t_in%invrank

 ABI_ALLOCATE(kptrank_t_out%multipl,(kptrank_t_out%npoints))
 kptrank_t_out%multipl = kptrank_t_in%multipl
 
end subroutine copy_kptrank
!!***

!----------------------------------------------------------------------

!!****f* m_kptrank/destroy_kptrank
!!
!! NAME
!! destroy_kptrank
!!
!! FUNCTION
!! This routine deallocates the arrays in a kptrank_type structure
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  kptrank_t = object containing ranking and inverse ranking, to be deallocated
!!
!! NOTES
!!
!! PARENTS
!!      defs_elphon,get_full_kgrid,m_gamma,m_tetrahedron,mkfskgrid,mknesting
!!      mkqptequiv,order_fs_kpts,outelph,read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_kptrank (kptrank_t)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_kptrank'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(kptrank_type), intent(out) :: kptrank_t

! *********************************************************************

 if (associated(kptrank_t%rank))  then
   ABI_DEALLOCATE(kptrank_t%rank)
 end if
 if (associated(kptrank_t%invrank))  then
   ABI_DEALLOCATE(kptrank_t%invrank)
 end if
 if (associated(kptrank_t%multipl))  then
   ABI_DEALLOCATE(kptrank_t%multipl)
 end if

end subroutine destroy_kptrank
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/nullify_kptrank
!!
!! NAME
!! nullify_kptrank
!!
!! FUNCTION
!! This routine nullifies the arrays in a kptrank_type structure
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  kptrank_t = object containing ranking and inverse ranking, to be nullified
!!
!! NOTES
!!
!! PARENTS
!!      defs_elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_kptrank (kptrank_t)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_kptrank'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(kptrank_type), intent(out) :: kptrank_t

!Local variables -------------------------

! *********************************************************************

 nullify(kptrank_t%invrank)
 nullify(kptrank_t%rank)
 nullify(kptrank_t%multipl)

end subroutine nullify_kptrank
!!***

end module m_kptrank
!!***
