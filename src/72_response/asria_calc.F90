!{\src2tex{textfont=tt}}
!!****f* ABINIT/asria_calc
!! NAME
!! asria_calc
!!
!! FUNCTION
!! Calculate the correction for the Acoustic sum rule on the InterAtomic Forces
!!  or on the dynamical matrix directly
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! asr=(0 => no ASR, 1 or 2=> the diagonal element is modified to give the ASR,
!!      5 => impose hermitian solution using lapack call)
!! d2cart=matrix of second derivatives of total energy, in cartesian coordinates
!! mpert =maximum number of ipert
!! natom=number of atom
!!
!! OUTPUT
!! d2asr=matrix used to store the correction needed to fulfill
!! the acoustic sum rule.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      anaddb,gath3
!!
!! CHILDREN
!!      wrtout,zgelss
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine asria_calc(asr,d2asr,d2cart,mpert,natom)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'asria_calc'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: asr,mpert,natom
!arrays
 real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
 real(dp),intent(out) :: d2asr(2,3,natom,3,natom)

!Local variables-------------------------------
!scalars
 integer :: idir1,idir2,ii,ipert1,ipert2
 character(len=500) :: message

 integer, allocatable :: packingindex(:,:,:,:)
 real(dp), allocatable :: constraints(:,:,:)
 real(dp), allocatable :: d2cart_packed(:,:)
 real(dp), allocatable :: singvals(:)
 real(dp), allocatable :: constr_rhs(:,:)
 real(dp), allocatable :: work(:,:),rwork(:)
 integer :: constrank, imatelem, iconst, nconst, nd2_packed, info

! *********************************************************************

 d2asr = zero

 if (asr==0) return

 write(message, '(a)' )&
& ' asria_calc : calculation of the correction to the ASR for the interatomic forces.'
 call wrtout(std_out,message,'COLL')
 do ipert1=1,natom
   do idir1=1,3
     do idir2=1,3

!      Compute d2asr
       do ipert2=1,natom
         d2asr(:,idir1,ipert1,idir2,ipert1)=&
&         d2asr(:,idir1,ipert1,idir2,ipert1)+&
&         d2cart(:,idir1,ipert1,idir2,ipert2)
       end do
     end do
   end do
 end do

!holistic method: overwrite d2asr with hermitian solution
 if (asr == 5) then
   nconst = 9*natom
   nd2_packed = 3*natom*(3*natom+1)/2
   ABI_ALLOCATE(constraints,(2,nconst, nd2_packed))
   ABI_ALLOCATE(d2cart_packed,(2,nd2_packed))
   ABI_ALLOCATE(constr_rhs,(2,nd2_packed))
   ABI_ALLOCATE(singvals,(nconst))
   ABI_ALLOCATE(work,(2,3*nd2_packed))
   ABI_ALLOCATE(rwork,(5*nd2_packed))
   ABI_ALLOCATE(packingindex,(3,natom,3,natom))
   ii=1
   packingindex=-1
   do ipert2=1,natom
     do idir2=1,3
       do ipert1=1,ipert2-1
         do idir1=1,3
           packingindex(idir1,ipert1,idir2,ipert2) = ii
           ii = ii+1
         end do
       end do
       do idir1=1,idir2
         packingindex(idir1,ipert2,idir2,ipert2) = ii
         ii = ii+1
       end do
     end do
   end do
!  setup constraint matrix
   constraints = zero
   do ipert1=1,natom
     do idir1=1,3
       do idir2=1,3
         iconst = idir2+3*(idir1-1 + 3*(ipert1-1))
!        set all atom forces, this component
         do ipert2=1,natom
           imatelem = packingindex(idir1,ipert1,idir2,ipert2)
           if (imatelem == -1) then
             imatelem = packingindex(idir2,ipert2,idir1,ipert1)
           end if
           constraints(1,iconst,imatelem) = one
         end do
       end do
     end do
   end do

   d2cart_packed = -999.0d0
   do ipert2=1,natom
     do idir2=1,3
       do ipert1=1,natom
         do idir1=1,3
           imatelem = packingindex(idir1,ipert1,idir2,ipert2)
           if (imatelem == -1) cycle
           d2cart_packed(:,imatelem) = d2cart(:,idir1,ipert1,idir2,ipert2)
         end do
       end do
     end do
   end do
   constr_rhs = zero
   constr_rhs(1,1:nconst) = matmul(constraints(1,:,:),d2cart_packed(1,:))
   constr_rhs(2,1:nconst) = matmul(constraints(1,:,:),d2cart_packed(2,:))

!  lwork = 3*nd2_packed
   call zgelss (nconst,nd2_packed,1,constraints,nconst,constr_rhs,nd2_packed,&
&   singvals,-one,constrank,work,3*nd2_packed,rwork,info)
   write(std_out,*) 'zgelss info ', info

!  unpack 
   do ipert2=1,natom
     do idir2=1,3
       do ipert1=1,natom
         do idir1=1,3
           imatelem = packingindex(idir1,ipert1,idir2,ipert2)
           if (imatelem == -1) then
             imatelem = packingindex(idir2,ipert2,idir1,ipert1)
!            NOTE: should complex conjugate the correction below.
           end if
           d2asr(:,idir1,ipert1,idir2,ipert2) = constr_rhs(:,imatelem)
         end do
       end do
     end do
   end do

   ABI_DEALLOCATE(constraints)
   ABI_DEALLOCATE(d2cart_packed)
   ABI_DEALLOCATE(singvals)
   ABI_DEALLOCATE(constr_rhs)
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(rwork)
   ABI_DEALLOCATE(packingindex)

 end if

end subroutine asria_calc
!!***
