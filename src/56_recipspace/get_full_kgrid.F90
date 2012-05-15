!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_full_kgrid
!! NAME
!! get_full_kgrid
!!
!! FUNCTION
!! create full grid of kpoints and find equivalent
!! irred ones. Duplicates work in getkgrid, but need all outputs
!! of klatt,kpt_fullbz, and indkpt
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (MVer,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  klatt(3,3)=reciprocal of lattice vectors for full kpoint grid
!!  kpt(3,nkpt)=irreducible kpoints
!!  kptrlatt(3,3)=lattice vectors for full kpoint grid
!!  nkpt=number of irreducible kpoints
!!  nkpt_fullbz=number of kpoints in full brillouin zone
!!  nshiftk=number of kpoint grid shifts
!!  nsym=number of symmetries
!!  shiftk(3,nshiftk)=kpoint shifts
!!  symrel(3,3,nsym)=symmetry matrices in real space
!!
!! OUTPUT
!!  indkpt(nkpt_fullbz)=non-symmetrized indices of the k-points (see symkpt.f)
!!  kpt_fullbz(3,nkpt_fullbz)=kpoints in full brillouin zone
!!
!! NOTES
!!  MG: The present inplementation always assumes kptopt==1 !!!!
!!
!! PARENTS
!!      m_bz_mesh,m_phdos,tetrahedron
!!
!! CHILDREN
!!      destroy_kptrank,get_rank_1kpt,mati3inv,mkkptrank,wrap2_pmhalf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine get_full_kgrid(indkpt,klatt,kpt,kpt_fullbz,kptrlatt,nkpt,&
& nkpt_fullbz,nshiftk,nsym,shiftk,symrel)

 use m_profiling

 use defs_basis
 use m_kptrank
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_full_kgrid'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nkpt_fullbz,nshiftk,nsym
!arrays
 integer,intent(in) :: kptrlatt(3,3),symrel(3,3,nsym)
 integer,intent(out) :: indkpt(nkpt_fullbz)
 real(dp),intent(in) :: klatt(3,3),kpt(3,nkpt),shiftk(3,nshiftk)
 real(dp),intent(out) :: kpt_fullbz(3,nkpt_fullbz)

!Local variables-------------------------------
!scalars
 integer :: ii,ikpt,ikshft,isym,itim,jj,kk,nn,timrev
 integer :: symrankkpt
 character(len=500) :: message
 type(kptrank_type) :: kptrank_t

!arrays
 integer :: boundmax(3),boundmin(3),inv_symrel(3,3,nsym)
 real(dp) :: k1(3),k2(3),shift(3)

! *********************************************************************

!Invert symrels => gives symrels for kpoints

 do isym=1,nsym
   call mati3inv (symrel(:,:,isym),inv_symrel(:,:,isym))
 end do

!Generate full grid with 1 shift
!Now, klatt contains the three primitive vectors of the k lattice,
!in reduced coordinates. One builds all k vectors that
!are contained in the first Brillouin zone, with coordinates
!in the interval [0,1[ . First generate boundaries of a big box.

 do jj=1,3
!  To accomodate the shifts, boundmin starts from -1
!  Well, this is not a complete solution ...
   boundmin(jj)=-1
   boundmax(jj)=0
   do ii=1,3
     if(kptrlatt(ii,jj)<0)boundmin(jj)=boundmin(jj)+kptrlatt(ii,jj)
     if(kptrlatt(ii,jj)>0)boundmax(jj)=boundmax(jj)+kptrlatt(ii,jj)
   end do
 end do

 nn=1
 do kk=boundmin(3),boundmax(3)
   do jj=boundmin(2),boundmax(2)
     do ii=boundmin(1),boundmax(1)
       do ikshft=1,nshiftk

!        Coordinates of the trial k point with respect to the k primitive lattice
         k1(1)=ii+shiftk(1,ikshft)
         k1(2)=jj+shiftk(2,ikshft)
         k1(3)=kk+shiftk(3,ikshft)

!        Reduced coordinates of the trial k point
         k2(:)=k1(1)*klatt(:,1)+k1(2)*klatt(:,2)+k1(3)*klatt(:,3)

!        Eliminate the point if outside [0,1[
         if(k2(1)<-tol10)cycle ; if(k2(1)>one-tol10)cycle
         if(k2(2)<-tol10)cycle ; if(k2(2)>one-tol10)cycle
         if(k2(3)<-tol10)cycle ; if(k2(3)>one-tol10)cycle

!        Wrap the trial values in the interval ]-1/2,1/2] .
         call wrap2_pmhalf(k2(1),k1(1),shift(1))
         call wrap2_pmhalf(k2(2),k1(2),shift(2))
         call wrap2_pmhalf(k2(3),k1(3),shift(3))
         if(nn > nkpt_fullbz) then
           write (message,'(a,i0)')' nkpt_fullbz mis-estimated, exceed nn=',nn
           MSG_BUG(message)
         end if
         kpt_fullbz(:,nn)=k1(:)
         nn=nn+1
       end do
     end do
   end do
 end do
 nn = nn-1

 if (nn /= nkpt_fullbz) then
   write (message,'(2(a,i0))')' nkpt_fullbz= ',nkpt_fullbz,' underestimated  nn=',nn
   MSG_BUG(message)
 end if

!make full k-point rank arrays
 call mkkptrank (kpt,nkpt,kptrank_t)

!
!find equivalence to irred kpoints in kpt
!
 indkpt(:) = 0
 timrev=1 ! includes the time inversion symmetry
 do ikpt=1,nkpt_fullbz
   do isym=1,nsym
     do itim=1,(1-2*timrev),-2

       k2(:) = itim*(inv_symrel(:,1,isym)*kpt_fullbz(1,ikpt) + &
&       inv_symrel(:,2,isym)*kpt_fullbz(2,ikpt) + &
&       inv_symrel(:,3,isym)*kpt_fullbz(3,ikpt))

       call get_rank_1kpt (k2,symrankkpt,kptrank_t)
       if (kptrank_t%invrank(symrankkpt) /= -1) &
&       indkpt(ikpt) = kptrank_t%invrank(symrankkpt)

     end do ! loop time reversal symmetry
   end do !  loop sym ops

   if (indkpt(ikpt) == 0) then
     write (message,'(a,i0)')' indkpt(ikpt) is still 0: no irred kpoint is equiv to ikpt ',ikpt
     MSG_BUG(message)
   end if
 end do !  loop full kpts

 call destroy_kptrank (kptrank_t)

end subroutine get_full_kgrid
!!***
