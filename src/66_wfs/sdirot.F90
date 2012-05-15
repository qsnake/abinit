!{\src2tex{textfont=tt}}
!!****f* ABINIT/sdirot
!!
!! NAME
!! sdirot
!!
!! FUNCTION
!! Rotate band states according to subspace eigenvectors to diagonalize
!! Hamiltonian in subspace.  Eigenvectors assumed normalized when
!! input.
!!    $cg_{rotated}(G,n)=Sum(m) [ cg(G,m) * evec(m,n) ]$  for all G s.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcg)=original wavefunctions (to be modified).
!!  evec(2*ndim,num)=normalized eigenvectors of subspace Hamlitonian.
!!  icg=shift to be applied on the location of data in the array cg
!!  mcg=second dimension of the cg array
!!  ndim=dimension of arrays as indicated (usually nband)
!!  num=number of bands under consideration.
!!  npw=number of planewaves in basis at this k point.
!!  work(2*ndim)=work space.
!!
!! OUTPUT
!!  cg(2*npw,num)=rotated wavefunctions (as described above).
!!
!! NOTES
!! xg : this routine should be faster, by blocking !
!! RM : Indeed, it is!
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

subroutine sdirot(cg,evec,icg,mcg,ndim,num,npw)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sdirot'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! Blocksize 8 should be ok for most machines. Can be savely changed.
!scalars
 integer,parameter :: blksz=8
 integer,intent(in) :: icg,mcg,ndim,npw,num
!arrays
 real(dp),intent(in) :: evec(2*ndim,num)
 real(dp),intent(inout) :: cg(2,mcg)

!Local variables-------------------------------
!scalars
 integer :: ig,index1,index2,n2,nb,nblk,nim,nn,nre
!arrays
 real(dp),allocatable :: work(:,:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' sdirot : enter  '
!if(.true.)stop
!ENDDEBUG

!
!Outer loop must run over G vectors because not enough
!storage for keeping intermediate wf results:
!$OMP PARALLEL PRIVATE(ig,index1,index2,nim,nn,nre,n2,nb,nblk) &
!$OMP PRIVATE(work) &
!$OMP SHARED(cg,evec,icg,ndim,npw,num)
 ABI_ALLOCATE(work,(2, blksz, ndim))
!$OMP DO
 do ig=1,npw,blksz
!  Loop over bands:
   index2=ig+icg-1

!  calculate actual block size
   nblk = min(blksz, npw-ig+1)

   if( nblk.eq.blksz ) then
     do nn=1,num
       do nb=1,blksz
         work(1,nb,nn)=0.0
         work(2,nb,nn)=0.0
       end do
     end do
     do n2=1,num
       nim=2*n2
       nre=nim-1
       do nn=1,num
         do nb=1,blksz
           work(1,nb,nn)=work(1,nb,nn)+cg(1,index2+nb)*evec(nre,nn)-&
&           cg(2,index2+nb)*evec(nim,nn)
           work(2,nb,nn)=work(2,nb,nn)+cg(2,index2+nb)*evec(nre,nn)+&
&           cg(1,index2+nb)*evec(nim,nn)
         end do
       end do
       index2=index2+npw
!      Define re and im of C(G,nn) at given G for each nn:
     end do
!    Copy re and im of C(G,nn) into CG over all bands (single G):
     index1=ig+icg-1
     do nn=1,num
       do nb=1,blksz
         cg(1,index1+nb)=work(1,nb,nn)
         cg(2,index1+nb)=work(2,nb,nn)
       end do
       index1=index1+npw
     end do
   else
     do nn=1,num
       do nb=1,nblk
         work(1,nb,nn)=0.0
         work(2,nb,nn)=0.0
       end do
     end do
     do n2=1,num
       nim=2*n2
       nre=nim-1
       do nn=1,num
         do nb=1,nblk
           work(1,nb,nn)=work(1,nb,nn)+cg(1,index2+nb)*evec(nre,nn)-&
&           cg(2,index2+nb)*evec(nim,nn)
           work(2,nb,nn)=work(2,nb,nn)+cg(2,index2+nb)*evec(nre,nn)+&
&           cg(1,index2+nb)*evec(nim,nn)
         end do
       end do
       index2=index2+npw
!      Define re and im of C(G,nn) at given G for each nn:
     end do
!    Copy re and im of C(G,nn) into CG over all bands (single G):
     index1=ig+icg-1
     do nn=1,num
       do nb=1,nblk
         cg(1,index1+nb)=work(1,nb,nn)
         cg(2,index1+nb)=work(2,nb,nn)
       end do
       index1=index1+npw
     end do
   end if
 end do
!$OMP END DO

 ABI_DEALLOCATE(work)
!$OMP END PARALLEL

!DEBUG
!write(std_out,*)' sdirot : exit  '
!ENDDEBUG

end subroutine sdirot
!!***
