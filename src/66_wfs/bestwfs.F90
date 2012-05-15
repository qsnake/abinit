!{\src2tex{textfont=tt}}
!!****f* ABINIT/bestwfs
!!
!! NAME
!! bestwfs
!!
!! FUNCTION
!! From an Hilbert space with nvectin vectors (gcc_block),
!! and the knowledge of the result of the application
!! of the Hamiltonian to these vectors (ghc_block), select a set
!! of nvectout vectors (nvectout<=nvectin), put it in gcc_block,
!! and update ghc_block, as well as gvnlv_block and (eventually)
!! gscc_block (results of the application of other operators).
!! The criterion might depend upon wfoptalg,
!! but one aims at constructing the eigenvectors of this
!! Hamiltonian (with the lowest eigenvalues if nvectout<nvectin).
!! The algorithm wfoptalg=1 involves simply diagonalizing
!! the Hamiltonian, taking into account the possible non-zero
!! overlap between vectors.
!!
!! The nvectin vectors are expected to be "close"
!! to normalized, but not orthogonal. Especially
!! the overlap for a pair of vectors (1,2), or (3,4) ...
!! might be very close to one, and this special possibility
!! is handled by freezing the first of these vectors.
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gscc_calc=
!!       0 : gscc_block unchanged in this routine
!!       1 : gscc_block has to be updated (usually in PAW calculations)
!!!  gvnlc_calc=
!!       0 : gvnlc_block unchanged in this routine (usually in PAW calculations)
!!       1 : gvnlc_block has to be updated
!  istwf_k=option parameter that describes the storage of wfs
!!  mpi_enreg=informations about MPI parallelization
!!  nbdblock=needed to dimension the set of vectors (nvectin<=nbdblock)
!!  npw_k=number of plane waves
!!  nspinor=number of spinorial components of the wavefunctions
!!  nvectin=number of input vectors, defining the Hilbert space
!!  nvectout=number of output vectors
!!  wfoptalg=algorithm to select the wavefunctions (only one presently)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  gcc_block(2,npw_k*nspinor,nbdblock)=set of nvectin (<=nbdblock)
!!   vectors defining the Hilbert space at input,
!!   of which nvectout vectors will be constructed at output
!!  ghc_block(2,npw_k*nspinor,nbdblock)=Hamiltonian applied
!!   to the vectors in gcc_block
!!  if(gvnlc_calc==1):
!!   gvnlc_block(2,npw_k*nspinor,nbdblock)=Vnl operator
!!    applied to the vectors in gcc_block
!!  if(gscc_calc==1):
!!   gscc_block(2,npw_k*nspinor,nbdblock*gscc_calc)=
!!    Overlap matrix S applied to the vectors in gcc_block
!!
!! WARNING
!!
!! TODO
!!  This routine should be optimized
!!
!! PARENTS
!!      cgwf
!!
!! CHILDREN
!!      leave_new,timab,wrtout,xcomm_init,xsum_mpi,zheev,zhegst,zpotrf,ztrsm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine bestwfs(gcc_block,ghc_block,gscc_block,gscc_calc,&
& gvnlc_block,gvnlc_calc,istwf_k,mpi_enreg,nbdblock,npw_k,nspinor,nvectin,nvectout,wfoptalg)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bestwfs'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none
!Following is the hand tuned interface to lapack routines
!the tuned part are written in lower case
! interface
!  SUBROUTINE ZPOTRF( UPLO, N, A, LDA, INFO )
!   CHARACTER          UPLO
!   INTEGER            INFO, LDA, N
!   real*8         A(2, LDA, * )
!  end subroutine
!  SUBROUTINE ZHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!   CHARACTER          UPLO
!   INTEGER            INFO, ITYPE, LDA, LDB, N
!   real*8         A(2, LDA, * ), B(2, LDB, * )
!  end subroutine
!  SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO)
!   CHARACTER          JOBZ, UPLO
!   INTEGER            INFO, LDA, LWORK, N
!   DOUBLE PRECISION   RWORK( * ), W( * )
!   real*8         A(2, LDA, * ), WORK(2, * )
!  end subroutine
!  SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!   DOUBLE COMPLEX ALPHA
!   INTEGER LDA,LDB,M,N
!   CHARACTER DIAG,SIDE,TRANSA,UPLO
!   double precision A(2, LDA,*),B(2, LDB,*)
!  end subroutine
! end interface

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: gscc_calc,gvnlc_calc,istwf_k,nbdblock,npw_k,nspinor
 integer,intent(in) :: nvectin,nvectout,wfoptalg
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: gcc_block(2,npw_k*nspinor,nbdblock)
 real(dp),intent(inout) :: ghc_block(2,npw_k*nspinor,nbdblock)
 real(dp),intent(inout) :: gscc_block(2,npw_k*nspinor,nbdblock*gscc_calc)
 real(dp),intent(inout) :: gvnlc_block(2,npw_k*nspinor,nbdblock)

!Local variables-------------------------------
 character(len=500) :: message
!scalars
 integer :: iband,ierr,info,ipw,ipw1,isp,jband,old_paral_level,spaceComm
 real(dp) :: ai,ar,bi,br,cosine2,gcim,gcre,norm1,norm2
! arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: eigenvalues(:),hamiltonian(:,:,:),overlap(:,:,:)
 real(dp),allocatable :: rwork(:),tmpgcc(:,:),tmpghc(:,:),tmpgscc(:,:)
 real(dp),allocatable :: tmpgvnlc(:,:),work(:,:)

! *************************************************************************

 if(mod(wfoptalg,10)==1) then

!  Compute the Hamiltonian and overlap matrices
   ABI_ALLOCATE(hamiltonian,(2,nvectin,nvectin))
   ABI_ALLOCATE(overlap,(2,nvectin,nvectin))

   if(istwf_k==1)then

     do jband=1,nvectin
       do iband=1,jband
         ar=0.0_dp ; ai=0.0_dp
         br=0.0_dp ; bi=0.0_dp
         if (gscc_calc==0) then
!          $OMP PARALLEL DO PRIVATE(ipw,gcre,gcim) REDUCTION(+:ai,ar,bi,br) &
!          $OMP&SHARED(gcc_block,ghc_block,iband,jband,npw_k,nspinor)
           do ipw=1,npw_k*nspinor
             gcre=gcc_block(1,ipw,iband); gcim=gcc_block(2,ipw,iband)
             ar=ar+gcre*ghc_block(1,ipw,jband) &
&             +gcim*ghc_block(2,ipw,jband)
             ai=ai+gcre*ghc_block(2,ipw,jband) &
&             -gcim*ghc_block(1,ipw,jband)
             br=br+gcre*gcc_block(1,ipw,jband) &
&             +gcim*gcc_block(2,ipw,jband)
             bi=bi+gcre*gcc_block(2,ipw,jband) &
&             -gcim*gcc_block(1,ipw,jband)
           end do
!          $OMP END PARALLEL DO
         else
!          $OMP PARALLEL DO PRIVATE(ipw,gcre,gcim) REDUCTION(+:ai,ar,bi,br) &
!          $OMP&SHARED(gcc_block,gscc_block,ghc_block,iband,jband,npw_k,nspinor)
           do ipw=1,npw_k*nspinor
             gcre=gcc_block(1,ipw,iband); gcim=gcc_block(2,ipw,iband)
             ar=ar+gcre*ghc_block(1,ipw,jband) &
&             +gcim*ghc_block(2,ipw,jband)
             ai=ai+gcre*ghc_block(2,ipw,jband) &
&             -gcim*ghc_block(1,ipw,jband)
             br=br+gcre*gscc_block(1,ipw,jband) &
&             +gcim*gscc_block(2,ipw,jband)
             bi=bi+gcre*gscc_block(2,ipw,jband) &
&             -gcim*gscc_block(1,ipw,jband)
           end do
!          $OMP END PARALLEL DO
         end if
         hamiltonian(1,iband,jband)=ar
         hamiltonian(2,iband,jband)=ai
         overlap(1,iband,jband)=br
         overlap(2,iband,jband)=bi
       end do
     end do

   else

     do jband=1,nvectin
       do iband=1,jband
         if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
           ar=0.5_dp*gcc_block(1,1,iband)*ghc_block(1,1,jband)
           br=0.5_dp*gcc_block(1,1,iband)*gcc_block(1,1,jband)
           ipw1=2
         else
           ar=0.0_dp ; br=0.0_dp ; ipw1=1
         end if
         if (gscc_calc==0) then
!          $OMP PARALLEL DO PRIVATE(ipw,gcre,gcim) REDUCTION(+:ar,br) &
!          $OMP&SHARED(gcc_block,ghc_block,iband,ipw1,jband,npw_k,nspinor)
           do isp=1,nspinor
             do ipw=ipw1+(isp-1)*npw_k,npw_k*isp
               gcre=gcc_block(1,ipw,iband); gcim=gcc_block(2,ipw,iband)
               ar=ar+gcre*ghc_block(1,ipw,jband) &
&               +gcim*ghc_block(2,ipw,jband)
               br=br+gcre*gcc_block(1,ipw,jband) &
&               +gcim*gcc_block(2,ipw,jband)
             end do
           end do
!          $OMP END PARALLEL DO
         else
!          $OMP PARALLEL DO PRIVATE(ipw,gcre,gcim) REDUCTION(+:ar,br) &
!          $OMP&SHARED(gcc_block,gscc_block,ghc_block,iband,ipw1,jband,npw_k,nspinor)
           do isp=1,nspinor
             do ipw=ipw1+(isp-1)*npw_k,npw_k*isp
               gcre=gcc_block(1,ipw,iband); gcim=gcc_block(2,ipw,iband)
               ar=ar+gcre*ghc_block(1,ipw,jband) &
&               +gcim*ghc_block(2,ipw,jband)
               br=br+gcre*gscc_block(1,ipw,jband) &
&               +gcim*gscc_block(2,ipw,jband)
             end do
           end do
!          $OMP END PARALLEL DO
         end if
         hamiltonian(1,iband,jband)=2.0_dp*ar
         hamiltonian(2,iband,jband)=0.0_dp
         overlap(1,iband,jband)=2.0_dp*br
         overlap(2,iband,jband)=0.0_dp
       end do
     end do

   end if ! istwf_k==1

!  XG030513 : MPIWF reduction on hamiltonian and overlap is needed here
!  Init mpi_comm
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   call timab(48,1,tsec)
   call xsum_mpi(overlap,spaceComm ,ierr)
   call xsum_mpi(hamiltonian,spaceComm ,ierr)
   call timab(48,2,tsec)
   mpi_enreg%paral_level=old_paral_level

!  Complete the matrices thanks to their hermiticity
   do jband=1,nvectin-1
     do iband=jband+1,nvectin
       hamiltonian(1,iband,jband)= hamiltonian(1,jband,iband)
       hamiltonian(2,iband,jband)=-hamiltonian(2,jband,iband)
       overlap(1,iband,jband)= overlap(1,jband,iband)
       overlap(2,iband,jband)=-overlap(2,jband,iband)
     end do
   end do

!  DEBUG
!  write(std_out,*)' bestwfs : overlap and hamiltonian matrices '
!  do jband=1,nvectin
!  do iband=1,nvectin
!  write(std_out,'(2i4,4es16.6)' )iband,jband,&
!  &   overlap(1:2,iband,jband),hamiltonian(1:2,iband,jband)
!  end do
!  end do
!  ENDDEBUG

!  Some pairs of vectors might be linearly dependent,
!  due to the special algorithm in which bestwfs
!  is called. The first vector of the pair is now eliminated,
!  if this happens.
   if(nvectin==2*nvectout)then
     do iband=1,nvectout
!      Compute the cosine between the two vectors
       norm1=0._dp ; norm2=0._dp ; ar=0._dp ; ai=0._dp
       do jband=1,nvectin
         norm1=norm1+overlap(1,jband,iband*2-1)**2 &
&         +overlap(2,jband,iband*2-1)**2
         norm2=norm2+overlap(1,jband,iband*2  )**2 &
&         +overlap(2,jband,iband*2  )**2
         ar=ar+overlap(1,jband,iband*2-1)*overlap(1,jband,iband*2)&
&         +overlap(2,jband,iband*2-1)*overlap(2,jband,iband*2)
         ai=ai+overlap(1,jband,iband*2-1)*overlap(2,jband,iband*2)&
&         -overlap(2,jband,iband*2-1)*overlap(1,jband,iband*2)
       end do
       cosine2=(ar**2+ai**2)/norm1/norm2
       if(abs(one-cosine2)<tol12)then
!        Isolate the first state of the pair from all the others,
!        with a large expectation value
         overlap(:,:,iband*2-1)=zero
         overlap(:,iband*2-1,:)=zero
         overlap(1,iband*2-1,iband*2-1)=one
         hamiltonian(:,:,iband*2-1)=zero
         hamiltonian(:,iband*2-1,:)=zero
         hamiltonian(1,iband*2-1,iband*2-1)=huge(1.0_dp)/1.d10/iband
       end if
     end do
   end if

!  The Hamiltonian and overlap matrices are known.
!  Determine the eigenvalues and eigenvectors of this
!  generalized eigenvalues problem.
   ABI_ALLOCATE(eigenvalues,(nvectin))
   ABI_ALLOCATE(rwork,(3*nvectin))
   ABI_ALLOCATE(work,(2,3*nvectin))

!  The overlap matrix is destroyed, to give its Cholesky factorisation
   call zpotrf('U',nvectin,overlap,nvectin,info)
!  The hamiltonian matrix is destroyed
   call zhegst(1,'U',nvectin,hamiltonian,nvectin,overlap,nvectin,info)
!  The eigenvalues and eigenvectors of the generalized eigenvalue
!  problem are found
   call zheev('V','U',nvectin,hamiltonian,nvectin,eigenvalues,&
&   work,3*nvectin,rwork,info)
!  The Cholesky factor is eliminated
   call ztrsm('L','U','N','N',nvectin,nvectin,(1.0_dp,0.0_dp),&
&   overlap,nvectin,hamiltonian,nvectin)

!  DEBUG
!  write(std_out,*)' bestwfs : overlap and hamiltonian matrices '
!  do iband=1,nvectin
!  write(std_out,'(i4,es16.6)' )iband,eigenvalues(iband)
!  end do
!  ENDDEBUG

   ABI_DEALLOCATE(eigenvalues)
   ABI_DEALLOCATE(rwork)
   ABI_DEALLOCATE(work)

!  Now, must transform the vectors and related quantities
   ABI_ALLOCATE(tmpgcc,(2,nvectout))
   ABI_ALLOCATE(tmpghc,(2,nvectout))
   do ipw=1,npw_k*nspinor
!    Compute the component ipw, for each vector
     do iband=1,nvectout
       ar=0.0_dp ; ai=0.0_dp
       br=0.0_dp ; bi=0.0_dp
       do jband=1,nvectin
         ar=ar+gcc_block(1,ipw,jband)*hamiltonian(1,jband,iband)&
&         -gcc_block(2,ipw,jband)*hamiltonian(2,jband,iband)
         ai=ai+gcc_block(1,ipw,jband)*hamiltonian(2,jband,iband)&
&         +gcc_block(2,ipw,jband)*hamiltonian(1,jband,iband)
         br=br+ghc_block(1,ipw,jband)*hamiltonian(1,jband,iband)&
&         -ghc_block(2,ipw,jband)*hamiltonian(2,jband,iband)
         bi=bi+ghc_block(1,ipw,jband)*hamiltonian(2,jband,iband)&
&         +ghc_block(2,ipw,jband)*hamiltonian(1,jband,iband)
       end do ! jband
       tmpgcc(1,iband)=ar
       tmpgcc(2,iband)=ai
       tmpghc(1,iband)=br
       tmpghc(2,iband)=bi
     end do ! iband
!    Store the component ipw, for each vector
     do iband=1,nvectout
       gcc_block(1,ipw,iband)=tmpgcc(1,iband)
       gcc_block(2,ipw,iband)=tmpgcc(2,iband)
       ghc_block(1,ipw,iband)=tmpghc(1,iband)
       ghc_block(2,ipw,iband)=tmpghc(2,iband)
     end do ! iband

   end do ! ipw
   ABI_DEALLOCATE(tmpgcc)
   ABI_DEALLOCATE(tmpghc)

!  Eventually transform gvnlc_block (according to gvnlc_calc)
   ABI_ALLOCATE(tmpgvnlc,(2,nvectout))
   if (gvnlc_calc==1) then
     do ipw=1,npw_k*nspinor
!      Compute the component ipw, for each vector
       do iband=1,nvectout
         ar=0.0_dp ; ai=0.0_dp
         do jband=1,nvectin
           ar=ar+gvnlc_block(1,ipw,jband)*hamiltonian(1,jband,iband)&
&           -gvnlc_block(2,ipw,jband)*hamiltonian(2,jband,iband)
           ai=ai+gvnlc_block(1,ipw,jband)*hamiltonian(2,jband,iband)&
&           +gvnlc_block(2,ipw,jband)*hamiltonian(1,jband,iband)
         end do ! jband
         tmpgvnlc(1,iband)=ar
         tmpgvnlc(2,iband)=ai
       end do ! iband
!      Store the component ipw, for each vector
       do iband=1,nvectout
         gvnlc_block(1,ipw,iband)=tmpgvnlc(1,iband)
         gvnlc_block(2,ipw,iband)=tmpgvnlc(2,iband)
       end do ! iband
     end do ! ipw
     ABI_DEALLOCATE(tmpgvnlc)
   end if

!  Eventually transform gscc_block (according to gscc_calc)
   if (gscc_calc==1) then
     ABI_ALLOCATE(tmpgscc,(2,nvectout))
     do ipw=1,npw_k*nspinor
!      Compute the component ipw, for each vector
       do iband=1,nvectout
         ar=0.0_dp ; ai=0.0_dp
         do jband=1,nvectin
           ar=ar+gscc_block(1,ipw,jband)*hamiltonian(1,jband,iband)&
&           -gscc_block(2,ipw,jband)*hamiltonian(2,jband,iband)
           ai=ai+gscc_block(1,ipw,jband)*hamiltonian(2,jband,iband)&
&           +gscc_block(2,ipw,jband)*hamiltonian(1,jband,iband)
         end do
         tmpgscc(1,iband)=ar
         tmpgscc(2,iband)=ai
       end do
!      Store the component ipw, for each vector
       do iband=1,nvectout
         gscc_block(1,ipw,iband)=tmpgscc(1,iband)
         gscc_block(2,ipw,iband)=tmpgscc(2,iband)
       end do
     end do
     ABI_DEALLOCATE(tmpgscc)
   end if

   ABI_DEALLOCATE(hamiltonian)
   ABI_DEALLOCATE(overlap)

 else !if(wfoptalg==1)
   write(message,'(a)')' bestwfs wfoptalg/=1 (or 11) is not yet implemented: exit '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL') ! Unsure of the correct arg value for leave_new here : PERS or COLL ?
 end if !if(wfoptalg==1)

end subroutine bestwfs
!!***
