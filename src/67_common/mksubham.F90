!{\src2tex{textfont=tt}}
!!****f* ABINIT/mksubham
!! NAME
!! mksubham
!!
!! FUNCTION
!! This routine build the Hamiltonian matrix in the eigenfunctions subspace,
!! for one given band (or for one given block of bands)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mcg)=wavefunctions
!!  ghc_block(2,npw_k*nspinor,nbdblock)=<G|H|C band,k> for the current block of bands
!!  gsc(2,mgsc)=<g|S|c> matrix elements (S=overlap)
!!  gvnlc_block(2,npw_k*nspinor,nbdblock*use_vnl)=<G|Vnl|C band,k> for the current block of bands
!!  iblock=index of block of bands
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array cg
!!  ikpt=number of the k-point
!!  isppol=spin polarization currently treated
!!  istwf_k=input parameter that describes the storage of wfs
!!  mcg=second dimension of the cg array
!!  mgsc=second dimension of the gsc array
!!  mpi_enreg=informations about MPI parallelization
!!  nband_k=number of bands at this k point for that spin polarization
!!  nbdblock=number of bands in a block
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  use_subovl=1 if the overlap matrix is not identity in WFs subspace
!!  use_vnl= 1 if <C band,k|H|C band_prime,k> has to be computed
!!  wfoptalg=govern the choice of algorithm for wf optimisation
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ghc(2,npw_k*nspinor)=<G|H|C band,k> for the current state
!!                       This is an input in non-blocked algorithm
!!                               an output in blocked algorithm
!!  gvnlc(2,npw_k*nspinor)=<G|Vnl|C band,k> for the current state
!!                       This is an input in non-blocked algorithm
!!                               an output in blocked algorithm
!!  isubh=index of current state in array subham
!!  isubo=index of current state in array subovl
!!  subham(nband_k*(nband_k+1))=Hamiltonian expressed in sthe WFs subspace
!!  subovl(nband_k*(nband_k+1)*use_subovl)=overlap matrix expressed in sthe WFs subspace
!!  subvnl(nband_k*(nband_k+1)*use_vnl)=non-local Hamiltonian expressed in sthe WFs subspace
!!
!! PARENTS
!!      cgwf
!!
!! CHILDREN
!!      xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mksubham(cg,ghc,ghc_block,gsc,gvnlc,gvnlc_block,iblock,icg,igsc,ikpt,isppol,istwf_k,&
&                    isubh,isubo,mcg,mgsc,mpi_enreg,nband_k,nbdblock,npw_k,&
&                    nspinor,subham,subovl,subvnl,use_subovl,use_vnl,wfoptalg)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mksubham'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iblock,icg,igsc,ikpt,isppol,istwf_k,mcg,mgsc,nband_k
 integer,intent(in) :: nbdblock,npw_k,nspinor,use_subovl,use_vnl
 integer,intent(in) :: wfoptalg
 integer,intent(inout) :: isubh,isubo
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: cg(2,mcg),ghc_block(2,npw_k*nspinor,nbdblock)
 real(dp),intent(in) :: gsc(2,mgsc)
 real(dp),intent(in) :: gvnlc_block(2,npw_k*nspinor,nbdblock*use_vnl)
 real(dp),intent(inout) :: ghc(2,npw_k*nspinor),gvnlc(2,npw_k*nspinor)
 real(dp),intent(inout) :: subham(nband_k*(nband_k+1))
 real(dp),intent(inout) :: subovl(nband_k*(nband_k+1)*use_subovl)
 real(dp),intent(inout) :: subvnl(nband_k*(nband_k+1)*use_vnl)

!Local variables-------------------------------
!scalars
 integer :: iband,ibdblock,ii,ipw,ipw1,isp,iwavef,jwavef,me
 real(dp) :: cgimipw,cgreipw,chcim,chcre,cscim,cscre,cvcim,cvcre

! *********************************************************************
!DEBUG
!write(std_out,*)' mksubham: debug, enter.'
!ENDDEBUG

 call xme_init(mpi_enreg,me)

!Loop over bands in a block
!This loop can be parallelized
 do iband=1+(iblock-1)*nbdblock,min(iblock*nbdblock,nband_k)
   ibdblock=iband-(iblock-1)*nbdblock

!  If blocked algorithm, update ghc and gvnlc
   if(mod(wfoptalg,10)==1)then
     if (mpi_enreg%paralbd >= 1) then
       if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/= me) then
         isubh=isubh+2*iband
         isubo=isubo+2*iband
         cycle
       end if
     end if
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(ghc,ghc_block,ibdblock,npw_k,nspinor)
     do ipw=1,npw_k*nspinor
!      If MPI-parallel, ghc and gvnlc should be sent from the master of this k point
!      to the particular processor that takes care of this ibdblock
!      As an alternative, ghc_block and gvnlc_block might be sent to
!      all processors of this k point, during the previous transmission of gcc_block
!      to all processors of this k point
       ghc(1,ipw)=ghc_block(1,ipw,ibdblock)
       ghc(2,ipw)=ghc_block(2,ipw,ibdblock)
     end do
!    $OMP END PARALLEL DO
     if (use_vnl==1) then
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(gvnlc,gvnlc_block,ibdblock,npw_k,nspinor)
       do ipw=1,npw_k*nspinor
         gvnlc(1,ipw)=gvnlc_block(1,ipw,ibdblock)
         gvnlc(2,ipw)=gvnlc_block(2,ipw,ibdblock)
       end do
!      $OMP END PARALLEL DO
     end if
   end if !(end if wfoptalg==1)

!  Compute elements of subspace Hamiltonian <C(i)|H|C(n)> and <C(i)|Vnl|C(n)>
   if(istwf_k==1)then
     do ii=1,iband
       iwavef=(ii-1)*npw_k*nspinor+icg
       chcre=zero ; chcim=zero
       if (use_vnl==0) then
         do ipw=1,npw_k*nspinor
           cgreipw=cg(1,ipw+iwavef)
           cgimipw=cg(2,ipw+iwavef)
           chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
           chcim=chcim+cgreipw*ghc(2,ipw)-cgimipw*ghc(1,ipw)
         end do
       else
         cvcre=zero ; cvcim=zero
         do ipw=1,npw_k*nspinor
           cgreipw=cg(1,ipw+iwavef)
           cgimipw=cg(2,ipw+iwavef)
           chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
           chcim=chcim+cgreipw*ghc(2,ipw)-cgimipw*ghc(1,ipw)
           cvcre=cvcre+cgreipw*gvnlc(1,ipw)+cgimipw*gvnlc(2,ipw)
           cvcim=cvcim+cgreipw*gvnlc(2,ipw)-cgimipw*gvnlc(1,ipw)
         end do
!        Store real and imag parts in Hermitian storage mode:
         subvnl(isubh  )=cvcre
         subvnl(isubh+1)=cvcim
       end if
       subham(isubh  )=chcre
       subham(isubh+1)=chcim
       isubh=isubh+2
     end do
   else if(istwf_k>=2)then
     do ii=1,iband
       iwavef=(ii-1)*npw_k+icg
!      Use the time-reversal symmetry, but should not double-count G=0
       if(istwf_k==2 .and. mpi_enreg%me_g0==1) then
         chcre=0.5_dp*cg(1,1+iwavef)*ghc(1,1)
         if (use_vnl==1) cvcre=0.5_dp*cg(1,1+iwavef)*gvnlc(1,1)
         ipw1=2
       else
         chcre=zero ; ipw1=1
         if (use_vnl==1) cvcre=zero
       end if
       if (use_vnl==0) then
         do isp=1,nspinor
           do ipw=ipw1+(isp-1)*npw_k,npw_k*isp
             cgreipw=cg(1,ipw+iwavef)
             cgimipw=cg(2,ipw+iwavef)
             chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
           end do
         end do
         chcre=2.0_dp*chcre
       else
         do isp=1,nspinor
           do ipw=ipw1+(isp-1)*npw_k,npw_k*isp
             cgreipw=cg(1,ipw+iwavef)
             cgimipw=cg(2,ipw+iwavef)
             chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
             cvcre=cvcre+cgreipw*gvnlc(1,ipw)+cgimipw*gvnlc(2,ipw)
           end do
         end do
         chcre=2.0_dp*chcre
         cvcre=2.0_dp*cvcre
!        Store real and imag parts in Hermitian storage mode:
         subvnl(isubh  )=cvcre
         subvnl(isubh+1)=zero
       end if
       subham(isubh  )=chcre
       subham(isubh+1)=zero
       isubh=isubh+2
     end do
   end if

!  Compute elements of subspace <C(i)|S|C(n)> (S=overlap matrix)
!  <C(i)|S|C(n)> should be closed to Identity.
   if (use_subovl==1) then
     jwavef=(iband-1)*npw_k*nspinor+igsc
     if(istwf_k==1)then
       do ii=1,iband
         iwavef=(ii-1)*npw_k*nspinor+icg
         cscre=zero ; cscim=zero
         do ipw=1,npw_k*nspinor
           cgreipw=cg(1,ipw+iwavef)
           cgimipw=cg(2,ipw+iwavef)
           cscre=cscre+cgreipw*gsc(1,ipw+jwavef)+cgimipw*gsc(2,ipw+jwavef)
           cscim=cscim+cgreipw*gsc(2,ipw+jwavef)-cgimipw*gsc(1,ipw+jwavef)
         end do
!        Store real and imag parts in Hermitian storage mode:
         subovl(isubo  )=cscre
         subovl(isubo+1)=cscim
         isubo=isubo+2
       end do
     else if(istwf_k>=2)then
       do ii=1,iband
         iwavef=(ii-1)*npw_k*nspinor+icg
         if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
           cscre=0.5_dp*cg(1,1+iwavef)*gsc(1,1+jwavef)
           ipw1=2
         else
           cscre=zero ; ipw1=1
         end if
         do isp=1,nspinor
           do ipw=ipw1+(isp-1)*npw_k,npw_k*isp
             cgreipw=cg(1,ipw+iwavef)
             cgimipw=cg(2,ipw+iwavef)
             cscre=cscre+cg(1,ipw+iwavef)*gsc(1,ipw+jwavef)&
&             +cg(2,ipw+iwavef)*gsc(2,ipw+jwavef)
           end do
         end do
         cscre=2.0_dp*cscre
!        Store real and imag parts in Hermitian storage mode:
         subovl(isubo  )=cscre
         subovl(isubo+1)=zero
         isubo=isubo+2
       end do
     end if
   end if

 end do ! iband in a block
!DEBUG
!write(std_out,*)' mksubham: debug, exit.'
!ENDDEBUG

end subroutine mksubham
!!***
