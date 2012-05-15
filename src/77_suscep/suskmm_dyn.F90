!{\src2tex{textfont=tt}}
!!****f* ABINIT/suskmm_dyn
!! NAME
!! suskmm_dyn
!!
!! FUNCTION
!! Compute the contribution of one k point to the imaginary-frequency
!! susceptibility matrix from input wavefunctions, band occupations,
!! and k point wts. Include the usual sum-over-state terms, but also the
!! corrections due to the change of the Fermi level in the metallic
!! case, as well as implicit sum over higher lying conduction
!! states, thanks to the closure relation (referred to as an extrapolation).
!!
!! This routine is similar to susk_dyn, but use blocking on wavefunctions
!! to decrease memory requirements, at the expense of CPU time.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2012 ABINIT group (MF, XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  bdtot_index=index for the number of the band
!!  cg(2,mcg)=wf in G space
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt
!!           the energy for each band and k point
!!    (used in commented code ... removed from arguments !)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  extrap: if==1, the closure relation (an extrapolation) must be used
!!  freq(nfreq)=array for frequencies (hartree)
!!  gbound(2*mgfftdiel+8,2)=G sphere boundary for going from WF sphere to
!!      medium size FFT grid
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for going from medium size
!!      FFT grid to small sphere.
!!  icg=index for cg
!!  ikpt=number of the k point
!!  isp=number of the current spin
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere.
!!  mband=maximum number of bands
!!  mcg=dimension of cg
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of
!!     the dielectric matrix
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum allowed value for npw
!!  nband_k=number of bands at this k point for that spin polarization
!!  ndiel4,ndiel5,ndiel6= FFT dimensions, modified to avoid cache trashing
!!  nfreq=size of frequency grid
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwdiel=third and fifth dimension of the susmat array.
!!  npw_k=number of plane waves at this k point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  occopt=option for occupancies
!!  occ_deavg(mband)=factor for extrapolation (occup. divided by an energy gap)
!!  occ_freq(2,mband,nfreq)=array holding average weights for extrapolation
!!  susopt=option for susceptibility matrix
!!        =0 static case (not used)
!!        =1 dynamical case
!!        =2 susmat_dyn holds square modulus of the one-particle density matrix
!!  ucvol=unit cell volume (Bohr**3)
!!  wtk(nkpt)=k point weights (they sum to 1.0)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! These quantities are accumulated in this routine:
!! drhode(2,npwdiel,nsppol)=weighted density, needed to compute the
!!   effect of change of fermi energy
!!    (used in commented code ... removed from arguments !)
!!  rhoextrap_dyn(2,ndiel4,ndiel5,ndiel6,ifreq)=density-like array,
!!   needed for extrapolation procedure.
!!  sumdocc=sum of weighted occupation numbers, needed to compute the
!!   effect of change of fermi energy
!!    (used in commented code ... removed from arguments !)
!!  susmat_dyn(2,npwdiel,nsppol,npwdiel,nsppol,nfreq)=the frequency
!!   dependent susceptibility matrix in reciprocal space
!!
!! NOTES
!! There is still room for optimization !!
!!
!! WARNINGS
!! a - Argument occopt>=3 not allowed, since contributions due to Fermi
!!     level changes are not implemented.
!! b - Only susopt==0 or 1 should be used. susopt==2 may be handled
!!     through routine susk_dyn instead.
!!
!! TODO
!! a - resolve Warnings
!! b - time and try to optimize wfprod in real space for large fft meshes,
!!     possibly save time by precomputing the wavefunction products for
!!     all bands, might require to put/take them to/from disk
!! c - why not save memory by making susmat_dyn single precision?
!!
!! PARENTS
!!      suscep_dyn,suscep_kxc_dyn
!!
!! CHILDREN
!!      fourwf,leave_new,timab,wrtout,zgerc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine suskmm_dyn(bdtot_index,cg,eigen,extrap,freq,&
&  gbound,gbound_diel,icg,ikpt,isp,istwfk,kg_diel,kg_k,&
&  mband,mcg,mgfftdiel,mpi_enreg,mpw,&
&  nband_k,ndiel4,ndiel5,ndiel6,nfreq,ngfftdiel,nkpt,&
&  npwdiel,npw_k,nsppol,occ,occopt,&
&  occ_deavg,occ_freq,paral_kgb,rhoextrap_dyn,susmat_dyn,susopt,ucvol,wtk)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'suskmm_dyn'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: bdtot_index,extrap,icg,ikpt,isp,mband,mcg,mgfftdiel
 integer,intent(in) :: mpw,nband_k,ndiel4,ndiel5,ndiel6,nfreq
 integer,intent(in) :: nkpt,npw_k,npwdiel,nsppol,occopt
 integer,intent(in) :: paral_kgb,susopt
 real(dp),intent(in) :: ucvol
! real(dp),intent(inout) :: sumdocc
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: gbound(2*mgfftdiel+8,2),gbound_diel(2*mgfftdiel+8,2)
 integer,intent(in) :: istwfk(nkpt),kg_diel(3,npwdiel),kg_k(3,npw_k)
 integer,intent(in) :: ngfftdiel(18)
 real(dp),intent(in) :: cg(2,mcg)
! real(dp),intent(in) :: doccde(mband*nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),freq(nfreq)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),occ_deavg(mband)
 real(dp),intent(in) :: occ_freq(2,mband,nfreq),wtk(nkpt)
! real(dp),intent(inout) :: drhode(2,npwdiel,nsppol)
 real(dp),intent(inout) :: rhoextrap_dyn(2,ndiel4,ndiel5,ndiel6,nfreq)
 real(dp),intent(inout) :: susmat_dyn(2,npwdiel,nsppol,npwdiel,nsppol,nfreq)

!Local variables-------------------------------
! real(dp), allocatable :: cg_disk(:,:)
!scalars
 integer :: i1,i2,i3,iband,iband_shift,iband_shift2,ibd1,ibd2,ibdshft1,ibdshft2
 integer :: iblk1,iblk2,ifreq,ipw1,ipw2,istwf_k,mblk,nblk
 integer :: nbnd_current,nbnd_in_blk,nbnd_in_blk1,ndiel1,ndiel2,ndiel3,testocc
 integer :: tim_fourwf
! integer :: ipw
 real(dp) :: eigdiff,occdiff,tolocc,weight,wght1
! real(dp) :: wght2
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),dummy(:,:),rhoaug(:,:,:),wfprod(:,:)
 real(dp),allocatable :: wfraug(:,:,:,:),wfrspa1(:,:,:,:,:),wfrspa2(:,:,:,:,:)
 real(dp),allocatable :: wght_dyn(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' suskmm_dyn : ENTER '
!write(std_out,*)' extrap:', extrap
!write(std_out,*)' nfreq :', nfreq
!write(std_out,*)' mpw   :', mpw
!write(std_out,*)' npw_k :', npw_k
!write(std_out,*)' occopt:', occopt
!write(std_out,*)' susopt:', susopt
!if(.true.)stop
!ENDDEBUG

 if(occopt>=3) then
   write(std_out,*) ' %suskmm_dyn: occopt>=3 not implemented: stop '
   stop
 end if

 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)
 istwf_k=istwfk(1)

 testocc=1

 ABI_ALLOCATE(cwavef,(2,mpw))
 ABI_ALLOCATE(dummy,(2,1))
 ABI_ALLOCATE(rhoaug,(ndiel4,ndiel5,ndiel6))
 ABI_ALLOCATE(wfraug,(2,ndiel4,ndiel5,ndiel6))
 ABI_ALLOCATE(wfprod,(2,npwdiel))
 ABI_ALLOCATE(wght_dyn,(2,nfreq))

 rhoaug(:,:,:)=0._dp
 wfraug(:,:,:,:)=0._dp
 wfprod(:,:)=0._dp
 wght_dyn(:,:)=0.0_dp

!Prepare the blocking : compute the number of blocks,
!the number of bands in each normal block,
!and the number in the first one, usually smaller.

!Consider that if the number of bands is large, there are at most 8 blocks
 nbnd_in_blk=0
 if(nband_k>=48)then
   mblk=8
   nbnd_in_blk=(nband_k-1)/mblk+1
!  If the number of bands is medium, place 6 bands per block
 else if(nband_k>=12)then
   nbnd_in_blk=6
!  Otherwise, must have at least 2 blocks
 else if(nband_k>=2)then
   mblk=2
   nbnd_in_blk=(nband_k-1)/mblk+1
 else
   write(message, '(a,a,a,a,a,a,i2,a,a,a)') ch10,&
&   ' suskmm_dyn : ERROR -',ch10,&
&   '  The number of bands must be larger or equal to 2, in suskmm_dyn.',ch10,&
&   '  It is equal to ',nband_k,'.',ch10,&
&   '  Action : choose another preconditioner.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

!Compute the effective number of blocks, and the number of bands in
!the first block.
 nblk=(nband_k-1)/nbnd_in_blk+1
 nbnd_in_blk1=nband_k-(nblk-1)*nbnd_in_blk

!DEBUG
!write(std_out,*)' suskmm_dyn : nband_k,nblk,nbnd_in_blk,nbnd_in_blk1 '
!write(std_out,*) nband_k,nblk,nbnd_in_blk,nbnd_in_blk1
!stop
!ENDDEBUG

!wfrspa1 will contain the wavefunctions of the slow sampling (iblk1)
 ABI_ALLOCATE(wfrspa1,(2,ndiel4,ndiel5,ndiel6,nbnd_in_blk))
!wfrspa2 will contain the wavefunctions of the rapid sampling (iblk2)
 ABI_ALLOCATE(wfrspa2,(2,ndiel4,ndiel5,ndiel6,nbnd_in_blk))

 wfrspa1(:,:,:,:,:)=0._dp
 wfrspa2(:,:,:,:,:)=0._dp

!First loop over blocks
 do iblk1=1,nblk

   call timab(87,1,tsec)

!  Initialisation
   if(iblk1==1)then

     nbnd_current=nbnd_in_blk1
     iband_shift=0
!    Loop over bands to fft and store Fourier transform of wavefunction
     do iband=1,nbnd_current
!      Obtain Fourier transform in fft box
       cwavef(:,1:npw_k)=cg(:,1+(iband-1)*npw_k+icg:iband*npw_k+icg)
       tim_fourwf=12
       call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&       istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg,1,ngfftdiel,npw_k,1,ndiel4,ndiel5,ndiel6,&
&       0,paral_kgb,tim_fourwf,weight,weight)
       wfrspa1(:,:,:,:,iband)=wfraug(:,:,:,:)
     end do

   else

!    The Fourier transform of wavefunctions have already been obtained
     nbnd_current=nbnd_in_blk
     iband_shift=nbnd_in_blk1+(iblk1-2)*nbnd_in_blk

   end if

!  Loop over bands of this block, to generate band-diagonal
!  contributions to sumdocc, drhode, rhoextrap, and susmat.

!  DEBUG
!  write(std_out,*)' suskmm_dyn : 1'
!  ENDDEBUG

   do iband=1,nbnd_current

     if( (occopt>=3 .and. testocc==1) .or. extrap==1 )then
!      In the case of metallic occupation, or if the extrapolation
!      over higher bands is included, must compute the
!      Fourier transform of the density of each band, then
!      generate the part of the susceptibility matrix due
!      varying occupation numbers.

!      MF Not needed as occopt>=3 not implemented
!      weight=-2.0_dp*occ_deavg(iband+iband_shift)*wtk(ikpt)/ucvol

!      Case where the band density does not need to be accumulated
       if(extrap==0 .or. susopt==2) then
         do i3=1,ndiel3
           do i2=1,ndiel2
             do i1=1,ndiel1
               wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,iband)**2&
&               +wfrspa1(2,i1,i2,i3,iband)**2
               wfraug(2,i1,i2,i3)=0.0_dp
             end do
           end do
         end do

!        Accumulate density in real space
       else if(extrap==1) then
         do ifreq=1,nfreq
           wght_dyn(:,ifreq)=2.0_dp*occ_freq(:,iband+iband_shift,ifreq)*wtk(ikpt)/ucvol

!          DEBUG
!          write(std_out,*) ' %suskmm_dyn: show on-diagonal weights:'
!          write(std_out,*) '  bands: iband    :',iband+iband_shift
!          write(std_out,*) '  occ_freq(1)(2)  :',occ_freq(:,iband+iband_shift,ifreq)
!          write(std_out,*) '  wght_dyn(1)(2)  :',wght_dyn(:,ifreq)
!          ENDDEBUG

           do i3=1,ndiel3
             do i2=1,ndiel2
               do i1=1,ndiel1
                 wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,iband)**2&
&                 +wfrspa1(2,i1,i2,i3,iband)**2
                 wfraug(2,i1,i2,i3)=0.0_dp
                 rhoextrap_dyn(1,i1,i2,i3,ifreq)=rhoextrap_dyn(1,i1,i2,i3,ifreq) &
&                 +wght_dyn(1,ifreq)*wfraug(1,i1,i2,i3)
               end do
             end do
           end do

!          End frequency loop
         end do

!        End condition extrap or susopt
       end if

!      Performs the Fourier Transform of the density of the band,
!      and store it in wfprod
       tim_fourwf=13

       call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&       istwf_k,kg_diel,kg_diel,&
&       mgfftdiel,mpi_enreg,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,paral_kgb,tim_fourwf,weight,weight)

!      For density matrix related case
       if(susopt==2) then
         ifreq=1
         weight=2._dp*wtk(ikpt)/ucvol
         do ipw2=1,npwdiel
           do ipw1=ipw2,npwdiel
             susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)+&
&             weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
             susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)+&
&             weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
           end do
         end do

!        For susceptibility matrix case, perform now the summation of terms related
!        to direct change of eigenvalues or extrapolation over higher bands
       else
!        Loop over frequencies
         do ifreq=1,nfreq

!          MF   Not needed as occopt>=3 not implemented
!          wght1=0.0_dp ; wght2=0.0_dp
!          if(occopt>=3 .and. testocc==1)then
!          wght1=doccde(iband+iband_shift+bdtot_index)*wtk(ikpt)/ucvol
!          end if
!          if(extrap==1) wght2=2.0_dp*occ_deavg(iband+iband_shift)*wtk(ikpt)/ucvol
!          weight=wght1+wght2
!          MF

           if(extrap==1) wght_dyn(:,ifreq)=-wght_dyn(:,ifreq)

           call ZGERC(npwdiel,npwdiel,wght_dyn(1,ifreq),&
&           wfprod(1,1),1,&
&           wfprod(1,1),1,&
&           susmat_dyn(1,1,isp,1,isp,ifreq),npwdiel)

!          MF   Not needed as occopt>=3 not implemented
!          if( occopt>=3 .and. testocc==1) then
!          Accumulate product of band densities by their doccde, for the
!          computation of the effect of change of Fermi level.
!          do ipw=1,npwdiel
!          drhode(1,ipw,isp)=drhode(1,ipw,isp)+wfprod(1,ipw)*wght1
!          drhode(2,ipw,isp)=drhode(2,ipw,isp)+wfprod(2,ipw)*wght1
!          end do
!          Also accumulate weighted sum of doccde
!          sumdocc=sumdocc+wght1
!          end if
!          MF

!          End frequency loop
         end do

!        End condition susopt
       end if

!      End condition of metallic occupancies or extrapolation
     end if

!    End loop on iband
   end do

!  DEBUG
!  write(std_out,*)' suskmm_dyn : 2'
!  ENDDEBUG

   call timab(87,2,tsec)

!  -- Compute now off-band-diagonal terms ------------------------------------

   call timab(88,1,tsec)

   tolocc=1.0d-3

!  Compute product of wavefunctions for different bands, inside the blok
   if(nbnd_current>1)then
     do ibd1=1,nbnd_current-1
       ibdshft1=ibd1+iband_shift
       do ibd2=ibd1+1,nbnd_current
         ibdshft2=ibd2+iband_shift

!        If the occupation numbers are sufficiently different, or
!        if extrapolation is used and the corresponding factor is not zero,
!        then there is a contribution
         occdiff=occ(ibdshft1+bdtot_index)-occ(ibdshft2+bdtot_index)
         if( abs(occdiff)>tolocc      .or. &
&         ( extrap==1 .and.            &
&         ( abs(occ_deavg(ibdshft1)) + abs(occ_deavg(ibdshft2)) ) >tolocc ) &
&         .or. susopt==2 ) then

           eigdiff=eigen(ibdshft1+bdtot_index) - eigen(ibdshft2+bdtot_index)

!          Store the contribution in wfraug
           do i3=1,ndiel3
             do i2=1,ndiel2
               do i1=1,ndiel1
                 wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,ibd1)*wfrspa1(1,i1,i2,i3,ibd2)&
&                 +wfrspa1(2,i1,i2,i3,ibd1)*wfrspa1(2,i1,i2,i3,ibd2)
                 wfraug(2,i1,i2,i3)=wfrspa1(2,i1,i2,i3,ibd1)*wfrspa1(1,i1,i2,i3,ibd2)&
&                 -wfrspa1(1,i1,i2,i3,ibd1)*wfrspa1(2,i1,i2,i3,ibd2)
               end do
             end do
           end do

!          Performs the Fourier Transform of the product, and store it in wfprod
           tim_fourwf=13
           call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&           istwf_k,kg_diel,kg_diel,&
&           mgfftdiel,mpi_enreg,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,paral_kgb,tim_fourwf,weight,weight)

!          For density matrix related case
           if(susopt==2) then

             ifreq=1
             weight=4._dp*wtk(ikpt)/ucvol
             do ipw2=1,npwdiel
               do ipw1=ipw2,npwdiel
                 susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)+&
&                 weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
                 susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)+&
&                 weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
               end do
             end do

!            For susceptibility matrix case
           else

!            Loop over frequencies
             do ifreq=1,nfreq

!              Step 1 - determine the weight for each frequency
               wght_dyn(:,ifreq)=0.0_dp
               if(abs(occdiff)>tolocc) then
                 wght1=occdiff/max( 1.e-20_dp,freq(ifreq)**2+eigdiff**2 ) * 2.0_dp*wtk(ikpt)/ucvol
!                MF note: eigdiff is *minus* the difference (e_unoccupied - e_occupied)
                 wght_dyn(1,ifreq)=eigdiff*wght1
                 wght_dyn(2,ifreq)=-freq(ifreq)*wght1
               else
                 wght1=0.0_dp
               end if

               if(extrap==1)then
                 wght_dyn(1,ifreq)=wght_dyn(1,ifreq)- 2.0_dp*wtk(ikpt)/ucvol &
&                 *(occ_freq(1,ibdshft1,ifreq)+occ_freq(1,ibdshft2,ifreq))
                 wght_dyn(2,ifreq)=wght_dyn(2,ifreq)- 2.0_dp*wtk(ikpt)/ucvol &
&                 *(occ_freq(2,ibdshft1,ifreq)+occ_freq(2,ibdshft2,ifreq))
               end if

!              Step 2 - sum contribution for each frequency
               call ZGERC(npwdiel,npwdiel,wght_dyn(1,ifreq),&
&               wfprod(1,1),1,&
&               wfprod(1,1),1,&
&               susmat_dyn(1,1,isp,1,isp,ifreq),npwdiel)

!              End frequency loop
             end do
!            End condition susopt
           end if
!          End condition of different occupation numbers or extrapolation
         end if
!        End internal loop over bands
       end do
!      End external loop over bands
     end do
!    End condition of having more than one band
   end if

!  DEBUG
!  write(std_out,*)' suskmm_dyn : 3'
!  ENDDEBUG

!  Loop on secondary block, with fast varying index, in decreasing order.
   if(iblk1/=nblk)then
     do iblk2=nblk,iblk1+1,-1

!      Loop over bands to fft and store Fourier transform of wavefunction
       iband_shift2=nbnd_in_blk1+(iblk2-2)*nbnd_in_blk
       do iband=1,nbnd_in_blk
!        Obtain Fourier transform in fft box
         cwavef(:,1:npw_k)= &
&         cg(:,1+(iband+iband_shift2-1)*npw_k+icg:(iband+iband_shift2)*npw_k+icg)
         tim_fourwf=12
         call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&         istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg,1,ngfftdiel,npw_k,1,&
&         ndiel4,ndiel5,ndiel6,0,paral_kgb,tim_fourwf,weight,weight)
         wfrspa2(:,:,:,:,iband)=wfraug(:,:,:,:)
       end do

       do ibd1=1,nbnd_current
         ibdshft1=ibd1+iband_shift
         do ibd2=1,nbnd_in_blk
           ibdshft2=ibd2+iband_shift2

!          If the occupation numbers are sufficiently different, or
!          if extrapolation is used and the corresponding factor is not zero,
!          then there is a contribution
           occdiff=occ(ibdshft1+bdtot_index)-occ(ibdshft2+bdtot_index)
           if( abs(occdiff)>tolocc      .or. &
&           ( extrap==1 .and.            &
&           ( abs(occ_deavg(ibdshft1)) + abs(occ_deavg(ibdshft2)) ) >tolocc ) &
&           .or. susopt==2 ) then

             eigdiff=eigen(ibdshft1+bdtot_index) - eigen(ibdshft2+bdtot_index)

!            Store the contribution in wfraug
             do i3=1,ndiel3
               do i2=1,ndiel2
                 do i1=1,ndiel1
                   wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,ibd1)*wfrspa2(1,i1,i2,i3,ibd2)&
&                   +wfrspa1(2,i1,i2,i3,ibd1)*wfrspa2(2,i1,i2,i3,ibd2)
                   wfraug(2,i1,i2,i3)=wfrspa1(2,i1,i2,i3,ibd1)*wfrspa2(1,i1,i2,i3,ibd2)&
&                   -wfrspa1(1,i1,i2,i3,ibd1)*wfrspa2(2,i1,i2,i3,ibd2)
                 end do
               end do
             end do

!            Performs the Fourier Transform of the product, and store it in wfprod
             tim_fourwf=13
             call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&             istwf_k,kg_diel,kg_diel,&
&             mgfftdiel,mpi_enreg,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,paral_kgb,tim_fourwf,weight,weight)

!            For density matrix related case
             if(susopt==2) then

               ifreq=1
               weight=4._dp*wtk(ikpt)/ucvol
               do ipw2=1,npwdiel
                 do ipw1=ipw2,npwdiel
                   susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)+&
&                   weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
                   susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)+&
&                   weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
                 end do
               end do

!              For susceptibility matrix case
             else

!              Loop over frequencies
               do ifreq=1,nfreq

!                Step 1 - determine the weight for each frequency
                 wght_dyn(:,ifreq)=0.0_dp
                 if(abs(occdiff)>tolocc) then
                   wght1=occdiff/max( 1.e-20_dp,freq(ifreq)**2+eigdiff**2 ) * 2.0_dp*wtk(ikpt)/ucvol
!                  MF note: eigdiff is *minus* the difference (e_unoccupied - e_occupied)
                   wght_dyn(1,ifreq)=eigdiff*wght1
                   wght_dyn(2,ifreq)=-freq(ifreq)*wght1
                 else
                   wght1=0.0_dp
                 end if

                 if(extrap==1)then
                   wght_dyn(1,ifreq)=wght_dyn(1,ifreq)- 2.0_dp*wtk(ikpt)/ucvol &
&                   *(occ_freq(1,ibdshft1,ifreq)+occ_freq(1,ibdshft2,ifreq))
                   wght_dyn(2,ifreq)=wght_dyn(2,ifreq)- 2.0_dp*wtk(ikpt)/ucvol &
&                   *(occ_freq(2,ibdshft1,ifreq)+occ_freq(2,ibdshft2,ifreq))
                 end if

!                Step 2 - sum contribution for each frequency
                 call ZGERC(npwdiel,npwdiel,wght_dyn(1,ifreq),&
&                 wfprod(1,1),1,&
&                 wfprod(1,1),1,&
&                 susmat_dyn(1,1,isp,1,isp,ifreq),npwdiel)

!                End frequency loop
               end do
!              End condition susopt
             end if
!            End condition of different occupation numbers or extrapolation
           end if
!          End internal loop over bands
         end do
!        End external loop over bands
       end do
!      End loop on bloks
     end do

!    Finish the loop on blok with iblk2=iblk1+1, so can use the
!    FFTd wavefunctions for the next iblk1.
     do iband=1,nbnd_in_blk
       wfrspa1(:,:,:,:,iband)=wfrspa2(:,:,:,:,iband)
     end do

!    End condition of iblk1/=nblk
   end if

   call timab(88,2,tsec)

!  End loop on iblk1
 end do

!DEBUG
!write(std_out,*)' suskmm_dyn : DONE & EXIT '
!do ipw1=1,npwdiel
!write(std_out,*)ipw1,susmat(1,ipw1,1,ipw1,1),susmat(2,ipw1,1,ipw1,1)
!end do
!write(std_out,*)' suskmm_dyn : end of susmat '
!stop
!ENDDEBUG

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(dummy)
 ABI_DEALLOCATE(rhoaug)
 ABI_DEALLOCATE(wfprod)
 ABI_DEALLOCATE(wfraug)
 ABI_DEALLOCATE(wfrspa1)
 ABI_DEALLOCATE(wfrspa2)
 ABI_DEALLOCATE(wght_dyn)

end subroutine suskmm_dyn
!!***
