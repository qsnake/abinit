!{\src2tex{textfont=tt}}
!!****f* ABINIT/suscep_dyn
!! NAME
!! suscep_dyn
!!
!! FUNCTION
!! Compute the imaginary-frequency susceptibility matrix
!! from input wavefunctions, band occupations, and k point wts.
!! Include the usual sum-over-state terms, but also the
!! corrections due to the change of the Fermi level in the metallic
!! case, as well as implicit sum over higher lying conduction
!! states, thanks to the closure relation (referred to as an extrapolation).
!! In addition, density matrix related quantities can be obtained,
!! see input variable 'susopt'.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG,MF,AR,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dielar(7)=input parameters for dielectric matrix and susceptibility:
!!              diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  freq(nfreq)=array for frequencies (hartree)
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for the dielectric matrix
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  irrzondiel(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  mband=maximum number of bands
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
!!  mkmem=maximum number of k points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum allowed value for npw
!!  nband(nkpt*nsppol)=number of bands to be included in summation
!!   at each k point for each spin channel
!!  nband_mx=like nband (see NOTES)
!!  nfftdiel=number of fft grid points for the computation of the diel matrix
!!  nfreq=size of frequency grid
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwarr(nkpt)=number of planewaves and boundary planewaves
!!   at each k, for going from the WF sphere to the medium size FFT grid.
!!  npwdiel=third and fifth dimension of the susmat_dyn array.
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in group (at least 1 for identity)
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  occopt=option for occupancies
!!  phnonsdiel(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  rprimd(3,3)=dimensional real space primitive translations
!!  susopt=option for functionality, affecting susmat
!!        =0, static susceptibility
!!        =1, dynamic imaginary frequency susceptility
!!        =2, squared modulus of density matrix (for exchange)
!!        =3, density weighted squared modulus of density matrix (for PGG-X kernel)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!  tnons(3,nsym)=reduced nonsymmorphic translations
!!     (symrel and tnons are in terms of real space primitive translations)
!!  ucvol=unit cell volume (Bohr**3)
!!  unkg=unit number for susceptibility (k+G) sphere data file
!!  wtk(nkpt)=k point weights (they sum to 1.0)
!!
!! OUTPUT
!!  sus_diag_dyn(2,npwdiel,nspden,nfreq)= 1st component is the diagonal part
!!   of the frequency-dependent susceptibility matrix, 2nd component is like
!!   1st but without the density contribution in the extrapolation (related
!!   to the density)
!!  susmat_dyn(2,npwdiel,nspden,npwdiel,nspden,nfreq)= the frequency dependent
!!   susceptibility (or density-density response) matrix in reciprocal space
!!
!! SIDE EFFECTS
!!  wff1=structured info about the wavefunction file
!!
!! WARNINGS
!! Restrictions (MF):
!! A Argument occopt>=3 not allowed, since contributions due to Fermi
!!   level changes are not implemented.
!! B For susopt==2 the square of the one-particle density matrix is
!!   computed (and later used to get the exchange energy). The present
!!   implementation is meant for closed shell cases only. dielam or
!!   dielar(6) must be zero in this case.
!! C Spin-polarized calculation not yet possible.
!!
!! NOTES
!! A 'nband_mx' is used to force the maximum number of bands to a value
!!   <= 'nband' after reading them from disk, which requires that 'nband'
!!   as a technical parameter is consistent with the disk file.
!!
!! PARENTS
!!      suscep,xcacfd
!!
!! CHILDREN
!!      fftpac,hdr_skip,leave_test,rdnpw,rwwf,sphereboundary,susk_dyn
!!      susk_dyn_pgg,suskmm_dyn,symg,symrhg,timab,wrtout,xcomm_init
!!      xmaster_init,xme_init,xsum_mpi,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine suscep_dyn(dielar,dtset,&
&  eigen,freq,gbound_diel,gprimd,irrzondiel,istwfk,kg,kg_diel,&
&  mband,mgfftdiel,mkmem,mpi_enreg,mpw,nband,nband_mx,nfftdiel,nfreq,&
&  ngfftdiel,nkpt,npwarr,&
&  npwdiel,nspden,nspinor,nsppol,nsym,occ,occopt,phnonsdiel,rprimd,&
&  susopt,sus_diag_dyn,susmat_dyn,symafm,symrel,tnons,ucvol,unkg,wff1,wtk)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'suscep_dyn'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_67_common
 use interfaces_77_suscep, except_this_one => suscep_dyn
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfftdiel,mkmem,mpw,nband_mx,nfftdiel,nfreq
 integer,intent(in) :: nkpt,npwdiel,nspden,nsppol,nsym,occopt,susopt
 integer,intent(in) :: unkg
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(wffile_type),intent(inout) :: wff1
!arrays
 integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
!no_abirules
!nfftdiel**(1-1/nsym) is 1 if nsym==1, and nfftdiel otherwise
 integer,intent(in) :: irrzondiel(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4)),&
 & istwfk(nkpt)
 integer,intent(in) :: kg(3,mpw*mkmem),kg_diel(3,npwdiel),nband(nkpt*nsppol),&
 & ngfftdiel(18)
 integer,intent(in) :: npwarr(nkpt),symafm(nsym),symrel(3,3,nsym)
 real(dp),intent(in) :: dielar(7)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),&
 & freq(nfreq)
 real(dp),intent(in) :: gprimd(3,3),occ(mband*nkpt*nsppol)
!nfftdiel**(1-1/nsym) is 1 if nsym==1, and nfftdiel otherwise
 real(dp),intent(in) :: phnonsdiel(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: sus_diag_dyn(2,npwdiel,nspden,nfreq)
 real(dp),intent(out) :: susmat_dyn(2,npwdiel,nspden,npwdiel,nspden,nfreq)
 real(dp),intent(in) :: tnons(3,nsym),wtk(nkpt)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,diag,extrap,i1,i2,i3,iband,icg
 integer :: ier,ierr,ifft,ifreq,ii,ikg,ikpt,index,ipw1,ipw2,isp
! integer :: isp1,isp2
 integer :: istwf_k,isym,j1,j2,j3,jj,k1,k2,k3,master,mcg_disk,me
 integer :: nband_k,nband_k_use,ndiel1,ndiel2,ndiel3,ndiel4,ndiel5,ndiel6
 integer :: npw_k,npwsp,spaceComm=0,t1,t2,testocc,tim_rwwf
 real(dp) :: ai,ar,diegap,dielam,emax,invnsym
 real(dp) :: phi1,phi12,phi2,phr1,phr12,phr2,sumdocc,weight
 character(len=500) :: message
!arrays
 integer,allocatable :: gbound(:,:),kg_dum(:,:),kg_k(:,:),sym_g(:,:),tmrev_g(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cg_disk(:,:),drhode(:,:,:),eig_diel(:),eig_dum(:)
 real(dp),allocatable :: occ_deavg(:),occ_dum(:),occ_freq(:,:,:),phdiel(:,:,:)
 real(dp),allocatable :: rhoextrap(:,:,:),rhoextrap_dyn(:,:,:,:,:)
 real(dp),allocatable :: rhoextrg(:,:),rhoextrr_cmplx(:,:),sush(:)
 real(dp),allocatable :: susvec(:,:,:),suswk(:,:,:),zhpev1(:,:),zhpev2(:)

! *************************************************************************

!DEBUG
!write(std_out,*)'%suscep_dyn: enter'
!write(std_out,*)' extrap : ', extrap
!write(std_out,*)' dielam : ', dielar(6), '=dielar(6)'
!write(std_out,*)' freq   : ', freq(:)
!write(std_out,*)' mkmem  : ', mkmem
!write(std_out,*)' nfreq  : ', nfreq
!write(std_out,*)' occopt : ', occopt
!write(std_out,*)' susopt : ', susopt
!write(std_out,*)' unwfnew: ', unwfnew
!call flush(6)
!if(.true.)stop
!ENDDEBUG

!DEBUG
!if(susopt==0 .or. susopt> 3) then
!susopt=1
!write(std_out,*)' %suscep_dyn: WARNING: force susopt=',susopt
!else if(susopt==2 .and. dielar(6) /= 0._dp) then
!write(std_out,*)' %suscep_dyn: BUG: susopt==2 and dielar(6)=dielam==',dielar(6),'not allowed.'
!stop
!end if
!ENDDEBUG

 if( occopt>=3 )then
   write(std_out,*)' suscep_dyn : occopt>= 3 does not yet work, stopping'
   stop
 end if

 if(nspinor==2)then
   write(std_out,*)' suscep_dyn : does not yet work for nspinor=2'
   stop
 end if

 call timab(740,1,tsec)

!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
!Init me
 call xme_init(mpi_enreg,me)
!Init master
 call xmaster_init(mpi_enreg,master)

!mkmem==0 means wf and kg info on disk file
 if (mkmem==0) then

!  Read wff1 header
   call hdr_skip(wff1,ierr)
!  MF now allocated inside spin loop and deallocated before symmetrization
!  allocate(cg_disk(2,mpw*nspinor*mband))
!  Should use xdefineOff, for MPI I/O

 end if

!Initialize some scalar quantities
 bdtot_index=0 ; icg=0

 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)

!ndiel4,ndiel5,ndiel6 are FFT dimensions, modified to avoid cache trashing
 ndiel4=ngfftdiel(4) ; ndiel5=ngfftdiel(5) ; ndiel6=ngfftdiel(6)
 diegap=dielar(5) ; dielam=dielar(6)
 extrap=0
!If dielam is too small, there is no extrapolation.
 if(dielam>1.0d-6)extrap=1

!Perform allocations
 if(occopt>=3)ABI_ALLOCATE(drhode,(2,npwdiel,nspden))
 ABI_ALLOCATE(occ_deavg,(mband))
 ABI_ALLOCATE(occ_freq,(2,mband,nfreq))
 ABI_ALLOCATE(rhoextrap,(ndiel4,ndiel5,ndiel6))
 ABI_ALLOCATE(rhoextrap_dyn,(2,ndiel4,ndiel5,ndiel6,nfreq))

!Perform initializations
 occ_deavg(:)=0._dp
 occ_freq(:,:,:)=0._dp
 rhoextrap(:,:,:)=0._dp
 rhoextrap_dyn(:,:,:,:,:)=0._dp
 susmat_dyn(:,:,:,:,:,:)=0.0_dp
 if(occopt>=3)then
   drhode(:,:,:)=0.0_dp
   sumdocc=0.0_dp
 end if

!testocc to be taken away
 testocc=1
!DEBUG
!write(std_out,*)' suscep : set testocc to 0 '
!testocc=0
!ENDDEBUG

!--BIG loop over spins -----------------------------------------------------

 do isp=1,nsppol

   ikg=0

   if (mkmem==0) then
!    rewind the kpgsph data file on unit unkg
     rewind (unkg)
   end if

   if(extrap==1)rhoextrap(:,:,:)=0.0_dp

   if(mkmem==0) then
     mcg_disk=mpw*nspinor*mband
     ABI_ALLOCATE(cg_disk,(2,mcg_disk))
   end if

!  --BIG loop over k-points -------------------------------------------------

   do ikpt=1,nkpt

!    DEBUG
!    write(std_out,*)' %suscep_dyn : only one k point '
!    do ikpt=1,1
!    ENDDEBUG

     nband_k=nband(ikpt+(isp-1)*nkpt)
     istwf_k=istwfk(ikpt)
     npw_k=npwarr(ikpt)

     if(mpi_enreg%paral_compil_kpt==1)then
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isp) &
&       -mpi_enreg%me))/=0) then
         bdtot_index=bdtot_index+nband_k
         cycle
       end if
     end if

     ABI_ALLOCATE(gbound,(2*mgfftdiel+8,2))
     ABI_ALLOCATE(kg_k,(3,npw_k))

!    Do i/o as needed
     if (mkmem==0) then

       call rdnpw(ikpt,isp,nband_k,npw_k,nspinor,0,unkg)
!      Read k+g data
       read (unkg) kg_k(1:3,1:npw_k)

       call sphereboundary(gbound,istwf_k,kg_k,mgfftdiel,npw_k)

!      DEBUG
!      write(std_out,*) '%suscep_dyn: mark 0 : mband,nband_k,nband_mx=',mband,nband_k,nband_mx
!      call flush(6)
!      ENDDEBUG

!      Read the wavefunction block for ikpt,isp
       tim_rwwf=9
       ABI_ALLOCATE(eig_dum,(mband))
       ABI_ALLOCATE(kg_dum,(3,0))
       ABI_ALLOCATE(occ_dum,(mband))
       call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isp,kg_dum,mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
&       npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wff1)

!      DEBUG
!      write(std_out,*) '%suscep_dyn: after rwwf: mband,nband_k=',mband,nband_k
!      write(std_out,*) '%suscep_dyn: setting nband_k==nband_mx'
!      call flush(6)
!      ENDDEBUG

       ABI_DEALLOCATE(eig_dum)
       ABI_DEALLOCATE(kg_dum)
       ABI_DEALLOCATE(occ_dum)

     else

       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       call sphereboundary(gbound,istwf_k,kg_k,mgfftdiel,npw_k)

!      End test for mkmem==0
     end if

!    DEBUG
!    write(std_out,*) 'suscep_dyn: mark 1'; call flush(6)
!    ENDDEBUG

     if(extrap==1)then
       occ_deavg(:)=0.0_dp; occ_freq(:,:,:)=0._dp
       nband_k_use=min(nband_k,nband_mx)

!      Compute inverse of average dielectric gap for each band
!      and multiply by occupation factor
       emax=maxval(eigen(1+bdtot_index:nband_k_use+bdtot_index))
       do iband=1,nband_k_use
         occ_deavg(iband)= occ(iband+bdtot_index)*dielam &
&         /max( 1.e-20_dp,(emax-eigen(iband+bdtot_index)+diegap) )
       end do

!      frequency dependent occupation factors for extrapolation
       do ifreq=1,nfreq
         do iband=1,nband_k_use
           weight=occ(iband+bdtot_index)*dielam &
&           /max( 1.e-20_dp,(emax-eigen(iband+bdtot_index)+diegap)**2 + freq(ifreq)**2 )
           occ_freq(1,iband,ifreq)=-weight*(emax-eigen(iband+bdtot_index)+diegap)
           occ_freq(2,iband,ifreq)=-weight*freq(ifreq)

!          DEBUG
!          write(std_out,*) ' %suscep_dyn: show on-diagonal weights:'
!          write(std_out,*) '  bdtot_index,freq  :',bdtot_index,freq
!          write(std_out,*) '  iband+...,occ,freq:',iband+bdtot_index,occ(iband+bdtot_index)
!          write(std_out,*) '  emax,emax-eigen   :',emax  ,emax-eigen(iband+bdtot_index)
!          write(std_out,*) '  diegap,...+diegap :',diegap,emax-eigen(iband+bdtot_index)+diegap
!          write(std_out,*) '  occ_freq(1)(2)    :',occ_freq(:,iband,ifreq), 'done'
!          ENDDEBUG

         end do
       end do

     else
       occ_deavg(:)=0.0_dp; occ_freq(:,:,:)=0._dp
       nband_k_use=min(nband_k,nband_mx)

!      for exchange or PGG case use only occupied bands
       if(susopt>=2) then
         do iband=1,nband_k
           if(abs(occ(iband+bdtot_index)) <= tol8 ) exit
         end do
         nband_k_use=iband-1
       end if

     end if

!    DEBUG
!    write(std_out,*) ' %suscep_dyn: nband_k_use', nband_k_use
!    write(std_out,*) ' %suscep_dyn: eigen=  :'
!    write(std_out,*) eigen(1+bdtot_index:min(10,nband_k_use)+bdtot_index)
!    ENDDEBUG

!    Compute the contribution of each k-point to susmat_dyn, rhoextrap_dyn, drhode,
!    and sumdocc.
     if(.true.)then
!      Use either the simpler implementation

       if(mkmem/=0 .or. susopt==2)then

!        DEBUG
!        write(std_out,*)' %suscep_dyn: call susk_dyn, mkmem=',mkmem
!        ENDDEBUG

         call susk_dyn(bdtot_index,cg_disk,dtset,eigen,extrap,freq,&
&         gbound,gbound_diel,icg,ikpt,&
&         isp,istwfk,kg_diel,kg_k,mband,mcg_disk,mgfftdiel,mpi_enreg,mpw,&
&         nband_k_use,ndiel4,ndiel5,ndiel6,nfreq,ngfftdiel,nkpt,&
&         npwdiel,npw_k,nsppol,occ,occopt,occ_deavg,occ_freq,&
&         rhoextrap_dyn,susmat_dyn,susopt,ucvol,wtk)

       else if(susopt==3) then

!        DEBUG
!        write(std_out,*)' %suscep_dyn: call susk_dyn_pgg, mkmem=',mkmem
!        ENDDEBUG
         call susk_dyn_pgg(bdtot_index,cg_disk,eigen,extrap,freq,&
&         gbound,gbound_diel,icg,ikpt,&
&         isp,istwfk,kg_diel,kg_k,mband,mcg_disk,mgfftdiel,mpi_enreg,mpw,&
&         nband_k_use,ndiel4,ndiel5,ndiel6,nfreq,ngfftdiel,nkpt,&
&         npwdiel,npw_k,nsppol,occ,occopt,occ_deavg,occ_freq,dtset%paral_kgb,&
&         rhoextrap_dyn,susmat_dyn,2,ucvol,wtk)

       else

!        DEBUG
!        write(std_out,*)' %suscep_dyn: call suskmm_dyn, mkmem=',mkmem
!        ENDDEBUG

         call suskmm_dyn(bdtot_index,cg_disk,eigen,extrap,freq,&
&         gbound,gbound_diel,icg,ikpt,&
&         isp,istwfk,kg_diel,kg_k,mband,mcg_disk,mgfftdiel,mpi_enreg,mpw,&
&         nband_k_use,ndiel4,ndiel5,ndiel6,nfreq,ngfftdiel,nkpt,&
&         npwdiel,npw_k,nsppol,occ,occopt,occ_deavg,occ_freq,dtset%paral_kgb,&
&         rhoextrap_dyn,susmat_dyn,susopt,ucvol,wtk)
       end if

     else
!      Or the more sophisticated one, needed to save memory.

!      DEBUG
!      write(std_out,*) '%sucep_dyn: call to suskmm not implemented.'; stop
!      ENDDEBUG

!      neglect_pawhat=0
!      if(mkmem/=0)then
!      BEAUTIFICATION NOTE: input variable cg has been removed
!      call suskmm(bdtot_index,cg,doccde,drhode,eigen,extrap,gbound,&
!      &     gbound_diel,icg,ikpt,&
!      &     isp,istwfk,kg_diel,kg_k,mband,mcg,mgfftdiel,mpi_enreg,&
!      &     nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,nfftdiel,ngfftdiel,nkpt,&
!      &     npwdiel,npw_k,nspden,nspinor,nsppol,occ,occopt,occ_deavg,rhoextrap,sumdocc,&
!      &     susmat,ucvol,wtk)
!      else
!      call suskmm(bdtot_index,cg_disk,doccde,drhode,eigen,extrap,gbound,&
!      &     gbound_diel,icg,ikpt,&
!      &     isp,istwfk,kg_diel,kg_k,mband,mcg_disk,mgfftdiel,mpi_enreg,&
!      &     nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,nfftdiel,ngfftdiel,nkpt,&
!      &     npwdiel,npw_k,nspden,nspinor,nsppol,occ,occopt,occ_deavg,rhoextrap,sumdocc,&
!      &     susmat,ucvol,wtk)
!      end if

     end if

     ABI_DEALLOCATE(gbound)
     ABI_DEALLOCATE(kg_k)

     bdtot_index=bdtot_index+nband_k

     if (mkmem/=0) then
       icg=icg+npw_k*nband_k
       ikg=ikg+npw_k
     end if

!    End loop on ikpt:  --------------------------------------------------------
   end do

   if(mkmem==0) then
     ABI_DEALLOCATE(cg_disk)
   end if

!  Here include the contribution from the extrapolation to susmat_dyn, diagonal part
   if(extrap==1)then

!    DEBUG
!    write(std_out,*)' %suscep_dyn: extrapolation: show density part'; call flush(6)
!    ENDDEBUG

     call timab(89,1,tsec)

     ABI_ALLOCATE(rhoextrg,(2,nfftdiel))
     ABI_ALLOCATE(rhoextrr_cmplx,(2*nfftdiel,1))

!    Compute real and imaginary parts of the density extrapolation
     rhoextrr_cmplx(:,:)=0._dp
     do ifreq=1,nfreq

!      Transfer extrapolating density on augmented fft grid to
!      normal fft grid in real space. Warning : must treat only one spin
!      at a time.
!      Complex density: option 10, fill real part into rhoextrr_cmplx(even,:)
       rhoextrap(:,:,:)=rhoextrap_dyn(1,:,:,:,ifreq)
       call fftpac(1,1,2*ndiel1,ndiel2,ndiel3,ndiel4,ndiel5,ndiel6,&
&       ngfftdiel,rhoextrr_cmplx,rhoextrap,10)

!      Complex density: option 11, fill imaginary part into rhoextrr_cmplx(odd,:)
       rhoextrap(:,:,:)=rhoextrap_dyn(2,:,:,:,ifreq)
       call fftpac(1,1,2*ndiel1,ndiel2,ndiel3,ndiel4,ndiel5,ndiel6,&
&       ngfftdiel,rhoextrr_cmplx,rhoextrap,11)

!      DEBUG
!      write(std_out,*) '%suscep_dyn: symmetrization mark 4'; call flush(6)
!      ENDDEBUG

!      rhoextrr_cmplx(:,:) now contains the extrapolating density on augmented
!      fft grid.
!      Generate the density in reciprocal space, and symmetrize it
!      (note symrhg also make the reverse FFT, to get symmetrized density;
!      this is useless here, and should be made an option).
!      Note call to symrhg(2,... instead of symrhg(1,... for complex density.
       call symrhg(2,gprimd,irrzondiel,mpi_enreg,nfftdiel,nfftdiel,ngfftdiel,1,1,nsym,&
&       dtset%paral_kgb,phnonsdiel,rhoextrg,rhoextrr_cmplx,rprimd,symafm,symrel)

!      DEBUG
!      write(std_out,*) '%suscep_dyn: symmetrization mark 5'; call flush(6)
!      ENDDEBUG

!      Save the real part of the diagonal
       do ipw1=1,npwdiel
         sus_diag_dyn(1,ipw1,isp,ifreq)=susmat_dyn(1,ipw1,isp,ipw1,isp,ifreq)
       end do

       do ipw2=1,npwdiel
         j1=kg_diel(1,ipw2) ; j2=kg_diel(2,ipw2) ; j3=kg_diel(3,ipw2)
!        static:    Only fills lower half of the matrix (here, the susceptibility matrix)
!        dynamical: fill all, will not affect susopt==2 for which extrap==0
         do ipw1=1,npwdiel
           i1=kg_diel(1,ipw1) ; i2=kg_diel(2,ipw1) ; i3=kg_diel(3,ipw1)
           k1=i1-j1; k1=modulo(k1,ndiel1)
           k2=i2-j2; k2=modulo(k2,ndiel2)
           k3=i3-j3; k3=modulo(k3,ndiel3)
           ifft=k1+1+ndiel1*(k2+ndiel2*k3)

!          DEBUG
!          write(std_out,'(a,2i3,2es14.6)' ) &
!          &     ' i1,i2,susmat separ',ipw1,ipw2,susmat(1,ipw1,isp,ipw2,isp),&
!          &                                     susmat(2,ipw1,isp,ipw2,isp)
!          write(std_out,'(a,2i3,i5,2es14.6)' ) &
!          &     ' i1,i2,ifft,rhoextrg',ipw1,ipw2,ifft,rhoextrg(1,ifft),rhoextrg(2,ifft)
!          ENDDEBUG

           susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)=   &
&           susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)+rhoextrg(1,ifft)
           susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=   &
&           susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)+rhoextrg(2,ifft)
         end do
       end do

!      Save the real part of the diagonal = squared modulus of the density matrix in G space
       do ipw1=1,npwdiel
         sus_diag_dyn(2,ipw1,isp,ifreq)=susmat_dyn(1,ipw1,isp,ipw1,isp,ifreq)-sus_diag_dyn(1,ipw1,isp,ifreq)
       end do

!      End frequency loop
     end do

     call timab(89,2,tsec)

     ABI_DEALLOCATE(rhoextrg)
     ABI_DEALLOCATE(rhoextrr_cmplx)

!    End condition extrap=1
   end if

!  End loop over spins ---------------------------------------------------------
 end do

 if(mpi_enreg%paral_compil_kpt==1)then
   call timab(86,1,tsec)
!  BEGIN TF_CHANGES
   call leave_test()
!  END TF_CHANGES
   write(message,*) 'suscep: loop on k-points and spins done in parallel'
   call wrtout(std_out,message,'COLL')
   call timab(86,2,tsec)
 end if

!DEBUG
!write(std_out,*) '%suscep_dyn: symmetrization mark 6'; call flush(6)
!ENDDEBUG


 ABI_DEALLOCATE(occ_deavg)
 ABI_DEALLOCATE(occ_freq)
 ABI_DEALLOCATE(rhoextrap_dyn)
!modified by MM 20010515
 ABI_DEALLOCATE(rhoextrap)


!MF has been done inside spin loop
!if(mkmem==0)deallocate(cg_disk)

 if( mpi_enreg%paral_compil_kpt ==1)then
   call timab(85,1,tsec)
!  Recreate full susmat on all proc.
!  This should be coded more efficiently,
!  since half of the matrix is still empty, and
!  it is spin-diagonal.
!  Recreate full drhode on all proc.
   if(occopt>=3 .and. testocc==1)then
     call xsum_mpi(drhode,spaceComm,ierr)
!    Should use only one mpi-allreduce call instead of the three
     call xsum_mpi(sumdocc,spaceComm,ierr)
   end if
   call timab(85,2,tsec)
 end if

 call timab(89,1,tsec)

 if( occopt>=3 .and. testocc==1 )then

!  MF not implemented for frequency dependent case
!  weight=1.0_dp/sumdocc
!  do isp2=1,nsppol
!  do ipw2=1,npwdiel
!  Presently fills complete susceptibility matrix, not only lower half
!  do isp1=1,nsppol
!  do ipw1=1,npwdiel
!  susmat(1,ipw1,isp1,ipw2,isp2)=susmat(1,ipw1,isp1,ipw2,isp2)- &
!  &     weight*( drhode(1,ipw1,isp1)*drhode(1,ipw2,isp2)  &
!  &             +drhode(2,ipw1,isp1)*drhode(2,ipw2,isp2) )
!  susmat(2,ipw1,isp1,ipw2,isp2)=susmat(2,ipw1,isp1,ipw2,isp2)- &
!  &     weight*( drhode(2,ipw1,isp1)*drhode(1,ipw2,isp2)  &
!  &             -drhode(1,ipw1,isp1)*drhode(2,ipw2,isp2) )
!  end do
!  end do
!  end do
!  end do
!  MF

   ABI_DEALLOCATE(drhode)

 end if

!-The susceptibility matrix has been generated---------------------------
!-Symmetries : hermitian, time-reversal, spatial-------------------------

!Generate upper half of the matrix (still the susceptibility matrix)

!DEBUG
!write(std_out,*) '%suscep_dyn: begin symmetrizing:'
!ENDDEBUG

 do ifreq=1,nfreq

!  DEBUG
!  write(std_out,*)'%suscep_dyn: show susmat_dyn() diagonal elements, for ifreq:', ifreq
!  do ipw1=1,npwdiel,20
!  write(std_out,*) ipw1,susmat_dyn(1,ipw1,1,ipw1,1,ifreq),susmat_dyn(2,ipw1,1,ipw1,1,ifreq)
!  end do
!  ENDDEBUG

!  hermitian symmetry
!  MF not valid for non-zero imaginary frequencies
   if(susopt==2 .or. susopt==3) then
     do isp=1,nsppol
       do ipw2=2,npwdiel
         do ipw1=1,ipw2-1
           susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)= susmat_dyn(1,ipw2,isp,ipw1,isp,ifreq)
           susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=-susmat_dyn(2,ipw2,isp,ipw1,isp,ifreq)
         end do
       end do
     end do
   end if

!  End frequency loop
 end do

!DEBUG
!write(std_out,*) '%suscep_dyn: symmetrization mark 7'; call flush(6)
!ENDDEBUG

!Compute symmetric of G-vectors and eventual phases
!(either time-reversal or spatial symmetries)
!Note: have some overhead here if more than one frequency
 ABI_ALLOCATE(suswk,(2,npwdiel,npwdiel))
 ABI_ALLOCATE(tmrev_g,(npwdiel))
 ABI_ALLOCATE(sym_g,(npwdiel,nsym))
 ABI_ALLOCATE(phdiel,(2,npwdiel,nsym))
 call symg(kg_diel,npwdiel,nsym,phdiel,sym_g,symrel,tmrev_g,tnons)
 invnsym=1.0_dp/dble(nsym)

!DEBUG
!write(std_out,*) '%suscep_dyn: symmetrization mark 8'; call flush(6)
!ENDDEBUG

 do isp=1,nsppol

   do ifreq=1,nfreq

     do ipw2=1,npwdiel
       do ipw1=1,npwdiel
         suswk(1,ipw1,ipw2)=susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)
         suswk(2,ipw1,ipw2)=susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)
       end do
     end do

     if(nsym > 1) then

!      Impose spatial symmetries to the susceptibility matrix
       do ipw2=1,npwdiel
         do ipw1=1,npwdiel
           ar=susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)
           ai=susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)
           if(nsym>1)then
             do isym=2,nsym
               t1=sym_g(ipw1,isym) ; t2=sym_g(ipw2,isym)
!              Not all symmetries are non-symmorphic. Should save time here ...
               phr1=phdiel(1,ipw1,isym) ; phi1=phdiel(2,ipw1,isym)
               phr2=phdiel(1,ipw2,isym) ; phi2=phdiel(2,ipw2,isym)
               phr12= phr1*phr2+phi1*phi2 ; phi12=phi1*phr2-phr1*phi2
               ar=ar+susmat_dyn(1,t1,isp,t2,isp,ifreq)*phr12-susmat_dyn(2,t1,isp,t2,isp,ifreq)*phi12
               ai=ai+susmat_dyn(2,t1,isp,t2,isp,ifreq)*phr12+susmat_dyn(1,t1,isp,t2,isp,ifreq)*phi12
             end do
           end if
           suswk(1,ipw1,ipw2)=ar*invnsym
           suswk(2,ipw1,ipw2)=ai*invnsym
         end do
       end do

     end if

!    time reversal
     if(.true.) then
       do ipw2=1,npwdiel
         t2=tmrev_g(ipw2)
         do ipw1=1,npwdiel
           t1=tmrev_g(ipw1)
           susmat_dyn(1,ipw1,isp,ipw2,isp,ifreq)=(suswk(1,ipw1,ipw2)+suswk(1,t1,t2))*0.5_dp
           susmat_dyn(2,ipw1,isp,ipw2,isp,ifreq)=(suswk(2,ipw1,ipw2)-suswk(2,t1,t2))*0.5_dp
         end do
       end do
     end if

!    End frequency loop
   end do

!  DEBUG
!  write(std_out,*) '%suscep_dyn: symmetrization mark 8'; call flush(6)
!  ENDDEBUG

!  End spin loop
 end do

 ABI_DEALLOCATE(suswk)
 ABI_DEALLOCATE(phdiel)
 ABI_DEALLOCATE(sym_g)
 ABI_DEALLOCATE(tmrev_g)

!DEBUG
!write(std_out,*)' %suscep_dyn: end symmetrising'
!ENDDBUG

!DEBUG
!write(std_out,*) '%suscep_dyn: write spin-up susmat_dyn diagonal elements, ifreq=1'
!ifreq=1
!do ipw1=1,npwdiel
!write(std_out,*)ipw1,susmat(1,ipw1,1,ipw1,1,ifreq),susmat(2,ipw1,1,ipw1,1,ifreq)
!end do
!ENDDEBUG

!-The full susceptibility matrix is computed ------------------------------
!-Now, eventually diagonalize it and stop ---------------------------------

!MF diagonalisation not tested
!Must turn on this flag to make the diagonalisation
 diag=0
 if(diag==1)then

   npwsp=npwdiel*nspden
   ABI_ALLOCATE(sush,(npwsp*(npwsp+1)))
   ABI_ALLOCATE(susvec,(2,npwsp,npwsp))
   ABI_ALLOCATE(eig_diel,(npwsp))
   ABI_ALLOCATE(zhpev1,(2,2*npwsp-1))
   ABI_ALLOCATE(zhpev2,(3*npwsp-2))
   ier=0

   do ifreq=1,nfreq
!    Store the susceptibility matrix in proper mode before calling zhpev
     index=1
     do ii=1,npwdiel
       do jj=1,ii
         sush(index  )=susmat_dyn(1,jj,1,ii,1,ifreq)
         sush(index+1)=susmat_dyn(2,jj,1,ii,1,ifreq)
         index=index+2
       end do
     end do

!    If spin-polarized, need to store other parts of the matrix
     if(nsppol/=1)then
       do ii=1,npwdiel
!        Here, spin-flip contribution
         do jj=1,npwdiel
           sush(index  )=susmat_dyn(1,jj,1,ii,2,ifreq)
           sush(index+1)=susmat_dyn(2,jj,1,ii,2,ifreq)
           index=index+2
         end do
!        Here spin down-spin down upper matrix
         do jj=1,ii
           sush(index  )=susmat_dyn(1,jj,2,ii,2,ifreq)
           sush(index+1)=susmat_dyn(2,jj,2,ii,2,ifreq)
           index=index+2
         end do
       end do
     end if

     call ZHPEV ('V','U',npwsp,sush,eig_diel,susvec,npwdiel,zhpev1,&
&     zhpev2,ier)

     write(std_out,*)' suscep_dyn: print eigenvalues of the susceptibility matrix'
     write(std_out,*)' ifreq:',ifreq
     do ii=1,npwdiel
       write(std_out,'(i5,es16.6)' )ii,eig_diel(ii)
     end do

!    End frequency loop
   end do

   ABI_DEALLOCATE(sush)
   ABI_DEALLOCATE(susvec)
   ABI_DEALLOCATE(eig_diel)
   ABI_DEALLOCATE(zhpev1)
   ABI_DEALLOCATE(zhpev2)

 end if

 call timab(89,2,tsec)
 call timab(740,2,tsec)

end subroutine suscep_dyn
!!***
