!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkrho
!! NAME
!! mkrho
!!
!! FUNCTION
!! Depending on option argument value:
!! --Compute charge density rho(r) and rho(G) in electrons/bohr**3
!!   from input wavefunctions, band occupations, and k point wts.
!! --Compute kinetic energy density tau(r) and tau(G) in bohr**-5
!!   from input wavefunctions, band occupations, and k point wts.
!! --Compute a given element of the kinetic energy density tensor
!!   tau_{alpha,beta}(r) and tau_{alpha,beta}(G) in bohr**-5
!!   from input wavefunctions, band occupations, and k point wts.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, LSI, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcg)=wf in G space
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | istwfk(nkpt)=input option parameter that describes the storage of wfs
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem=maximum number of k points in core memory
!!   | mpw=maximum allowed value for npw
!!   | nband(nkpt*nsppol)=number of bands to be included in summation
!!   |  at each k point for each spin channel
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | ngfft(18)=contain all needed information about 3D FFT,
!!   |  see ~abinit/doc/input_variables/vargs.htm#ngfft
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in group (at least 1 for identity)
!!   | symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!   | symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!   | wtk(nkpt)=k point weights (they sum to 1.0)
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=informations about MPI parallelization
!!  npwarr(nkpt)=number of planewaves and boundary planewaves at each k
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  option if 0: compute rhor (electron density)
!!         if 1: compute taur (kinetic energy density)
!!               (i.e. Trace over the kinetic energy density tensor)
!!         if 2: compute taur_{alpha,beta} !!NOT YET IMPLEMENTED
!!               (a given element of the kinetic energy density tensor)
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  rprimd(3,3)=dimensional real space primitive translations
!!  tim_mkrho=timing code of the calling routine(can be set to 0 if not attributed)
!!  ucvol=unit cell volume (Bohr**3)
!!  unkg=unit number for (k+G) sphere data file
!!  wffnow=struct info for current wf disk file
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!! rhog(2,nfft)=total electron density in G space
!! rhor(nfft,nspden)=electron density in r space
!!   (if spin polarized, array contains total density in first half and
!!    spin-up density in second half)
!!   (for non-collinear magnetism, first element: total density, 3 next ones: mx,my,mz in units of hbar/2)
!!
!! PARENTS
!!      afterscfloop,energy,gstate,respfn,vtorho
!!
!! CHILDREN
!!      fftpac,fourwf,hdr_skip,leave_test,prep_fourwf,prtrhomxmn,rdnpw,rwwf
!!      sphereboundary,symrhg,timab,wrtout,wvl_mkrho,xcomm_init,xdefineoff
!!      xmaster_init,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&
&                rhog,rhor,rprimd,tim_mkrho,ucvol,unkg,wffnow,wfs,wvl,&
&                option) !optional

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_wffile
 use m_errors

 use m_paw_dmft, only: paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkrho'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_66_wfs
 use interfaces_67_common, except_this_one => mkrho
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,tim_mkrho,unkg
 integer,intent(in),optional :: option
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(wffile_type),intent(inout) :: wffnow
 type(wvl_wf_type),intent(inout) :: wfs
 type(wvl_internal_type), intent(in) :: wvl
!no_abirules
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,  &
&               (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)
 real(dp), intent(in) :: gprimd(3,3)
 real(dp), intent(in) :: cg(2,mcg)
 real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 real(dp), intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3))**(1-1/dtset%nsym),  &
&                                 (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(out) :: rhor(dtset%nfft,dtset%nspden),rhog(2,dtset%nfft)

!Local variables-------------------------------
!scalars
 integer,save :: nskip=0
 integer :: alpha,use_nondiag_occup_dmft,bdtot_index,beta,blocksize,formeig,iband,iband1,ibandc1,ib,iblock,icg,ierr
 integer :: ifft,ikg,ikpt,ikpt_this_proc,ioption,ipw,ipwsp,ishf,ispden,ispinor,ispinor_index
 integer :: isppol,istwf_k,jspinor_index,master
 integer :: mcg_disk,me,my_nspinor,n1,n2,n3,n4,n5,n6,nalpha,nband_k,nbandc1,nbdblock,nbeta
 integer :: nfftot,npw_k,nsp,spaceComm
 integer :: tim_fourwf,tim_rwwf
 real(dp) :: kpt_cart,kg_k_cart,gp2pi1,gp2pi2,gp2pi3,cwftmp
 real(dp) :: weight
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk_local
!arrays
 integer,allocatable :: gbound(:,:),kg_dum(:,:),kg_k(:,:)
 logical :: locc_test,nspinor1TreatedByThisProc,nspinor2TreatedByThisProc
 real(dp) :: dummy(2,1),tsec(2)
 real(dp),allocatable :: cg_disk(:,:),cwavef(:,:,:),cwavefb(:,:),cwavef_x(:,:)
 real(dp),allocatable :: cwavef_y(:,:),eig_dum(:),kg_k_cart_block(:)
 real(dp),allocatable :: occ_dum(:),occ_k(:),rhoaug(:,:,:),rhoaug_down(:,:,:)
 real(dp),allocatable :: rhoaug_mx(:,:,:),rhoaug_my(:,:,:),rhoaug_up(:,:,:)
 real(dp),allocatable :: taur_alphabeta(:,:,:,:),wfraug(:,:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(790+tim_mkrho,1,tsec)
 call timab(799,1,tsec)

 if(.not.(present(option))) then
   ioption=0
 else
   ioption=option
 end if

 if(ioption/=0.and.paw_dmft%use_dmft==1) then
   message = ' option argument value of this routines should be 0 if usedmft=1. '
   MSG_ERROR(message)
 end if
 if(paw_dmft%use_dmft/=0) then
   nbandc1=(paw_dmft%mbandc-1)*paw_dmft%use_dmft+1
 else
   nbandc1=1
 end if
 use_nondiag_occup_dmft=0

 if(dtset%nspinor==2.and.paw_dmft%use_dmft==1) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' mkrho : ERROR -',ch10,&
&   '  nspinor argument value of this routines should be 1 if usedmft=1. '
   call wrtout(std_out,message,'COLL')
!  call leave_new('COLL')
 end if

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
 if (mpi_enreg%paral_spin==0) then
   ispinor_index=1;jspinor_index=1
   nspinor1TreatedByThisProc=.true.
   nspinor2TreatedByThisProc=(dtset%nspinor==2)
 else
   ispinor_index=mpi_enreg%me_spin+1;jspinor_index=3-ispinor_index
   nspinor1TreatedByThisProc=(mpi_enreg%me_spin==0)
   nspinor2TreatedByThisProc=(mpi_enreg%me_spin==1)
 end if

!Set local variable which depend on option argument

!nalpha*nbeta is the number of element of the kinetic energy density tensor
!to be computed in the irreducible Brillouin Zone (BZ) to get the result in the full BZ.
!In case of electron density calculation, nalpha=nbeta=1
 select case (ioption)
   case (0)
     nalpha = 1
     nbeta = 1
   case (1)
     nalpha = 3
     nbeta = 1
     ABI_ALLOCATE(taur_alphabeta,(dtset%nfft,dtset%nspden,3,1))
   case (2)
     nalpha = 3
     nbeta = 3
     ABI_ALLOCATE(taur_alphabeta,(dtset%nfft,dtset%nspden,3,3))
     case default
     MSG_BUG(' ioption argument value should be 0,1 or 2.')
 end select

!Init me
 call xme_init(mpi_enreg,me)

!Init master
 call xmaster_init(mpi_enreg,master)

!zero the charge density array in real space
 do ispden=1,dtset%nspden
!  $OMP PARALLEL DO PRIVATE(ifft) SHARED(dtset%nfft,rhor,zero,ispden)
   do ifft=1,dtset%nfft
     rhor(ifft,ispden)=zero
   end do
!  $OMP END PARALLEL DO
 end do

!WVL - Branching with a separate mkrho procedure
!in wavelet.
 if (dtset%usewvl == 1) then
   select case(ioption)
     case (0)
       call wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wfs, wvl)
       return
     case (1)
!      call wvl_mkrho(dtset, mpi_enreg, occ, rhor, wfs)
       message = ' Sorry, kinetic energy density (taur) is not yet implemented in wavelet formalism.'
       MSG_ERROR(message)
     case (2)
!      call wvl_mkrho(dtset, mpi_enreg, occ, rhor, wfs)
       message = '  Sorry, kinetic energy density tensor (taur_(alpha,beta)) is not yet implemented in wavelet formalism.'
       MSG_BUG(message)
   end select
 end if
!WVL - Following is done in plane waves.

!start loop over alpha and beta

 do alpha=1,nalpha
   do beta=1,nbeta

!    dtset%mkmem==0 means wf and kg info on disk file
     if (dtset%mkmem==0) then
!      Skip header of wffnow
       call hdr_skip(wffnow,ierr)
!      Define offsets, in case of MPI I/O
       formeig=0
       call xdefineOff(formeig,wffnow,mpi_enreg,dtset%nband,npwarr,dtset%nspinor,dtset%nsppol,dtset%nkpt)

       mcg_disk=dtset%mpw*my_nspinor*dtset%mband
       ABI_ALLOCATE(cg_disk,(2,mcg_disk))
     end if

!    start loop over spin and k points
     bdtot_index=0
     icg=0

!    n4,n5,n6 are FFT dimensions, modified to avoir cache trashing
     n1 = dtset%ngfft(1) ; n2 = dtset%ngfft(2) ; n3 = dtset%ngfft(3)
     n4 = dtset%ngfft(4) ; n5 = dtset%ngfft(5) ; n6 = dtset%ngfft(6)
     ABI_ALLOCATE(cwavef,(2,dtset%mpw,my_nspinor))
     ABI_ALLOCATE(rhoaug,(n4,n5,n6))
     ABI_ALLOCATE(wfraug,(2,n4,n5,n6))
     ABI_ALLOCATE(cwavefb,(2,dtset%mpw*my_nspinor*paw_dmft%use_dmft))
     if(dtset%nspden==4) then
       ABI_ALLOCATE(rhoaug_up,(n4,n5,n6))
       ABI_ALLOCATE(rhoaug_down,(n4,n5,n6))
       ABI_ALLOCATE(rhoaug_mx,(n4,n5,n6))
       ABI_ALLOCATE(rhoaug_my,(n4,n5,n6))
       rhoaug_up(:,:,:)=zero
       rhoaug_down(:,:,:)=zero
       rhoaug_mx(:,:,:)=zero
       rhoaug_my(:,:,:)=zero
     end if

     do isppol=1,dtset%nsppol

!      Rewind the kpgsph data file on unit unkg
       if (dtset%mkmem==0) rewind (unkg)
       ikg=0

       rhoaug(:,:,:)=0.0_dp
       do ikpt=1,dtset%nkpt

         nband_k = dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
         npw_k=npwarr(ikpt)
         istwf_k = dtset%istwfk(ikpt)

         if(mpi_enreg%paral_compil_kpt==1)then
           if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) then
             bdtot_index=bdtot_index+nband_k
             cycle
           end if
         end if

         ABI_ALLOCATE(gbound,(2*dtset%mgfft+8,2))
         ABI_ALLOCATE(kg_k,(3,npw_k))

!        Do i/o as needed
         if (dtset%mkmem==0) then

           nsp=dtset%nspinor
           call rdnpw(ikpt,isppol,nband_k,npw_k,nsp,0,unkg)

!          Read k+g data
           read (unkg) kg_k(1:3,1:npw_k)

           call sphereboundary(gbound,istwf_k,kg_k,dtset%mgfft,npw_k)

!          Read the wavefunction block for ikpt,isppol
           if((mpi_enreg%paralbd==0) .or. (mpi_enreg%paralbd>1)) tim_rwwf=5
           if(mpi_enreg%paralbd==1)tim_rwwf=12
           ABI_ALLOCATE(eig_dum,(dtset%mband))
           ABI_ALLOCATE(kg_dum,(3,0))
           ABI_ALLOCATE(occ_dum,(dtset%mband))
           call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,dtset%mband,mcg_disk,&
&           mpi_enreg,nband_k,nband_k,npw_k,my_nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
           ABI_DEALLOCATE(eig_dum)
           ABI_DEALLOCATE(kg_dum)
           ABI_DEALLOCATE(occ_dum)

         else

           kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
           call sphereboundary(gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
         end if ! dtset%mkmem==0

!        Loop over bands to fft and square for rho(r)
!        Shoulb be changed to treat bands by batch always

!        DEBUG
!        write(std_out,*)' mkrho : mpi_enreg%mode_para=',mpi_enreg%mode_para
!        ENDDEBUG

         if(mpi_enreg%mode_para /= 'b') then  ! Not yet parallelized on spinors
           do iband=1,nband_k
!            if(paw_dmft%use_dmft==1) then
!            write(std_out,*) 'iband  ',iband,occ(iband+bdtot_index),paw_dmft%occnd(iband,iband,ikpt,isppol)
!            else
!            write(std_out,*) 'iband  ',iband,occ(iband+bdtot_index)
!            endif

             do ibandc1=1,nbandc1 ! in case of DMFT

               if(mpi_enreg%paral_compil_kpt==1)then
                 if(mpi_enreg%paralbd>=1)then
                   if(mpi_enreg%proc_distrb(ikpt, iband, isppol) /= me) cycle
                 end if
               end if

!              Check if DMFT and only treat occupied states (check on occ.)
               if(paw_dmft%use_dmft == 1) then
                 iband1 = paw_dmft%include_bands(ibandc1)
                 if(paw_dmft%band_in(iband)) then
                   if(.not. paw_dmft%band_in(iband1))  stop
                   use_nondiag_occup_dmft = 1
                   locc_test = abs(paw_dmft%occnd(iband,iband1,ikpt,isppol))>tol8
!                  write(std_out,*) "mkrho,ikpt,iband,use_occnd",ikpt,iband
                 else
                   use_nondiag_occup_dmft = 0
                   locc_test = abs(occ(iband+bdtot_index))>tol8
                   if(ibandc1 /=1 .and. .not. paw_dmft%band_in(iband)) cycle
                 end if
               else
                 use_nondiag_occup_dmft = 0
                 locc_test = abs(occ(iband+bdtot_index))>tol8
               end if

               if (locc_test) then
!                Obtain Fourier transform in fft box and accumulate the density or the kinetic energy density
!                Not yet parallise on nspinor if paral_kgb non equal to 1
                 if(dtset%mkmem/=0)then
                   ipwsp=(iband-1)*npw_k*my_nspinor +icg
                   cwavef(:,1:npw_k,1)=cg(:,1+ipwsp:ipwsp+npw_k)
                   if (my_nspinor==2) cwavef(:,1:npw_k,2)=cg(:,ipwsp+npw_k+1:ipwsp+2*npw_k)

                 else
                   ipwsp=(iband-1)*npw_k*my_nspinor
                   cwavef(:,1:npw_k,1)=cg_disk(:,1+ipwsp:ipwsp+npw_k)
                   if (my_nspinor==2) cwavef(:,1:npw_k,2)=cg_disk(:,ipwsp+npw_k+1:ipwsp+2*npw_k)
                 end if

                 if(ioption==1)then
!                  Multiplication by 2pi i (k+G)_alpha
                   gp2pi1=gprimd(alpha,1)*two_pi ; gp2pi2=gprimd(alpha,2)*two_pi ; gp2pi3=gprimd(alpha,3)*two_pi
                   kpt_cart=gp2pi1*dtset%kptns(1,ikpt)+gp2pi2*dtset%kptns(2,ikpt)+gp2pi3*dtset%kptns(3,ikpt)
                   do ispinor=1,my_nspinor
                     do ipw=1,npw_k
                       kg_k_cart=gp2pi1*kg_k(1,ipw)+gp2pi2*kg_k(2,ipw)+gp2pi3*kg_k(3,ipw)+kpt_cart
                       cwftmp=-cwavef(2,ipw,ispinor)*kg_k_cart
                       cwavef(2,ipw,ispinor)=cwavef(1,ipw,ispinor)*kg_k_cart
                       cwavef(1,ipw,ispinor)=cwftmp
                     end do
                   end do
                 else if(ioption==2)then
                   message = ' Sorry, kinetic energy density tensor (taur_(alpha,beta)) is not yet implemented.'
                   MSG_ERROR(message)
                 end if

!                Non diag occupation in DMFT.
                 if(use_nondiag_occup_dmft==1) then
                   if(dtset%mkmem/=0)then
                     do ipw=1,npw_k*my_nspinor
                       cwavefb(1,ipw)=cg(1,ipw+(iband1-1)*npw_k*my_nspinor+icg)
                       cwavefb(2,ipw)=cg(2,ipw+(iband1-1)*npw_k*my_nspinor+icg)
                     end do
                   else
                     do ipw=1,npw_k*my_nspinor
                       cwavefb(1,ipw)=cg_disk(1,ipw+(iband1-1)*npw_k*my_nspinor)
                       cwavefb(2,ipw)=cg_disk(2,ipw+(iband1-1)*npw_k*my_nspinor)
                     end do
                   end if
                   weight=paw_dmft%occnd(iband,iband1,ikpt,isppol)*dtset%wtk(ikpt)/ucvol
                 else
                   weight=occ(iband+bdtot_index)*dtset%wtk(ikpt)/ucvol
                 end if

                 if((mpi_enreg%paralbd==0) .or. (mpi_enreg%paralbd>1)) tim_fourwf=3
                 if(mpi_enreg%paralbd==1)tim_fourwf=6

!                The same section of code is also found in vtowfk.F90 : should be rationalized !

                 call fourwf(1,rhoaug,cwavef(:,:,1),dummy,wfraug,gbound,gbound,&
&                 istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                 npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight,weight,&
&                 use_ndo=use_nondiag_occup_dmft,fofginb=cwavefb,&
&                 use_gpu_cuda=dtset%use_gpu_cuda)


                 if(dtset%nspinor==2)then
!                  DEBUG GZ !To obtain a x-directed magnetization(test)
!                  cwavef1(1,1:npw_k)=-cwavef(2,1:npw_k)
!                  cwavef1(2,1:npw_k)= cwavef(1,1:npw_k)
!                  ENDDEBUG

                   if(dtset%nspden==1) then

!                    We need only the total density : accumulation continues on top of rhoaug

                     call fourwf(1,rhoaug,cwavef(:,:,2),dummy,wfraug,gbound,gbound,&
&                     istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                     npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight,weight,&
&                     use_gpu_cuda=dtset%use_gpu_cuda)


                   else if(dtset%nspden==4) then

!                    Build the four components of rho. We use only norm quantities and, so fourwf.
!                    $\sum_{n} f_n \Psi^{* \alpha}_n \Psi^{\alpha}_n =\rho^{\alpha \alpha}$
!                    $\sum_{n} f_n (\Psi^{1}+\Psi^{2})^*_n (\Psi^{1}+\Psi^{2})_n=rho+m_x$
!                    $\sum_{n} f_n (\Psi^{1}-i \Psi^{2})^*_n (\Psi^{1}-i \Psi^{2})_n=rho+m_y$
                     ABI_ALLOCATE(cwavef_x,(2,npw_k))
                     ABI_ALLOCATE(cwavef_y,(2,npw_k))
!                    $(\Psi^{1}+\Psi^{2})$
                     cwavef_x(:,:)=cwavef(:,1:npw_k,1)+cwavef(:,1:npw_k,2)
!                    $(\Psi^{1}-i \Psi^{2})$
                     cwavef_y(1,:)=cwavef(1,1:npw_k,1)+cwavef(2,1:npw_k,2)
                     cwavef_y(2,:)=cwavef(2,1:npw_k,1)-cwavef(1,1:npw_k,2)
                     rhoaug_up(:,:,:)=rhoaug(:,:,:) !Already computed
                     call fourwf(1,rhoaug_down,cwavef(:,:,2),dummy,wfraug,gbound,gbound,&
&                     istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                     npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight,weight,&
&                     use_gpu_cuda=dtset%use_gpu_cuda)
                     call fourwf(1,rhoaug_mx,cwavef_x,dummy,wfraug,gbound,gbound,&
&                     istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                     npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight,weight,&
&                     use_gpu_cuda=dtset%use_gpu_cuda)
                     call fourwf(1,rhoaug_my,cwavef_y,dummy,wfraug,gbound,gbound,&
&                     istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                     npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight,weight,&
&                     use_gpu_cuda=dtset%use_gpu_cuda)

                     ABI_DEALLOCATE(cwavef_x)
                     ABI_DEALLOCATE(cwavef_y)

                   end if ! dtset%nspden/=4


                 end if
!                DEBUG
!                write(std_out,*)' ikpt, iband, rhoaug',ikpt,iband,rhoaug(1,1,1)
!                ENDDEBUG

               else
!                Accumulate the number of one-way 3D ffts skipped
                 nskip=nskip+1
               end if ! abs(occ(iband+bdtot_index))>tol8
!              End loop on iband
             end do ! iband1=1,(nband_k-1)*paw_dmft%use_dmft+1
           end do ! iband=1,nband_k

         else !mode_para==b
           if(minval(abs(mpi_enreg%proc_distrb(ikpt,:,isppol)-me))/=0) then
             cycle
           end if
           ikpt_this_proc=mpi_enreg%tab_kpt_distrib(ikpt)
           ABI_ALLOCATE(gs_hamk_local%gbound,(2*dtset%mgfft+8,2))
           if ((mpi_enreg%bandfft_kpt(ikpt_this_proc)%flag1_is_allocated==1).and.&
&           (ikpt_this_proc <= dtset%mkmem).and.(ikpt_this_proc/=0)) then
             gs_hamk_local%gbound(:,:)=mpi_enreg%bandfft_kpt(ikpt_this_proc)%gbound(:,:)
           else
             MSG_BUG('the bandfft tabs are not allocated !')
           end if
           gs_hamk_local%ngfft(:)=dtset%ngfft(:)
           gs_hamk_local%ucvol=ucvol
           gs_hamk_local%use_gpu_cuda=dtset%use_gpu_cuda

           nbdblock=nband_k/(mpi_enreg%nproc_band * mpi_enreg%bandpp)
           blocksize=nband_k/nbdblock
           if(allocated(cwavef))  then
             ABI_DEALLOCATE(cwavef)
           end if
           ABI_ALLOCATE(cwavef,(2,npw_k*blocksize,dtset%nspinor))
           if(ioption==1)  then
             ABI_ALLOCATE(kg_k_cart_block,(npw_k))
           end if
           ABI_ALLOCATE(occ_k,(nband_k))
           occ_k(:)=occ(bdtot_index+1:bdtot_index+nband_k)

           do iblock=1,nbdblock
             if (dtset%nspinor==1) then
               cwavef(:,1:npw_k*blocksize,1)=cg(:,1+(iblock-1)*npw_k*blocksize+icg:iblock*npw_k*blocksize+icg)
             else
               if (mpi_enreg%paral_spin==0) then
                 ishf=(iblock-1)*npw_k*my_nspinor*blocksize+icg
                 do ib=1,blocksize
                   cwavef(:,(ib-1)*npw_k+1:ib*npw_k,1)=cg(:,1+(2*ib-2)*npw_k+ishf:(2*ib-1)*npw_k+ishf)
                   cwavef(:,(ib-1)*npw_k+1:ib*npw_k,2)=cg(:,1+(2*ib-1)*npw_k+ishf:ib*2*npw_k+ishf)
                 end do
               else
                 ishf=(iblock-1)*npw_k*my_nspinor*blocksize+icg
                 do ib=1,blocksize
                   cwavef(:,(ib-1)*npw_k+1:ib*npw_k,ispinor_index)=&
&                   cg(:,1+(ib-1)*npw_k+ishf:ib*npw_k+ishf)
                   cwavef(:,(ib-1)*npw_k+1:ib*npw_k,jspinor_index)=zero
                 end do
                 call xsum_mpi(cwavef,mpi_enreg%comm_spin,ierr)
               end if
             end if
             if(ioption==1)then
!              Multiplication by 2pi i (k+G)_alpha
               gp2pi1=gprimd(alpha,1)*two_pi ; gp2pi2=gprimd(alpha,2)*two_pi ; gp2pi3=gprimd(alpha,3)*two_pi
               kpt_cart=gp2pi1*dtset%kptns(1,ikpt)+gp2pi2*dtset%kptns(2,ikpt)+gp2pi3*dtset%kptns(3,ikpt)
               kg_k_cart_block(1:npw_k)=gp2pi1*kg_k(1,1:npw_k)+gp2pi2*kg_k(2,1:npw_k)+gp2pi3*kg_k(3,1:npw_k)+kpt_cart
               do ib=1,blocksize
                 do ipw=1,npw_k
                   cwftmp=-cwavef(2,ipw+(ib-1)*npw_k,1)*kg_k_cart_block(ipw)
                   cwavef(2,ipw,1)=cwavef(1,ipw+(ib-1)*npw_k,1)*kg_k_cart_block(ipw)
                   cwavef(1,ipw,1)=cwftmp
                   if (my_nspinor==2) then
                     cwftmp=-cwavef(2,ipw+(ib-1)*npw_k,2)*kg_k_cart_block(ipw)
                     cwavef(2,ipw,2)=cwavef(1,ipw+(ib-1)*npw_k,2)*kg_k_cart_block(ipw)
                     cwavef(1,ipw,2)=cwftmp
                   end if
                 end do
               end do
             else if(ioption==2)then
               message = '  Sorry, kinetic energy density tensor (taur_(alpha,beta)) is not yet implemented.'
               MSG_ERROR(message)
             end if

             call timab(538,1,tsec)
             if (nspinor1TreatedByThisProc) then
               call prep_fourwf(rhoaug,blocksize,cwavef(:,:,1),wfraug,&
&               gs_hamk_local,iblock,ikpt,istwf_k,dtset%mgfft,mpi_enreg,&
&               nband_k,npw_k,n4,n5,n6,occ_k,dtset%paral_kgb,dtset%wtk(ikpt))
             end if
             call timab(538,2,tsec)
             if(dtset%nspinor==2)then
               if (dtset%nspden==1) then
                 if (nspinor2TreatedByThisProc) then
                   call prep_fourwf(rhoaug,blocksize,cwavef(:,:,2),wfraug,&
&                   gs_hamk_local,iblock,ikpt,istwf_k,dtset%mgfft,mpi_enreg,&
&                   nband_k,npw_k,n4,n5,n6,occ_k,dtset%paral_kgb,dtset%wtk(ikpt))
                 end if
               else if(dtset%nspden==4 ) then
                 ABI_ALLOCATE(cwavef_x,(2,npw_k*blocksize))
                 ABI_ALLOCATE(cwavef_y,(2,npw_k*blocksize))
                 cwavef_x(:,:)=cwavef(:,:,1)+cwavef(:,:,2)
                 cwavef_y(1,:)=cwavef(1,:,1)+cwavef(2,:,2)
                 cwavef_y(2,:)=cwavef(2,:,1)-cwavef(1,:,2)
                 call timab(538,1,tsec)
                 if (nspinor1TreatedByThisProc) then
                   call prep_fourwf(rhoaug_down,blocksize,cwavef(:,:,2),wfraug,&
&                   gs_hamk_local,iblock,ikpt,istwf_k,dtset%mgfft,mpi_enreg,&
&                   nband_k,npw_k,n4,n5,n6,occ_k,dtset%paral_kgb,dtset%wtk(ikpt))
                 end if
                 if (nspinor2TreatedByThisProc) then
                   call prep_fourwf(rhoaug_mx,blocksize,cwavef_x,wfraug,&
&                   gs_hamk_local,iblock,ikpt,istwf_k,dtset%mgfft,mpi_enreg,&
&                   nband_k,npw_k,n4,n5,n6,occ_k,dtset%paral_kgb,dtset%wtk(ikpt))
                   call prep_fourwf(rhoaug_my,blocksize,cwavef_y,wfraug,&
&                   gs_hamk_local,iblock,ikpt,istwf_k,dtset%mgfft,mpi_enreg,&
&                   nband_k,npw_k,n4,n5,n6,occ_k,dtset%paral_kgb,dtset%wtk(ikpt))
                 end if
                 call timab(538,2,tsec)
                 ABI_DEALLOCATE(cwavef_x)
                 ABI_DEALLOCATE(cwavef_y)
               end if
             end if
           end do !iblock
           if(ioption==1)  then
             ABI_DEALLOCATE(kg_k_cart_block)
           end if
           if (allocated(cwavef))  then
             ABI_DEALLOCATE(cwavef)
           end if
           ABI_DEALLOCATE(occ_k)
           ABI_DEALLOCATE(gs_hamk_local%gbound)
         end if

         ABI_DEALLOCATE(gbound)
         ABI_DEALLOCATE(kg_k)

         bdtot_index=bdtot_index+nband_k

         if (dtset%mkmem/=0) then
           icg=icg+npw_k*my_nspinor*nband_k
           ikg=ikg+npw_k
         end if

!        End loop on ikpt:
       end do

       if(mpi_enreg%mode_para == 'b') then
         if (dtset%nspden==4) then
!          Sum the contribution of the band and of the FFT
           call xsum_mpi(rhoaug     ,mpi_enreg%commcart_3d, ierr)
           call xsum_mpi(rhoaug_down,mpi_enreg%commcart_3d, ierr)
           call xsum_mpi(rhoaug_mx ,mpi_enreg%commcart_3d, ierr)
           call xsum_mpi(rhoaug_my ,mpi_enreg%commcart_3d, ierr)
           rhoaug_up(:,:,:) = rhoaug(:,:,:)
         else
           call xsum_mpi(rhoaug,mpi_enreg%commcart_3d,ierr)
         end if
       end if

!      Write the number of one-way 3D ffts skipped until now
       if(mpi_enreg%paral_compil_kpt==0)then
         write(message, '(a,i8)' )' mkrho : number of one-way 3D ffts skipped in mkrho until now =',nskip
         call wrtout(std_out,message,'PERS')
       end if

!      Transfer density on augmented fft grid to normal fft grid in real space
!      Take also into account the spin, to place it correctly in rhor.
       if(dtset%nspden==1 .or. dtset%nspden==2) then
         call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug,1)
       else if(dtset%nspden==4) then
         ispden=1
         call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_up,1)
         ispden=2
         call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_mx,1)
         ispden=3
         call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_my,1)
         ispden=4
         call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_down,1)
         ABI_DEALLOCATE(rhoaug_up)
         ABI_DEALLOCATE(rhoaug_down)
         ABI_DEALLOCATE(rhoaug_mx)
         ABI_DEALLOCATE(rhoaug_my)
       end if

     end do !  isppol=1,dtset%nsppol

     if(mpi_enreg%paral_compil_kpt==1)then
       call timab(63,1,tsec)
       call leave_test()
       write(message,'(a)') 'mkrho: loop on k-points and spins done in parallel'
       call wrtout(std_out,message,'COLL')
       call timab(63,2,tsec)
     end if

     if(allocated(cwavef))  then
       ABI_DEALLOCATE(cwavef)
     end if
     ABI_DEALLOCATE(cwavefb)
     ABI_DEALLOCATE(rhoaug)
     ABI_DEALLOCATE(wfraug)
     if(dtset%mkmem==0) ABI_DEALLOCATE(cg_disk)

     if(mpi_enreg%paral_compil_kpt==1)then
!      Recreate full rhor on all proc.
       call timab(48,1,tsec)
       call timab(71,1,tsec)
       call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_kpt)
       call xsum_mpi(rhor,spaceComm,ierr)
       call timab(71,2,tsec)
       call timab(48,2,tsec)
     end if

!    DEBUG
!    write(std_out,*) 'mkrho : dtset%nfft,dtset%nsppol,dtset%nsym',dtset%nfft,dtset%nsppol,dtset%nsym
!    write(std_out,*) 'ngfft',ngfft
!    write(std_out,*) ' ir irrzon phnons '
!    do ipw=1,dtset%nfft,31
!    write(std_out,'(i5,2i5,2es16.8)' )ipw,irrzon(ipw,:,1),phnons(:,ipw,1)
!    end do
!    select case (ioption)
!    case(0)
!    write(std_out,*)' mkrho : density before symrhg'
!    case(1)
!    write(std_out,*)' mkrho(tau) : kinetic energy density before symrhg'
!    case(2)
!    write(std_out,*)' mkrho(tau_alphabeta) : element of kinetic energy density tensor before symtaug'
!    end select
!    if(ioption==1 .or. ioption==2) write(std_out,'(i5,i5)' )alpha,beta
!    do ipw=1,dtset%nfft,31
!    write(std_out,'(i5,es16.6)' )ipw,rhor(ipw,1)
!    end do
!    ENDDEBUG

     call timab(799,2,tsec)
     call timab(549,1,tsec)

     if(ioption==1 .or. ioption==2) then
       do ispden=1,dtset%nspden
!        $OMP PARALLEL DO PRIVATE(ifft) &
!        $OMP&SHARED(dtset%nfft,taur_alphabeta,rhor,alpha,beta)
         do ifft=1,dtset%nfft
           taur_alphabeta(ifft,ispden,alpha,beta) = rhor(ifft,ispden)
         end do
!        $OMP END PARALLEL DO
       end do
     end if

   end do !  beta=1,nbeta
 end do !  alpha=1,nalpha


!Compute the trace over the kinetic energy denisty tensor
!i.e. Sum of the 3 diagonal elements.
 if(ioption==1)then
!  zero rhor array in real space
   do ispden=1,dtset%nspden
!    $OMP PARALLEL DO PRIVATE(ifft) &
!    $OMP&SHARED(dtset%nfft,rhor,ispden,zero)
     do ifft=1,dtset%nfft
       rhor(ifft,ispden)=zero
     end do
!    $OMP END PARALLEL DO
   end do
   do alpha = 1, nalpha
     do ispden=1,dtset%nspden
!      $OMP PARALLEL DO PRIVATE(ifft) &
!      $OMP&SHARED(dtset%nfft,rhor,ispden,zero)
       do ifft=1,dtset%nfft
         rhor(ifft,ispden) = rhor(ifft,ispden) + taur_alphabeta(ifft,ispden,alpha,1)
       end do
!      $OMP END PARALLEL DO
     end do
   end do
 end if

 nfftot=dtset%ngfft(1) * dtset%ngfft(2) * dtset%ngfft(3)

 select case (ioption)
   case(0,1)
     call symrhg(1,gprimd,irrzon,mpi_enreg,dtset%nfft,nfftot,dtset%ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
     dtset%paral_kgb,phnons,rhog,rhor,rprimd,dtset%symafm,dtset%symrel)
     if(ioption==1)then
!      $OMP PARALLEL DO PRIVATE(ifft) &
!      $OMP&SHARED(dtset%nfft,dtset%nspden,ispden,rhog,rhor,zero)
       do ifft=1,dtset%nfft
         do ispden=1,dtset%nspden
           rhor(ifft,ispden)=1.0d0/2.0d0*rhor(ifft,ispden)
         end do
         rhog(:,ifft)=1.0d0/2.0d0*rhog(:,ifft)
       end do
!      $OMP END PARALLEL DO
     end if
   case(2)
     message = ' Sorry, kinetic energy density tensor (taur_(alpha,beta)) is not yet implemented.'
     MSG_BUG(message)

!    call symtaug(1,gprimd,irrzon,mpi_enreg,dtset%nfft,nfftot,dtset%ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
!    dtset%paral_kgb,phnons,rhog,rhor,rprimd,dtset%symafm,dtset%symrel)
 end select

 call timab(549,2,tsec)

!DEBUG
!select case (ioption)
!case(0)
!write(std_out,*)' mkrho : density after symrhg'
!case(1)
!write(std_out,*)' mkrho(tau) : kinetic energy density after symrhg'
!case(2)
!write(std_out,*)' mkrho(tau_alphabeta) : element of kinetic energy density tensor after symtaug'
!end select
!do ipw=1,dtset%nfft,31
!write(std_out,'(i5,es16.6)' )ipw,rhor(ipw,1)
!end do
!ENDDEBUG

!We now have both rho(r) and rho(G), symmetrized, and if dtset%nsppol=2
!we also have the spin-up density, symmetrized, in rhor(:,2).
!In case of non collinear magnetism, we have rho,mx,my,mz. No symmetry is applied

!Debugging output
!select case (ioption)
!case(0)
!write(std_out,*)' Debugging from mkrho: rhog values'
!case(1)
!write(std_out,*)' Debugging from mkrho(tau): taug values'
!case(1)
!write(std_out,*)' Debugging from mkrho(tau_alphabeta): taug_alphabeta values'
!end select
!do ipw=1,dtset%nfft
!if (abs(rhog(1,ipw))>1.d-09.or.abs(rhog(2,ipw))>1.d-09) then
!write(std_out,2000) ipw,rhog(1,ipw),rhog(2,ipw)
!end if
!2000  format(i10,1p,2e15.5)
!end do

!DEBUG
!write(std_out,*)' rhor after sym',rhor(1,:)
!write(std_out,*)'nsym',nsym
!ENDDEBUG

 call timab(799,1,tsec)

 if(ioption==1 .or. ioption==2)  then
   ABI_DEALLOCATE(taur_alphabeta)
 end if

!Find and print minimum and maximum total electron density
!(or total kinetic energy density, or total element of kinetic energy density tensor) and locations
 write(message,'(a)') 'mkrho: echo density (plane-wave part only) '
 call wrtout(std_out,message,'COLL')
 call prtrhomxmn(std_out,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%nspden,1,rhor,optrhor=ioption,ucvol=ucvol)

 call timab(799,2,tsec)
 call timab(790+tim_mkrho,2,tsec)

 DBG_EXIT("COLL")

end subroutine mkrho
!!***
