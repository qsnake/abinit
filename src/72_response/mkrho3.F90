!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkrho3
!! NAME
!! mkrho3
!!
!! FUNCTION
!! Compute RF charge density rho1(r) and rho1(G) in electrons/bohr**3
!! from input RF and GS wavefunctions, band occupations, and k point weights.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, LSI, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=wf in G space
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=first-order wf in G space
!!  cplex=1 if rhor1 is real, 2 if rhor1 is complex
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the
!!    storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates, GS data.
!!  kg1(3,mpw1*mkmem1)=reduced planewave coordinates, RF data.
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=maximum number of k points in core memory for GS data
!!  mk1mem=maximum number of k points in core memory for RF data
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum allowed value for npw (GS wfs)
!!  mpw1=maximum allowed value for npw1 (RF data)
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands to be included in summation
!!   at each k point for each spin channel.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt_rbz=number of k points in the reduced Brillouin zone
!!  npwarr(nkpt_rbz)=number of planewaves and boundary planewaves at k points
!!  npwar1(nkpt_rbz)=number of planewaves and boundary planewaves at k+q points
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in group (at least 1 for identity)
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation numbers for each band
!!   (usually 2.0) at each k point of the reduced Brillouin zone
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  rprimd(3,3)=dimensional real space primitive translations
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!  ucvol=unit cell volume (Bohr**3)
!!  unkg=unit number for (k+G) sphere data file for GS data
!!  unkg1=unit number for (k+G+q) sphere data file for RF data
!!  wffnow=struct info for current RF wf disk file
!!  wfftgs=struct info for GS wf disk file
!!  wtk_rbz(nkpt_rbz)=k point weights (they sum to 1.0).
!!
!! OUTPUT
!!  rhog1(2,nfft)=total electron density in G space
!!  rhor1(cplex*nfft,nspden)=electron density in r space
!!   (if spin polarized, array contains total density in first half and
!!    spin-up density in second half)
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      fftpac,fourwf,hdr_skip,leave_test,rdnpw,rwwf,sphereboundary,symrhg
!!      timab,wrtout,xcomm_world,xdefineoff,xmaster_init,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkrho3(cg,cg1,cplex,gprimd,irrzon,istwfk_rbz,&
& kg,kg1,mband,mgfft,mkmem,mk1mem,mpi_enreg,mpw,mpw1,nband_rbz,&
& nfft,ngfft,nkpt_rbz,npwarr,npwar1,nspden,nspinor,nsppol,nsym,&
& occ_rbz,paral_kgb,phnons,rhog1,rhor1,rprimd,symafm,symrel,&
& ucvol,unkg,unkg1,wffnow,wfftgs,wtk_rbz)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkrho3'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!scalars
 integer,intent(in) :: cplex,mband,mgfft,mk1mem,mkmem,mpw,mpw1,nfft,nkpt_rbz
 integer,intent(in) :: nspden,nspinor,nsppol,nsym,paral_kgb,unkg,unkg1
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wffile_type),intent(inout) :: wffnow,wfftgs
!arrays
 integer,intent(in) :: irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
 integer,intent(in) :: istwfk_rbz(nkpt_rbz),kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol),ngfft(18),npwar1(nkpt_rbz)
 integer,intent(in) :: npwarr(nkpt_rbz),symafm(nsym),symrel(3,3,nsym)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol),gprimd(3,3)
 real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: rprimd(3,3),wtk_rbz(nkpt_rbz)
 real(dp),intent(out) :: rhog1(2,nfft),rhor1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,save :: nskip=0
 integer :: bdtot_index,formeig,i1,i2,i3,iband,icg,icg1,ierr,ifft,ikg
 integer :: ikg1,ikpt,ipw,ispden,ispinor,isppol,istwf_k,master
 integer :: mcg1_disk,mcg_disk,me,n1,n2,n3,n4,n5,n6,nband_k,npw1_k
 integer :: npw_k,nsp,sender,spaceworld,tim_fourwf,tim_rwwf
 real(dp) :: im0,im1,re0,re1,weight
 character(len=500) :: message
!arrays
 integer,allocatable :: gbound(:,:),gbound1(:,:),kg1_k(:,:),kg_dum(:,:)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cg1_disk(:,:),cg_disk(:,:),cwavef(:,:),cwavef1(:,:)
 real(dp),allocatable :: dummy(:,:),eig_dum(:),occ_dum(:),rhoaug(:,:,:)
 real(dp),allocatable :: rhoaug1(:,:,:),wfraug(:,:,:,:),wfraug1(:,:,:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' mkrho3 : enter '
!if(.true.)stop
!ENDDEBUG

 if(nspden==4)then
!  NOTE : see mkrho for the modifications needed for non-collinear treatment
   write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&   ' mkrho3 : ERROR -',ch10,&
&   '  Linear-response calculations are not yet possible with nspden=4',ch10,&
&   ' Action : modify value of nspden in input file.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

!Init spaceworld
 call xcomm_world(mpi_enreg,spaceworld)
!Init me
 call xme_init(mpi_enreg,me)
!Init master
 call xmaster_init(mpi_enreg,master)
 sender=master

!zero the charge density array in real space
 do ispden=1,nspden
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(cplex,nfft,rhor1)
   do ifft=1,cplex*nfft
     rhor1(ifft,ispden)=zero
   end do
!  $OMP END PARALLEL DO
 end do

!mkmem==0 means wf and kg info on disk file
 if (mkmem==0) then

!  Read header of wfftgs
   call hdr_skip(wfftgs,ierr)

!  Define offsets, in case of MPI I/O
   formeig=0
   call xdefineOff(formeig,wfftgs,mpi_enreg,nband_rbz,npwarr,nspinor,nsppol,nkpt_rbz)

   mcg_disk=mpw*nspinor*mband
   ABI_ALLOCATE(cg_disk,(2,mcg_disk))

 end if

!mk1mem==0 means RF wf and kgq info on disk file
 if (mk1mem==0) then

!  Read header of wffnow
   call hdr_skip(wffnow,ierr)

!  Define offsets, in case of MPI I/O
   formeig=1
   call xdefineOff(formeig,wffnow,mpi_enreg,nband_rbz,npwar1,nspinor,nsppol,nkpt_rbz)

   mcg1_disk=mpw1*nspinor*mband
!  One should not allocate an array as large as cg1_disk here
   ABI_ALLOCATE(cg1_disk,(2,mcg1_disk))

 end if

!start loop over spin and k points
 bdtot_index=0 ; icg=0 ; icg1=0

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
!n4,n5,n6 are FFT dimensions, modified to avoid cache trashing
 n4=ngfft(4) ; n5=ngfft(5) ; n6=ngfft(6)

!Note that the dimensioning of cwavef and cwavef1 do not include nspinor
 ABI_ALLOCATE(cwavef,(2,mpw))
 ABI_ALLOCATE(cwavef1,(2,mpw1))
!Actually, rhoaug is not needed, except for strong dimensioning requirement
 ABI_ALLOCATE(dummy,(2,1))
 ABI_ALLOCATE(rhoaug,(n4,n5,n6))
 ABI_ALLOCATE(rhoaug1,(cplex*n4,n5,n6))
 ABI_ALLOCATE(wfraug,(2,n4,n5,n6))
 ABI_ALLOCATE(wfraug1,(2,n4,n5,n6))

 do isppol=1,nsppol

!  Eventually rewind the kpgsph data file on unit unkg and/or unkg1
   if (mkmem==0) rewind (unkg)
   if (mk1mem==0) rewind (unkg1)
   ikg=0 ; ikg1=0


   rhoaug1(:,:,:)=zero

   do ikpt=1,nkpt_rbz

     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt)
     npw1_k=npwar1(ikpt)

     if(mpi_enreg%paral_compil_kpt==1)then
!      BEGIN TF_CHANGES
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) &
&       -me))/=0) then
!        END TF_CHANGES
         bdtot_index=bdtot_index+nband_k
         cycle
       end if
     end if

     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     ABI_ALLOCATE(kg_k,(3,npw_k))
     ABI_ALLOCATE(gbound1,(2*mgfft+8,2))
     ABI_ALLOCATE(kg1_k,(3,npw1_k))

!    Do GS i/o as needed
     if (mkmem==0) then

       nsp=nspinor
       call rdnpw(ikpt,isppol,nband_k,npw_k,nsp,0,unkg)

!      Read k+g data
       read (unkg) kg_k(:,1:npw_k)

       call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

!      Read the wavefunction block for ikpt,isppol
       tim_rwwf=15
       ABI_ALLOCATE(eig_dum,(mband))
       ABI_ALLOCATE(kg_dum,(3,0))
       ABI_ALLOCATE(occ_dum,(mband))
       call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
&       npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wfftgs)
       ABI_DEALLOCATE(eig_dum)
       ABI_DEALLOCATE(kg_dum)
       ABI_DEALLOCATE(occ_dum)


     else

       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

!      End test for mkmem==0
     end if

!    Do RF i/o as needed
     if (mk1mem==0) then

       nsp=nspinor
       call rdnpw(ikpt,isppol,nband_k,npw1_k,nsp,0,unkg1)

!      Read k+g data
       read (unkg1) kg1_k(:,1:npw1_k)

       call sphereboundary(gbound1,istwf_k,kg1_k,mgfft,npw1_k)

!      Read the wavefunction block for ikpt,isppol
       tim_rwwf=15
       ABI_ALLOCATE(eig_dum,(2*mband*mband))
       ABI_ALLOCATE(kg_dum,(3,0))
       ABI_ALLOCATE(occ_dum,(mband))
       call rwwf(cg1_disk,eig_dum,1,0,0,ikpt,isppol,kg_dum,mband,mcg1_disk,mpi_enreg,nband_k,nband_k,&
&       npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
       ABI_DEALLOCATE(eig_dum)
       ABI_DEALLOCATE(kg_dum)
       ABI_DEALLOCATE(occ_dum)

     else

       kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
       call sphereboundary(gbound1,istwf_k,kg1_k,mgfft,npw1_k)

!      End test for mk1mem==0
     end if

!    Loop over bands to fft and square for rho(r)
     do iband=1,nband_k

       if(mpi_enreg%paral_compil_kpt==1)then
!        BEGIN TF_CHANGES
         if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me) then
!          END TF_CHANGES
           cycle
         end if
       end if

!      Only treat occupied states
       if (abs(occ_rbz(iband+bdtot_index))>tol8) then

!        Treat separately the two spinor components
         do ispinor=1,nspinor

!          Obtain Fourier transform in fft box and accumulate the density
           if(mkmem/=0)then
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(cg,cwavef,iband,icg,ispinor,npw_k,nspinor)
             do ipw=1,npw_k
               cwavef(1,ipw)=cg(1,ipw+(ispinor-1)*npw_k+(iband-1)*npw_k*nspinor+icg)
               cwavef(2,ipw)=cg(2,ipw+(ispinor-1)*npw_k+(iband-1)*npw_k*nspinor+icg)
             end do
!            $OMP END PARALLEL DO
           else
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(cg_disk,cwavef,iband,ispinor,npw_k,nspinor)
             do ipw=1,npw_k
               cwavef(1,ipw)=cg_disk(1,ipw+(ispinor-1)*npw_k+(iband-1)*npw_k*nspinor)
               cwavef(2,ipw)=cg_disk(2,ipw+(ispinor-1)*npw_k+(iband-1)*npw_k*nspinor)
             end do
!            $OMP END PARALLEL DO
           end if
           if(mk1mem/=0)then
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(cg1,cwavef1,iband,icg1,ispinor,npw1_k,nspinor)
             do ipw=1,npw1_k
               cwavef1(1,ipw)=cg1(1,ipw+(ispinor-1)*npw1_k+(iband-1)*npw1_k*nspinor+icg1)
               cwavef1(2,ipw)=cg1(2,ipw+(ispinor-1)*npw1_k+(iband-1)*npw1_k*nspinor+icg1)
             end do
!            $OMP END PARALLEL DO
           else
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(cg1_disk,cwavef1,iband,ispinor,npw1_k,nspinor)
             do ipw=1,npw1_k
               cwavef1(1,ipw)=cg1_disk(1,ipw+(ispinor-1)*npw1_k+(iband-1)*npw1_k*nspinor)
               cwavef1(2,ipw)=cg1_disk(2,ipw+(ispinor-1)*npw1_k+(iband-1)*npw1_k*nspinor)
             end do
!            $OMP END PARALLEL DO
           end if

!          In these two calls, rhoaug, rhoaug1 and weight are dummy variables, and are
!          not modified
           tim_fourwf=7
           call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&           istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,1,n4,n5,n6,0,paral_kgb,tim_fourwf,weight,weight)
           call fourwf(cplex,rhoaug1,cwavef1,dummy,wfraug1,gbound1,gbound1,&
&           istwf_k,kg1_k,kg1_k,mgfft,mpi_enreg,1,ngfft,npw1_k,1,n4,n5,n6,0,&
&           paral_kgb,tim_fourwf,weight,weight)

!          Compute the weight, note that the factor 2 is
!          not the spin factor (see Eq.44 of PRB55,10337 (1997))
           weight=two*occ_rbz(iband+bdtot_index)*wtk_rbz(ikpt)/ucvol

!          Accumulate density
           if(cplex==2)then
!            $OMP PARALLEL DO PRIVATE(i1,i2,i3,im0,im1,re0,re1) &
!            $OMP&SHARED(n1,n2,n3,rhoaug1,weight,wfraug,wfraug1)
             do i3=1,n3
               do i2=1,n2
                 do i1=1,n1
                   re0=wfraug(1,i1,i2,i3)  ; im0=wfraug(2,i1,i2,i3)
                   re1=wfraug1(1,i1,i2,i3) ; im1=wfraug1(2,i1,i2,i3)
                   rhoaug1(2*i1-1,i2,i3)=rhoaug1(2*i1-1,i2,i3)+weight*(re0*re1+im0*im1)
                   rhoaug1(2*i1  ,i2,i3)=rhoaug1(2*i1  ,i2,i3)+weight*(re0*im1-im0*re1)
                 end do
               end do
             end do
!            $OMP END PARALLEL DO
           else
!            $OMP PARALLEL DO PRIVATE(i1,i2,i3) &
!            $OMP&SHARED(n1,n2,n3,rhoaug1,weight,wfraug,wfraug1)
             do i3=1,n3
               do i2=1,n2
                 do i1=1,n1
                   rhoaug1(i1,i2,i3)=rhoaug1(i1,i2,i3)+&
&                   weight*( wfraug(1,i1,i2,i3)*wfraug1(1,i1,i2,i3) &
&                   +wfraug(2,i1,i2,i3)*wfraug1(2,i1,i2,i3)  )
                 end do
               end do
             end do
!            $OMP END PARALLEL DO
           end if

         end do ! ispinor=1,nspinor

       else  ! if the state is not occupied

!        Accumulate the number of one-way 3D ffts skipped
         nskip=nskip+1

       end if ! end of the state occupied or not

!      End loop on iband
     end do

     ABI_DEALLOCATE(gbound)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(gbound1)
     ABI_DEALLOCATE(kg1_k)

     bdtot_index=bdtot_index+nband_k

     if (mkmem/=0) then
       icg=icg+npw_k*nband_k
       ikg=ikg+npw_k
     end if
     if (mk1mem/=0) then
       icg1=icg1+npw1_k*nband_k
       ikg1=ikg1+npw1_k
     end if

!    End loop on ikpt:
   end do

!  Write the number of one-way 3D ffts skipped until now
   if(mpi_enreg%paral_compil_kpt==0)then
     write(message, '(a,i8)' )&
&     ' mkrho3 : number of one-way 3D ffts skipped in mkrho3 until now =',nskip
     call wrtout(std_out,message,'PERS')
   end if

!  DEBUG
!  write(std_out,*)' rhoaug ',rhoaug(1,1,1)
!  ENDDEBUG

!  Transfer density on augmented fft grid to normal fft grid in real space
!  Take also into account the spin, to place it correctly in rhor1.
!  Note the use of cplex
   call fftpac(isppol,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,ngfft,rhor1,rhoaug1,1)

!  DEBUG
!  write(std_out,*)' rhor1 ',rhor1(1,1)
!  ENDDEBUG

!  End loop over spins
 end do

 if(mpi_enreg%paral_compil_kpt==1)then
   call timab(63,1,tsec)
!  BEGIN TF_CHANGES
   call leave_test()
!  END TF_CHANGES
   write(message,*) 'mkrho3: loop on k-points and spins done in parallel'
   call wrtout(std_out,message,'COLL')
   call timab(63,2,tsec)
 end if

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(cwavef1)
 ABI_DEALLOCATE(dummy)
 ABI_DEALLOCATE(rhoaug)
 ABI_DEALLOCATE(rhoaug1)
 ABI_DEALLOCATE(wfraug)
 ABI_DEALLOCATE(wfraug1)
 if(mkmem==0)ABI_DEALLOCATE(cg_disk)
 if(mk1mem==0)ABI_DEALLOCATE(cg1_disk)

 if(mpi_enreg%paral_compil_kpt==1)then
!  Recreate full rhor1 on all proc.
   call timab(48,1,tsec)
   call timab(71,1,tsec)
   call xsum_mpi(rhor1,spaceworld,ierr)
   call timab(71,2,tsec)
   call timab(48,2,tsec)
 end if

!DEBUG
!write(std_out,*) 'mkrho3 : nfft,nspden,nsym',nfft,nspden,nsym
!write(std_out,*) 'ngfft',ngfft
!write(std_out,*) ' ir irrzon phnons '
!do ipw=1,nfft,31
!write(std_out,'(i5,2i5,2es16.8)' )ipw,irrzon(ipw,:,1),phnons(:,ipw,1)
!end do
!ENDDEBUG

 call symrhg(cplex,gprimd,irrzon,mpi_enreg,nfft,nfft,ngfft,nspden,nsppol,nsym,paral_kgb,phnons,&
& rhog1,rhor1,rprimd,symafm,symrel)

!We now have both rho(r) and rho(G), symmetrized, and if nsppol=2
!we also have the spin-up density, symmetrized, in rhor1(:,2).

!Debugging output
!write(std_out,*)' Debugging from mkrho3: rhog1 values'
!do ipw=1,nfft
!if (abs(rhog1(1,ipw))>1.d-09.or.abs(rhog1(2,ipw))>1.d-09)
!& then
!write(std_out,2000) ipw,rhog1(1,ipw),rhog1(2,ipw)
!end if
!2000  format(i10,1p,2e15.5)
!end do

!DEBUG
!write(std_out,*)' rhor1 ',rhor1(1,1)
!ENDDEBUG

!DEBUG
!write(std_out,*)' mkrho3 : exit '
!if(.true.)stop
!ENDDEBUG

end subroutine mkrho3
!!***
