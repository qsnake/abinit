!{\src2tex{textfont=tt}}
!!****f* ABINIT/resp3dte
!! NAME
!! resp3dte
!!
!! FUNCTION
!! Compute the linear response part to the 3dte
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see loop3dte.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol) = array for planewave
!!                                          coefficients of wavefunctions
!!  cg1 = first-order wavefunction relative to the perturbations i1pert
!!  cg3 = first-order wavefunction relative to the perturbations i3pert
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!          if 2, COMPLEX
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1)
!!  i1dir,i2dir,i3dir=directions of the corresponding perturbations
!!  i1pert,i2pert,i3pert = type of perturbation that has to be computed
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  mband = maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem = maximum number of k points which can fit in core memory
!!  mk1mem = maximum number of k points for first-order WF
!!           which can fit in core memory
!!  mpert =maximum number of ipert
!!  mpi_enreg=MPI-parallelisation information
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  natom = number of atoms in unit cell
!!  nfft  = (effective) number of FFT grid points (for this processor)
!!  nkpt  = number of k points
!!  nspden = number of spin-density components
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  npwarr(nkpt) = array holding npw for each k point
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  ucvol = unit cell volume (bohr^3)
!!  vtrial1(cplex*nfft,nspden)=firs-order local potential
!!  xred(3,natom) = reduced atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= spherical harmonics for
!!       each G and k point
!!
!! OUTPUT
!!  d3lo(2,3,mpert,3,mpert,3,mpert) = matrix of the 3DTEs
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! PARENTS
!!      loop3dte
!!
!! CHILDREN
!!      fftpac,fourwf,mkffnl,mkkpg,nonlop,sphereboundary,status,xcomm_world
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine resp3dte(atindx,atindx1,cg,cg1,cg3,cplex,dtfil,dtset,d3lo,&
& gmet,gprimd,i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,&
& kg,mband,mgfft,mkmem,mk1mem,&
& mpert,mpi_enreg,mpsang,mpw,natom,nfft,nkpt,nspden,nspinor,nsppol,&
& npwarr,occ,ph1d,psps,rmet,ucvol,vtrial1,xred,ylm)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
!***********************************************************************

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'resp3dte'
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,mband,mgfft
 integer,intent(in) :: mk1mem,mkmem,mpert,mpsang,mpw,natom,nfft,nkpt,nspden
 integer,intent(in) :: nspinor,nsppol
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),kg(3,mpw*mkmem)
 integer,intent(in) :: npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol),gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),rmet(3,3)
 real(dp),intent(in) :: xred(3,natom),ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(inout) :: vtrial1(cplex*nfft,nspden)
 real(dp),intent(out) :: d3lo(2,3,mpert,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=52
 integer :: bantot,choice,counter,cpopt,dimffnl,iband,icg0,ider,ierr,iexit
 integer :: ii,ikg,ikpt,ilm,ilmn,ipw,isppol,istwf_k,itypat,jband,jj
 integer :: matblk_der,me,n1,n2,n3,n4,n5,n6,natom_der,nband_k,nkpg,nnlout
 integer :: npw1_k,npw_k,nqpg,ntypat_der,option,paw_opt
 integer :: shift1,shift2,shift3,signs,spaceComm,tim_fourwf,tim_nonlop
 real(dp) :: arg,dot1i,dot1r,dot2i,dot2r,doti,dotr,dum,lagi,lagr,sumi,sumr
 real(dp) :: weight
!arrays
 integer :: nattyp_der(1),nloalg_der(5)
 integer,allocatable :: gbound(:,:),indlmn_typ(:,:,:),kg_k(:,:)
 real(dp) :: buffer(2),enlout(3),kpq(3),kpt(3)
 real(dp) :: dum_sij(1,1),dum_svectout(1,1)
 real(dp) :: phkxredin(2,1),phkxredout(2,1),tmp(2),xred_der(3),ylmgr_dum(1)
 real(dp),allocatable :: cwave0(:,:),cwavef3(:,:)
 real(dp),allocatable :: ekb_typ(:,:,:),ffnl(:,:,:,:),ffnlk(:,:,:,:)
 real(dp),allocatable :: ffnlkq(:,:,:,:),gh0(:,:),gh1(:,:),gvnl(:,:),kpg_k(:,:)
 real(dp),allocatable :: ph1d_der(:,:),ph3din(:,:,:),ph3dout(:,:,:),qpg_k(:,:)
 real(dp),allocatable :: vlocal1(:,:,:),wfraug(:,:,:,:),ylm_k(:,:)
 type(cprj_type) :: cprj_dum(1,1)

!***********************************************************************

!DEBUG
!write(std_out,*)'resp3dte : enter'
!write(std_out,*)xred(:,:)
!call leave_new('COLL')
!ENDDEBUG

 me = mpi_enreg%me
 if (mpi_enreg%paral_compil_kpt == 1) then
!  BEGIN TF_CHANGES
   call xcomm_world(mpi_enreg,spaceComm)
!  END TF_CHANGES
 end if

 call status(0,dtfil%filstat,iexit,level,'enter         ')

 bantot = 0
 icg0 = 0
 option = 2
 xred_der(:) = zero
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)

 ABI_ALLOCATE(vlocal1,(cplex*n4,n5,n6))
 ABI_ALLOCATE(wfraug,(2,n4,n5,n6))
 ABI_ALLOCATE(ph1d_der,(2,(2*n1+1)+(2*n2+1)+(2*n3+1)))

 if (i2pert <= natom) then
!  Store at the right place the 1d phases
   shift1=(atindx(i2pert)-1)*(2*n1+1)
   ph1d_der(:,1:2*n1+1)=ph1d(:,1+shift1:2*n1+1+shift1)
   shift2=(atindx(i2pert)-1)*(2*n2+1)+natom*(2*n1+1)
   ph1d_der(:,1+2*n1+1:2*n2+1+2*n1+1)=ph1d(:,1+shift2:2*n2+1+shift2)
   shift3=(atindx(i2pert)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
   ph1d_der(:,1+2*n1+1+2*n2+1:2*n3+1+2*n2+1+2*n1+1)=&
&   ph1d(:,1+shift3:2*n3+1+shift3)
   xred_der(:) = xred(:,i2pert)
 end if

 nloalg_der(:)=dtset%nloalg(:)
 nloalg_der(1)=-abs(dtset%nloalg(1))
 nloalg_der(4)=1


 sumr = 0._dp ; sumi = 0._dp


!Loop over spins

 tmp(:) = 0._dp
 do isppol = 1, nsppol

   call status(0,dtfil%filstat,iexit,level,'call fftpac   ')
   call fftpac(isppol,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,vtrial1,vlocal1,option)
!  Loop over k-points

   ikg = 0
   do ikpt = 1, nkpt

     if (mpi_enreg%paral_compil_kpt == 1) then
       if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:mband,1:dtset%nsppol) &
&       - mpi_enreg%me)) /= 0) cycle
     end if

     counter = 100*ikpt

     nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
     npw_k = npwarr(ikpt)
     npw1_k = npw_k              ! this has to be changed in case q /= 0
     istwf_k = dtset%istwfk(ikpt)

     ABI_ALLOCATE(cwave0,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gh0,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gvnl,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(cwavef3,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gh1,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gbound,(2*dtset%mgfft+8,2))
     ABI_ALLOCATE(kg_k,(3,npw_k))
     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     dimffnl=1
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
     ABI_ALLOCATE(ffnlk,(npw_k,dimffnl,psps%lmnmax,1))
     ABI_ALLOCATE(ffnlkq,(npw_k,dimffnl,psps%lmnmax,1))
     ABI_ALLOCATE(ph3din,(2,npw_k,1))
     ABI_ALLOCATE(ph3dout,(2,npw1_k,1))
     ABI_ALLOCATE(indlmn_typ,(6,psps%lmnmax,1))
     ABI_ALLOCATE(ekb_typ,(psps%dimekb,1,nspinor**2))

     kg_k(:,1:npw_k) = kg(:,1+ikg:npw_k+ikg)

     kpt(:) = dtset%kptns(:,ikpt)
     kpq(:) = dtset%kptns(:,ikpt)     ! In case of non zero q, kpt = kpt + q

     call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

     if (psps%useylm==1) then
       do ilm=1,mpsang*mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if


!    Compute (k+G) and (k+q+G) vectors (only if useylm=1)
     nkpg=0;if (i2pert<natom+1) nkpg=3*dtset%nloalg(5)
     nqpg=0;if (i2pert<natom+1) nqpg=3*dtset%nloalg(5)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     ABI_ALLOCATE(qpg_k,(npw1_k,nqpg))
     if (nkpg>0) call mkkpg(kg_k,kpg_k,kpt,nkpg,npw_k)
     if (nqpg>0) call mkkpg(kg_k,qpg_k,kpq,nqpg,npw1_k)

!    Compute nonlocal form factors ffnl at (k+G), for all atoms
     if (i2pert < natom + 1) then
       ider=0
       if (me==0) call status(counter,dtfil%filstat,iexit,level,'call mkffnl  ')
       call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,gmet,gprimd,ider,ider,&
&       psps%indlmn,kg_k,kpg_k,kpt,psps%lmnmax,&
&       psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&       npw_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&       psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
     end if


!    Compute here phkxred for kpt and kpq
     arg=two_pi*(kpt(1)*xred_der(1)+&
&     kpt(2)*xred_der(2)+kpt(3)*xred_der(3))
     phkxredin(1,1)=cos(arg)
     phkxredin(2,1)=sin(arg)
     arg=two_pi*(kpq(1)*xred_der(1)+kpq(2)*xred_der(2)+kpq(3)*xred_der(3))
     phkxredout(1,1)=cos(arg)
     phkxredout(2,1)=sin(arg)


!    Loop over bands

     do iband = 1,nband_k

       cwave0(:,:)=cg(:,1+(iband - 1)*npw_k*dtset%nspinor+icg0:&
&       iband*npw_k*dtset%nspinor+icg0)
       cwavef3(:,:)=cg3(:,1+(iband-1)*npw_k*dtset%nspinor+icg0:&
&       iband*npw_k*dtset%nspinor+icg0)

!      Compute vtrial1 | cwafef3 >
       gh1(:,:) = 0._dp
       tim_fourwf = 0 ; weight = 1._dp
       if (me==0) call status(counter,dtfil%filstat,iexit,level,'call fourwf  ')
       call fourwf(cplex,vlocal1,cwavef3,gh1,wfraug,gbound,gbound,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,dtset%ngfft,npw_k,npw_k,n4,n5,n6,option,&
&       dtset%paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=dtset%use_gpu_cuda)

!      In case i2pert = phonon-type perturbation
!      add first-order change in the nonlocal potential

       if (i2pert < natom + 1) then

         gvnl(:,:) = 0._dp
         itypat = dtset%typat(i2pert)
         signs=2 ; choice=2 ; nnlout=3 ; natom_der=1 ; nattyp_der(1)=1
         ntypat_der=1 ; matblk_der=1 ; tim_nonlop = 0 ; paw_opt=0 ; cpopt=-1
         ekb_typ(:,1,1)=psps%ekb(:,itypat)
         if (nspinor==2) then
           ekb_typ(:,1,2)=psps%ekb(:,itypat)
           ekb_typ(:,1,3:4)=zero
         end if
         indlmn_typ(:,:,1)=psps%indlmn(:,:,itypat)

         do ilmn=1,psps%lmnmax
           ffnlkq(:,:,ilmn,1)=ffnl(:,:,ilmn,itypat)
           ffnlk(:,:,ilmn,1)=ffnl(:,:,ilmn,itypat)
         end do

         if (me == 0) call status(counter,dtfil%filstat,iexit,level,'call nonlop  ')
         call nonlop(atindx1,choice,cpopt,cprj_dum,psps%dimekb,ntypat_der,dimffnl,dimffnl,ekb_typ,&
&         enlout,ffnlk,ffnlkq,gmet,gprimd,&
&         i2dir,indlmn_typ,istwf_k,kg_k,kg_k,kpg_k,qpg_k,&
&         kpt,kpq,dum,psps%lmnmax,&
&         matblk_der,mgfft,mpi_enreg,psps%mpsang,psps%mpssoang,&
&         natom_der,nattyp_der,dtset%ngfft,nkpg,nqpg,&
&         nloalg_der,nnlout,npw_k,npw1_k,dtset%nspinor,dtset%nspinor,&
&         ntypat_der,0,paw_opt,phkxredin,phkxredout,&
&         ph1d_der,ph3din,ph3dout,&
&         signs,dum_sij,dum_svectout,&
&         tim_nonlop,ucvol,psps%useylm,cwavef3,gvnl,&
&         use_gpu_cuda=dtset%use_gpu_cuda)

         gh1(:,:) = gh1(:,:) + gvnl(:,:)

       end if

       ii = (iband - 1)*npw_k*dtset%nspinor + icg0
       dotr = 0._dp ; doti = 0._dp
       do ipw = 1,npw_k
         ii = ii + 1
         dotr = dotr + cg1(1,ii)*gh1(1,ipw) + cg1(2,ii)*gh1(2,ipw)
         doti = doti + cg1(1,ii)*gh1(2,ipw) - cg1(2,ii)*gh1(1,ipw)
       end do

!      Compute vtrial1 | cwave0 >
       gh0(:,:) = 0._dp
       tim_fourwf = 0 ; weight = 1._dp
       if (me==0) call status(counter,dtfil%filstat,iexit,level,'call fourwf  ')
       call fourwf(cplex,vlocal1,cwave0,gh0,wfraug,gbound,gbound,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,dtset%ngfft,npw_k,npw_k,n4,n5,n6,option,&
&       dtset%paral_kgb,tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

!      In case i2pert = phonon-type perturbation
!      add first-order change in the nonlocal potential

       if (i2pert < natom + 1) then

         gvnl(:,:) = 0._dp
         itypat = dtset%typat(i2pert)
         signs=2 ; choice=2 ; nnlout=3 ; natom_der=1 ; nattyp_der(1)=1
         ntypat_der=1 ; matblk_der=1 ; tim_nonlop = 0 ; paw_opt=0 ; cpopt=-1
         ekb_typ(:,1,1)=psps%ekb(:,itypat)
         if (nspinor==2) then
           ekb_typ(:,1,2)=psps%ekb(:,itypat)
           ekb_typ(:,1,3:4)=zero
         end if
         indlmn_typ(:,:,1)=psps%indlmn(:,:,itypat)


         do ilmn=1,psps%lmnmax
           ffnlkq(:,:,ilmn,1)=ffnl(:,:,ilmn,itypat)
           ffnlk(:,:,ilmn,1)=ffnl(:,:,ilmn,itypat)
         end do

         if (me==0) call status(counter,dtfil%filstat,iexit,level,'call nonlop  ')
         call nonlop(atindx1,choice,cpopt,cprj_dum,psps%dimekb,ntypat_der,dimffnl,dimffnl,ekb_typ,&
&         enlout,ffnlk,ffnlkq,gmet,gprimd,&
&         i2dir,indlmn_typ,istwf_k,kg_k,kg_k,kpg_k,qpg_k,&
&         kpt,kpq,dum,psps%lmnmax,&
&         matblk_der,mgfft,mpi_enreg,psps%mpsang,psps%mpssoang,&
&         natom_der,nattyp_der,dtset%ngfft,nkpg,nqpg,&
&         nloalg_der,nnlout,npw_k,npw1_k,dtset%nspinor,dtset%nspinor,&
&         ntypat_der,0,paw_opt,phkxredin,phkxredout,&
&         ph1d_der,ph3din,ph3dout,&
&         signs,dum_sij,dum_svectout,&
&         tim_nonlop,ucvol,psps%useylm,cwave0,gvnl,&
&         use_gpu_cuda=dtset%use_gpu_cuda)


         gh0(:,:) = gh0(:,:) + gvnl(:,:)

       end if



!      Compute the dft contribution to the Lagrange multiplier
!      cwavef3 and cwave0 have been transferred to gh1 and gh0
!      these vectors will be used to store the wavefunctions of band iband
!      cg1 and gh0 contain the wavefunctions of band jband

       lagr = 0._dp ; lagi = 0._dp
       do jband = 1, nband_k

         ii = (jband - 1)*npw_k*dtset%nspinor + icg0
         jj = (iband - 1)*npw_k*dtset%nspinor + icg0

!        dot1r and dot1i contain < u_mk | v^(1) | u_nk >
!        dot2r and dot2i contain < u_nk^(1) | u_mk^(1) >
!        m -> jband and n -> iband

         dot1r = 0._dp ; dot1i = 0._dp
         dot2r = 0._dp ; dot2i = 0._dp
         do ipw = 1, npw_k

           ii = ii + 1 ; jj = jj + 1

           dot1r = dot1r + cg(1,ii)*gh0(1,ipw) + cg(2,ii)*gh0(2,ipw)
           dot1i = dot1i + cg(1,ii)*gh0(2,ipw) - cg(2,ii)*gh0(1,ipw)

           dot2r = dot2r + cg1(1,jj)*cg3(1,ii) + &
&           cg1(2,jj)*cg3(2,ii)
           dot2i = dot2i + cg1(1,jj)*cg3(2,ii) - &
&           cg1(2,jj)*cg3(1,ii)

         end do  !  ipw

         lagr = lagr + dot1r*dot2r - dot1i*dot2i
         lagi = lagi + dot1r*dot2i + dot1i*dot2r

       end do    ! jband

       sumr = sumr + &
&       dtset%wtk(ikpt)*occ(bantot+iband)*(dotr-lagr)

       sumi = sumi + &
&       dtset%wtk(ikpt)*occ(bantot+iband)*(doti-lagi)

!      DEBUG
!      tmp(1) = tmp(1) + dtset%wtk(ikpt)*occ(bantot+iband)*lagr
!      tmp(2) = tmp(2) + dtset%wtk(ikpt)*occ(bantot+iband)*lagi
!      ENDDEBUG


     end do   ! end loop over bands

     bantot = bantot + nband_k
     icg0 = icg0 + npw_k*dtset%nspinor*nband_k
     ikg = ikg + npw_k

     ABI_DEALLOCATE(cwave0)
     ABI_DEALLOCATE(cwavef3)
     ABI_DEALLOCATE(gh0)
     ABI_DEALLOCATE(gh1)
     ABI_DEALLOCATE(gvnl)
     ABI_DEALLOCATE(gbound)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(ffnlk)
     ABI_DEALLOCATE(ffnlkq)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(qpg_k)
     ABI_DEALLOCATE(ph3din)
     ABI_DEALLOCATE(ph3dout)
     ABI_DEALLOCATE(ekb_typ)
     ABI_DEALLOCATE(indlmn_typ)

   end do   ! end loop over k-points

 end do   ! end loop over spins


 if (mpi_enreg%paral_compil_kpt == 1) then
   buffer(1) = sumr ; buffer(2) = sumi
   call xsum_mpi(buffer,spaceComm,ierr)
   sumr = buffer(1) ; sumi = buffer(2)
 end if


 d3lo(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumr
!d3lo(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumi

!In some cases, the imaginary part is /= 0 because of the
!use of time reversal symmetry

 d3lo(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 0._dp

 ABI_DEALLOCATE(vlocal1)
 ABI_DEALLOCATE(wfraug)
 ABI_DEALLOCATE(ph1d_der)

 call status(0,dtfil%filstat,iexit,level,'exit          ')

end subroutine resp3dte
!!***
