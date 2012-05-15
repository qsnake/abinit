!{\src2tex{textfont=tt}}
!!****f* ABINIT/getghc
!!
!! NAME
!! getghc
!!
!! FUNCTION
!! Compute <G|H|C> for input vector |C> expressed in reciprocal space
!! Result is put in array ghc.
!! <G|Vnonlocal|C> is also returned in gvnlc.
!! if required, <G|S|C> is returned in gsc (S=overlap - PAW only)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, LSI, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cpopt=flag defining the status of cwaveprj%cp(:)=<Proj_i|Cnk> scalars (PAW only)
!!       (same meaning as in nonlop.F90 routine)
!!       if cpopt=-1, <p_lmn|in> (and derivatives) are computed here (and not saved)
!!       if cpopt= 0, <p_lmn|in> are computed here and saved
!!       if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!       if cpopt= 2  <p_lmn|in> are already in memory;
!!       if cpopt= 3  <p_lmn|in> are already in memory; first derivatives are computed here and saved
!!       if cpopt= 4  <p_lmn|in> and first derivatives are already in memory;
!! cwavef(2,npw*nspinor*ndat)=planewave coefficients of wavefunction.
!! dimffnl=second dimension of ffnl (1+number of derivatives)
!! ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!! filstat=name of the status file
!! gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!! kg_k(3,npw)=G vec coordinates wrt recip lattice transl.
!! kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!! lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!        Typically lambda is the eigenvalue (or its guess)
!! lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!       =if useylm=0, max number of (l,n)   comp. over all type of psps
!! matblk=dimension of the array ph3d
!! mgfft=maximum size for 1D FFTs
!! mpi_enreg=informations about MPI parallelization
!! mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!! mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!! natom=number of atoms in unit cell.
!! ndat=number of FFT to do in //
!! npw=number of planewaves in basis for given k point.
!! nspinor=number of spinorial components of the wavefunctions (on current proc)
!! ntypat=number of types of atoms in cell.
!! nvloc=final dimension of vlocal (usually 1, but 4 for non-collinear)
!! n4,n5,n6 used for dimensionning of vlocal
!! ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!! prtvol=control print volume and debugging output
!! sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!    (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                      if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!! tim_getghc=timing code of the calling subroutine(can be set to 0 if not attributed)
!! type_calc= option governing which part of Hamitonian is to be applied:
!             0: whole Hamiltonian
!!            1: local part only
!!            2: non-local+kinetic only (added to the exixting Hamiltonian)
!!            3: local + kinetic only (added to the existing Hamiltonian)
!! vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!! vxctaulocal(n4,n5,n6,nvloc,4)= local potential corresponding to the derivative of XC energy with respect to
!!  kinetic energy density, in real space, on the augmented fft grid. (optional argument)
!!  This array contains also the gradient of vxctaulocal (gvxctaulocal) in vxctaulocal(:,:,:,:,2:4).
!!
!! OUTPUT
!!  ghc(2,npw*nspinor*ndat)=matrix elements <G|H|C> (if sij_opt>=0)
!!                                       or <G|H-lambda.S|C> (if sij_opt=-1)
!!  gvnlc(2,npw*nspinor*ndat)=matrix elements <G|Vnonlocal|C> (if sij_opt>=0)
!!                                         or <G|Vnonlocal-lambda.S|C> (if sij_opt=-1)
!!    (sometimes desired for computing nonlocal part of total energy, but can be ignored).
!!  if (sij_opt=1)
!!    gsc(2,npw*nspinor*ndat)=matrix elements <G|S|C> (S=overlap).
!!
!! SIDE EFFECTS
!!  ====== if gs_ham%usepaw==1 and
!!  cwaveprj(natom,nspinor*(1+cpopt))= wave functions at k projected with nl projectors
!!
!! PARENTS
!!      cgwf,cgwf3,ks_ddiago,lobpcgIIwf,lobpcgccIIIwf,lobpcgccIIwf,lobpcgwf
!!      m_lobpcgIIIwf,mkresi,prep_getghc,update_mmat
!!
!! CHILDREN
!!      fourwf,leave_new,nonlop,status,timab,wrtout,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getghc(cpopt,cwavef,cwaveprj,dimffnl,ffnl,filstat,ghc,gsc,gs_ham,&
&  gvnlc,kg_k,kinpw,lambda,lmnmax,&
&  matblk,mgfft,mpi_enreg,mpsang,mpssoang,&
&  natom,ndat,npw,nspinor,ntypat,nvloc,n4,n5,n6,&
&  paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,type_calc,vlocal,&
&  vxctaulocal) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
  use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getghc'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_53_ffts
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpopt,dimffnl,lmnmax,matblk,mgfft,mpsang,mpssoang,n4,n5
 integer,intent(in) :: n6,natom,ndat,npw,nspinor,ntypat,nvloc,paral_kgb, prtvol
 integer,intent(in) :: sij_opt,tim_getghc,type_calc
 real(dp) :: lambda
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: gs_ham
!arrays
 integer,intent(in) :: kg_k(3,npw)
 real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat),kinpw(npw)
 real(dp),intent(inout) :: cwavef(2,npw*nspinor*ndat)
 real(dp),intent(inout) :: ghc   (2,npw*nspinor*ndat)
 real(dp),intent(inout) :: gvnlc (2,npw*nspinor*ndat)
 real(dp),intent(inout) :: ph3d(2,npw,matblk),vlocal(n4,n5,n6,nvloc)
 real(dp),intent(out) :: gsc(2,npw*nspinor*ndat*((sij_opt+1)/2))
 real(dp),intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)
 type(cprj_type),intent(inout) :: cwaveprj(natom,nspinor*((cpopt+5)/5)*gs_ham%usepaw)
!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,level=114,re=1
 integer :: choice,cplex,cpopt_here,i1,i2,i3,idat,idir,ierr,iexit
 integer :: ig,igspinor,ii,iispinor,ipw,ispinor,nkpg,nnlout,nspinortot
 integer :: paw_opt,shift,signs,tim_fourwf,tim_nonlop
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc
 real(dp) :: ghcim,ghcre,gp2pi1,gp2pi2,gp2pi3,kpt_cart,kg_k_cart,weight
 character(len=500) :: message
!arrays
 real(dp) :: enlout(1),nonlop_dum(1,1),nonlop_dum2(1,1),tsec(2)
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:),cwavef_nonlop(:,:)
 real(dp),allocatable :: gcwavef(:,:,:),gcwavef1(:,:,:),gcwavef2(:,:,:)
 real(dp),allocatable :: ghc1(:,:),ghc2(:,:),ghc3(:,:),ghc4(:,:),gsc_nonlop(:,:)
 real(dp),allocatable :: gvnlc_nonlop(:,:),kpg_dum(:,:)
 real(dp),allocatable :: lcwavef(:,:),lcwavef1(:,:),lcwavef2(:,:)
 real(dp),allocatable :: vlocal_tmp(:,:,:),work(:,:,:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' enter getghc '
!write(std_out,*)' getghc : cwavef(:,1)=',cwavef(:,1)
!stop
!ENDDEBUG

!Keep track of total time spent in getghc:
 call timab(200+tim_getghc,1,tsec)

 if(prtvol<0)then
   call status(0,filstat,iexit,level,'enter         ')
 end if

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
   write(message,'(80a,a,a)') ('=',ii=1,80),ch10,&
&   ' getghc : enter, debugging '
   call wrtout(std_out,message,'PERS')
 end if

!Parallelization over spinors management
 nspinortot=min(2,(1+mpi_enreg%paral_spin)*nspinor)
 if (mpi_enreg%paral_spin==0) then
   shift=npw
   nspinor1TreatedByThisProc=.true.
   nspinor2TreatedByThisProc=(nspinortot==2)
 else
   shift=0
   nspinor1TreatedByThisProc=(mpi_enreg%me_spin==0)
   nspinor2TreatedByThisProc=(mpi_enreg%me_spin==1)
 end if

 if ((type_calc==0).or.(type_calc==1).or.(type_calc==3)) then

   ABI_ALLOCATE(work,(2,n4,n5,n6*ndat))

!  Apply the local potential to the wavefunction
!  Start from wavefunction in reciprocal space cwavef
!  End with function ghc in reciprocal space also.
   weight=1.0_dp

!  DEBUG
!  write(std_out,*)' getghc : will call fourwf '
!  write(std_out,*)' vlocal='
!  do i1=1,n4
!  do i2=1,n5
!  do i3=1,n6
!  write(std_out,*)i1,i2,i3,vlocal(i1,i2,i3,1)
!  end do
!  end do
!  end do
!  do i1 =1, nspinor
!  do ii=1,npw
!  write(std_out,*)ii,i1,cwavef(1,ii),cwavef(2,ii)
!  end do
!  end do
!  ENDDEBUG

!  Application of the local potential
   tim_fourwf=1

   if (nspinortot==2) then
     ABI_ALLOCATE(cwavef1,(2,npw*ndat))
     ABI_ALLOCATE(cwavef2,(2,npw*ndat))
     do idat=1,ndat
       do ipw=1,npw
         cwavef1(1:2,ipw+(idat-1)*npw)=cwavef(1:2,ipw+(idat-1)*nspinor*npw)
         cwavef2(1:2,ipw+(idat-1)*npw)=cwavef(1:2,ipw+(idat-1)*nspinor*npw+shift)
       end do
     end do
   end if

!  Treat scalar local potentials
   if (nvloc==1) then

     if (nspinortot==1) then

       call fourwf(1,vlocal,cwavef,ghc,work,gs_ham%gbound,gs_ham%gbound,&
&       gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)

       if(present(vxctaulocal))then  !metaGGA
!        Do it in 3 STEPs:
!        STEP1: Compute grad of cwavef and Laplacian of cwavef
         ABI_ALLOCATE(gcwavef,(2,npw*ndat,3))
         ABI_ALLOCATE(lcwavef,(2,npw*ndat))
!        $OMP PARALLEL DO PRIVATE(idat,ipw) &
!        $OMP&SHARED(npw,ndat,gcwavef,lcwavef)
         do idat=1,ndat
           do ipw=1,npw
             gcwavef(:,ipw+(idat-1)*npw,1:3)=zero
             lcwavef(:,ipw+(idat-1)*npw)  =zero
           end do
         end do
!        $OMP END PARALLEL DO
         do idir=1,3
           gp2pi1=gs_ham%gprimd(idir,1)*two_pi
           gp2pi2=gs_ham%gprimd(idir,2)*two_pi
           gp2pi3=gs_ham%gprimd(idir,3)*two_pi
           kpt_cart=gp2pi1*gs_ham%kpoint(1)+gp2pi2*gs_ham%kpoint(2)+gp2pi3*gs_ham%kpoint(3)
!          Multiplication by 2pi i (G+k)_idir for gradient
!          Multiplication by -(2pi (G+k)_idir )**2 for Laplacian
           do idat=1,ndat
             do ipw=1,npw
               kg_k_cart=gp2pi1*kg_k(1,ipw)+gp2pi2*kg_k(2,ipw)+gp2pi3*kg_k(3,ipw)+kpt_cart
               gcwavef(1,ipw+(idat-1)*npw,idir)= cwavef(2,ipw+(idat-1)*npw)*kg_k_cart
               gcwavef(2,ipw+(idat-1)*npw,idir)=-cwavef(1,ipw+(idat-1)*npw)*kg_k_cart
               lcwavef(1,ipw+(idat-1)*npw)=lcwavef(1,ipw+(idat-1)*npw)-cwavef(1,ipw+(idat-1)*npw)*kg_k_cart**2
               lcwavef(2,ipw+(idat-1)*npw)=lcwavef(2,ipw+(idat-1)*npw)-cwavef(2,ipw+(idat-1)*npw)*kg_k_cart**2
             end do
           end do
         end do ! idir
!        STEP2: Compute (vxctaulocal)*(Laplacian of cwavef) and add it to ghc
         ABI_ALLOCATE(ghc1,(2,npw*ndat))
         call fourwf(1,vxctaulocal(:,:,:,:,1),lcwavef,ghc1,work,gs_ham%gbound,gs_ham%gbound,&
&         gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)
!        $OMP PARALLEL DO PRIVATE(idat,ipw) &
!        $OMP&SHARED(npw,ndat,ghc,ghc1)
         do idat=1,ndat
           do ipw=1,npw
             ghc(:,ipw+(idat-1)*npw)=ghc(:,ipw+(idat-1)*npw)-half*ghc1(:,ipw+(idat-1)*npw)
           end do
         end do
!        $OMP END PARALLEL DO
         ABI_DEALLOCATE(ghc1)
         ABI_DEALLOCATE(lcwavef)
!        STEP3: Compute sum of (grad components of vxctaulocal)*(grad components of cwavef)
         ABI_ALLOCATE(ghc1,(2,npw*ndat))
         do idir=1,3
           call fourwf(1,vxctaulocal(:,:,:,:,1+idir),gcwavef(:,:,idir),ghc1,work,gs_ham%gbound,gs_ham%gbound,&
&           gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&           npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&           use_gpu_cuda=gs_ham%use_gpu_cuda)
!          $OMP PARALLEL DO PRIVATE(idat,ipw) &
!          $OMP&SHARED(ndat,npw,ghc,ghc1)
           do idat=1,ndat
             do ipw=1,npw
               ghc(:,ipw+(idat-1)*npw)=ghc(:,ipw+(idat-1)*npw)-half*ghc1(:,ipw+(idat-1)*npw)
             end do
           end do
!          $OMP END PARALLEL DO
         end do ! idir
         ABI_DEALLOCATE(ghc1)
         ABI_DEALLOCATE(gcwavef)

       end if ! if present(vxctaulocal) i.e. metaGGA

     else ! nspinortot==2

       if (nspinor1TreatedByThisProc) then

         ABI_ALLOCATE(ghc1,(2,npw*ndat))
         ghc1(:,:)=zero
         call fourwf(1,vlocal,cwavef1,ghc1,work,gs_ham%gbound,gs_ham%gbound,&
&         gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)
         do idat=1,ndat
           do ipw =1, npw
             ghc(1:2,ipw+(idat-1)*nspinor*npw)    =ghc1(1:2,ipw+(idat-1)*npw)
           end do
         end do
         ABI_DEALLOCATE(ghc1)

         if(present(vxctaulocal))then  !metaGGA
!          Do it in 3 STEPs:
!          STEP1: Compute grad of cwavef and Laplacian of cwavef
           ABI_ALLOCATE(gcwavef1,(2,npw*ndat,3))
           ABI_ALLOCATE(lcwavef1,(2,npw*ndat))
!          $OMP PARALLEL DO PRIVATE(idat,ipw) &
!          $OMP&SHARED(npw,ndat,gcwavef1,gcwavef2,lcwavef1,lcwavef2)
           do idat=1,ndat
             do ipw=1,npw
               gcwavef1(:,ipw+(idat-1)*npw,1:3)=zero
               lcwavef1(:,ipw+(idat-1)*npw)  =zero
             end do
           end do
!          $OMP END PARALLEL DO
           do idir=1,3
             gp2pi1=gs_ham%gprimd(idir,1)*two_pi
             gp2pi2=gs_ham%gprimd(idir,2)*two_pi
             gp2pi3=gs_ham%gprimd(idir,3)*two_pi
             kpt_cart=gp2pi1*gs_ham%kpoint(1)+gp2pi2*gs_ham%kpoint(2)+gp2pi3*gs_ham%kpoint(3)
!            Multiplication by 2pi i (G+k)_idir for gradient
!            Multiplication by -(2pi (G+k)_idir )**2 for Laplacian
             do idat=1,ndat
               do ipw=1,npw
                 kg_k_cart=gp2pi1*kg_k(1,ipw)+gp2pi2*kg_k(2,ipw)+gp2pi3*kg_k(3,ipw)+kpt_cart
                 gcwavef1(1,ipw+(idat-1)*npw,idir)= cwavef1(2,ipw+(idat-1)*npw)*kg_k_cart
                 gcwavef1(2,ipw+(idat-1)*npw,idir)=-cwavef1(1,ipw+(idat-1)*npw)*kg_k_cart
                 lcwavef1(1,ipw+(idat-1)*npw)=lcwavef1(1,ipw+(idat-1)*npw)-cwavef1(1,ipw+(idat-1)*npw)*kg_k_cart**2
                 lcwavef1(2,ipw+(idat-1)*npw)=lcwavef1(2,ipw+(idat-1)*npw)-cwavef1(2,ipw+(idat-1)*npw)*kg_k_cart**2
               end do
             end do
           end do ! idir
!          STEP2: Compute (vxctaulocal)*(Laplacian of cwavef) and add it to ghc
           ABI_ALLOCATE(ghc1,(2,npw*ndat))
           call fourwf(1,vxctaulocal(:,:,:,:,1),lcwavef1,ghc1,work,gs_ham%gbound,gs_ham%gbound,&
&           gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&           npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&           use_gpu_cuda=gs_ham%use_gpu_cuda)
!          $OMP PARALLEL DO PRIVATE(idat,ipw) &
!          $OMP&SHARED(npw,ndat,ghc,ghc1,ghc2)
           do idat=1,ndat
             do ipw=1,npw
               ghc(:,ipw+(idat-1)*nspinor*npw)    =ghc(:,ipw+(idat-1)*nspinor*npw)&
&               -half*ghc1(:,ipw+(idat-1)*npw)
             end do
           end do
!          $OMP END PARALLEL DO
           ABI_DEALLOCATE(ghc1)
           ABI_DEALLOCATE(lcwavef1)
!          STEP3: Compute (grad components of vxctaulocal)*(grad components of cwavef)
           ABI_ALLOCATE(ghc1,(2,npw*ndat))
           do idir=1,3
             call fourwf(1,vxctaulocal(:,:,:,:,1+idir),gcwavef1(:,:,idir),ghc1,work,gs_ham%gbound,&
&             gs_ham%gbound,gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&             npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&             use_gpu_cuda=gs_ham%use_gpu_cuda)
!            $OMP PARALLEL DO PRIVATE(idat,ipw) &
!            $OMP&SHARED(ndat,npw,ghc,ghc1,ghc2)
             do idat=1,ndat
               do ipw=1,npw
                 ghc(:,ipw+(idat-1)*nspinor*npw)    =ghc(:,ipw+(idat-1)*nspinor*npw)&
&                 -half*ghc1(:,ipw+(idat-1)*npw)
               end do
             end do
!            $OMP END PARALLEL DO
           end do ! idir
           ABI_DEALLOCATE(ghc1)
           ABI_DEALLOCATE(gcwavef1)
         end if ! if present(vxctaulocal) i.e. metaGGA

       end if ! spin 1 treated by this proc

       if (nspinor2TreatedByThisProc) then

         ABI_ALLOCATE(ghc2,(2,npw*ndat))
         ghc2(:,:)=zero
         
         call fourwf(1,vlocal,cwavef2,ghc2,work,gs_ham%gbound,gs_ham%gbound,&
&         gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)
         do idat=1,ndat
           do ipw=1,npw
             ghc(1:2,ipw+(idat-1)*nspinor*npw+shift)=ghc2(1:2,ipw+(idat-1)*npw)
           end do
         end do
         ABI_DEALLOCATE(ghc2)

         if(present(vxctaulocal))then  !metaGGA
!          Do it in 3 STEPs:
!          STEP1: Compute grad of cwavef and Laplacian of cwavef
           ABI_ALLOCATE(gcwavef1,(2,npw*ndat,3))
           ABI_ALLOCATE(gcwavef2,(2,npw*ndat,3))
           ABI_ALLOCATE(lcwavef1,(2,npw*ndat))
           ABI_ALLOCATE(lcwavef2,(2,npw*ndat))
!          $OMP PARALLEL DO PRIVATE(idat,ipw) &
!          $OMP&SHARED(npw,ndat,gcwavef1,gcwavef2,lcwavef1,lcwavef2)
           do idat=1,ndat
             do ipw=1,npw
               gcwavef2(:,ipw+(idat-1)*npw,1:3)=zero
               lcwavef2(:,ipw+(idat-1)*npw)  =zero
             end do
           end do
!          $OMP END PARALLEL DO
           do idir=1,3
             gp2pi1=gs_ham%gprimd(idir,1)*two_pi
             gp2pi2=gs_ham%gprimd(idir,2)*two_pi
             gp2pi3=gs_ham%gprimd(idir,3)*two_pi
             kpt_cart=gp2pi1*gs_ham%kpoint(1)+gp2pi2*gs_ham%kpoint(2)+gp2pi3*gs_ham%kpoint(3)
!            Multiplication by 2pi i (G+k)_idir for gradient
!            Multiplication by -(2pi (G+k)_idir )**2 for Laplacian
             do idat=1,ndat
               do ipw=1,npw
                 kg_k_cart=gp2pi1*kg_k(1,ipw)+gp2pi2*kg_k(2,ipw)+gp2pi3*kg_k(3,ipw)+kpt_cart
                 gcwavef2(1,ipw+(idat-1)*npw,idir)= cwavef2(2,ipw+(idat-1)*npw)*kg_k_cart
                 gcwavef2(2,ipw+(idat-1)*npw,idir)=-cwavef2(1,ipw+(idat-1)*npw)*kg_k_cart
                 lcwavef2(1,ipw+(idat-1)*npw)=lcwavef2(1,ipw+(idat-1)*npw)-cwavef2(1,ipw+(idat-1)*npw)*kg_k_cart**2
                 lcwavef2(2,ipw+(idat-1)*npw)=lcwavef2(2,ipw+(idat-1)*npw)-cwavef2(2,ipw+(idat-1)*npw)*kg_k_cart**2
               end do
             end do
           end do ! idir
!          STEP2: Compute (vxctaulocal)*(Laplacian of cwavef) and add it to ghc
           ABI_ALLOCATE(ghc2,(2,npw*ndat))
           call fourwf(1,vxctaulocal(:,:,:,:,1),lcwavef2,ghc2,work,gs_ham%gbound,gs_ham%gbound,&
&           gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&           npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&           use_gpu_cuda=gs_ham%use_gpu_cuda)
!          $OMP PARALLEL DO PRIVATE(idat,ipw) &
!          $OMP&SHARED(npw,ndat,ghc,ghc1,ghc2)
           do idat=1,ndat
             do ipw=1,npw
               ghc(:,ipw+(idat-1)*nspinor*npw)=ghc(:,ipw+(idat-1)*nspinor*npw)&
&               -half*ghc2(:,ipw+(idat-1)*npw)
             end do
           end do
!          $OMP END PARALLEL DO
           ABI_DEALLOCATE(ghc2)
           ABI_DEALLOCATE(lcwavef2)
!          STEP3: Compute sum of (grad components of vxctaulocal)*(grad components of cwavef)
           ABI_ALLOCATE(ghc2,(2,npw*ndat))
           do idir=1,3
             call fourwf(1,vxctaulocal(:,:,:,:,1+idir),gcwavef2(:,:,idir),ghc2,work,gs_ham%gbound,gs_ham%gbound,&
&             gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&             npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&             use_gpu_cuda=gs_ham%use_gpu_cuda)
!            $OMP PARALLEL DO PRIVATE(idat,ipw) &
!            $OMP&SHARED(ndat,npw,ghc,ghc1,ghc2)
             do idat=1,ndat
               do ipw=1,npw
                 ghc(:,ipw+(idat-1)*nspinor*npw)    =ghc(:,ipw+(idat-1)*nspinor*npw)&
&                 -half*ghc2(:,ipw+(idat-1)*npw)
               end do
             end do
!            $OMP END PARALLEL DO
           end do ! idir
           ABI_DEALLOCATE(ghc2)
           ABI_DEALLOCATE(gcwavef2)
         end if ! if present(vxctaulocal) i.e. metaGGA

       end if ! spin 2 treated by this proc
     end if ! npsinortot

!    Treat non-collinear local potentials
   else if (nvloc==4) then
     ABI_ALLOCATE(ghc1,(2,npw*ndat))
     ABI_ALLOCATE(ghc2,(2,npw*ndat))
     ABI_ALLOCATE(ghc3,(2,npw*ndat))
     ABI_ALLOCATE(ghc4,(2,npw*ndat))
     ghc1(:,:)=zero; ghc2(:,:)=zero; ghc3(:,:)=zero ;  ghc4(:,:)=zero
     ABI_ALLOCATE(vlocal_tmp,(n4,n5,n6))
!    ghc1=v11*phi1
     vlocal_tmp(:,:,:)=vlocal(:,:,:,1)
     if (nspinor1TreatedByThisProc) then
       call fourwf(1,vlocal_tmp,cwavef1,ghc1,work,gs_ham%gbound,gs_ham%gbound,&
&       gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
!    ghc2=v22*phi2
     vlocal_tmp(:,:,:)=vlocal(:,:,:,2)
     if (nspinor2TreatedByThisProc) then
       call fourwf(1,vlocal_tmp,cwavef2,ghc2,work,gs_ham%gbound,gs_ham%gbound,&
&       gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
     ABI_DEALLOCATE(vlocal_tmp)
     cplex=2
     ABI_ALLOCATE(vlocal_tmp,(cplex*n4,n5,n6))
!    ghc3=(re(v12)-im(v12))*phi1
     do i3=1,n6
       do i2=1,n5
         do i1=1,n4
           vlocal_tmp(2*i1-1,i2,i3)= vlocal(i1,i2,i3,3)
           vlocal_tmp(2*i1  ,i2,i3)=-vlocal(i1,i2,i3,4)
         end do
       end do
     end do
     if (nspinor1TreatedByThisProc) then
       call fourwf(cplex,vlocal_tmp,cwavef1,ghc3,work,gs_ham%gbound,gs_ham%gbound,&
&       gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
!    ghc4=(re(v12)+im(v12))*phi2
     if (nspinor2TreatedByThisProc) then
       do i3=1,n6
         do i2=1,n5
           do i1=1,n4
             vlocal_tmp(2*i1,i2,i3)=-vlocal_tmp(2*i1,i2,i3)
           end do
         end do
       end do
       call fourwf(cplex,vlocal_tmp,cwavef2,ghc4,work,gs_ham%gbound,gs_ham%gbound,&
&       gs_ham%istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw,npw,n4,n5,n6,2,paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
       ABI_DEALLOCATE(vlocal_tmp)
     end if
!    Build ghc from pieces
!    (v11,v22,Re(v12)+iIm(v12);Re(v12)-iIm(v12))(psi1;psi2): matrix product
     if (mpi_enreg%paral_spin==0) then
       do idat=1,ndat
         do ipw=1,npw
           ghc(1:2,ipw+(idat-1)*nspinor*npw)    =ghc1(1:2,ipw+(idat-1)*npw)+ghc4(1:2,ipw+(idat-1)*npw)
           ghc(1:2,ipw+(idat-1)*nspinor*npw+npw)=ghc3(1:2,ipw+(idat-1)*npw)+ghc2(1:2,ipw+(idat-1)*npw)
         end do
       end do
     else
       call xsum_mpi(ghc4,mpi_enreg%comm_spin,ierr)
       call xsum_mpi(ghc3,mpi_enreg%comm_spin,ierr)
       if (nspinor1TreatedByThisProc) then
         do idat=1,ndat
           do ipw=1,npw
             ghc(1:2,ipw+(idat-1)*nspinor*npw)    =ghc1(1:2,ipw+(idat-1)*npw)+ghc4(1:2,ipw+(idat-1)*npw)
           end do
         end do
       else if (nspinor2TreatedByThisProc) then
         do idat=1,ndat
           do ipw=1,npw
             ghc(1:2,ipw+(idat-1)*nspinor*npw)=ghc3(1:2,ipw+(idat-1)*npw)+ghc2(1:2,ipw+(idat-1)*npw)
           end do
         end do
       end if
     end if
     ABI_DEALLOCATE(ghc1)
     ABI_DEALLOCATE(ghc2)
     ABI_DEALLOCATE(ghc3)
     ABI_DEALLOCATE(ghc4)
   end if ! nvloc

   if (nspinortot==2)  then
     ABI_DEALLOCATE(cwavef1)
     ABI_DEALLOCATE(cwavef2)
   end if
   ABI_DEALLOCATE(work)
   if(prtvol<0)then
     call status(0,filstat,iexit,level,'call nonlop   ')
   end if

 end if ! type_calc

 if ((type_calc==0).or.(type_calc==2).or.(type_calc==3)) then

   if ((type_calc==0).or.(type_calc==2)) then
     signs=2 ; choice=1 ; nnlout=1 ; idir=0 ; tim_nonlop=1 ; nkpg=0
     cpopt_here=-1;if (gs_ham%usepaw==1) cpopt_here=cpopt
     paw_opt=gs_ham%usepaw ; if (sij_opt/=0) paw_opt=sij_opt+3
     if(ndat==1)then
       if(gs_ham%usepaw==0)then
         call nonlop(gs_ham%atindx1,choice,cpopt_here,cwaveprj,gs_ham%dimekb1,gs_ham%dimekb2,dimffnl,dimffnl,&
&         gs_ham%ekb,enlout,ffnl,ffnl,gs_ham%gmet,gs_ham%gprimd,idir,&
&         gs_ham%indlmn,gs_ham%istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_ham%kpoint,gs_ham%kpoint,&
&         lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_ham%nattyp,&
&         gs_ham%ngfft,nkpg,nkpg,gs_ham%nloalg,nnlout,npw,npw,nspinor,nspinortot,ntypat,0,paw_opt,&
&         gs_ham%phkxred,gs_ham%phkxred,gs_ham%ph1d,ph3d,ph3d,&
&         signs,nonlop_dum,nonlop_dum2,tim_nonlop,gs_ham%ucvol,gs_ham%useylm,cwavef,gvnlc,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)
       else
         call nonlop(gs_ham%atindx1,choice,cpopt_here,cwaveprj,gs_ham%dimekb1,gs_ham%dimekb2,dimffnl,dimffnl,&
&         gs_ham%ekb,enlout,ffnl,ffnl,gs_ham%gmet,gs_ham%gprimd,idir,&
&         gs_ham%indlmn,gs_ham%istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_ham%kpoint,gs_ham%kpoint,&
&         lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_ham%nattyp,&
&         gs_ham%ngfft,nkpg,nkpg,gs_ham%nloalg,nnlout,npw,npw,nspinor,nspinortot,ntypat,0,paw_opt,&
&         gs_ham%phkxred,gs_ham%phkxred,gs_ham%ph1d,ph3d,ph3d,&
&         signs,gs_ham%sij,gsc,tim_nonlop,gs_ham%ucvol,gs_ham%useylm,cwavef,gvnlc,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)
       end if
     else
       ABI_ALLOCATE(cwavef_nonlop,(2,npw*nspinor))
       ABI_ALLOCATE(gvnlc_nonlop,(2,npw*nspinor))
       if (sij_opt==1)  then
         ABI_ALLOCATE(gsc_nonlop,(2,npw*nspinor))
       end if
       do idat=1,ndat
         cwavef_nonlop(re,1:npw*nspinor)=cwavef(re,1+npw*nspinor*(idat-1):npw*nspinor*idat)
         cwavef_nonlop(im,1:npw*nspinor)=cwavef(im,1+npw*nspinor*(idat-1):npw*nspinor*idat)
         if(gs_ham%usepaw==0)then
           call nonlop(gs_ham%atindx1,choice,cpopt_here,cwaveprj,gs_ham%dimekb1,gs_ham%dimekb2,dimffnl,dimffnl,&
&           gs_ham%ekb,enlout,ffnl,ffnl,gs_ham%gmet,gs_ham%gprimd,idir,&
&           gs_ham%indlmn,gs_ham%istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_ham%kpoint,gs_ham%kpoint,&
&           lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_ham%nattyp,&
&           gs_ham%ngfft,nkpg,nkpg,gs_ham%nloalg,nnlout,npw,npw,nspinor,nspinortot,ntypat,0,paw_opt,&
&           gs_ham%phkxred,gs_ham%phkxred,gs_ham%ph1d,ph3d,ph3d,&
&           signs,nonlop_dum,nonlop_dum2,tim_nonlop,gs_ham%ucvol,gs_ham%useylm,cwavef_nonlop,gvnlc_nonlop,&
&           use_gpu_cuda=gs_ham%use_gpu_cuda)
         else
           call nonlop(gs_ham%atindx1,choice,cpopt_here,cwaveprj,gs_ham%dimekb1,gs_ham%dimekb2,dimffnl,dimffnl,&
&           gs_ham%ekb,enlout,ffnl,ffnl,gs_ham%gmet,gs_ham%gprimd,idir,&
&           gs_ham%indlmn,gs_ham%istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_ham%kpoint,gs_ham%kpoint,&
&           lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_ham%nattyp,&
&           gs_ham%ngfft,nkpg,nkpg,gs_ham%nloalg,nnlout,npw,npw,nspinor,nspinortot,ntypat,0,paw_opt,&
&           gs_ham%phkxred,gs_ham%phkxred,gs_ham%ph1d,ph3d,ph3d,&
&           signs,gs_ham%sij,gsc_nonlop,tim_nonlop,gs_ham%ucvol,gs_ham%useylm,cwavef_nonlop,gvnlc_nonlop,&
&           use_gpu_cuda=gs_ham%use_gpu_cuda)
         end if
         gvnlc(re,1+npw*nspinor*(idat-1):npw*nspinor*idat)=gvnlc_nonlop(re,1:npw*nspinor)
         gvnlc(im,1+npw*nspinor*(idat-1):npw*nspinor*idat)=gvnlc_nonlop(im,1:npw*nspinor)
         if (sij_opt==1) then
           gsc(re,1+npw*nspinor*(idat-1):npw*nspinor*idat)=gsc_nonlop(re,1:npw*nspinor)
           gsc(im,1+npw*nspinor*(idat-1):npw*nspinor*idat)=gsc_nonlop(im,1:npw*nspinor)
         end if
       end do
       ABI_DEALLOCATE(cwavef_nonlop)
       ABI_DEALLOCATE(gvnlc_nonlop)
       if (sij_opt==1)  then
         ABI_DEALLOCATE(gsc_nonlop)
       end if
     end if
   end if ! end type_calc 0 or 2 for nonlop application

   if(prtvol<0)then
     call status(0,filstat,iexit,level,'assemble      ')
   end if

   if(prtvol==-level)then
     write(message,'(a)')' getghc : components of ghc '
     call wrtout(std_out,message,'PERS')
     write(message,'(a,a)')'icp ig ispinor igspinor re/im     ghc     ',&
&     '   kinpw         cwavef      glocc        gvnlc  gsc'
     call wrtout(std_out,message,'PERS')
   end if

!  Assemble modified kinetic, local and nonlocal contributions
!  to <G|H|C(n,k)>. Take also into account build-in debugging.
   if(prtvol/=-level)then
     do idat=1,ndat
       do ispinor=1,nspinor
!        $OMP PARALLEL DO PRIVATE(idat,ig,igspinor) &
!        $OMP&SHARED(cwavef,ghc,gvnlc,kinpw,ndat,npw)
         do ig=1,npw
           igspinor=ig+npw*(ispinor-1)+npw*nspinor*(idat-1)
           if(kinpw(ig)<huge(zero)*1.d-11)then
             ghc(re,igspinor)=kinpw(ig)*cwavef(re,igspinor)+ghc(re,igspinor)+gvnlc(re,igspinor)
             ghc(im,igspinor)=kinpw(ig)*cwavef(im,igspinor)+ghc(im,igspinor)+gvnlc(im,igspinor)
           else
             ghc(re,igspinor)=zero
             ghc(im,igspinor)=zero
             if (sij_opt==1) then
               gsc(re,igspinor)=zero
               gsc(im,igspinor)=zero
             end if
           end if
         end do ! ig
       end do ! ispinor
!      $OMP END PARALLEL DO
     end do ! idat
   else
!    Here, debugging section
     do idat=1,ndat
       do ispinor=1,nspinor
!        $OMP PARALLEL DO PRIVATE(ghcre,ghcim,idat,ig,igspinor) &
!        $OMP&SHARED(cwavef,ghc,gvnlc,kinpw,ndat,npw)
         do ig=1,npw
           igspinor=ig+npw*(ispinor-1)+npw*nspinor*(idat-1)
           if(kinpw(ig)<huge(zero)*1.d-11)then
             ghcre=kinpw(ig)*cwavef(re,igspinor)+ghc(re,igspinor)+gvnlc(re,igspinor)
             ghcim=kinpw(ig)*cwavef(im,igspinor)+ghc(im,igspinor)+gvnlc(im,igspinor)
           else
             ghcre=zero
             ghcim=zero
             if (sij_opt==1) then
               gsc(re,igspinor)=zero
               gsc(im,igspinor)=zero
             end if
           end if
           iispinor=ispinor;if (mpi_enreg%paral_spin==1) iispinor=mpi_enreg%me_spin+1
           if (sij_opt == 1) then
             write(message,'(a,3(1x,i5),6(1x,es13.6))') '  1 ', ig, iispinor, igspinor,ghcre,&
&             kinpw(ig),cwavef(re,igspinor),ghc(re,igspinor),gvnlc(re,igspinor), gsc(re,igspinor)
             call wrtout(std_out,message,'PERS')
             write(message,'(a,3(1x,i5),6(1x,es13.6))') '  2 ', ig, iispinor, igspinor,ghcim,&
&             kinpw(ig),cwavef(im,igspinor),ghc(im,igspinor),gvnlc(im,igspinor), gsc(im,igspinor)
             call wrtout(std_out,message,'PERS')
           else
             write(message,'(a,3(1x,i5),6(1x,es13.6))') '  1 ', ig, iispinor, igspinor,ghcre,&
&             kinpw(ig),cwavef(re,igspinor),ghc(re,igspinor),gvnlc(re,igspinor)
             call wrtout(std_out,message,'PERS')
             write(message,'(a,3(1x,i5),6(1x,es13.6))') '  2 ', ig, iispinor, igspinor,ghcim,&
&             kinpw(ig),cwavef(im,igspinor),ghc(im,igspinor),gvnlc(im,igspinor)
             call wrtout(std_out,message,'PERS')
           end if
           ghc(re,igspinor)=ghcre
           ghc(im,igspinor)=ghcim
         end do ! ig
       end do ! ispinor
!      $OMP END PARALLEL DO
     end do ! idat
   end if

!  Structured debugging : if prtvol=-level, stop here.
   if(prtvol==-level)then
     write(message,'(a,a,a,a,i2,a)') ch10,&
&     ' getghc : exit ',&
&     ch10,'  prtvol=-',level,', debugging mode => stop '
     call wrtout(std_out,message,'PERS')
     call status(0,filstat,iexit,level,'debug => stop ')
     call leave_new('PERS')
   end if

   if(prtvol<0)then
     call status(0,filstat,iexit,level,'exit          ')
   end if

 end if ! type_calc

 call timab(200+tim_getghc,2,tsec)

!DEBUG
!write(std_out,*)' getghc : exit'
!ENDDEBUG

end subroutine getghc
!!***
