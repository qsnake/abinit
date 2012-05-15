!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_getghc
!! NAME
!! prep_getghc
!!
!! FUNCTION
!! this routine prepares the data to the call of getghc.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FBottin,MT,GZ,MD,FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blocksize= size of block for FFT
!!  cwavef(2,npw*nspinor*ndat)=planewave coefficients of wavefunction.
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  gvnlc=matrix elements <G|Vnonlocal|C>
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (hartree)
!!  lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!         Typically lambda is the eigenvalue (or its guess)
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mgfft=maximum size of 1d ffts
!!  mpi_enreg=informations about mpi parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  ntypat=number of types of atoms in unit cell.
!!  nvloc=final dimension of vlocal (usually 1, but 4 for non-collinear)
!!  n4,n5,n6 used for dimensionning of vlocal
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!     (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                       if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!  vxctaulocal(n4,n5,n6,nvloc,4)= local potential corresponding to the derivative of XC energy with respect to
!!   kinetic energy density, in real space, on the augmented fft grid. (optional argument)
!!   This array contains also the gradient of vxctaulocal (gvxctaulocal) in vxctaulocal(:,:,:,:,2:4).
!!
!! OUTPUT
!!  gwavef=(2,npw*nspinor*ndat)=matrix elements <G|H|C> (if sij_opt>=0)
!!                                  or <G|H-lambda.S|C> (if sij_opt=-1).
!!  swavef=(2,npw*nspinor*ndat)=matrix elements <G|S|C>.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!      getghc,prep_index_wavef_bandpp,prep_sort_wavef_spin,prep_wavef_sym_do
!!      prep_wavef_sym_undo,timab,xalltoallv_mpi,xcomm_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_getghc(cwavef,dimffnl,dtfil,gs_hamk,gvnlc,gwavef,swavef,ikpt,istwf_k,&
& lambda,lmnmax,matblk,blocksize,mgfft,mpi_enreg,mpsang,mpssoang,natom,npw_k,&
& nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,prtvol,sij_opt,vlocal, &
& vxctaulocal) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_getghc'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_66_wfs, except_this_one => prep_getghc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: blocksize,dimffnl,ikpt,istwf_k,lmnmax,matblk,mgfft
 integer :: mpsang,mpssoang,n4,n5,n6,natom,npw_k,nspinor
 integer :: ntypat,nvloc,paral_kgb,prtvol,sij_opt
 real(dp) :: lambda
 type(datafiles_type) :: dtfil
 type(gs_hamiltonian_type) :: gs_hamk
 type(mpi_type) :: mpi_enreg
!arrays
 real(dp) :: cwavef(2,npw_k*nspinor*blocksize)
 real(dp) :: gvnlc (2,npw_k*nspinor*blocksize)
 real(dp) :: gwavef(2,npw_k*nspinor*blocksize)
 real(dp) :: swavef(2,npw_k*nspinor*blocksize)
 real(dp) :: vlocal(n4,n5,n6,nvloc)
 real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)

!Local variables-------------------------------
!local variables for mpialltoallv
!local variable for bandpp and inversion by symetry of time
!scalars
 integer :: bandpp,bandpp_sym,binf,bsup,cpopt,ier,iibandpp,ikpt_this_proc
 integer :: iscalc,nbval,nproc_band,nproc_fft
 integer :: old_me_g0,old_paral_level,spaceComm=0,tim_getghc
 logical :: flag_inv_sym
!arrays
 integer,pointer :: idatarecv0,ndatarecv,ndatarecv_tot,ndatasend_sym
 integer,pointer :: kg_k_gather(:,:),kg_k_gather_sym(:,:)
 integer,pointer :: index_wavef_band(:),index_wavef_send(:),index_wavef_spband(:)
 integer,pointer :: rdispls(:),rdispls_sym(:)
 integer,pointer :: recvcounts(:),recvcounts_sym(:),recvcounts_sym_tot(:)
 integer,pointer :: sdispls(:),sdispls_sym(:)
 integer,pointer :: sendcounts(:),sendcounts_sym(:),sendcounts_sym_all(:)
 integer,pointer :: tab_proc(:)
 integer,allocatable :: rdisplsloc(:)
 integer,allocatable :: recvcountsloc(:),sdisplsloc(:)
 integer,allocatable :: sendcountsloc(:)
 real(dp) :: tsec(2)
 real(dp),pointer :: ffnl_gather(:,:,:,:),kinpw_gather(:),ph3d_gather(:,:,:)
 real(dp),allocatable :: cwavef_alltoall(:,:),gvnlc_alltoall(:,:)
 real(dp),allocatable :: gwavef_alltoall(:,:)
 real(dp),allocatable :: swavef_alltoall(:,:)
 real(dp),allocatable :: tmp_ffnl_gather(:,:,:,:),tmp_kinpw_gather(:)
 real(dp),allocatable :: tmp_ph3d_gather(:,:,:)
 real(dp),pointer :: ewavef_alltoall_sym(:,:),gvnlc_alltoall_sym(:,:)
 real(dp),pointer :: gwavef_alltoall_sym(:,:),swavef_alltoall_sym(:,:)
 real(dp),pointer :: swavef_alltoall_sym_tmp(:,:)
 type(cprj_type) :: cwaveprj_alltoall_dum(1,1)

! *************************************************************************

 call timab(630,1,tsec)
 call timab(631,3,tsec)

 nproc_band = mpi_enreg%nproc_band
 nproc_fft  = mpi_enreg%nproc_fft
 bandpp     = mpi_enreg%bandpp

 flag_inv_sym = ((gs_hamk%ngfft(7)==401) .and. (gs_hamk%istwf_k==2))

 if (flag_inv_sym) then
   istwf_k         = 1
   gs_hamk%istwf_k = 1
   if (modulo(bandpp,2)==0) then
     bandpp_sym   = bandpp/2
   else
     bandpp_sym   = bandpp
   end if
 end if
!====================================================================================

 tim_getghc=6 ; cpopt=-1
 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_band)
 ikpt_this_proc=mpi_enreg%tab_kpt_distrib(ikpt)

 ABI_ALLOCATE(sendcountsloc,(nproc_band))
 ABI_ALLOCATE(sdisplsloc   ,(nproc_band))
 ABI_ALLOCATE(recvcountsloc,(nproc_band))
 ABI_ALLOCATE(rdisplsloc   ,(nproc_band))

 recvcounts   =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts(:)
 sendcounts   =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts(:)
 rdispls      =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%rdispls   (:)
 sdispls      =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%sdispls   (:)
 ndatarecv    =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%ndatarecv

 kg_k_gather           =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather(:,:)
 kinpw_gather          =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%kinpw_gather(:)
 ffnl_gather           =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%ffnl_gather(:,:,:,:)
 ph3d_gather           =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%ph3d_gather(:,:,:)
 gs_hamk%gbound(:,:)   = mpi_enreg%bandfft_kpt(ikpt_this_proc)%gbound(:,:)

 if (flag_inv_sym ) then
   idatarecv0           =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%idatarecv0
   ndatarecv_tot        =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%ndatarecv_tot
   ndatasend_sym        =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%ndatasend_sym
   kg_k_gather_sym      =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather_sym(:,:)
   rdispls_sym          =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%rdispls_sym(:)
   recvcounts_sym       =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym(:)
   recvcounts_sym_tot   =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym_tot(:)
   sdispls_sym          =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%sdispls_sym(:)
   sendcounts_sym       =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym(:)
   sendcounts_sym_all   =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym_all(:)
   tab_proc             =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%tab_proc(:)
 end if
 iscalc=(sij_opt+1)/2  ! 0 if S not calculated, 1 otherwise
 nbval=(ndatarecv*nspinor*bandpp)*iscalc

 ABI_ALLOCATE(cwavef_alltoall,(2,ndatarecv*nspinor*bandpp))
 ABI_ALLOCATE(gwavef_alltoall,(2,ndatarecv*nspinor*bandpp))
 ABI_ALLOCATE(swavef_alltoall,(2,ndatarecv*nspinor*bandpp))
 ABI_ALLOCATE(gvnlc_alltoall,(2,ndatarecv*nspinor*bandpp))
 swavef_alltoall(:,:)=zero
 gvnlc_alltoall(:,:)=zero

 recvcountsloc(:)=recvcounts(:)*2*nspinor*bandpp
 rdisplsloc(:)=rdispls(:)*2*nspinor*bandpp
 sendcountsloc(:)=sendcounts(:)*2*nspinor
 sdisplsloc(:)=sdispls(:)*2*nspinor
 call timab(631,2,tsec)

 call timab(545,3,tsec)
 call xalltoallv_mpi(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall,&
& recvcountsloc,rdisplsloc,spaceComm,ier)
 call timab(545,2,tsec)

 if(gs_hamk%istwf_k==2) then
   old_me_g0=mpi_enreg%me_g0
   if (mpi_enreg%me_fft==0) then
     mpi_enreg%me_g0=1
   else
     mpi_enreg%me_g0=0
   end if
 end if
!====================================================================
 if ((.not.(flag_inv_sym)) .and. (bandpp==1)) then

   if (mpi_enreg%paral_spin==0.and.nspinor==2)then
     call timab(632,3,tsec)
!    Sort to have all nspinor=1 first, then all nspinor=2
     call prep_sort_wavef_spin(nproc_band,nspinor,ndatarecv,recvcounts,rdispls,index_wavef_spband)
     cwavef_alltoall(:,:)=cwavef_alltoall(:,index_wavef_spband)
     call timab(632,2,tsec)
   end if

   call timab(635,3,tsec)
   if(present(vxctaulocal))then
     call getghc(cpopt,cwavef_alltoall,cwaveprj_alltoall_dum,dimffnl,ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall,swavef_alltoall(:,1:nbval),gs_hamk,gvnlc_alltoall,kg_k_gather,&
&     kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,&
&     ndatarecv,nspinor,ntypat,&
&     nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
   else
     call getghc(cpopt,cwavef_alltoall,cwaveprj_alltoall_dum,dimffnl,ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall,swavef_alltoall(:,1:nbval),gs_hamk,gvnlc_alltoall,kg_k_gather,&
&     kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,&
&     ndatarecv,nspinor,ntypat,&
&     nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal)
   end if
   call timab(635,2,tsec)

   if (mpi_enreg%paral_spin==0.and.nspinor==2)then
     call timab(634,3,tsec)
     gwavef_alltoall(:,index_wavef_spband)=gwavef_alltoall(:,:)
     if (sij_opt==1) swavef_alltoall(:,index_wavef_spband)=swavef_alltoall(:,:)
     gvnlc_alltoall(:,index_wavef_spband)=gvnlc_alltoall(:,:)
     ABI_DEALLOCATE(index_wavef_spband)
     call timab(634,2,tsec)
   end if

 else if ((.not.(flag_inv_sym)) .and. (bandpp>1)) then

!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------

   call timab(632,3,tsec)
   call prep_index_wavef_bandpp(nproc_band,bandpp,&
&   nspinor,ndatarecv, recvcounts,rdispls, index_wavef_band)

!  -------------------------------------------------------
!  Sorting of the waves functions below bandpp
!  -------------------------------------------------------
   cwavef_alltoall(:,:) = cwavef_alltoall(:,index_wavef_band)
   call timab(632,2,tsec)

!  ----------------------
!  Fourier transformation
!  ----------------------
   call timab(636,3,tsec)
   if (gs_hamk%ngfft(7)/=401.or.gs_hamk%use_gpu_cuda==1) then
     if(present(vxctaulocal))then
       call getghc(cpopt,cwavef_alltoall,cwaveprj_alltoall_dum,dimffnl,ffnl_gather,dtfil%filstat,&
&       gwavef_alltoall,swavef_alltoall,gs_hamk,gvnlc_alltoall,kg_k_gather,&
&       kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,bandpp,&
&       ndatarecv,nspinor,ntypat,&
&       nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
     else
       call getghc(cpopt,cwavef_alltoall,cwaveprj_alltoall_dum,dimffnl,ffnl_gather,dtfil%filstat,&
&       gwavef_alltoall,swavef_alltoall,gs_hamk,gvnlc_alltoall,kg_k_gather,&
&       kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,bandpp,&
&       ndatarecv,nspinor,ntypat,&
&       nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal)
     end if
   else
     do iibandpp=1,bandpp
       binf = nspinor* ndatarecv*(iibandpp-1)+1
       bsup = nspinor* ndatarecv * iibandpp
       if(present(vxctaulocal))then
         call getghc(cpopt,cwavef_alltoall(:,binf:bsup),cwaveprj_alltoall_dum,dimffnl,&
&         ffnl_gather,dtfil%filstat,gwavef_alltoall(:,binf:bsup),&
&         swavef_alltoall(:,binf*iscalc:bsup*iscalc),gs_hamk,&
&         gvnlc_alltoall(:,binf: bsup),kg_k_gather,kinpw_gather,lambda,lmnmax,&
&         matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,ndatarecv,nspinor,ntypat,&
&         nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
       else
         call getghc(cpopt,cwavef_alltoall(:,binf:bsup),cwaveprj_alltoall_dum,dimffnl,&
&         ffnl_gather,dtfil%filstat,gwavef_alltoall(:,binf:bsup),&
&         swavef_alltoall(:,binf*iscalc:bsup*iscalc),gs_hamk,&
&         gvnlc_alltoall(:,binf: bsup),kg_k_gather,kinpw_gather,lambda,lmnmax,&
&         matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,ndatarecv,nspinor,ntypat,&
&         nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal)
       end if
     end do
   end if
   call timab(636,2,tsec)

!  -----------------------------------------------------
!  Sorting of waves functions below the prossecors
!  -----------------------------------------------------
   call timab(634,3,tsec)
   cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)
   gwavef_alltoall(:,index_wavef_band) = gwavef_alltoall(:,:)
   if (sij_opt==1) swavef_alltoall(:,index_wavef_band) = swavef_alltoall(:,:)
   gvnlc_alltoall(:,index_wavef_band)  = gvnlc_alltoall(:,:)

   ABI_DEALLOCATE(index_wavef_band)
   call timab(634,2,tsec)


 else if (flag_inv_sym) then

!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------
   call timab(632,3,tsec)
   call prep_index_wavef_bandpp(nproc_band,bandpp,&
&   nspinor,ndatarecv,&
&   recvcounts,rdispls,&
&   index_wavef_band)

!  -------------------------------------------------------
!  Sorting of de the waves functions below bandpp
!  -------------------------------------------------------
   cwavef_alltoall(:,:) = cwavef_alltoall(:,index_wavef_band)

!  ------------------------------------------------------------
!  We associate the waves functions by two
!  ------------------------------------------------------------
   call prep_wavef_sym_do(mpi_enreg,bandpp,nspinor,&
   ndatarecv,&
   ndatarecv_tot,ndatasend_sym,tab_proc,&
   cwavef_alltoall,&
   sendcounts_sym,sdispls_sym,&
   recvcounts_sym,rdispls_sym,&
   ewavef_alltoall_sym,&
   index_wavef_send)

!  ------------------------------------------------------------
!  Allocation
!  ------------------------------------------------------------
   ABI_ALLOCATE(gwavef_alltoall_sym,(2,ndatarecv_tot*bandpp_sym))
   ABI_ALLOCATE(swavef_alltoall_sym,(2,(ndatarecv_tot*bandpp_sym)*iscalc))
   ABI_ALLOCATE(gvnlc_alltoall_sym ,(2,ndatarecv_tot*bandpp_sym))

   gwavef_alltoall_sym(:,:)=zero
   swavef_alltoall_sym(:,:)=zero
   gvnlc_alltoall_sym(:,:)=zero

   ABI_ALLOCATE(tmp_ffnl_gather ,(ndatarecv_tot,dimffnl,lmnmax,ntypat))
   ABI_ALLOCATE(tmp_kinpw_gather,(ndatarecv_tot))
   ABI_ALLOCATE(tmp_ph3d_gather ,(2,ndatarecv_tot,matblk))
   call timab(632,2,tsec)

!  ------------------------------------------------------------
!  Fourier calculcation
!  ------------------------------------------------------------
   call timab(637,3,tsec)
   if(gs_hamk%use_gpu_cuda==1) then
     if(present(vxctaulocal))then
       call getghc(cpopt,ewavef_alltoall_sym,&
&       cwaveprj_alltoall_dum,dimffnl,tmp_ffnl_gather,dtfil%filstat,&
&       gwavef_alltoall_sym, swavef_alltoall_sym,gs_hamk,gvnlc_alltoall_sym,&
&       kg_k_gather_sym,tmp_kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,&
&       bandpp_sym,ndatarecv_tot,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,tmp_ph3d_gather,prtvol,&
&       sij_opt,tim_getghc,1,vlocal,vxctaulocal=vxctaulocal)
     else
       call getghc(cpopt,ewavef_alltoall_sym,&
&       cwaveprj_alltoall_dum,dimffnl,tmp_ffnl_gather,dtfil%filstat,&
&       gwavef_alltoall_sym, swavef_alltoall_sym,gs_hamk,gvnlc_alltoall_sym,&
&       kg_k_gather_sym,tmp_kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,&
&       bandpp_sym,ndatarecv_tot,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,tmp_ph3d_gather,prtvol,&
&       sij_opt,tim_getghc,1,vlocal)
     end if
   else
     do iibandpp=1,bandpp_sym
       if (iscalc==1) then
         swavef_alltoall_sym_tmp => &
&         swavef_alltoall_sym(:,(ndatarecv_tot*(iibandpp-1)+1):ndatarecv_tot* iibandpp)
       else
         swavef_alltoall_sym_tmp => swavef_alltoall_sym
       end if
       if(present(vxctaulocal))then
         call getghc(cpopt,ewavef_alltoall_sym(:,(ndatarecv_tot*(iibandpp-1))+1:(ndatarecv_tot*iibandpp)),&
&         cwaveprj_alltoall_dum,dimffnl,tmp_ffnl_gather,dtfil%filstat,&
&         gwavef_alltoall_sym(:,(ndatarecv_tot*(iibandpp-1))+1:(ndatarecv_tot*iibandpp)),&
&         swavef_alltoall_sym_tmp,&
&         gs_hamk, gvnlc_alltoall_sym(:,(ndatarecv_tot*(iibandpp-1))+1:(ndatarecv_tot*iibandpp)),&
&         kg_k_gather_sym,tmp_kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,&
&         1,ndatarecv_tot,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,tmp_ph3d_gather,prtvol,&
&         sij_opt,tim_getghc,1,vlocal,vxctaulocal=vxctaulocal)
       else
         call getghc(cpopt,ewavef_alltoall_sym(:,(ndatarecv_tot*(iibandpp-1))+1:(ndatarecv_tot*iibandpp)),&
&         cwaveprj_alltoall_dum,dimffnl,tmp_ffnl_gather,dtfil%filstat,&
&         gwavef_alltoall_sym(:,(ndatarecv_tot*(iibandpp-1))+1:(ndatarecv_tot*iibandpp)),&
&         swavef_alltoall_sym_tmp,&
&         gs_hamk, gvnlc_alltoall_sym(:,(ndatarecv_tot*(iibandpp-1))+1:(ndatarecv_tot*iibandpp)),&
&         kg_k_gather_sym,tmp_kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,&
&         1,ndatarecv_tot,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,tmp_ph3d_gather,prtvol,&
&         sij_opt,tim_getghc,1,vlocal)
       end if
     end do
   end if
   call timab(637,2,tsec)

   call timab(633,3,tsec)
   ABI_DEALLOCATE(tmp_ffnl_gather)
   ABI_DEALLOCATE(tmp_kinpw_gather)
   ABI_DEALLOCATE(tmp_ph3d_gather)

!  ------------------------------------------------------------
!  We dissociate each wave function in two waves functions
!  gwavef is classed below of bandpp
!  ------------------------------------------------------------
   call prep_wavef_sym_undo(mpi_enreg,bandpp,nspinor,&
   ndatarecv,&
   ndatarecv_tot,ndatasend_sym,idatarecv0,&
   gwavef_alltoall,&
   sendcounts_sym,sdispls_sym,&
   recvcounts_sym,rdispls_sym,&
   gwavef_alltoall_sym,&
   index_wavef_send)

   ABI_DEALLOCATE(ewavef_alltoall_sym)
   ABI_DEALLOCATE(index_wavef_send)
   ABI_DEALLOCATE(gwavef_alltoall_sym)
   ABI_DEALLOCATE(swavef_alltoall_sym)
   ABI_DEALLOCATE(gvnlc_alltoall_sym)

!  -------------------------------------------
!  We call getghc for the calcul of ffnl,...
!  --------------------------------------------
   gs_hamk%istwf_k=2

   old_me_g0=mpi_enreg%me_g0
   if (mpi_enreg%me_fft==0) then
     mpi_enreg%me_g0=1
   else
     mpi_enreg%me_g0=0
   end if
   call timab(633,2,tsec)

   call timab(638,3,tsec)
   if(present(vxctaulocal))then
     call getghc(cpopt,cwavef_alltoall,cwaveprj_alltoall_dum,dimffnl,ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall,swavef_alltoall,gs_hamk,gvnlc_alltoall,kg_k_gather,&
&     kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,bandpp,&
&     ndatarecv,nspinor,ntypat,&
&     nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,2,vlocal,vxctaulocal=vxctaulocal)
   else
     call getghc(cpopt,cwavef_alltoall,cwaveprj_alltoall_dum,dimffnl,ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall,swavef_alltoall,gs_hamk,gvnlc_alltoall,kg_k_gather,&
&     kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,bandpp,&
&     ndatarecv,nspinor,ntypat,&
&     nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,2,vlocal)
   end if

   call timab(638,2,tsec)
   call timab(634,3,tsec)
   mpi_enreg%me_g0=old_me_g0

   gs_hamk%istwf_k=1

!  -------------------------------------------------------
!  Sorting of waves functions below the processors
!  -------------------------------------------------------
   cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)
   gwavef_alltoall(:,index_wavef_band) = gwavef_alltoall(:,:)
   if (sij_opt==1) swavef_alltoall(:,index_wavef_band) = swavef_alltoall(:,:)
   gvnlc_alltoall(:,index_wavef_band)  = gvnlc_alltoall(:,:)
   ABI_DEALLOCATE(index_wavef_band)
   call timab(634,2,tsec)

 end if
!====================================================================

 if (gs_hamk%istwf_k==2) mpi_enreg%me_g0=old_me_g0
 call timab(545,3,tsec)
 if (sij_opt==1) then
   call xalltoallv_mpi(swavef_alltoall,recvcountsloc,rdisplsloc,swavef,&
&   sendcountsloc,sdisplsloc,spaceComm,ier)
 end if

 call xalltoallv_mpi(gvnlc_alltoall,recvcountsloc,rdisplsloc,gvnlc,&
& sendcountsloc,sdisplsloc,spaceComm,ier)

 call xalltoallv_mpi(gwavef_alltoall,recvcountsloc,rdisplsloc,gwavef,&
& sendcountsloc,sdisplsloc,spaceComm,ier)

 call timab(545,2,tsec)

!====================================================================
 if (flag_inv_sym) then
   istwf_k         = 2
   gs_hamk%istwf_k = 2
 end if
!====================================================================

 mpi_enreg%paral_level= old_paral_level
 ABI_DEALLOCATE(sendcountsloc)
 ABI_DEALLOCATE(sdisplsloc)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)
 ABI_DEALLOCATE(cwavef_alltoall)
 ABI_DEALLOCATE(gwavef_alltoall)
 ABI_DEALLOCATE(gvnlc_alltoall)
 ABI_DEALLOCATE(swavef_alltoall)
 call timab(630,2,tsec)

end subroutine prep_getghc
!!***
