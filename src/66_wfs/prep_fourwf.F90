!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_fourwf
!! NAME
!! prep_fourwf
!!
!! FUNCTION
!! this routine prepares the data to the call of fourwf.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FBottin,MT,GZ,FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blocksize= size of block for FFT
!!  cwavef(2,npw*ndat)=planewave coefficients of wavefunction (one spinorial component?).
!!  dtfil <type(datafiles_type)>=variables related to files
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  gvnlc=matrix elements <G|Vnonlocal|C>
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  mgfft=maximum size of 1d ffts
!!  mpi_enreg=informations about mpi parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in unit cell.
!!  nvloc=final dimension of vlocal (usually 1, but 4 for non-collinear)
!!  n4,n5,n6 used for dimensionning of vlocal
!!  prtvol=control print volume and debugging output
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!
!! OUTPUT
!!  gwavef=(2,npw*ndat)=matrix elements <G|H|C>.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      mkrho,vtowfk
!!
!! CHILDREN
!!      fourwf,gpu_fourwf,prep_index_wavef_bandpp,prep_wavef_sym_do
!!      prep_wavef_sym_undo,timab,xalltoallv_mpi,xcomm_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_fourwf(rhoaug,blocksize,cwavef,wfraug,gs_hamk,iblock,ikpt,istwf_k,&
& mgfft,mpi_enreg,nband_k,npw_k,n4,n5,n6,occ_k,paral_kgb,wtk)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_fourwf'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_66_wfs, except_this_one => prep_fourwf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: blocksize,iblock,ikpt,istwf_k,mgfft,n4,n5,n6,nband_k
 integer :: npw_k,paral_kgb
 real(dp) :: wtk
 type(gs_hamiltonian_type) :: gs_hamk
 type(mpi_type) :: mpi_enreg
!arrays
 real(dp) :: cwavef(2,npw_k*blocksize),occ_k(nband_k),rhoaug(n4,n5,n6)
 real(dp) :: wfraug(2,n4,n5,n6)

!Local variables-------------------------------
!local variables for mpialltoallv
!local variable for bandpp and inversion by symetry of time
!scalars
 integer :: bandpp,bandpp_sym,ier,iibandpp,ikpt_this_proc,ind_occ,ind_occ1,ind_occ2
 integer :: nproc_band,nproc_fft
 integer :: old_me_g0=0,old_paral_level,spaceComm=0,tim_fourwf
 real(dp) :: weight,weight1,weight2
 logical :: flag_inv_sym
!arrays
 integer,pointer :: idatarecv0,ndatarecv,ndatarecv_tot,ndatasend_sym
 integer,pointer :: kg_k_gather(:,:),kg_k_gather_sym(:,:)
 integer,pointer :: rdispls(:),rdispls_sym(:)
 integer,pointer :: recvcounts(:),recvcounts_sym(:),recvcounts_sym_tot(:)
 integer,pointer :: sdispls(:),sdispls_sym(:)
 integer,pointer :: sendcounts(:),sendcounts_sym(:),sendcounts_sym_all(:)
 integer,pointer :: tab_proc(:)
 integer,allocatable :: rdisplsloc(:)
 integer,allocatable :: recvcountsloc(:),sdisplsloc(:)
 integer,allocatable :: sendcountsloc(:)
 integer,pointer :: index_wavef_band(:),index_wavef_send(:)
 real(dp) :: dummy(2,1),tsec(2)
 real(dp),allocatable :: cwavef_alltoall(:,:)
 real(dp),allocatable :: weight_t(:),weight1_t(:),weight2_t(:)
 real(dp),pointer :: ewavef_alltoall_sym(:,:)
!no_abirules
!correspondence with abinit. here for real wf but in complex mode
!this is the index of a given band

! *************************************************************************

!====================================================================================
 nproc_band = mpi_enreg%nproc_band
 nproc_fft  = mpi_enreg%nproc_fft
 bandpp     = mpi_enreg%bandpp

 flag_inv_sym = ((gs_hamk%ngfft(7)==401) .and. (istwf_k==2))

 if (flag_inv_sym) then
   istwf_k         = 1
   if (modulo(bandpp,2)==0) then
     bandpp_sym   = bandpp/2
   else
     bandpp_sym   = bandpp
   end if
 end if
!====================================================================================

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

 kg_k_gather          =>mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather(:,:)
 gs_hamk%gbound(:,:)  = mpi_enreg%bandfft_kpt(ikpt_this_proc)%gbound(:,:)

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

 ABI_ALLOCATE(cwavef_alltoall,(2,ndatarecv*bandpp))

 recvcountsloc(:)=recvcounts(:)*2*bandpp
 rdisplsloc(:)=rdispls(:)*2*bandpp
 sendcountsloc(:)=sendcounts(:)*2
 sdisplsloc(:)=sdispls(:)*2

 call timab(547,1,tsec)
 call xalltoallv_mpi(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall,&
& recvcountsloc,rdisplsloc,spaceComm,ier)
 call timab(547,2,tsec)

!If me_fft==0, I have the G=0 vector, but keep for the record the old value
 if (mpi_enreg%me_fft==0) then
   old_me_g0=mpi_enreg%me_g0
   mpi_enreg%me_g0=1
 end if

 tim_fourwf=16

!====================================================================
 if ((.not.(flag_inv_sym)) .and. (bandpp==1)) then

!  Compute the index of the band
   ind_occ = (iblock-1)*blocksize + mpi_enreg%me_band + 1

   if(abs(occ_k(ind_occ)) >=tol8) then

!    Compute the weight of the band
     weight=occ_k(ind_occ)*wtk/gs_hamk%ucvol

     call fourwf(1,rhoaug,cwavef_alltoall,dummy,wfraug,&
&     gs_hamk%gbound,gs_hamk%gbound,&
&     istwf_k,kg_k_gather,kg_k_gather,mgfft,mpi_enreg,1,&
&     gs_hamk%ngfft,ndatarecv,1,n4,n5,n6,1,paral_kgb,tim_fourwf,weight,weight,&
&     use_gpu_cuda=gs_hamk%use_gpu_cuda)

   end if

 else if ((.not.(flag_inv_sym)) .and. (bandpp>1) ) then

!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------
   call prep_index_wavef_bandpp(nproc_band,bandpp,&
&   1,ndatarecv,&
&   recvcounts,rdispls,&
&   index_wavef_band)

!  -------------------------------------------------------
!  Sorting of the waves functions below bandpp
!  -------------------------------------------------------
   cwavef_alltoall(:,:) = cwavef_alltoall(:,index_wavef_band)

!  -------------------
!  Fourier calculation
!  -------------------

   if(gs_hamk%use_gpu_cuda==1) then
     ABI_ALLOCATE(weight_t,(bandpp))
     do iibandpp=1,bandpp
!      Compute the index of the band
       ind_occ = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + iibandpp
!      Compute the weight of the band
       weight_t(iibandpp)=occ_k(ind_occ)*wtk/gs_hamk%ucvol
       if(abs(occ_k(ind_occ)) < tol8) then
         weight_t(iibandpp) = 0
       end if
     end do
!    Accumulate time because it is not done in gpu_fourwf
     call timab(240+tim_fourwf,1,tsec)
#if defined HAVE_GPU_CUDA
     call gpu_fourwf(1,rhoaug,&
&     cwavef_alltoall,&
&     dummy,wfraug,&
&     gs_hamk%gbound,gs_hamk%gbound,&
&     istwf_k,kg_k_gather,kg_k_gather,mgfft,mpi_enreg,bandpp,&
&     gs_hamk%ngfft,ndatarecv,1,n4,n5,n6,1,paral_kgb,&
&     tim_fourwf,weight_t,weight_t)
#endif
     call timab(240+tim_fourwf,2,tsec)
     ABI_DEALLOCATE(weight_t)
   else

     do iibandpp=1,bandpp

!      Compute the index of the band
       ind_occ = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + iibandpp
!      Compute the weight of the band
       weight=occ_k(ind_occ)*wtk/gs_hamk%ucvol

       if(abs(occ_k(ind_occ)) >=tol8) then

         call fourwf(1,rhoaug,&
&         cwavef_alltoall(:,(ndatarecv*(iibandpp-1))+1:(ndatarecv*iibandpp)),&
&         dummy,wfraug,&
&         gs_hamk%gbound,gs_hamk%gbound,&
&         istwf_k,kg_k_gather,kg_k_gather,mgfft,mpi_enreg,1,&
&         gs_hamk%ngfft,ndatarecv,1,n4,n5,n6,1,paral_kgb,&
&         tim_fourwf,weight,weight)
       end if
     end do

   end if ! (use_gpu_cuda==1)

!  -----------------------------------------------------
!  Sorting of waves functions below the prossecors
!  -----------------------------------------------------
   cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)
   ABI_DEALLOCATE(index_wavef_band)


 else if (flag_inv_sym) then

!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------
   call prep_index_wavef_bandpp(nproc_band,bandpp,&
&   1,ndatarecv,&
&   recvcounts,rdispls,&
&   index_wavef_band)

!  -------------------------------------------------------
!  Sorting of de the waves functions below bandpp
!  -------------------------------------------------------
   cwavef_alltoall(:,:) = cwavef_alltoall(:,index_wavef_band)

!  ------------------------------------------------------------
!  We associate the waves functions by two
!  ------------------------------------------------------------
   call prep_wavef_sym_do(mpi_enreg,bandpp,1,&
   ndatarecv,&
   ndatarecv_tot,ndatasend_sym,tab_proc,&
   cwavef_alltoall,&
   sendcounts_sym,sdispls_sym,&
   recvcounts_sym,rdispls_sym,&
   ewavef_alltoall_sym,&
   index_wavef_send)

!  ------------------------------------------------------------
!  Fourier calculcation
!  ------------------------------------------------------------
   if(gs_hamk%use_gpu_cuda==1) then
     ABI_ALLOCATE(weight1_t,(bandpp_sym))
     ABI_ALLOCATE(weight2_t,(bandpp_sym))
     do iibandpp=1,bandpp_sym
!      Compute the index of the bands
       if (bandpp/=1) then
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp-1)
         ind_occ2 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp  )
       else
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + 1
         ind_occ2 = ind_occ1
       end if
!      Compute the weight of the band
       weight1_t(iibandpp) = occ_k(ind_occ1)*wtk/gs_hamk%ucvol
       weight2_t(iibandpp) = occ_k(ind_occ2)*wtk/gs_hamk%ucvol
     end do
     call timab(240+tim_fourwf,1,tsec)
#if defined HAVE_GPU_CUDA
     call gpu_fourwf(1,rhoaug,&
&     ewavef_alltoall_sym,&
&     dummy,wfraug,&
&     gs_hamk%gbound,gs_hamk%gbound,&
&     istwf_k,kg_k_gather_sym,kg_k_gather_sym,mgfft,mpi_enreg,bandpp_sym,&
&     gs_hamk%ngfft,ndatarecv_tot,1,n4,n5,n6,1,paral_kgb,&
&     tim_fourwf,weight1_t,weight2_t)
#endif
     call timab(240+tim_fourwf,2,tsec)
     ABI_DEALLOCATE(weight1_t)
     ABI_DEALLOCATE(weight2_t)
   else

     do iibandpp=1,bandpp_sym

!      Compute the index of the bands
       if (bandpp/=1) then
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp-1)
         ind_occ2 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp  )
       else
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + 1
         ind_occ2 = ind_occ1
       end if

       weight1 = occ_k(ind_occ1)*wtk/gs_hamk%ucvol
       weight2 = occ_k(ind_occ2)*wtk/gs_hamk%ucvol

       call fourwf(1,rhoaug,&
&       ewavef_alltoall_sym(:,(ndatarecv_tot*(iibandpp-1))+1:(ndatarecv_tot*iibandpp)),&
&       dummy,wfraug,&
&       gs_hamk%gbound,gs_hamk%gbound,&
&       istwf_k,kg_k_gather_sym,kg_k_gather_sym,mgfft,mpi_enreg,1,&
&       gs_hamk%ngfft,ndatarecv_tot,1,n4,n5,n6,1,paral_kgb,&
&       tim_fourwf,weight1,weight2)

     end do

   end if ! (use_gpu_cuda==1)

!  ------------------------------------------------------------
!  We dissociate each wave function in two waves functions
!  gwavef is classed below of bandpp
!  ------------------------------------------------------------
   call prep_wavef_sym_undo(mpi_enreg,bandpp,1,&
   ndatarecv,&
   ndatarecv_tot,ndatasend_sym,idatarecv0,&
   cwavef_alltoall,&
   sendcounts_sym,sdispls_sym,&
   recvcounts_sym,rdispls_sym,&
   ewavef_alltoall_sym,&
   index_wavef_send)

   ABI_DEALLOCATE(ewavef_alltoall_sym)
   ABI_DEALLOCATE(index_wavef_send)

!  -------------------------------------------------------
!  Soritng of waves functions below the processors
!  -------------------------------------------------------
   cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)

   ABI_DEALLOCATE(index_wavef_band)

 end if

!====================================================================
 if (flag_inv_sym) then
   istwf_k         = 2
 end if
!====================================================================

 if (mpi_enreg%me_fft==0) mpi_enreg%me_g0=old_me_g0
 mpi_enreg%paral_level= old_paral_level
 ABI_DEALLOCATE(sendcountsloc)
 ABI_DEALLOCATE(sdisplsloc)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)
 ABI_DEALLOCATE(cwavef_alltoall)

end subroutine prep_fourwf
!!***
