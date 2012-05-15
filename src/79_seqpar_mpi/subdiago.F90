!{\src2tex{textfont=tt}}
!!****f* ABINIT/subdiago
!! NAME
!! subdiago
!!
!! FUNCTION
!! This routine diagonalizes the Hamiltonian in the eigenfunction subspace
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  istwf_k=input parameter that describes the storage of wfs
!!  mcg=second dimension of the cg array
!!  mgsc=second dimension of the gsc array
!!  mpi_enreg=informations about MPI parallelization
!!  nband_k=number of bands at this k point for that spin polarization
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  subham(nband_k*(nband_k+1))=Hamiltonian expressed in sthe WFs subspace
!!  subovl(nband_k*(nband_k+1)*use_subovl)=overlap matrix expressed in sthe WFs subspace
!!  use_subovl=1 if the overlap matrix is not identity in WFs subspace
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  eig_k(nband_k)=array for holding eigenvalues (hartree)
!!  evec(2*nband_k,nband_k)=array for holding eigenvectors
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=wavefunctions
!!  gsc(2,mgsc)=<g|S|c> matrix elements (S=overlap)
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!      compute_eigen_problem,compute_generalized_eigen_problem,dcopy
!!      destruction_matrix_scalapack,dgemm,dspev,dspgv,end_scalapack,hermit
!!      init_matrix_scalapack,init_scalapack,leave_new,matrix_from_global
!!      matrix_to_global,matrix_to_reference,mpi_allreduce,mpi_bcast,normev
!!      timab,wrtout,xcomm_world,zgemm,zhpev,zhpgv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


 subroutine subdiago(cg,eig_k,evec,gsc,icg,igsc,istwf_k,&
&                    mcg,mgsc,mpi_enreg,nband_k,npw_k,nspinor,paral_kgb,&
&                    subham,subovl,use_subovl,usepaw)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_abitypes
 use defs_scalapack
#if defined HAVE_MPI2
 use mpi
#endif
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subdiago'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: icg,igsc,istwf_k,mcg,mgsc,nband_k,npw_k
 integer,intent(in) :: nspinor,paral_kgb,use_subovl,usepaw
 type(MPI_type),intent(inout) :: mpi_enreg
 real(dp),intent(inout) :: subham(nband_k*(nband_k+1)),subovl(nband_k*(nband_k+1)*use_subovl)
 real(dp),intent(out) :: eig_k(nband_k),evec(2*nband_k,nband_k)
 real(dp),intent(inout) :: cg(2,mcg),gsc(2,mgsc)

!Local variables-------------------------------
 integer :: iband,ii,ierr
 character(len=500) :: message
 real(dp) :: tsec(2)
 real(dp),allocatable :: work(:,:),zhpev1(:,:),zhpev2(:),dspev1(:)
 integer :: rvectsize,vectsize
 real(dp),allocatable :: blockvectora(:,:),blockvectorb(:,:),blockvectorc(:,:)
 integer :: cgindex,gscindex
!integer :: iwavef

#if defined HAVE_LINALG_SCALAPACK
!Scalapack variables----------------------------
 TYPE(matrix_scalapack)    :: sca_subham,sca_subovl,sca_evec
 TYPE(processor_scalapack) :: processor

 real(dp),dimension(:,:),allocatable ::tmp_evec
 INTEGER         :: communicator
#endif

!no_abirules
!Function definitions
 cgindex(iband) =npw_k*nspinor*(iband-1)+icg+1
 gscindex(iband)=npw_k*nspinor*(iband-1)+igsc+1


! *********************************************************************

!DEBUG
!write(std_out,*)' subdiago : enter '
!ENDDEBUG

 if(paral_kgb<0)then
   write(message,'(4a)' )ch10,&
&   ' cgwf : BUG ',ch10,&
&   '   paral_kgb should be positive '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 rvectsize = npw_k*nspinor
 if (mpi_enreg%me_g0 == 1) then
   vectsize=2*npw_k*nspinor-1
 else
   vectsize=2*npw_k*nspinor
 end if

!Impose Hermiticity on diagonal elements of subham (and subovl, if needed)
 call hermit(subham,subham,ierr,nband_k)
 if (use_subovl==1) then
   call hermit(subovl,subovl,ierr,nband_k)
 end if

!Diagonalize the Hamitonian matrix

#if defined HAVE_MPI
!============================================

#if defined HAVE_LINALG_SCALAPACK
!============================================
 if(paral_kgb == 1) then

!  call timab(570,1,tsec)

!  ===============================
!  INITIALISATION WORK VARIABLE
!  ===============================
   ABI_ALLOCATE(tmp_evec,(2*nband_k,nband_k))
   tmp_evec(:,:)=0._DP

!  ============================
!  INITIALISATION COMMUNICATOR
!  ===========================
   if (mpi_enreg%paralbd <=1) then

     if (mpi_enreg%paral_compil_kpt==1) then
       communicator  = mpi_enreg%commcart_3d
     else
       call xcomm_world(mpi_enreg,communicator)
     end if

   else
     communicator  = mpi_enreg%kpt_comm(mpi_enreg%num_group)
   end if

!  ========================
!  INITIALISATION SCALAPACK
!  ========================
   call init_scalapack(processor,communicator)

!  ================================
!  INITIALISATION SCALAPACK MATRIX
!  ================================
   call init_matrix_scalapack(sca_subham,nband_k,nband_k,processor,istwf_k,10)
   call init_matrix_scalapack(sca_evec,nband_k,nband_k,processor,istwf_k,10)

!  ==============================
!  FILLING SCALAPACK MATRIX
!  ==============================
   call matrix_from_global(sca_subham,subham,istwf_k)

   if (use_subovl==1) then
     call init_matrix_scalapack(sca_subovl,nband_k,nband_k,processor,istwf_k,10)
     call matrix_from_global(sca_subovl,subovl,istwf_k)
   end if


   if (use_subovl==1) then
!    ================================
!    COMPUTE EIGEN VALUES AND VECTORS
!    FOR THE EIGEN PROBLEM :
!    A * X = lambda * B * X
!    ================================
!    write(std_out,*) 'I am using Scalapack : compute eigen values and vectors for generalized eigen problem'
     call compute_generalized_eigen_problem(processor,sca_subham,sca_subovl,&
&     sca_evec,eig_k,&
&     communicator,istwf_k)

     call matrix_to_global(sca_subham,subham,istwf_k)
     call matrix_to_global(sca_subovl,subovl,istwf_k)

   else
!    ================================
!    COMPUTE EIGEN VALUES AND VECTORS
!    ================================
!    write(std_out,*) 'I am using Scalapack : compute eigen values and vectors for eigen problem'
     call compute_eigen_problem(processor,sca_subham,&
&     sca_evec,eig_k,&
&     communicator,istwf_k)

     call matrix_to_global(sca_subham,subham,istwf_k)

   end if

!  ==============================
!  CONCATENATE EIGEN VECTORS
!  ==============================
   call matrix_to_reference(sca_evec,tmp_evec,istwf_k)

   CALL MPI_ALLREDUCE(tmp_evec, evec, 2*nband_k*nband_k, MPI_DOUBLE_PRECISION, &
   MPI_SUM, communicator,ierr)

!  ====================================
!  DESTRUCTION SCALAPACK AND TMP MATRICES
!  ====================================
   CALL destruction_matrix_scalapack(sca_subham)
   CALL destruction_matrix_scalapack(sca_evec)

   if (use_subovl==1) then
     CALL destruction_matrix_scalapack(sca_subovl)
   end if

!  ===========================
!  CLOSE SCALAPACK
!  ===========================
   CALL end_scalapack(processor)

   ABI_DEALLOCATE(tmp_evec)
!  call timab(570,2,tsec)

 else ! paral_kgb=0

#endif
!  END HAVE_LINALG_SCALAPACK
!  =====================

   if ((mpi_enreg%paralbd <=1) .or. ((mpi_enreg%paralbd >1) .and. &
&   (mpi_enreg%me_group==0))) then

#endif
!    END MPI
!    ===============

!    ------------------
     if (istwf_k==2) then
!      ------------------

       ABI_ALLOCATE(dspev1,(3*nband_k))
       evec(:,:)=0._dp

       if (use_subovl==1) then

!        write(std_out,*) ' I am using DSPGV'
         call DSPGV(1,'V','U',nband_k,&
&         subham(1:nband_k*(nband_k+1):2),&
&         subovl(1:nband_k*(nband_k+1):2),&
&         eig_k,&
&         evec(1:2*nband_k:2,:),&
&         nband_k,&
&         dspev1,ierr)
       else

!        write(std_out,*) ' I am using DSPEV'
         call DSPEV ('V','U',nband_k,&
&         subham(1:nband_k*(nband_k+1):2),&
&         eig_k,&
&         evec(1:2*nband_k:2,:),&
&         nband_k,&
&         dspev1,ierr)
       end if

       ABI_DEALLOCATE(dspev1)

!      -----------------------
     else
!      -----------------------

       ABI_ALLOCATE(zhpev1,(2,2*nband_k-1))
       ABI_ALLOCATE(zhpev2,(3*nband_k-2))

       if (use_subovl==1) then
!        write(std_out,*) 'I am using ZHPGV'
         call ZHPGV(1,'V','U',nband_k,subham,subovl,eig_k,evec,nband_k,&
&         zhpev1,zhpev2,ierr)
       else
!        write(std_out,*) 'I am using ZHPEV'
         call ZHPEV ('V','U',nband_k,subham,eig_k,evec,nband_k,zhpev1,&
&         zhpev2,ierr)
       end if

       ABI_DEALLOCATE(zhpev1)
       ABI_DEALLOCATE(zhpev2)

!      ----------------
     end if
!    ----------------

#if defined HAVE_MPI
!    ===============
   end if
   if (mpi_enreg%paralbd >1) then
     call timab(48,1,tsec)
     call MPI_BCAST(evec,2*nband_k*nband_k, &
&     MPI_DOUBLE_PRECISION,0,mpi_enreg%kpt_comm(mpi_enreg%num_group),ierr)
     call timab(48,2,tsec)
   end if
#endif
!  END MPI
!  ===============


#if defined HAVE_LINALG_SCALAPACK
 end if ! paral_kgb
#endif


!DEBUG
!write(std_out,*)' subdiago : after zhpev '
!stop
!ENDDEBUG

!Normalize each eigenvector and set phase:
 call normev(evec,nband_k,nband_k)

!if(prtvol==-level)then
!write(message,'(a)')&
!&  ' subdiago : iband band  evec(re:im)'
!call wrtout(std_out,message,'PERS')
!do iband=1,nband_k
!do ii=1,nband_k
!write(message,'(2i5,2es16.6)')&
!&    iband,ii,evec(2*ii-1,iband),evec(2*ii,iband)
!call wrtout(std_out,message,'PERS')
!end do
!end do
!end if

 if(istwf_k==2)then
   do iband=1,nband_k
     do ii=1,nband_k
       if(abs(evec(2*ii,iband))>1.0d-10)then
         write(message,'(a,a,a,a,2i5,2es16.6,a,a)')ch10,&
&         ' subdiago : BUG ',&
&         '  For istwf_k=2, observed the following element of evec :',ch10,&
&         iband,ii,evec(2*ii-1,iband),evec(2*ii,iband),ch10,&
&         '  with a non-negligible imaginary part.'
         call wrtout(std_out,message,'PERS')
         call leave_new('PERS')
       end if
     end do
   end do
 end if

!Carry out rotation of bands C(G,n) according to evecs:

!==============================
!SDIROT --> ZGEMM if istwfk==1
!--> DGEMM if istwfk==2
!==============================

!-----------------
 if (istwf_k==2)then
!  -----------------

   ABI_ALLOCATE(blockvectora,(vectsize,nband_k))
   ABI_ALLOCATE(blockvectorb,(nband_k,nband_k))
   ABI_ALLOCATE(blockvectorc,(vectsize,nband_k))

   do iband=1,nband_k
     if (mpi_enreg%me_g0 == 1) then
       call dcopy(1          ,cg(1,cgindex(iband))                     ,1,blockvectora(1                   ,iband),1)
       call dcopy(rvectsize-1,cg(1,cgindex(iband)+1:cgindex(iband+1)-1),1,blockvectora(2:rvectsize         ,iband),1)
       call dcopy(rvectsize-1,cg(2,cgindex(iband)+1:cgindex(iband+1)-1),1,blockvectora(rvectsize+1:vectsize,iband),1)
     else
       call dcopy(rvectsize,cg(1,cgindex(iband):cgindex(iband+1)-1),1,blockvectora(1:rvectsize         ,iband),1)
       call dcopy(rvectsize,cg(2,cgindex(iband):cgindex(iband+1)-1),1,blockvectora(rvectsize+1:vectsize,iband),1)
     end if

     call dcopy(nband_k,evec(2*iband-1,1:nband_k),1,blockvectorb(iband,1:nband_k),1)
   end do

!  write(std_out,*) 'I am using DGEMM'
   call dgemm('N','N',vectsize,nband_k,nband_k,&
&   one, &
&   blockvectora,vectsize, &
&   blockvectorb,nband_k,&
&   zero,blockvectorc,vectsize)

   do iband=1,nband_k
     if (mpi_enreg%me_g0 == 1) then
       call dcopy(1        ,blockvectorc(1                   ,iband),1,cg(1,cgindex(iband))                     ,1)
       call dcopy(rvectsize-1,blockvectorc(2:rvectsize         ,iband),1,cg(1,cgindex(iband)+1:cgindex(iband+1)-1),1)
       call dcopy(rvectsize-1,blockvectorc(rvectsize+1:vectsize,iband),1,cg(2,cgindex(iband)+1:cgindex(iband+1)-1),1)
     else
       call dcopy(rvectsize,blockvectorc(1:rvectsize         ,iband),1,cg(1,cgindex(iband):cgindex(iband+1)-1),1)
       call dcopy(rvectsize,blockvectorc(rvectsize+1:vectsize,iband),1,cg(2,cgindex(iband):cgindex(iband+1)-1),1)
     end if
   end do

!  If paw, musb also rotate S.C(G,n):
   if (usepaw==1) then

     do iband=1,nband_k
       if (mpi_enreg%me_g0 == 1) then
         call dcopy(1          ,gsc(1,gscindex(iband))                      ,1,blockvectora(1                   ,iband),1)
         call dcopy(rvectsize-1,gsc(1,gscindex(iband)+1:gscindex(iband+1)-1),1,blockvectora(2:rvectsize         ,iband),1)
         call dcopy(rvectsize-1,gsc(2,gscindex(iband)+1:gscindex(iband+1)-1),1,blockvectora(rvectsize+1:vectsize,iband),1)
       else
         call dcopy(rvectsize,gsc(1,gscindex(iband):gscindex(iband+1)-1),1,blockvectora(1:rvectsize         ,iband),1)
         call dcopy(rvectsize,gsc(2,gscindex(iband):gscindex(iband+1)-1),1,blockvectora(rvectsize+1:vectsize,iband),1)
       end if
       call dcopy(nband_k,evec(2*iband-1,1:nband_k),1,blockvectorb(iband,1:nband_k),1)
     end do

     call dgemm('N','N',vectsize,nband_k,nband_k,&
&     one, &
&     blockvectora,vectsize, &
&     blockvectorb,nband_k, &
&     zero,blockvectorc,vectsize)

     do iband=1,nband_k
       if (mpi_enreg%me_g0 == 1) then
         call dcopy(1        ,blockvectorc(1                   ,iband),1,gsc(1,gscindex(iband))                      ,1)
         call dcopy(rvectsize-1,blockvectorc(2:rvectsize         ,iband),1,gsc(1,gscindex(iband)+1:gscindex(iband+1)-1),1)
         call dcopy(rvectsize-1,blockvectorc(rvectsize+1:vectsize,iband),1,gsc(2,gscindex(iband)+1:gscindex(iband+1)-1),1)
       else
         call dcopy(rvectsize,blockvectorc(1:rvectsize         ,iband),1,gsc(1,gscindex(iband):gscindex(iband+1)-1),1)
         call dcopy(rvectsize,blockvectorc(rvectsize+1:vectsize,iband),1,gsc(2,gscindex(iband):gscindex(iband+1)-1),1)
       end if
     end do

   end if

   ABI_DEALLOCATE(blockvectora)
   ABI_DEALLOCATE(blockvectorb)
   ABI_DEALLOCATE(blockvectorc)

!  -----------------------
 else
!  -----------------------

   ABI_ALLOCATE(work,(2,npw_k*nspinor*nband_k))
   work=zero

!  call sdirot(cg,evec,icg,mcg,nband_k,nband_k,npw_k*nspinor)
!  write(std_out,*) 'I am using ZGEMM'
   call zgemm('N','N',npw_k*nspinor,nband_k,nband_k,&
&   dcmplx(1._dp), &
&   cg(1,icg+1),npw_k*nspinor, &
&   evec,nband_k,&
&   dcmplx(0._dp), &
&   work,npw_k*nspinor)
   cg(:,1+icg:npw_k*nspinor*nband_k+icg)=work(:,:)

!  If paw, must also rotate S.C(G,n):
!  if (usepaw==1) then
!  call sdirot(gsc,evec,icg,mcg,nband_k,nband_k,npw_k*nspinor)
!  endif
   if (usepaw==1) then
     call zgemm('N','N',npw_k*nspinor,nband_k,nband_k,&
&     dcmplx(1._dp), &
&     gsc(1,igsc+1),npw_k*nspinor, &
&     evec,nband_k, &
&     dcmplx(zero), &
&     work,npw_k*nspinor)
     gsc(:,1+igsc:npw_k*nspinor*nband_k+igsc)=work(:,:)
   end if

   ABI_DEALLOCATE(work)

!  ----------------
 end if
!----------------

!DEBUG
!write(std_out,*)' subdiago : cg(1:2) for different bands (3) '
!do iband=1,nband_k
!iwavef=(iband-1)*npw_k+icg
!write(std_out,'(4es16.6)' )cg(1:2,1+iwavef:2+iwavef)
!end do
!ENDDEBUG

!DEBUG
!write(std_out,*)' subdiago : exit '
!stop
!ENDDEBUG

end subroutine subdiago
!!***
