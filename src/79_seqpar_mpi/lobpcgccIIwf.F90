!{\src2tex{textfont=tt}}
!!****f* abinit/lobpcgccIIwf
!! NAME
!! lobpcgccIIwf
!!
!! FUNCTION
!! this routine updates the whole wave functions at a given k-point,
!! using the lobpcg method
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvaluesg at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (GZ,AR,MT)
!! this file is distributed under the terms of the
!! gnu general public license, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/infos/contributors .
!!
!! INPUTS
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variales for this dataset
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (hartree)
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mcg=second dimension of the cg array
!!  mgfft=maximum size of 1d ffts
!!  mgsc=second dimension of the gsc array
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  nbdblock : number of blocks
!!  npw_k=number of plane waves at this k point
!!  ntypat=number of types of atoms in unit cell.
!!  nvloc=final dimension of vlocal (usually 1, but 4 for non-collinear)
!!  n4,n5,n6 used for dimensionning of vlocal
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  use_subovl= 1 if "subovl" array is computed (see below)
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!  vxctaulocal(n4,n5,n6,nvloc,4)= local potential corresponding to the derivative of XC energy with respect to
!!    kinetic energy density, in real space, on the augmented fft grid. (optional argument)
!!    This array contains also the gradient of vxctaulocal (gvxctaulocal) in vxctaulocal(:,:,:,:,2:4).
!!
!! OUTPUT
!!  resid_k(nband_k)=residuals for each states
!!  subham(nband_k*(nband_k+1))=the matrix elements of h
!!  If gs_hamk%usepaw==0:
!!    gsc(2,mgsc)=<g|s|c> matrix elements (s=overlap)
!!    subvnl(nband_k*(nband_k+1)*(1-gs_hamk%usepaw))=the matrix elements of vnl
!!  If use_subovl==0:
!!    subovl(nband_k*(nband_k+1)*use_subovl)=the matrix elements of s
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=updated wavefunctions
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!      getghc,leave_new,lobpcgcciiiwf,nonlop,timab,wrtout,xcomm_init,xsum_mpi
!!      zgemm,zheev,zhegv,zorthonormalize,ztrsm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine lobpcgccIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
     &           kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&
     &           nband_k,nbdblock,npw_k,ntypat,nvloc,n4,n5,n6,ph3d,prtvol,&
     &           resid_k,subham,subovl,subvnl,use_subovl,vlocal,&
     &           vxctaulocal) ! optional argument

 use m_profiling

  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use m_xmpi

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lobpcgccIIwf'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_abiutil
 use interfaces_65_nonlocal
 use interfaces_66_wfs
 use interfaces_79_seqpar_mpi, except_this_one => lobpcgccIIwf
!End of the abilint section

  implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif
  !Arguments ------------------------------------
  integer :: dimffnl,icg,igsc,lmnmax,matblk,mcg,mgsc,mgfft,mpsang,mpssoang,n4,n5,n6
  integer :: natom,nband_k,nbdblock,npw_k,ntypat,nvloc,prtvol,use_subovl
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  integer :: kg_k(3,npw_k)
  real(dp) :: cg(2,mcg),ffnl(npw_k,dimffnl,lmnmax,ntypat),gsc(2,mgsc)
  real(dp) :: kinpw(npw_k),ph3d(2,npw_k,matblk),resid_k(nband_k)
  real(dp) :: vlocal(n4,n5,n6,nvloc)
  real(dp) :: subham(nband_k*(nband_k+1)),subvnl(nband_k*(nband_k+1)*(1-gs_hamk%usepaw))
  real(dp) :: subovl(nband_k*(nband_k+1)*use_subovl)
  real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)

  !Local variables-------------------------------
  integer :: spacecomm=0
  integer :: bblocksize,blocksize,cgindex,choice,cpopt
  integer :: gscindex,iblocksize
  integer :: iblock,iband,idir,ierr,ii,info
  integer :: ipw,istwf_k,isubo,isubh,iterationnumber,iwavef,jblocksize,littleblocksize,lwork
  integer :: maxiterations,my_nspinor,nkpg,nnlout,old_paral_level,optekin,optpcon,paw_opt
  integer :: signs,sij_opt,tim_getghc,tim_nonlop,vectsize
  logical :: gen_eigenpb
  real(dp) :: cgreipw,cgimipw,cscre,cscim,chcre,chcim,cvcre,cvcim,dum,lambda_i,sq2
  character(len=500) :: message
  logical, allocatable :: pflag(:)
  real(dp) :: tsec(2)
  real(dp), allocatable :: gwavef(:,:),cwavef(:,:),gvnlc(:,:)
  real(dp), allocatable :: residualnorms(:),eigen(:),rwork(:),lambda(:,:),kpg_dum(:,:)
  real(dp),allocatable :: pcon(:,:)
  complex(dp), allocatable :: blockvectorx(:,:),blockvectorax(:,:),blockvectorbx(:,:),&
       & blockvectorr(:,:),blockvectorar(:,:),blockvectorbr(:,:),&
       & blockvectorp(:,:),blockvectorap(:,:),blockvectorbp(:,:),blockvectordumm(:,:),&
       & blockvectory(:,:),blockvectorby(:,:),&
       & gramxax(:,:),gramxar(:,:),gramxap(:,:),gramrar(:,:),gramrap(:,:),&
       & grampap(:,:),&
       & gramxbx(:,:),gramxbr(:,:),gramxbp(:,:),gramrbr(:,:),gramrbp(:,:),&
       & grampbp(:,:),&
       & coordx(:,:),diagcoordx(:,:),&!lambda(:,:),&
       & gramyx(:,:),&
       & transf(:,:,:),work(:),&
       & blockvectorxc(:,:)
 real(dp), allocatable :: dummy1(:),dummy2(:,:),dummy3(:,:,:)
!This is for the call to lobpiii, transfer of information from blockvector to vector
! following variables appear below but are commented out
! complex(dp),allocatable :: vectorap(:),vectorax(:),vectorbp(:),vectorbx(:),vectorby(:)
! complex(dp),allocatable :: vectorp(:),vectorr(:),vectorx(:),vectory(:)
 type(cprj_type) :: cprj_dum(1,1)

  !no_abirules
  !correspondence with abinit. here for real wf but in complex mode
  !this is the index of a given band
  cgindex(iblocksize)=npw_k*my_nspinor*(iblocksize-1)+icg+1
  gscindex(iblocksize)=npw_k*my_nspinor*(iblocksize-1)+igsc+1

! *************************************************************************

 gen_eigenpb=(gs_hamk%usepaw==1)
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
 optekin=0;if (dtset%wfoptalg>10) optekin=0
 optpcon=1;if (dtset%wfoptalg>10) optpcon=0
 resid_k=zero

 if(mod(nband_k,nbdblock)/=0) then
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' vtowfk : ERROR -',ch10,&
&   '  For the moment, nband must be a multiple of nbdblock with wfoptalg=5 !',ch10,&
&   '  Action : raise nband or change nbdblock '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 sq2=sqrt(two)
 vectsize=npw_k*my_nspinor
 istwf_k=gs_hamk%istwf_k
 maxiterations=dtset%nline

!Big loop bands inside blocks
 do iblock=1,nbdblock

   blocksize=(nband_k-1)/nbdblock+1
   bblocksize=(iblock-1)*blocksize

!  allocations
   ABI_ALLOCATE(pcon,(npw_k,blocksize))
   ABI_ALLOCATE(blockvectorx,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectorax,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectorbx,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectorr,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectorar,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectorbr,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectorp,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectorap,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectorbp,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectordumm,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectorxc,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectory,(vectsize,bblocksize))
   ABI_ALLOCATE(blockvectorby,(vectsize,bblocksize))
   ABI_ALLOCATE(gramyx,(bblocksize,blocksize))
   ABI_ALLOCATE(gramxax,(blocksize,blocksize))
   ABI_ALLOCATE(gramxar,(blocksize,blocksize))
   ABI_ALLOCATE(gramxap,(blocksize,blocksize))
   ABI_ALLOCATE(gramrar,(blocksize,blocksize))
   ABI_ALLOCATE(gramrap,(blocksize,blocksize))
   ABI_ALLOCATE(grampap,(blocksize,blocksize))
   ABI_ALLOCATE(gramxbx,(blocksize,blocksize))
   ABI_ALLOCATE(gramxbr,(blocksize,blocksize))
   ABI_ALLOCATE(gramxbp,(blocksize,blocksize))
   ABI_ALLOCATE(gramrbr,(blocksize,blocksize))
   ABI_ALLOCATE(gramrbp,(blocksize,blocksize))
   ABI_ALLOCATE(grampbp,(blocksize,blocksize))
   ABI_ALLOCATE(lambda,(blocksize,blocksize))
   ABI_ALLOCATE(transf,(blocksize,blocksize,3))
   ABI_ALLOCATE(residualnorms,(blocksize))
   ABI_ALLOCATE(pflag,(blocksize))

   pflag = .false.

!  transfer array of wf coeff in block to blockvectorx (complex to complex)
   do iblocksize=1,blocksize
     iband=iblocksize+bblocksize
     blockvectorx(1:vectsize,iblocksize)=dcmplx(cg(1,cgindex(iband):cgindex(iband+1)-1),&
&     cg(2,cgindex(iband):cgindex(iband+1)-1))
   end do

!  transfer array of wf coeff less than iblock to blockvectory (not done)
   if(iblock /=1) then
!    transfer cg to blockvectory, for the previous band index
     do iblocksize=1,bblocksize
       iband=iblocksize
       blockvectory(1:vectsize,iblocksize)=dcmplx(cg(1,cgindex(iband):cgindex(iband+1)-1),&
&       cg(2,cgindex(iband):cgindex(iband+1)-1))
     end do
!    call operators(blockvectory,blockvectorby)
     if(gen_eigenpb) then
       ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
       ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor))
       do iblocksize=1,bblocksize
         cwavef(1,1:npw_k*my_nspinor)=real (blockvectory(1:npw_k*my_nspinor,iblocksize))
         cwavef(2,1:npw_k*my_nspinor)=aimag(blockvectory(1:npw_k*my_nspinor,iblocksize))
!        Call to nonlop: compute <g|S|c>
         choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
         call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
&         dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
&         istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
&         mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
&         gs_hamk%nloalg,nnlout,npw_k,npw_k,my_nspinor,dtset%nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
&         gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,signs,gs_hamk%sij,&
&         gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef,&
&         use_gpu_cuda=dtset%use_gpu_cuda)
         blockvectorby(1:npw_k*my_nspinor,iblocksize)=dcmplx(gwavef(1,1:npw_k*my_nspinor),gwavef(2,1:npw_k*my_nspinor))
       end do
       ABI_DEALLOCATE(cwavef)
       ABI_DEALLOCATE(gwavef)
     else
       blockvectorby(:,:)=blockvectory(:,:)
     end if

!    orthogonalize x to the constraint y(supposed orthonormal)
!    blockvectorx=blockvectorx-&
!    &matmul(blockvectory,matmul(transpose(blockvectorby),blockvectorx))
     call zgemm('c','n',bblocksize,blocksize,vectsize,cone,blockvectorby,&
&     vectsize,blockvectorx,vectsize,czero,gramyx,bblocksize)
     old_paral_level= mpi_enreg%paral_level
     mpi_enreg%paral_level=3
     call xcomm_init(mpi_enreg,spaceComm)
     call timab(48,1,tsec)
     call xsum_mpi(gramyx,spaceComm,ierr)
     call timab(48,2,tsec)
     mpi_enreg%paral_level= old_paral_level
     call zgemm('n','n',vectsize,blocksize,bblocksize,cone,blockvectory,&
&     vectsize,gramyx,bblocksize,czero,blockvectordumm,vectsize)
     blockvectorx=blockvectorx-blockvectordumm
   end if
!  compute right hand side
!  call operators(blockvectorx,blockvectorbx)
   if(gen_eigenpb) then
     ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
     ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor))
     do iblocksize=1,blocksize
       cwavef(1,1:npw_k*my_nspinor)=real (blockvectorx(1:npw_k*my_nspinor,iblocksize))
       cwavef(2,1:npw_k*my_nspinor)=aimag(blockvectorx(1:npw_k*my_nspinor,iblocksize))
!      Call to nonlop: compute <g|S|c>
       choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
       call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
&       dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
&       istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
&       mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
&       gs_hamk%nloalg,nnlout,npw_k,npw_k,my_nspinor,dtset%nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
&       gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,signs,gs_hamk%sij,&
&       gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef,use_gpu_cuda=dtset%use_gpu_cuda)
       blockvectorbx(1:npw_k*my_nspinor,iblocksize)=dcmplx(gwavef(1,1:npw_k*my_nspinor),gwavef(2,1:npw_k*my_nspinor))
     end do
     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(gwavef)
   else
     blockvectorbx(:,:)=blockvectorx(:,:)
   end if

!  orthogonalize x
!  call zorthonormalize(blockvectorx,blockvectorbx)
   call zorthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,gramxbx,vectsize)
   call ztrsm('r','u','n','n',vectsize,blocksize,cone,gramxbx,blocksize,&
&   blockvectorbx,vectsize)
!  call operatorh(blockvectorx,blockvectorax)
   ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(gvnlc,(2,npw_k*my_nspinor*blocksize))
   do iblocksize=1,blocksize
     cwavef(1,npw_k*my_nspinor*(iblocksize-1)+1:npw_k*my_nspinor*iblocksize)=real(blockvectorx(1:npw_k*my_nspinor,iblocksize))
     cwavef(2,npw_k*my_nspinor*(iblocksize-1)+1:npw_k*my_nspinor*iblocksize)=aimag(blockvectorx(1:npw_k*my_nspinor,iblocksize))
   end do
   tim_getghc=7 ; sij_opt=0
   if(present(vxctaulocal))then
     call getghc(-1,cwavef,cprj_dum,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
&     kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,my_nspinor,ntypat,&
&     nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
   else
     call getghc(-1,cwavef,cprj_dum,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
&     kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,my_nspinor,ntypat,&
&     nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
   end if
   do iblocksize=1,blocksize
     blockvectorax(1:npw_k*my_nspinor,iblocksize)=&
&     dcmplx(gwavef(1,npw_k*my_nspinor*(iblocksize-1)+1:npw_k*my_nspinor*iblocksize),&
&     gwavef(2,npw_k*my_nspinor*(iblocksize-1)+1:npw_k*my_nspinor*iblocksize))
   end do
   ABI_DEALLOCATE(cwavef)
   ABI_DEALLOCATE(gwavef)
   ABI_DEALLOCATE(gvnlc)

!  do rayleigh ritz on a in space x
!  gramxax=matmul(transpose(blockvectorx),blockvectorax)
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
&   vectsize,blockvectorax,vectsize,czero,gramxax,blocksize)
   old_paral_level= mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   call timab(48,1,tsec)
   call xsum_mpi(gramxax,spaceComm,ierr)
   call timab(48,2,tsec)
   mpi_enreg%paral_level= old_paral_level
   ABI_ALLOCATE(eigen,(blocksize))

!  call la_syev(gramxax,eigen,jobz='v')
   lwork=3*blocksize-2
   ABI_ALLOCATE(work,(lwork))
   ABI_ALLOCATE(rwork,(lwork))
!  write(std_out,*)'gramxax bef',gramxax(:,:)
   do iblocksize=1,blocksize
     do jblocksize=1,blocksize
       if(abs(gramxax(iblocksize,jblocksize)) < 1.e-14) then
         gramxax(iblocksize,jblocksize)=czero
       else
!        write(std_out,*)'gramxax non nul',gramxax(iblocksize,jblocksize)
       end if
     end do
   end do
   call zheev('v','u',blocksize,gramxax,blocksize,eigen,work,lwork,rwork,info)
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(rwork)
!  blockvectorx=matmul(blockvectorx,gramxax)
   call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
&   vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
   blockvectorx=blockvectordumm
!  blockvectorax=matmul(blockvectorax,gramxax)
   call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
&   vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
   blockvectorax=blockvectordumm
!  blockvectorbx=matmul(blockvectorbx,gramxax)
   call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
&   vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
   blockvectorbx=blockvectordumm
   do iblocksize=1,blocksize
     lambda(iblocksize,iblocksize)=eigen(iblocksize)
   end do
   write(std_out,*)'lambda',eigen

!  now the main alogrithm
!  !!     allocate(vectorx(vectsize),vectorbx(vectsize),vectorax(vectsize),vectory(vectsize),&
!  !!          &  vectorby(vectsize),vectorp(vectsize),vectorbp(vectsize),vectorap(vectsize),&
!  !!          &  vectorr(vectsize))
   iter: do iterationnumber=1,maxiterations
!    write(std_out,*)'bvx',blockvectorx(10,7)
     write(std_out,*)'iterationnumber',iterationnumber
     do iblocksize=1,blocksize
!      vectorx(:)=blockvectorx(:,iblocksize)
!      vectorbx(:)=blockvectorbx(:,iblocksize)
!      vectorax(:)=blockvectorax(:,iblocksize)
!      vectory(:)=blockvectory(:,iblocksize)! attention, purement blanc si pas de y
!      vectorby(:)=blockvectorby(:,iblocksize)! a remplacer par by etc....
!      vectorp(:)=blockvectorp(:,iblocksize)
!      vectorbp(:)=blockvectorbp(:,iblocksize)
!      vectorap(:)=blockvectorap(:,iblocksize)
       lambda_i=lambda(iblocksize,iblocksize)
!      if(iblock > 1) stop('huh')

       littleblocksize=1
       call lobpcgccIIIwf(dimffnl,dtfil,dtset,ffnl,gs_hamk,iterationnumber,&
&       kg_k,kinpw,lmnmax,matblk,mgfft,mpi_enreg,mpsang,&
&       mpssoang,natom,npw_k,ntypat,&
&       nvloc,n4,n5,n6,pcon,ph3d,prtvol,vlocal,&
&       littleblocksize,bblocksize,vectsize,pflag(iblocksize), &
&       blockvectorx (:,iblocksize:iblocksize),&
&       blockvectorbx(:,iblocksize:iblocksize),&
&       blockvectorax(:,iblocksize:iblocksize),&
&       blockvectory,blockvectorby,lambda_i,&
&       blockvectorp (:,iblocksize:iblocksize),&
&       blockvectorbp(:,iblocksize:iblocksize),&
&       blockvectorap(:,iblocksize:iblocksize)&
&       )

!      blockvectorx(:,iblocksize)=vectorx(:)
!      blockvectorbx(:,iblocksize)=vectorbx(:)
!      blockvectorax(:,iblocksize)=vectorax(:)
!      blockvectorp(:,iblocksize)=vectorp(:)
!      blockvectorbp(:,iblocksize)=vectorbp(:)
!      blockvectorap(:,iblocksize)=vectorap(:)
     end do
!    gramxax
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
&     vectsize,blockvectorax,vectsize,czero,gramxax,blocksize)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
&     vectsize,blockvectorx,vectsize,czero,gramxbx,blocksize)
!    write(std_out,*)'in iii,xax'!,gramxax(7,7)
!    write(std_out,*)'in iii xbx'!,gramxbx(7,7)
     lwork=3*blocksize-2
     ABI_ALLOCATE(work,(lwork))
     ABI_ALLOCATE(rwork,(lwork))
     call zhegv(1,'v','u',blocksize,gramxax,blocksize,gramxbx,blocksize,eigen,&
&     work,lwork,rwork,info)
     ABI_DEALLOCATE(work)
     ABI_DEALLOCATE(rwork)
     write(std_out,*)'lambda',eigen
     write(std_out,*)' '
     lambda(:,:)=zero
     do iblocksize=1,blocksize
       lambda(iblocksize,iblocksize)=eigen(iblocksize)
     end do
     ABI_ALLOCATE(coordx,(blocksize,blocksize))
     ABI_ALLOCATE(diagcoordx,(blocksize,blocksize))
     coordx=gramxax
!    rotate all the vectors according to coordx

!    choix de p
     diagcoordx=czero
     do iblocksize=1,blocksize
       diagcoordx(iblocksize,iblocksize) = coordx(iblocksize,iblocksize)
       coordx(iblocksize,iblocksize) = czero
     end do
!    blockvectorxc = matmul(blockvectorx,coordx)
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
&     vectsize,coordx,blocksize,czero,blockvectorxc,vectsize)
!    blockvectorx = matmul(blockvectorx,diagcoordx) + blockvectorxc
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
&     vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
     blockvectorx = blockvectordumm + blockvectorxc
!    blockvectorp = matmul(blockvectorp,diagcoordx) + blockvectorxc
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorp,&
&     vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
     blockvectorp = blockvectordumm + blockvectorxc
!    blockvectorxc = matmul(blockvectorbx,coordx)
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
&     vectsize,coordx,blocksize,czero,blockvectorxc,vectsize)
!    blockvectorbx = matmul(blockvectorbx,diagcoordx) + blockvectorxc
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
&     vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
     blockvectorbx = blockvectordumm + blockvectorxc
!    blockvectorbp = matmul(blockvectorbp,diagcoordx) + blockvectorxc
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbp,&
&     vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
     blockvectorbp = blockvectordumm + blockvectorxc
!    blockvectorxc = matmul(blockvectorax,coordx)
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
&     vectsize,coordx,blocksize,czero,blockvectorxc,vectsize)
!    blockvectorax = matmul(blockvectorax,diagcoordx) + blockvectorxc
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
&     vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
     blockvectorax = blockvectordumm + blockvectorxc
!    blockvectorap = matmul(blockvectorap,diagcoordx) + blockvectorxc
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorap,&
&     vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
     blockvectorap = blockvectordumm + blockvectorxc

!    !!    !autre choix possible
!    !!    !blockvectorx = matmul(blockvectorx,coordx)
!    !!    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
!    !!         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!    !!    blockvectorx = blockvectordumm
!    !!    !blockvectorbx = matmul(blockvectorbx,coordx)
!    !!    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
!    !!         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!    !!    blockvectorbx = blockvectordumm
!    !!    !blockvectorax = matmul(blockvectorax,coordx)
!    !!    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
!    !!         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!    !!    blockvectorax = blockvectordumm
!    !!    !blockvectorp = matmul(blockvectorp,coordx)
!    !!    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorp,&
!    !!         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!    !!    blockvectorp = blockvectordumm
!    !!    !blockvectorbp = matmul(blockvectorbp,coordx)
!    !!    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbp,&
!    !!         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!    !!    blockvectorbp = blockvectordumm
!    !!    !blockvectorap = matmul(blockvectorap,coordx)
!    !!    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorap,&
!    !!         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!    !!    blockvectorap = blockvectordumm
!    gramxax
!    call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
!    &               vectsize,blockvectorax,vectsize,czero,gramxax,blocksize)
!    call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
!    &               vectsize,blockvectorx,vectsize,czero,gramxbx,blocksize)
!    write(std_out,*)'in iii,xax after',gramxax
!    write(std_out,*)'in iii xbx after',gramxbx
     ABI_DEALLOCATE(coordx)
     ABI_DEALLOCATE(diagcoordx)
   end do iter
   ABI_DEALLOCATE(eigen)
!  epilogue
!  gramxbx=matmul(transpose(blockvectorx),blockvectorx)
   if(.true.) then !epilogue
!    call operators(blockvectorx,blockvectorbx)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
&     vectsize,blockvectorbx,vectsize,czero,gramxbx,blocksize)
!    blockvectorax=matmul(operatora,blockvectorx)
!    call operatorh(blockvectorx,blockvectorax)
!    gramxax=matmul(transpose(blockvectorax),blockvectorx)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
&     vectsize,blockvectorax,vectsize,czero,gramxax,blocksize)
     ABI_ALLOCATE(eigen,(blocksize))
!    call la_sygv(gramxax,gramxbx,eigen,itype=1,jobz='v')
     lwork=3*blocksize-2
     ABI_ALLOCATE(work,(lwork))
     ABI_ALLOCATE(rwork,(lwork))
     call zhegv(1,'v','u',blocksize,gramxax,blocksize,gramxbx,blocksize,eigen,&
&     work,lwork,rwork,info)
     ABI_DEALLOCATE(work)
     ABI_DEALLOCATE(rwork)
     lambda=czero
     do iblocksize=1,blocksize
       lambda(iblocksize,iblocksize)=eigen(iblocksize)
     end do
!    write(std_out,*)'gramxax'
!    write(std_out,*)gramxax
     write(std_out,*)'eigen at the end',eigen
!    debug
!    blockvectorr=blockvectorax-matmul(blockvectorx,lambda)
!    residualnorms=sqrt(sum(blockvectorr**2,dim=1))
!    write(std_out,*)'residualnorm at the end bef orth',residualnorms
!    debug
!    blockvectorx=matmul(blockvectorx,gramxax)
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
&     vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
     blockvectorx=blockvectordumm
!    blockvectorax=matmul(blockvectorax,gramxax)
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
&     vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
     blockvectorax=blockvectordumm
!    blockvectorbx=matmul(blockvectorbx,gramxax)
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
&     vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
     blockvectorbx=blockvectordumm
!    blockvectorr=blockvectorax-matmul(blockvectorbx,lambda)
!    call dgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
!    &               vectsize,lambda,blocksize,czero,blockvectordumm,vectsize)
!    blockvectorr=blockvectorax-blockvectordumm
     do iblocksize=1,blocksize
       blockvectorr(:,iblocksize)=blockvectorax(:,iblocksize)-eigen(iblocksize)*blockvectorbx(:,iblocksize)
     end do
     ABI_DEALLOCATE(eigen)
   end if !epilogue



   residualnorms=sum(abs(blockvectorr)**2,dim=1)
   old_paral_level= mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   call timab(48,1,tsec)
   call xsum_mpi(residualnorms,spaceComm,ierr)
   call timab(48,2,tsec)
   mpi_enreg%paral_level= old_paral_level
   residualnorms=sqrt(residualnorms)
   do iblocksize=1,blocksize
     iband=iblocksize+(iblock-1)*blocksize
     cg(1,cgindex(iband):cgindex(iband+1)-1)=real(blockvectorx(1:vectsize,iblocksize))
     cg(2,cgindex(iband):cgindex(iband+1)-1)=aimag(blockvectorx(1:vectsize,iblocksize))
   end do
   if(gen_eigenpb) then
     do iblocksize=1,blocksize
       iband=iblocksize+(iblock-1)*blocksize
       gsc(1,gscindex(iband):gscindex(iband+1)-1)=real(blockvectorbx(1:vectsize,iblocksize))
       gsc(2,gscindex(iband):gscindex(iband+1)-1)=aimag(blockvectorbx(1:vectsize,iblocksize))
     end do
   end if
!  this should not exist,since this induce one too much getghc.lazy programming....
!  call operatorh(blockvectorx,blockvectorax,subham,subvnl)!fill also subham, subvnl
   ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(gvnlc,(2,npw_k*my_nspinor*blocksize))
   isubh=1+2*(iblock-1)*blocksize*((iblock-1)*blocksize+1)/2
   do iblocksize=1,blocksize
     cwavef(1,npw_k*my_nspinor*(iblocksize-1)+1:npw_k*my_nspinor*iblocksize)=real(blockvectorx(1:npw_k*my_nspinor,iblocksize))
     cwavef(2,npw_k*my_nspinor*(iblocksize-1)+1:npw_k*my_nspinor*iblocksize)=aimag(blockvectorx(1:npw_k*my_nspinor,iblocksize))
   end do
   tim_getghc=7 ; sij_opt=0
   if(present(vxctaulocal))then
     call getghc(-1,cwavef,cprj_dum,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
&     kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,my_nspinor,ntypat,&
&     nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
   else
     call getghc(-1,cwavef,cprj_dum,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
&     kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,my_nspinor,ntypat,&
&     nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
   end if
   do iblocksize=1,blocksize
     blockvectorax(1:npw_k*my_nspinor,iblocksize)=&
&     dcmplx(gwavef(1,npw_k*my_nspinor*(iblocksize-1)+1:npw_k*my_nspinor*iblocksize),&
&     gwavef(2,npw_k*my_nspinor*(iblocksize-1)+1:npw_k*my_nspinor*iblocksize))
     do ii=1,(iblock-1)*blocksize+iblocksize
       iwavef=(ii-1)*npw_k*my_nspinor+icg
       chcre=zero ; chcim=zero
       if (gs_hamk%usepaw==1) then
         do ipw=1,npw_k*my_nspinor
           cgreipw=cg(1,ipw+iwavef);cgimipw=cg(2,ipw+iwavef)
           chcre=chcre+cgreipw*gwavef(1,ipw+(iblocksize-1)*npw_k*my_nspinor)+cgimipw*gwavef(2,ipw+(iblocksize-1)*npw_k*my_nspinor)
           chcim=chcim+cgreipw*gwavef(2,ipw+(iblocksize-1)*npw_k*my_nspinor)-cgimipw*gwavef(1,ipw+(iblocksize-1)*npw_k*my_nspinor)
         end do
       else
         cvcre=zero ; cvcim=zero
         do ipw=1,npw_k*my_nspinor
           cgreipw=cg(1,ipw+iwavef);cgimipw=cg(2,ipw+iwavef)
           chcre=chcre+cgreipw*gwavef(1,ipw+(iblocksize-1)*npw_k*my_nspinor)+cgimipw*gwavef(2,ipw+(iblocksize-1)*npw_k*my_nspinor)
           chcim=chcim+cgreipw*gwavef(2,ipw+(iblocksize-1)*npw_k*my_nspinor)-cgimipw*gwavef(1,ipw+(iblocksize-1)*npw_k*my_nspinor)
           cvcre=cvcre+cgreipw*gvnlc(1,ipw+(iblocksize-1)*npw_k*my_nspinor)+cgimipw*gvnlc(2,ipw+(iblocksize-1)*npw_k*my_nspinor)
           cvcim=cvcim+cgreipw*gvnlc(2,ipw+(iblocksize-1)*npw_k*my_nspinor)-cgimipw*gvnlc(1,ipw+(iblocksize-1)*npw_k*my_nspinor)
         end do
!        Store real and imag parts in hermitian storage mode:
         subvnl(isubh)=cvcre ; subvnl(isubh+1)=cvcim
       end if
!      Store real and imag parts in hermitian storage mode:
       subham(isubh)=chcre ; subham(isubh+1)=chcim
       isubh=isubh+2
     end do
   end do
!  comm for subham and subvnl is made in vtowfk

   ABI_DEALLOCATE(cwavef)
   ABI_DEALLOCATE(gwavef)
   ABI_DEALLOCATE(gvnlc)
!  call operators(blockvectorx,blockvectorbx,subovl)!fill also  subovl
   if((gen_eigenpb).and.(use_subovl==1)) then
     ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
     ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor))
     isubo=1+2*(iblock-1)*blocksize*((iblock-1)*blocksize+1)/2
     do iblocksize=1,blocksize
       cwavef(1,1:npw_k*my_nspinor)=real (blockvectorx(1:npw_k*my_nspinor,iblocksize))
       cwavef(2,1:npw_k*my_nspinor)=aimag(blockvectorx(1:npw_k*my_nspinor,iblocksize))
!      Call to nonlop: compute <g|S|c>
       choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
       call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
&       dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
&       istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
&       mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
&       gs_hamk%nloalg,nnlout,npw_k,npw_k,my_nspinor,dtset%nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
&       gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,signs,gs_hamk%sij,&
&       gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef,&
&       use_gpu_cuda=dtset%use_gpu_cuda)
       blockvectorbx(1:npw_k*my_nspinor,iblocksize)=dcmplx(gwavef(1,1:npw_k*my_nspinor),gwavef(2,1:npw_k*my_nspinor))
       do ii=1,(iblock-1)*blocksize+iblocksize
         iwavef=(ii-1)*npw_k*my_nspinor+icg
         cscre=zero;cscim=zero
         do ipw=1,npw_k*my_nspinor
           cgreipw=cg(1,ipw+iwavef);cgimipw=cg(2,ipw+iwavef)
           cscre=cscre+cgreipw*gwavef(1,ipw)+cgimipw*gwavef(2,ipw)
           cscim=cscim+cgreipw*gwavef(2,ipw)-cgimipw*gwavef(1,ipw)
         end do
!        Store real and imag parts in hermitian storage mode:
         subovl(isubo)=cscre ; subovl(isubo+1)=cscim
         isubo=isubo+2
       end do
     end do
     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(gwavef)
   end if
   write(std_out,*) "mytest"
!  stop

   write(std_out,*)'residualnorm at the end end',residualnorms
   ABI_DEALLOCATE(pcon)
   ABI_DEALLOCATE(blockvectorx)
   ABI_DEALLOCATE(blockvectorax)
   ABI_DEALLOCATE(blockvectorbx)
   ABI_DEALLOCATE(blockvectorr)
   ABI_DEALLOCATE(blockvectorar)
   ABI_DEALLOCATE(blockvectorbr)
   ABI_DEALLOCATE(blockvectorp)
   ABI_DEALLOCATE(blockvectorap)
   ABI_DEALLOCATE(blockvectorbp)
   ABI_DEALLOCATE(blockvectory)
   ABI_DEALLOCATE(blockvectorby)
   ABI_DEALLOCATE(gramyx)
   ABI_DEALLOCATE(transf)
   ABI_DEALLOCATE(blockvectordumm)
   ABI_DEALLOCATE(blockvectorxc)
   ABI_DEALLOCATE(gramxax)
   ABI_DEALLOCATE(gramxar)
   ABI_DEALLOCATE(gramxap)
   ABI_DEALLOCATE(gramrar)
   ABI_DEALLOCATE(gramrap)
   ABI_DEALLOCATE(grampap)
   ABI_DEALLOCATE(gramxbx)
   ABI_DEALLOCATE(gramxbr)
   ABI_DEALLOCATE(gramxbp)
   ABI_DEALLOCATE(gramrbr)
   ABI_DEALLOCATE(gramrbp)
   ABI_DEALLOCATE(grampbp)
   ABI_DEALLOCATE(lambda)
   ABI_DEALLOCATE(residualnorms)
   ABI_DEALLOCATE(pflag)
!  End big loop over bands inside blocks
 end do

end subroutine lobpcgccIIwf
!!***
