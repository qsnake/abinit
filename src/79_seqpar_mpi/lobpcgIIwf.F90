!{\src2tex{textfont=tt}}
!!****f* abinit/lobpcgIIwf
!! NAME
!! lobpcgIIwf
!!
!! FUNCTION
!! this routine updates the whole wave functions at a given k-point,
!! using the lobpcg method
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
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
!!      dgemm,dsyev,dsygv,dtrsm,getghc,leave_new,lobpcgiiiwf,nonlop
!!      orthonormalize,precon2,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine lobpcgIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
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
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lobpcgIIwf'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_53_abiutil
 use interfaces_65_nonlocal
 use interfaces_66_wfs
 use interfaces_79_seqpar_mpi, except_this_one => lobpcgIIwf
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
! integer, save :: frozen_count = 0 ! commented out below
integer :: iblock,ii,ipw,ipw1,isubo,isubh,iterationnumber,iwavef,maxiterations
integer :: info,lwork,tim_getghc,choice,idir,tim_nonlop
integer :: cpopt,paw_opt, nkpg,nnlout,optekin=0,optpcon=1,signs,sij_opt
integer :: iblocksize
integer :: rvectsize,vectsize,blocksize,bblocksize,iband,istwf_k
integer :: cgindex,gscindex,littleblocksize,my_nspinor
logical :: gen_eigenpb
real(dp) :: cgreipw,cgimipw,cscre,chcre,cvcre,dum,sq2
real(dp) :: lambda_i(1,1)
character(len=500) :: message
logical, allocatable :: pflag(:)
real(dp),allocatable :: pcon(:,:)
real(dp), allocatable :: blockvectorx(:,:),blockvectorax(:,:),blockvectorbx(:,:),&
& blockvectorr(:,:),blockvectorar(:,:),blockvectorbr(:,:),&
& blockvectorz(:,:),blockvectoraz(:,:),blockvectorbz(:,:),&
& blockvectordumm(:,:),&
& blockvectory(:,:),blockvectorby(:,:),&
& gramxax(:,:),gramxar(:,:),gramxap(:,:),gramrar(:,:),gramrap(:,:),&
& grampap(:,:),&
& gramxbx(:,:),gramxbr(:,:),gramxbp(:,:),gramrbr(:,:),gramrbp(:,:),&
& grampbp(:,:),&
& coordx(:,:),diagcoordx(:,:),blockvectorxc(:,:),lambda(:,:),&
& gramyx(:,:),&
& kpg_dum(:,:),work(:),dummy1(:),dummy2(:,:),dummy3(:,:,:)
real(dp), allocatable ::blockvectorp(:,:),blockvectorap(:,:),blockvectorbp(:,:)
real(dp), allocatable :: gwavef(:,:),cwavef(:,:),gvnlc(:,:)
real(dp), allocatable :: residualnorms(:),eigen(:)
!this is for the call to lobpiii, transfer of information from blockvector to vector
real(dp), allocatable ::vectorx(:,:),vectorbx(:,:),vectorax(:,:),vectory(:,:),vectorby(:,:),&
&                             vectorp(:,:),vectorbp(:,:),vectorap(:,:),vectorr(:,:)
 type(cprj_type) :: cprj_dum(1,1)

!NO_ABIRULES
!correspondence with abinit. here for real wf but in complex mode
!this is the index of a given band
cgindex(iblocksize)=npw_k*my_nspinor*(iblocksize-1)+icg+1
gscindex(iblocksize)=npw_k*my_nspinor*(iblocksize-1)+igsc+1

! *********************************************************************

!debug
!write(std_out,*) npw_k*my_nspinor
!write(std_out,*) cgindex
!write(std_out,*) size(cg,1),size(cg,2)
!enddebug

 resid_k=zero

 if(mod(nband_k,nbdblock)/=0) then
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' vtowfk : ERROR -',ch10,&
&   '  For the moment, nband must be a multiple of nbdblock with wfoptalg=5 !',ch10,&
&   '  Action : raise nband or change nbdblock '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 gen_eigenpb=(gs_hamk%usepaw==1)
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
 sq2=sqrt(2.0_dp)
 if (mpi_enreg%me_g0 == 1) then
   vectsize=2*npw_k*my_nspinor-1
 else
   vectsize=2*npw_k*my_nspinor
 end if
 rvectsize=npw_k*my_nspinor

 istwf_k=gs_hamk%istwf_k
 maxiterations=dtset%nline

!Big loop bands inside blocks
 do iblock=1,nbdblock

   blocksize=(nband_k-1)/nbdblock+1
   bblocksize=(iblock-1)*blocksize

   if (bblocksize > 0) then
     write(std_out,*) 'BUG - current version of LOBPCG II algorithm does not hold for bblocksize /= 0'
     call leave_new('PERS')
   end if

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
   ABI_ALLOCATE(residualnorms,(blocksize))
   ABI_ALLOCATE(pflag,(blocksize))
   pflag = .false.

!  transfer array of wf coeff in block to blockvectorx
   blockvectorr(:,:)=zero
   do iblocksize=1,blocksize
     iband=iblocksize+bblocksize
     if (mpi_enreg%me_g0 == 1) then
       blockvectorx(1,iblocksize)=cg(1,cgindex(iband))
       blockvectorx(2:rvectsize,iblocksize)=cg(1,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
       blockvectorx(rvectsize+1:vectsize,iblocksize)=cg(2,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
     else
       blockvectorx(1:rvectsize,iblocksize)=cg(1,cgindex(iband):cgindex(iband+1)-1)*sq2
       blockvectorx(rvectsize+1:vectsize,iblocksize)=cg(2,cgindex(iband):cgindex(iband+1)-1)*sq2
     end if
   end do
!  if(mpi_enreg%me_group==0)write(std_out,*) 'bvx entree',blockvectorx(10,7)

   if (iblock /=1) then
!    transfer array of wf coeff less than iblock to blockvectory (not done)
!    transfer cg to blockvectory, for the previous band index
     do iblocksize=1,bblocksize
       iband=iblocksize
       if (mpi_enreg%me_g0 == 1) then
         blockvectory(1,iblocksize)=cg(1,cgindex(iband))
         blockvectory(2:rvectsize,iblocksize)=cg(1,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
         blockvectory(rvectsize+1:vectsize,iblocksize)=cg(2,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
       else
         blockvectory(1:rvectsize,iblocksize)=cg(1,cgindex(iband):cgindex(iband+1)-1)*sq2
         blockvectory(rvectsize+1:vectsize,iblocksize)=cg(2,cgindex(iband):cgindex(iband+1)-1)*sq2
       end if
     end do
!    b-orthogonalize x to the constraint y(supposed b-orthonormal)
!    blockvectorx=blockvectorx-&
!    &matmul(blockvectory,matmul(transpose(blockvectory),blockvectorx))
!    call operators(blockvectory,blockvectorby)
     if(gen_eigenpb) then
       ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
       ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor))
       do iblocksize=1,bblocksize
         if (mpi_enreg%me_g0 == 1) then
           cwavef(1,2:npw_k*my_nspinor)=blockvectory(2:npw_k*my_nspinor,iblocksize)/sq2
           cwavef(2,2:npw_k*my_nspinor)=blockvectory(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)/sq2
           cwavef(1,1)=blockvectory(1,iblocksize)
           cwavef(2,1)=zero
         else
           cwavef(1,1:npw_k*my_nspinor)=blockvectory(1:npw_k*my_nspinor,iblocksize)/sq2
           cwavef(2,1:npw_k*my_nspinor)=blockvectory(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)/sq2      !a verifier
         end if
!        call to nonlop: compute <g|s|c>
         choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; paw_opt=3 ; cpopt=-1 ; nnlout=0 ; nkpg=0
         call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
&         dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
&         istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
&         mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
&         gs_hamk%nloalg,nnlout,npw_k,npw_k,my_nspinor,dtset%nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
&         gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,signs,gs_hamk%sij,&
&         gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef,use_gpu_cuda=dtset%use_gpu_cuda)
         if (mpi_enreg%me_g0 == 1) then
           blockvectorby(2:npw_k*my_nspinor,iblocksize)=gwavef(1,2:npw_k*my_nspinor)*sq2
           blockvectorby(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)=gwavef(2,2:npw_k*my_nspinor)*sq2
           blockvectorby(1,iblocksize)=gwavef(1,1)
         else
           blockvectorby(1:npw_k*my_nspinor,iblocksize)=gwavef(1,1:npw_k*my_nspinor)*sq2
           blockvectorby(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)=gwavef(2,1:npw_k*my_nspinor)*sq2
         end if
       end do
       ABI_DEALLOCATE(cwavef)
       ABI_DEALLOCATE(gwavef)
     else
       blockvectorby(:,:)=blockvectory(:,:)
     end if

!    orthogonalize x to the constraint y(supposed orthonormal)
!    blockvectorx=blockvectorx-&
!    &matmul(blockvectory,matmul(transpose(blockvectorby),blockvectorx))
     call dgemm('t','n',bblocksize,blocksize,vectsize,one,blockvectorby,&
&     vectsize,blockvectorx,vectsize,zero,gramyx,bblocksize)
     call dgemm('n','n',vectsize,blocksize,bblocksize,one,blockvectory,&
&     vectsize,gramyx,bblocksize,zero,blockvectordumm,vectsize)
     blockvectorx=blockvectorx-blockvectordumm
   end if

!  compute right hand side
!  call operators(blockvectorx,blockvectorbx)
   if(gen_eigenpb) then
     ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
     ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor))
     do iblocksize=1,blocksize
       if (mpi_enreg%me_g0 == 1) then
         cwavef(1,2:npw_k*my_nspinor)=blockvectorx(2:npw_k*my_nspinor,iblocksize)/sq2
         cwavef(2,2:npw_k*my_nspinor)=blockvectorx(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)/sq2
         cwavef(1,1)=blockvectorx(1,iblocksize)
         cwavef(2,1)=zero
       else
         cwavef(1,1:npw_k*my_nspinor)=blockvectorx(1:npw_k*my_nspinor,iblocksize)/sq2
         cwavef(2,1:npw_k*my_nspinor)=blockvectorx(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)/sq2      !a verifier
       end if
!      call to nonlop: compute <g|s|c>
       choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
       call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
&       dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
&       istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
&       mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
&       gs_hamk%nloalg,nnlout,npw_k,npw_k,my_nspinor,dtset%nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
&       gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,signs,gs_hamk%sij,&
&       gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef,use_gpu_cuda=dtset%use_gpu_cuda)
       if (mpi_enreg%me_g0 == 1) then
         blockvectorbx(2:npw_k*my_nspinor,iblocksize)=gwavef(1,2:npw_k*my_nspinor)*sq2
         blockvectorbx(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)=gwavef(2,2:npw_k*my_nspinor)*sq2
         blockvectorbx(1,iblocksize)=gwavef(1,1)
       else
         blockvectorbx(1:npw_k*my_nspinor,iblocksize)=gwavef(1,1:npw_k*my_nspinor)*sq2
         blockvectorbx(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)=gwavef(2,1:npw_k*my_nspinor)*sq2
       end if
     end do
     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(gwavef)
   else
     blockvectorbx(:,:)=blockvectorx(:,:)
   end if

!  orthogonalize x
!  call zorthonormalize(blockvectorx,blockvectorbx)
   call orthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,gramxbx,vectsize)
!  if(mpi_enreg%me_group==0)write(std_out,*) 'bvx ortho',blockvectorx(10,7)
   call dtrsm('r','u','n','n',vectsize,blocksize,one,gramxbx,blocksize,&
&   blockvectorbx,vectsize)
!  write(std_out,*) 'gramxbx 1'
!  do ii = 1,blocksize
!  write(std_out,*) gramxbx(ii,:)
!  enddo

!  call operatorh(blockvectorx,blockvectorax)
   ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(gvnlc,(2,npw_k*my_nspinor*blocksize))
   do iblocksize=1,blocksize
     iband=iblocksize
     if (mpi_enreg%me_g0 == 1) then
       cwavef(1,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(2:npw_k*my_nspinor,iblocksize)/sq2
       cwavef(2,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)/sq2
       cwavef(2,cgindex(iband))=zero
       cwavef(1,cgindex(iband))=blockvectorx(1,iblocksize)
     else
       cwavef(1,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(1:npw_k*my_nspinor,iblocksize)/sq2
       cwavef(2,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)/sq2
     end if
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
     iband=iblocksize
     if (mpi_enreg%me_g0 == 1) then
       blockvectorax(2:npw_k*my_nspinor,iblocksize)=gwavef(1,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
       blockvectorax(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)=gwavef(2,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
       blockvectorax(1,iblocksize)=gwavef(1,cgindex(iband))
     else
       blockvectorax(1:npw_k*my_nspinor,iblocksize)=gwavef(1,cgindex(iband):cgindex(iband+1)-1)*sq2
       blockvectorax(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)=gwavef(2,cgindex(iband):cgindex(iband+1)-1)*sq2
     end if
   end do
   ABI_DEALLOCATE(cwavef)
   ABI_DEALLOCATE(gwavef)
   ABI_DEALLOCATE(gvnlc)

!  do rayleigh ritz on a in space x
!  gramxax=matmul(transpose(blockvectorx),blockvectorax)
   call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorx,&
&   vectsize,blockvectorax,vectsize,zero,gramxax,blocksize)

   ABI_ALLOCATE(eigen,(blocksize))

!  write(std_out,*) 'gramxax 2'
!  do ii = 1,blocksize
!  write(std_out,*) gramxax(ii,:)
!  enddo


!  write(std_out,*)'gra',gramxax
!  write(std_out,*)'bx',blockvectorx
!  write(std_out,*)'bax',blockvectorax
!  write(std_out,*)'bbx',blockvectorbx
!  call la_syev(gramxax,eigen,jobz='v')
   lwork=3*blocksize
   ABI_ALLOCATE(work,(lwork))
   call dsyev('v','u',blocksize,gramxax,blocksize,eigen,work,lwork,info)
!  if(mpi_enreg%me_group==0)write(std_out,*)'gramxax apres zheev',gramxax(:,:)
   ABI_DEALLOCATE(work)
!  blockvectorx=matmul(blockvectorx,gramxax)
   call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
&   vectsize,gramxax,blocksize,zero,blockvectordumm,vectsize)
   blockvectorx=blockvectordumm
!  blockvectorax=matmul(blockvectorax,gramxax)
   call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorax,&
&   vectsize,gramxax,blocksize,zero,blockvectordumm,vectsize)
   blockvectorax=blockvectordumm
!  blockvectorbx=matmul(blockvectorbx,gramxax)
   call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbx,&
&   vectsize,gramxax,blocksize,zero,blockvectordumm,vectsize)
   blockvectorbx=blockvectordumm
   do iblocksize=1,blocksize
     lambda(iblocksize,iblocksize)=eigen(iblocksize)
   end do
!  DEBUG
!  write(std_out,*)'lambda',eigen

!  now the main alogrithm
   ABI_ALLOCATE(vectorx,(vectsize,1))
   ABI_ALLOCATE(vectorbx,(vectsize,1))
   ABI_ALLOCATE(vectorax,(vectsize,1))
   ABI_ALLOCATE(vectory,(vectsize,1))
   ABI_ALLOCATE(vectorby,(vectsize,1))
   ABI_ALLOCATE(vectorp,(vectsize,1))
   ABI_ALLOCATE(vectorbp,(vectsize,1))
   ABI_ALLOCATE(vectorap,(vectsize,1))
   ABI_ALLOCATE(vectorr,(vectsize,1))

   iter: do iterationnumber=1,maxiterations
!    DEBUG
     write(std_out,*)'iterationnumber',iterationnumber
!    if(mpi_enreg%me_group==0)write(std_out,*) 'bvx',blockvectorx(10,7)

!    passing x into z
     ABI_ALLOCATE(blockvectorz,(vectsize,blocksize))
     ABI_ALLOCATE(blockvectoraz,(vectsize,blocksize))
     ABI_ALLOCATE(blockvectorbz,(vectsize,blocksize))
     blockvectorz = blockvectorx
     blockvectoraz = blockvectorax
     blockvectorbz = blockvectorbx

     do iblocksize=1,blocksize
!      DEBUG
!      write(std_out,*)'eig number',iblocksize
       vectorx(:,1)=blockvectorx(:,iblocksize)
       vectorbx(:,1)=blockvectorbx(:,iblocksize)
       vectorax(:,1)=blockvectorax(:,iblocksize)
       vectorp(:,1)=blockvectorp(:,iblocksize)
       vectorbp(:,1)=blockvectorbp(:,iblocksize)
       vectorap(:,1)=blockvectorap(:,iblocksize)
       lambda_i=lambda(iblocksize,iblocksize)
       littleblocksize=1
       write(std_out,*) 'eigenvalue number in ',iblocksize, pflag(iblocksize)
       if(present(vxctaulocal))then
         call lobpcgIIIwf(dimffnl,dtfil,dtset,&
&         ffnl,gs_hamk,iterationnumber,&
&         kg_k,kinpw,lmnmax,matblk,mgfft,mpi_enreg,mpsang,&
&         mpssoang,natom,npw_k,ntypat,&
&         nvloc,n4,n5,n6,pcon,ph3d,prtvol,vlocal,&
&         littleblocksize,iblocksize,vectsize,& !bblocksize,vectsize,&
&         pflag(iblocksize),vectorx,vectorbx,vectorax,&
&         blockvectorbz(:,1:iblocksize),& !& blockvectory,blockvectorby,&
&         lambda_i,vectorp,vectorbp,vectorap,vxctaulocal=vxctaulocal)
       else
         call lobpcgIIIwf(dimffnl,dtfil,dtset,&
&         ffnl,gs_hamk,iterationnumber,&
&         kg_k,kinpw,lmnmax,matblk,mgfft,mpi_enreg,mpsang,&
&         mpssoang,natom,npw_k,ntypat,&
&         nvloc,n4,n5,n6,pcon,ph3d,prtvol,vlocal,&
&         littleblocksize,iblocksize,vectsize,& !bblocksize,vectsize,&
&         pflag(iblocksize),vectorx,vectorbx,vectorax,&
&         blockvectorbz(:,1:iblocksize),& !& blockvectory,blockvectorby,&
&         lambda_i,vectorp,vectorbp,vectorap)
       end if
       write(std_out,*) 'eigenvalue number out ',iblocksize, pflag(iblocksize)
       blockvectorx(:,iblocksize)=vectorx(:,1)
       blockvectorbx(:,iblocksize)=vectorbx(:,1)
       blockvectorax(:,iblocksize)=vectorax(:,1)
       blockvectorp(:,iblocksize)=vectorp(:,1)
       blockvectorbp(:,iblocksize)=vectorbp(:,1)
       blockvectorap(:,iblocksize)=vectorap(:,1)
     end do

     ABI_DEALLOCATE(blockvectorz)
     ABI_DEALLOCATE(blockvectoraz)
     ABI_DEALLOCATE(blockvectorbz)

!    gramxax
     call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorx,&
&     vectsize,blockvectorax,vectsize,zero,gramxax,blocksize)
     call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
&     vectsize,blockvectorx,vectsize,zero,gramxbx,blocksize)
!    write(std_out,*)'in iii,xax',gramxax(7,7)
!    write(std_out,*)'in iii xbx',gramxbx(7,7)

!    write(std_out,*) 'gramxax 3'
!    do ii = 1,blocksize
!    write(std_out,*) gramxax(ii,:)
!    enddo
!    write(std_out,*) 'gramxbx 4'
!    do ii = 1,blocksize
!    write(std_out,*) gramxbx(ii,:)
!    enddo


     lwork=3*blocksize
     ABI_ALLOCATE(work,(lwork))
     call dsygv(1,'v','u',blocksize,gramxax,blocksize,gramxbx,blocksize,eigen,&
&     work,lwork,info)
     ABI_DEALLOCATE(work)
!    DEBUG
!    write(std_out,*)'lambda',eigen
!    write(std_out,*)' '
     lambda(:,:)=zero
     do iblocksize=1,blocksize
       lambda(iblocksize,iblocksize)=eigen(iblocksize)
     end do
     ABI_ALLOCATE(coordx,(blocksize,blocksize))
     ABI_ALLOCATE(diagcoordx,(blocksize,blocksize))
     coordx=gramxax

!    rotate all the vectors according to coordx

!    here is a choice for p
     diagcoordx=zero
     do iblocksize=1,blocksize
       diagcoordx(iblocksize,iblocksize) = coordx(iblocksize,iblocksize)
       coordx(iblocksize,iblocksize) = zero
     end do
!    blockvectorxc = matmul(blockvectorx,coordx)
     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
&     vectsize,coordx,blocksize,zero,blockvectorxc,vectsize)
!    blockvectorx = matmul(blockvectorx,diagcoordx) + blockvectorxc
     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
&     vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
     blockvectorx = blockvectordumm + blockvectorxc
!    blockvectorp = matmul(blockvectorp,diagcoordx) + blockvectorxc
     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorp,&
&     vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
     blockvectorp = blockvectordumm + blockvectorxc

!    blockvectorxc = matmul(blockvectorbx,coordx)
     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbx,&
&     vectsize,coordx,blocksize,zero,blockvectorxc,vectsize)
!    blockvectorbx = matmul(blockvectorbx,diagcoordx) + blockvectorxc
     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbx,&
&     vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
     blockvectorbx = blockvectordumm + blockvectorxc
!    blockvectorbp = matmul(blockvectorbp,diagcoordx) + blockvectorxc
     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbp,&
&     vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
     blockvectorbp = blockvectordumm + blockvectorxc

!    blockvectorxc = matmul(blockvectorax,coordx)
     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorax,&
&     vectsize,coordx,blocksize,zero,blockvectorxc,vectsize)
!    blockvectorax = matmul(blockvectorax,diagcoordx) + blockvectorxc
     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorax,&
&     vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
     blockvectorax = blockvectordumm + blockvectorxc
!    blockvectorap = matmul(blockvectorap,diagcoordx) + blockvectorxc
     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorap,&
&     vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
     blockvectorap = blockvectordumm + blockvectorxc

!    DEBUG another choice for p
!    !blockvectorx = matmul(blockvectorx,coordx)
!    call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
!    &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
!    blockvectorx = blockvectordumm
!    !blockvectorbx = matmul(blockvectorbx,coordx)
!    call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbx,&
!    &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
!    blockvectorbx = blockvectordumm
!    !blockvectorax = matmul(blockvectorax,coordx)
!    call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorax,&
!    &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
!    blockvectorax = blockvectordumm
!    !blockvectorp = matmul(blockvectorp,coordx)
!    call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorp,&
!    &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
!    blockvectorp = blockvectordumm
!    !blockvectorbp = matmul(blockvectorbp,coordx)
!    call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbp,&
!    &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
!    blockvectorbp = blockvectordumm
!    !blockvectorap = matmul(blockvectorap,coordx)
!    call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorap,&
!    &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
!    blockvectorap = blockvectordumm
!    ENDDEBUG

!    DEBUG
!    gramxax
!    call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorx,&
!    &               vectsize,blockvectorax,vectsize,zero,gramxax,blocksize)
!    call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
!    &               vectsize,blockvectorx,vectsize,zero,gramxbx,blocksize)
!    write(std_out,*)'in iii,xax after',gramxax
!    write(std_out,*)'in iii xbx after',gramxbx
!    ENDDEBUG
     ABI_DEALLOCATE(coordx)
     ABI_DEALLOCATE(diagcoordx)

   end do iter

   call precon2(blockvectorbx,lambda,blocksize,&
&   iterationnumber,kinpw,mpi_enreg,npw_k,my_nspinor,&
&   optekin,optpcon,pcon,blockvectorax,blockvectorr,vectsize)

   residualnorms=sum(blockvectorr**2,dim=1)
   resid_k(bblocksize+1:bblocksize+blocksize)=residualnorms(1:blocksize)
   residualnorms=sqrt(residualnorms)
   write(std_out,*) 'residualnorm at the end',residualnorms

   ABI_DEALLOCATE(eigen)
!  write(std_out,*)'residualnorm at the end',residualnorms

!  epilogue
!  residualnorms=sqrt(sum(abs(blockvectorr)**2,dim=1))
!  write(std_out,*)'residualnorm at the end',residualnorms
   do iblocksize=1,blocksize
     iband=iblocksize+(iblock-1)*blocksize
     if (mpi_enreg%me_g0 == 1) then
       cg(1,cgindex(iband))=blockvectorx(1,iblocksize)
       cg(2,cgindex(iband))=zero
       cg(1,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(2:rvectsize,iblocksize)/sq2
       cg(2,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(rvectsize+1:vectsize,iblocksize)/sq2
     else
       cg(1,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(1:rvectsize,iblocksize)/sq2
       cg(2,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(rvectsize+1:vectsize,iblocksize)/sq2
     end if
   end do
   if(gen_eigenpb) then
     do iblocksize=1,blocksize
       iband=iblocksize+(iblock-1)*blocksize
       if (mpi_enreg%me_g0 == 1) then
         gsc(1,gscindex(iband))=blockvectorbx(1,iblocksize)
         gsc(2,gscindex(iband))=zero
         gsc(1,gscindex(iband)+1:gscindex(iband+1)-1)=blockvectorbx(2:rvectsize,iblocksize)/sq2
         gsc(2,gscindex(iband)+1:gscindex(iband+1)-1)=blockvectorbx(rvectsize+1:vectsize,iblocksize)/sq2
       else
         gsc(1,gscindex(iband):gscindex(iband+1)-1)=blockvectorbx(1:rvectsize,iblocksize)/sq2
         gsc(2,gscindex(iband):gscindex(iband+1)-1)=blockvectorbx(rvectsize+1:vectsize,iblocksize)/sq2
       end if
     end do
   end if

!  this should not exist,since this induce one too much getghc.lazy programming....
!  call operatorh(blockvectorx,blockvectorax,subham,subvnl)!fill also subham, subvnl
   ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(gvnlc,(2,npw_k*my_nspinor*blocksize))
   isubh=1+2*(iblock-1)*blocksize*((iblock-1)*blocksize+1)/2
   do iblocksize=1,blocksize
     iband=iblocksize
     if (mpi_enreg%me_g0 == 1) then
       cwavef(1,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(2:npw_k*my_nspinor,iblocksize)/sq2
       cwavef(2,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)/sq2
       cwavef(2,cgindex(iband))=zero
       cwavef(1,cgindex(iband))=blockvectorx(1,iblocksize)
     else
       cwavef(1,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(1:npw_k*my_nspinor,iblocksize)/sq2
       cwavef(2,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)/sq2
     end if
   end do
   tim_getghc=7; sij_opt=0
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
     iband=iblocksize
     if (mpi_enreg%me_g0 == 1) then
       blockvectorax(2:npw_k*my_nspinor,iblocksize)=gwavef(1,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
       blockvectorax(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)=gwavef(2,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
       blockvectorax(1,iblocksize)=gwavef(1,cgindex(iband))
     else
       blockvectorax(1:npw_k*my_nspinor,iblocksize)=gwavef(1,cgindex(iband):cgindex(iband+1)-1)*sq2
       blockvectorax(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)=gwavef(2,cgindex(iband):cgindex(iband+1)-1)*sq2
     end if
   end do
   do iblocksize=1,blocksize
     do ii=1,(iblock-1)*blocksize+iblocksize
       iwavef=(ii-1)*npw_k*my_nspinor+icg
       if (mpi_enreg%me_g0 == 1) then
         ipw1=2;chcre=0.5_dp*cg(1,1+iwavef)*gwavef(1,cgindex(iblocksize))
       else
         ipw1=1;chcre=zero
       end if
       if (gs_hamk%usepaw==1) then
         do ipw=ipw1,npw_k*my_nspinor
           cgreipw=cg(1,ipw+iwavef);cgimipw=cg(2,ipw+iwavef)
           chcre=chcre+cgreipw*gwavef(1,ipw+(iblocksize-1)*npw_k*my_nspinor)+cgimipw*gwavef(2,ipw+(iblocksize-1)*npw_k*my_nspinor)
         end do
       else
         if (mpi_enreg%me_g0 == 1) then
           cvcre=0.5_dp*cg(1,1+iwavef)*gvnlc(1,cgindex(iblocksize))
         else
           cvcre=zero
         end if
         do ipw=ipw1,npw_k*my_nspinor
           cgreipw=cg(1,ipw+iwavef);cgimipw=cg(2,ipw+iwavef)
           chcre=chcre+cgreipw*gwavef(1,ipw+(iblocksize-1)*npw_k*my_nspinor)+cgimipw*gwavef(2,ipw+(iblocksize-1)*npw_k*my_nspinor)
           cvcre=cvcre+cgreipw*gvnlc(1,ipw+(iblocksize-1)*npw_k*my_nspinor)+cgimipw*gvnlc(2,ipw+(iblocksize-1)*npw_k*my_nspinor)
         end do
!        store real and imag parts in hermitian storage mode:
         subvnl(isubh)=2.0_dp*cvcre ; subvnl(isubh+1)=zero
       end if
!      store real and imag parts in hermitian storage mode:
       subham(isubh)=2.0_dp*chcre ; subham(isubh+1)=zero
       isubh=isubh+2
     end do
   end do
   ABI_DEALLOCATE(cwavef)
   ABI_DEALLOCATE(gwavef)
   ABI_DEALLOCATE(gvnlc)

!  call operators(blockvectorx,blockvectorbx,subovl)!fill also  subovl
   if((gen_eigenpb).and.(use_subovl==1)) then
     ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
     ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor))
     isubo=1+2*(iblock-1)*blocksize*((iblock-1)*blocksize+1)/2
     do iblocksize=1,blocksize
       if (mpi_enreg%me_g0 == 1) then
         cwavef(1,2:npw_k*my_nspinor)=blockvectorx(2:npw_k*my_nspinor,iblocksize)/sq2
         cwavef(2,2:npw_k*my_nspinor)=blockvectorx(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)/sq2
         cwavef(1,1)=blockvectorx(1,iblocksize)
         cwavef(2,1)=zero
       else
         cwavef(1,1:npw_k*my_nspinor)=blockvectorx(1:npw_k*my_nspinor,iblocksize)/sq2
         cwavef(2,1:npw_k*my_nspinor)=blockvectorx(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)/sq2
       end if
!      call to nonlop: compute <g|s|c>
       choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
       call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
&       dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
&       istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
&       mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
&       gs_hamk%nloalg,nnlout,npw_k,npw_k,my_nspinor,dtset%nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
&       gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,signs,gs_hamk%sij,&
&       gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef,use_gpu_cuda=dtset%use_gpu_cuda)
       if (mpi_enreg%me_g0 == 1) then
         blockvectorbx(2:npw_k*my_nspinor,iblocksize)=gwavef(1,2:npw_k*my_nspinor)*sq2
         blockvectorbx(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)=gwavef(2,2:npw_k*my_nspinor)*sq2
         blockvectorbx(1,iblocksize)=gwavef(1,1)
       else
         blockvectorbx(1:npw_k*my_nspinor,iblocksize)=gwavef(1,1:npw_k*my_nspinor)*sq2
         blockvectorbx(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)=gwavef(2,1:npw_k*my_nspinor)*sq2
       end if
       do ii=1,(iblock-1)*blocksize+iblocksize
         iwavef=(ii-1)*npw_k*my_nspinor+icg
         if (istwf_k==2 .and. mpi_enreg%me_g0 == 1) then
           ipw1=2;cscre=0.5_dp*cg(1,1+iwavef)*gwavef(1,1)
         else
           ipw1=1;cscre=zero
         end if
         do ipw=ipw1,npw_k*my_nspinor
           cscre=cscre+cg(1,ipw+iwavef)*gwavef(1,ipw)+cg(2,ipw+iwavef)*gwavef(2,ipw)
         end do
         cscre=2.0_dp*cscre
!        store real and imag parts in hermitian storage mode:
         subovl(isubo)=cscre ; subovl(isubo+1)=zero
         isubo=isubo+2
       end do
     end do
     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(gwavef)
   end if
!  write(std_out,*)'residualnorm at the end',residualnorms

!  DEBUG
!  write(std_out,*)'frozen_count',frozen_count,'restart_count',restart_count
!  ENDDEBUG
   ABI_DEALLOCATE(vectorx)
   ABI_DEALLOCATE(vectorbx)
   ABI_DEALLOCATE(vectorax)
   ABI_DEALLOCATE(vectory)
   ABI_DEALLOCATE(vectorby)
   ABI_DEALLOCATE(vectorp)
   ABI_DEALLOCATE(vectorbp)
   ABI_DEALLOCATE(vectorap)
   ABI_DEALLOCATE(vectorr)
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

!write(std_out,*)'mpi_enreg%me,icount_ghc',mpi_enreg%me,icount_ghc
!DEBUG
!write(std_out,*)'end lobpcg'
!stop
!ENDDEBUG

end subroutine lobpcgIIwf
!!***
