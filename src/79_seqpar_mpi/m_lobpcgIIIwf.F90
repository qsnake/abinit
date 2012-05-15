!{\src2tex{textfont=tt}}
!!****f* ABINIT/lobpcgIIIwf
!! NAME
!! lobpcgIIIwf
!!
!! FUNCTION
!! TO BE DESCRIBED 090830
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      lobpcgIIwf
!!
!! CHILDREN
!!      dgemm,dtrsm,getghc,nonlop,orthonormalize
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine lobpcgIIIwf(dimffnl,dtfil,dtset,&
& ffnl,gs_hamk,iterationnumber,&
& kg_k,kinpw,lmnmax,matblk,mgfft,mpi_enreg,mpsang,&
& mpssoang,natom,npw_k,ntypat,&
& nvloc,n4,n5,n6,&
& pcon,ph3d,prtvol,vlocal,&
& blocksize,bblocksize,vectsize,pflag,&
& blockvectorx,blockvectorbx,blockvectorax,blockvectorby,lambda,&
& blockvectorp,blockvectorbp,blockvectorap,&
& vxctaulocal) ! optional argument

 use m_profiling

use defs_basis
use defs_datatypes
 use defs_abitypes

#if defined HAVE_MPI2
 use mpi
#endif
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lobpcgIIIwf'
 use interfaces_53_abiutil
 use interfaces_66_wfs
 use interfaces_79_seqpar_mpi, except_this_one => lobpcgIIIwf
!End of the abilint section

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif
!Arguments -------------------------------
type(gs_hamiltonian_type) :: gs_hamk
integer :: lmnmax,matblk,mgfft,mpsang,mpssoang,n4,n5
integer :: n6,natom,npw_k,ntypat,nvloc,prtvol
integer :: dimffnl
integer :: iterationnumber
type(datafiles_type) :: dtfil
type(dataset_type) :: dtset
type(mpi_type) :: mpi_enreg
integer :: kg_k(3,npw_k)
real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
real(dp) :: kinpw(npw_k),ph3d(2,npw_k,matblk)
real(dp) :: vlocal(n4,n5,n6,nvloc)
integer:: blocksize,vectsize,bblocksize
logical :: pflag(blocksize)
real(dp)    :: lambda(blocksize,blocksize)
real(dp) :: pcon(npw_k,blocksize)
real(dp) :: blockvectorx(vectsize,blocksize),blockvectorax(vectsize,blocksize),&
&                 blockvectorbx(vectsize,blocksize),&
&                 blockvectorby(vectsize,bblocksize)
real(dp) :: blockvectorp(vectsize,blocksize),blockvectorap(vectsize,blocksize),&
&                 blockvectorbp(vectsize,blocksize)
real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)
!Local variables -------------------------
real(dp) :: sq2
integer:: ii,maxiterations
integer:: bigorder,info,lwork,my_nspinor,tim_getghc,sij_opt,iminresid
integer:: iblocksize,optekin=0, optpcon=1, restart, cond_try, i1, i2, i3, i4,nkpg
integer:: iband,istwf_k
integer :: choice, cpopt, signs, idir, tim_nonlop, paw_opt, nnlout

logical::gen_eigenpb
real(dp) :: dum
real(dp), parameter :: tolerance1 = 1.e-13
real(dp) :: tolerance2=1.e2
 type(cprj_type) :: cprj_dum(1,1)
real(dp) :: blockvectorz(vectsize,blocksize),blockvectoraz(vectsize,blocksize),&
&                 blockvectorbz(vectsize,blocksize)
real(dp) :: blockvectorp_old(vectsize,blocksize),blockvectorap_old(vectsize,blocksize),&
&                 blockvectorbp_old(vectsize,blocksize)
real(dp), allocatable ::  blockvectorr(:,:),blockvectorar(:,:),blockvectorbr(:,:),&
& blockvectorr1(:,:),&
& blockvectordumm(:,:),&
& blockvectorritz(:,:),&
& blockvectoraritz(:,:),&
& blockvectorbritz(:,:),&
& blockvectorrritz(:,:),&
& blockvectorritzdumm(:,:),&
& gramxax(:,:),gramxar(:,:),gramxap(:,:),gramrar(:,:),gramrap(:,:),&
& grampap(:,:),&
& gramxbx(:,:),gramxbr(:,:),gramxbp(:,:),gramrbr(:,:),gramrbp(:,:),&
& grampbp(:,:),&
& coordx(:,:),&
& grama(:,:),gramb(:,:),gramyx(:,:),&
      & kpg_dum(:,:),work(:),dummy1(:),dummy2(:,:)
real(dp), allocatable :: gwavef(:,:),cwavef(:,:),gvnlc(:,:)
!  real(dp) :: zero,one
real(dp), allocatable :: residualnorms(:),eigen(:),residualnorms1(:)
real(dp), allocatable :: residualnormritz(:),normritz(:)

! *********************************************************************

if (dtset%useria == 0) then
 tolerance2 = 100.d0
else
 tolerance2 = real(dtset%useria,dp)
endif

!correspondence with abinit. here for real wf but in complex mode
!this is the index of a given band
!  cgindex(iblocksize)=npw_k*my_nspinor*(iblocksize-1)+icg+1
gen_eigenpb=(gs_hamk%usepaw==1)
my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
!  czero=dcmplx(zero,zero)
!  cone=dcmplx(one,zero)
sq2=sqrt(2.0_dp)
!vectsize=npw_k*my_nspinor
!blocksize=(nband_k-1)/nbdblock+1
!bblocksize=(iblock-1)*blocksize
istwf_k=gs_hamk%istwf_k
maxiterations=dtset%nline

!passing x into z and p into p_old for restart
blockvectorz = blockvectorx
blockvectoraz = blockvectorax
blockvectorbz = blockvectorbx

!allocations
ABI_ALLOCATE(blockvectorr,(vectsize,blocksize))
ABI_ALLOCATE(blockvectorar,(vectsize,blocksize))
ABI_ALLOCATE(blockvectorbr,(vectsize,blocksize))
ABI_ALLOCATE(blockvectorr1,(vectsize,blocksize))
ABI_ALLOCATE(blockvectordumm,(vectsize,blocksize))
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
ABI_ALLOCATE(residualnorms,(blocksize))
ABI_ALLOCATE(residualnorms1,(blocksize))

!construct residual
!blockvectorr=blockvectorax-matmul(blockvectorx,lambda)

 call precon2(blockvectorbx,lambda,blocksize,&
&             iterationnumber,kinpw,mpi_enreg,npw_k,my_nspinor,&
&             optekin,optpcon,pcon,blockvectorax,blockvectorr,vectsize)
!  blockvectorr(:,iblocksize)=blockvectorax(:,iblocksize)-lambda(iblocksize,iblocksize)*blockvectorbx(:,iblocksize)

residualnorms=sqrt(sum(abs(blockvectorr)**2,dim=1))

!DEBUG
! resid_k(bblocksize+1:bblocksize+blocksize)=residualnorms(1:blocksize)
write(std_out,*)'residualnorm before lobpcgiii',residualnorms
!ENDEBUG
if(residualnorms(1) > tolerance1)then  !DEBUG this is the wrong condition if blocksize /= 1

 call start_lobpcg ! r orthonormal to x and compute ar
 !call orthonormalize(blockvectorp,blockvectorbp)
 if(pflag(1)) then   !DEBUG this is the wrong condition if blocksize /= 1
  !DEBUG
  !write(std_out,*)'blockvectorp,blockvectorbp,blockvectorap'
  !write(std_out,*)blockvectorp
  !write(std_out,*)blockvectorbp
  !write(std_out,*)blockvectorap
  !ENDDEBUG
  !call zorthonormalize(blockvectorp,blockvectorbp,blockvectorap)
  call orthonormalize(blockvectorp,blockvectorbp,blocksize,mpi_enreg,grampbp,vectsize)
  call dtrsm('r','u','n','n',vectsize,blocksize,one,grampbp,blocksize,&
  &              blockvectorbp,vectsize)
  !blockvectorap=matmul(blockvectorap,grampbp)
  call dtrsm('r','u','n','n',vectsize,blocksize,one,grampbp,blocksize,&
  &              blockvectorap,vectsize)
 end if
 if (.not.pflag(1)) then  !DEBUG  this is the wrong condition if blocksize /= 1
  restart=1
 else
  restart=0
  blockvectorp_old = blockvectorp
  blockvectorap_old = blockvectorap
  blockvectorbp_old = blockvectorbp
 end if
 !gramxar=matmul(transpose(blockvectorax),blockvectorr)
 !gramrar=matmul(transpose(blockvectorar),blockvectorr)
 !gramxax=matmul(transpose(blockvectorax),blockvectorx)=lambda
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorax,&
 &               vectsize,blockvectorr,vectsize,zero,gramxar,blocksize)
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorar,&
 &               vectsize,blockvectorr,vectsize,zero,gramrar,blocksize)
 !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorax,&
 !     &               vectsize,blockvectorx,vectsize,zero,gramxax,blocksize)
 gramxax = lambda

 !gramxbx=matmul(transpose(blockvectorbx),blockvectorx)=identity
 !gramrbr=matmul(transpose(blockvectorbr),blockvectorr)=identity
 !gramxbr=matmul(transpose(blockvectorbx),blockvectorr)=zero
 !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
 !     &               vectsize,blockvectorx,vectsize,zero,gramxbx,blocksize)
 !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbr,&
 !     &               vectsize,blockvectorr,vectsize,zero,gramrbr,blocksize)
 !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
 !     &               vectsize,blockvectorr,vectsize,zero,gramxbr,blocksize)
 gramxbx = zero
 gramrbr = zero
 gramxbr = zero
 do iblocksize = 1,blocksize
  gramxbx(iblocksize,iblocksize)=one
  gramrbr(iblocksize,iblocksize)=one
 end do

 i1=0;i2=blocksize;i3=2*blocksize;i4=3*blocksize
 cond: do cond_try=1,2 !2 when restart
  call construct_gram
  !write(std_out,*) 'grama 5'
  !do ii = 1,bigorder
  !write(std_out,*) grama(ii,:)
  !enddo

  !write(std_out,*) 'gramb 6'
  !do ii = 1,bigorder
  !write(std_out,*) gramb(ii,:)
  !enddo


  !DEBUG
  !call la_sygv(grama,gramb,eigen,itype=1,jobz='v')
  !ENDDEBUG
  lwork=3*bigorder
  ABI_ALLOCATE(work,(lwork))

  call dsygv(1,'v','u',bigorder,grama,bigorder,gramb,bigorder,eigen,&
  &               work,lwork,info)
  ABI_DEALLOCATE(work)
  do iblocksize=1,blocksize
   lambda(iblocksize,iblocksize)=eigen(iblocksize)
  end do
  !DEBUG
  !write(std_out,*)'eigen',eigen(1:blocksize)
  !ENDDEBUG
  coordx=grama(:,1:blocksize)
  call rotate_vectors


  call precon2(blockvectorbx,lambda,blocksize,&
&      iterationnumber,kinpw,mpi_enreg,npw_k,my_nspinor,&
&      optekin,optpcon,pcon,blockvectorax,blockvectorr1,vectsize)

  residualnorms1=sqrt(sum(abs(blockvectorr1)**2,dim=1))
  !DEBUG
  write(std_out,*)'residualnorm after lobpcgiii',residualnorms1
  !ENDDEBUG

  !DEBUG
  !print out the B-scalar product of new approximation vector and former iteration smallest eigenvector
  !we use gramb matrix which is of no use anymore
  ABI_DEALLOCATE(gramb)
  ABI_ALLOCATE(gramb,(bblocksize,1))
  call dgemm('t','n',bblocksize,1,vectsize,one,blockvectorby,&
  &               vectsize,blockvectorx,vectsize,zero,gramb,bblocksize)
  write(std_out,*) 'B-scalar product with former lesser-order iterate',  gramb
  blockvectorr1(:,1) = blockvectorx(:,1)
  if (bblocksize>1)then
   call dgemm('n','n',vectsize,1,bblocksize-1,one,blockvectorby(:,1:bblocksize-1),&
   &               vectsize,gramb(1:bblocksize-1,1),bblocksize-1,zero,blockvectorr1(:,1),vectsize)
   write(std_out,*) 'B-distance XY, sqrt(XX)',sqrt(sum(abs(blockvectorx(:,1)-blockvectorr1(:,1))**2)),sqrt(abs(gramb(bblocksize,1)))
  else
   write(std_out,*) 'B-distance XY, sqrt(XX)',' NA',sqrt(abs(gramb(bblocksize,1)))
  endif

  !ENDDEBUG

  do iblocksize=1,blocksize    !DEBUG this do will work only if blocksize = 1
   if (residualnorms1(iblocksize) > tolerance2*residualnorms(iblocksize)) then
    write(std_out,*) 'restart apply',restart
    !the eigenvector we seek is one of the other Ritz vector
    if (restart==0) then
     call  apply_flip_flop
     exit cond
     ABI_DEALLOCATE(blockvectorritz)
     ABI_DEALLOCATE(blockvectoraritz)
     ABI_DEALLOCATE(blockvectorbritz)
     ABI_DEALLOCATE(blockvectorritzdumm)
     ABI_DEALLOCATE(blockvectorrritz)
     ABI_DEALLOCATE(residualnormritz)
     ABI_DEALLOCATE(normritz)
    end if
   else
    pflag = .true.
    write(std_out,*) 'set pftrue'
    exit cond
   end if
  end do
  ABI_DEALLOCATE(grama)
  ABI_DEALLOCATE(gramb)
  ABI_DEALLOCATE(eigen)
 end do cond
else
 blockvectorp = zero
 blockvectorap = zero
 blockvectorbp = zero
 pflag = .false.
 write(std_out,*) 'set pffalse'
end if

!write(std_out,*)'blockvectorr',blockvectorr
!write(std_out,*)'blockvectorx',blockvectorx
!write(std_out,*)'blockvectorax',blockvectorax
ABI_DEALLOCATE(blockvectorr)
ABI_DEALLOCATE(blockvectorar)
ABI_DEALLOCATE(blockvectorbr)
ABI_DEALLOCATE(blockvectorr1)
ABI_DEALLOCATE(gramyx)
ABI_DEALLOCATE(blockvectordumm)
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
ABI_DEALLOCATE(residualnorms)
ABI_DEALLOCATE(residualnorms1)

contains

!!***

!!****f* ABINIT/construct_gram
!! NAME
!! construct_gram
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_lobpcgIIIwf
!!
!! CHILDREN
!!      dgemm,dtrsm,getghc,nonlop,orthonormalize
!!
!! SOURCE
subroutine construct_gram

! *********************************************************************
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'construct_gram'
!End of the abilint section

if (restart==0) then
 !gramxap=matmul(transpose(blockvectorax),blockvectorp)
 !gramrap=matmul(transpose(blockvectorar),blockvectorp)
 !grampap=matmul(transpose(blockvectorap),blockvectorp)
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorax,&
 &                    vectsize,blockvectorp,vectsize,zero,gramxap,blocksize)
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorar,&
 &                    vectsize,blockvectorp,vectsize,zero,gramrap,blocksize)
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorap,&
 &                    vectsize,blockvectorp,vectsize,zero,grampap,blocksize)
 bigorder=i4
 ABI_ALLOCATE(grama,(i4,i4))
 ABI_ALLOCATE(gramb,(i4,i4))
 ABI_ALLOCATE(eigen,(i4))
 ABI_ALLOCATE(coordx,(i4,blocksize))
 grama(i1+1:i2,i1+1:i2)=gramxax
 grama(i1+1:i2,i2+1:i3)=gramxar
 grama(i1+1:i2,i3+1:i4)=gramxap
 !grama(i2+1:i3,i1+1:i2)=transpos(gramxar)
 grama(i2+1:i3,i2+1:i3)=gramrar
 grama(i2+1:i3,i3+1:i4)=gramrap
 !grama(i3+1:i4,i1+1:i2)=transpos(gramxap)
 !grama(i3+1:i4,i2+1:i3)=transpos(gramrap)
 grama(i3+1:i4,i3+1:i4)=grampap

 !gramxbp=matmul(transpose(blockvectorbx),blockvectorp)
 !gramrbp=matmul(transpose(blockvectorbr),blockvectorp)
 !grampbp=matmul(transpose(blockvectorbp),blockvectorp)=identity
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
 &                    vectsize,blockvectorp,vectsize,zero,gramxbp,blocksize)
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbr,&
 &                    vectsize,blockvectorp,vectsize,zero,gramrbp,blocksize)
 !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbp,&
 !     &                    vectsize,blockvectorp,vectsize,zero,grampbp,blocksize)
 grampbp = zero
 do iblocksize = 1,blocksize
  grampbp(iblocksize,iblocksize)=one
 end do

 gramb(i1+1:i2,i1+1:i2)=gramxbx
 gramb(i1+1:i2,i2+1:i3)=gramxbr
 gramb(i1+1:i2,i3+1:i4)=gramxbp
 !gramb(i2+1:i3,i1+1:i2)=transpos(gramxbr)
 gramb(i2+1:i3,i2+1:i3)=gramrbr
 gramb(i2+1:i3,i3+1:i4)=gramrbp
 !gramb(i3+1:i4,i1+1:i2)=transpos(gramxbp)
 !gramb(i3+1:i4,i2+1:i3)=transpos(gramrbp)
 gramb(i3+1:i4,i3+1:i4)=grampbp
else
 bigorder=i3
 ABI_ALLOCATE(grama,(i3,i3))
 ABI_ALLOCATE(gramb,(i3,i3))
 ABI_ALLOCATE(eigen,(i3))
 ABI_ALLOCATE(coordx,(i3,blocksize))
 grama(i1+1:i2,i1+1:i2)=gramxax
 grama(i1+1:i2,i2+1:i3)=gramxar
 !grama(i2+1:i3,i1+1:i2)=transpos(gramxar)
 grama(i2+1:i3,i2+1:i3)=gramrar
 gramb(i1+1:i2,i1+1:i2)=gramxbx
 gramb(i1+1:i2,i2+1:i3)=gramxbr
 !gramb(i2+1:i3,i1+1:i2)=transpos(gramxbr)
 gramb(i2+1:i3,i2+1:i3)=gramrbr
end if

end  subroutine construct_gram
!!***

!!****f* ABINIT/rotate_vectors
!! NAME
!! rotate_vectors
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_lobpcgIIIwf
!!
!! CHILDREN
!!      dgemm,dtrsm,getghc,nonlop,orthonormalize
!!
!! SOURCE
subroutine rotate_vectors

! *********************************************************************
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rotate_vectors'
!End of the abilint section

if (restart==0) then
 !    blockvectorp=matmul(blockvectorr,coordx(i2+1:i3,:))+&
 !          &matmul(blockvectorp,coordx(i3+1:i4,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorr,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectordumm,vectsize)
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorp,&
 &               vectsize,coordx(i3+1:i4,:),blocksize,one,blockvectordumm,vectsize)
 blockvectorp=blockvectordumm
 !    blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))+&
 !          &matmul(blockvectorap,coordx(i3+1:i4,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorar,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectordumm,vectsize)
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorap,&
 &               vectsize,coordx(i3+1:i4,:),blocksize,one,blockvectordumm,vectsize)
 blockvectorap=blockvectordumm
 !    blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))+&
 !          &matmul(blockvectorbp,coordx(i3+1:i4,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbr,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectordumm,vectsize)
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbp,&
 &               vectsize,coordx(i3+1:i4,:),blocksize,one,blockvectordumm,vectsize)
 blockvectorbp=blockvectordumm
else
 !blockvectorp =matmul(blockvectorr,coordx(i2+1:i3,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorr,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectorp,vectsize)
 !blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorar,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectorap,vectsize)
 !    blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbr,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectorbp,vectsize)
end if

!blockvectorx = matmul(blockvectorx,coordx(i1+1:i2,:))+blockvectorp
call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
&               vectsize,coordx(i1+1:i2,:),blocksize,zero,blockvectordumm,vectsize)
blockvectorx = blockvectordumm+blockvectorp
!blockvectorax= matmul(blockvectorax,coordx(i1+1:i2,:))+blockvectorap
call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorax,&
&               vectsize,coordx(i1+1:i2,:),blocksize,zero,blockvectordumm,vectsize)
blockvectorax = blockvectordumm+blockvectorap
!blockvectorbx= matmul(blockvectorbx,coordx(i1+1:i2,:))+blockvectorbp
call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbx,&
&               vectsize,coordx(i1+1:i2,:),blocksize,zero,blockvectordumm,vectsize)
blockvectorbx = blockvectordumm+blockvectorbp
ABI_DEALLOCATE(coordx)
end subroutine rotate_vectors
!!***

!!****f* ABINIT/apply_flip_flop
!! NAME
!! apply_flip_flop
!!
!! FUNCTION
!!restore former vector
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_lobpcgIIIwf
!!
!! CHILDREN
!!      dgemm,dtrsm,getghc,nonlop,orthonormalize
!!
!! SOURCE
subroutine apply_flip_flop

 use m_linalg_interfaces
!Arguments -------------------------------

!Local variables -------------------------

! *********************************************************************

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'apply_flip_flop'
 use interfaces_66_wfs
!End of the abilint section

 blockvectorp = blockvectorp_old
 blockvectorap = blockvectorap_old
 blockvectorbp = blockvectorbp_old

 blockvectorx = blockvectorz
 blockvectorax = blockvectoraz
 blockvectorbx = blockvectorbz
!compute all Ritz vectors
 ABI_ALLOCATE(blockvectorritz,(vectsize,bigorder))
 ABI_ALLOCATE(blockvectoraritz,(vectsize,bigorder))
 ABI_ALLOCATE(blockvectorbritz,(vectsize,bigorder))
 ABI_ALLOCATE(blockvectorritzdumm,(vectsize,bigorder))
 ABI_ALLOCATE(blockvectorrritz,(vectsize,bigorder))
 ABI_ALLOCATE(residualnormritz,(bigorder))
 ABI_ALLOCATE(normritz,(bigorder))
 blockvectorritz(:,i1+1:i2) = blockvectorx(:,:)
 blockvectorritz(:,i2+1:i3) = blockvectorr(:,:)
 blockvectorritz(:,i3+1:i4) = blockvectorp(:,:)
 blockvectoraritz(:,i1+1:i2) = blockvectorax(:,:)
 blockvectoraritz(:,i2+1:i3) = blockvectorar(:,:)
 blockvectoraritz(:,i3+1:i4) = blockvectorap(:,:)
 blockvectorbritz(:,i1+1:i2) = blockvectorbx(:,:)
 blockvectorbritz(:,i2+1:i3) = blockvectorbr(:,:)
 blockvectorbritz(:,i3+1:i4) = blockvectorbp(:,:)
!blockvectorritz=matmul(blockvectorritz,grama)
 call dgemm('n','n',vectsize,bigorder,bigorder,one,blockvectorritz,&
&               vectsize,grama,bigorder,zero,blockvectorritzdumm,vectsize)
 blockvectorritz(:,:)=blockvectorritzdumm(:,:)
!write(std_out,*) 'diff3'
!write(std_out,*) blockvectorritz(1:10,1)-(grama(1,1)*blockvectorx(1:10,1)+grama(2,1)*blockvectorr(1:10,1)+grama(3,1)*blockvectorp(1:10,1))

!blockvectoraritz=matmul(blockvectoraritz,grama)
 call dgemm('n','n',vectsize,bigorder,bigorder,one,blockvectoraritz,&
&               vectsize,grama,bigorder,zero,blockvectorritzdumm,vectsize)
 blockvectoraritz(:,:)=blockvectorritzdumm(:,:)
!blockvectorbritz=matmul(blockvectorritz,grama)
 call dgemm('n','n',vectsize,bigorder,bigorder,one,blockvectorbritz,&
&               vectsize,grama,bigorder,zero,blockvectorritzdumm,vectsize)
 blockvectorbritz(:,:)=blockvectorritzdumm(:,:)
!     !blockvectorp=matmul(blockvectorz,coordx(i1+1:i2,:))+blockvectorp
!     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorz,&
!          &               vectsize,coordx(i1+1:i2,:),blocksize,one,blockvectorp,vectsize)
!     !blockvectorap=matmul(blockvectoraz,coordx(i1+1:i2,:))+blockvectorap
!     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectoraz,&
!          &               vectsize,coordx(i1+1:i2,:),blocksize,one,blockvectorap,vectsize)
!     !blockvectorbp=matmul(blockvectorbz,coordx(i1+1:i2,:))+blockvectorbp
!     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbz,&
!          &               vectsize,coordx(i1+1:i2,:),blocksize,one,blockvectorbp,vectsize)

!compute the residuals for the Ritz vectors other than the first
!     residualnormritz(i1+1:i2) = residualnorm1(:)

 do ii=i1+1,bigorder
 call precon2(blockvectorbritz(:,ii),eigen(ii),1,&
 &        iterationnumber,kinpw,mpi_enreg,npw_k,my_nspinor,&
 &        optekin,optpcon,pcon,&
 &        blockvectoraritz(:,ii),blockvectorrritz(:,ii),vectsize)
 end do
 residualnormritz=sqrt(sum(abs(blockvectorrritz)**2,dim=1))
 normritz=sqrt(sum(abs(blockvectorritz)**2,dim=1))

 write(std_out,*) 'residualnorms1',residualnorms1
 write(std_out,*) 'residualnormritz',residualnormritz
 write(std_out,*) 'normritz',normritz
 write(std_out,*) 'quotient',residualnormritz/normritz
 iminresid = 1
 do ii = 2,bigorder
 if(residualnormritz(ii)/normritz(ii)<residualnormritz(iminresid)/normritz(iminresid))then
  iminresid = ii
 endif
 enddo
 write(std_out,*) 'iminresid',iminresid
!following lines will not work with blocksize > 1
 blockvectorp(:,1)  = grama(2,iminresid)*blockvectorr(:,1)+grama(3,iminresid)*blockvectorp(:,1)
 blockvectorap(:,1) = grama(2,iminresid)*blockvectorar(:,1)+grama(3,iminresid)*blockvectorap(:,1)
 blockvectorbp(:,1) = grama(2,iminresid)*blockvectorbr(:,1)+grama(3,iminresid)*blockvectorbp(:,1)
 blockvectorx(:,1)  = blockvectorritz(:,iminresid)
 blockvectorax(:,1) = blockvectoraritz(:,iminresid)
 blockvectorbx(:,1) = blockvectorbritz(:,iminresid)
 pflag = .true.
 write(std_out,*) 'set pftrue2'

end subroutine apply_flip_flop
!!***

!!****f* ABINIT/start_lobpcg
!! NAME
!! start_lobpcg
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_lobpcgIIIwf
!!
!! CHILDREN
!!      dgemm,dtrsm,getghc,nonlop,orthonormalize
!!
!! SOURCE
 subroutine start_lobpcg

!if(abs(sum(residualnorms)) < 1.d-10) exit

 use m_profiling
!!$if (.false.)then!(bbblocksize>0) then !(iblock /=1) then !residuals orthogonal to blockvectorby
!!! !   blockvectorr=blockvectorr-&
!!! !           &matmul(blockvectory,matmul(transpose(blockvectorby),blockvectorr))
!!! call dgemm('t','n',bblocksize,blocksize,vectsize,one,blockvectorby,&
!!! &               vectsize,blockvectorr,vectsize,zero,gramyx,bblocksize)
!!! call dgemm('n','n',vectsize,blocksize,bblocksize,one,blockvectory,&
!!! &               vectsize,gramyx,bblocksize,zero,blockvectordumm,vectsize)
!!! blockvectorr=blockvectorr-blockvectordumm
!!$end if
!residuals orthogonal to blockvectorx
!  blockvectorr=blockvectorr-&
!          &matmul(blockvectorx,matmul(transpose(blockvectorbx),blockvectorr))
! *********************************************************************
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'start_lobpcg'
 use interfaces_53_abiutil
 use interfaces_65_nonlocal
 use interfaces_66_wfs
!End of the abilint section

 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
&               vectsize,blockvectorr,vectsize,zero,gramxax,blocksize)
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
&               vectsize,gramxax,blocksize,zero,blockvectordumm,vectsize)
 blockvectorr=blockvectorr-blockvectordumm
!and now (b)orthornormalize r
!call operators(blockvectorr,blockvectorbr)
 if (gen_eigenpb) then
   ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
   ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor))
   do iblocksize=1,blocksize
     if (mpi_enreg%me_g0 == 1) then
       cwavef(1,2:npw_k*my_nspinor)=blockvectorr(2:npw_k*my_nspinor,iblocksize)/sq2
       cwavef(2,2:npw_k*my_nspinor)=blockvectorr(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)/sq2
       cwavef(1,1)=blockvectorr(1,iblocksize)
       cwavef(2,1)=zero
     else
       cwavef(1,1:npw_k*my_nspinor)=blockvectorr(1:npw_k*my_nspinor,iblocksize)/sq2
       cwavef(2,1:npw_k*my_nspinor)=blockvectorr(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)/sq2
     end if
  !   call to nonlop: compute <g|s|c>
     choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
     call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy2,&
  &               dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
  &               istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
  &               mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
  &               gs_hamk%nloalg,nnlout,npw_k,npw_k,my_nspinor,dtset%nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
  &               gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,signs,gs_hamk%sij,&
  &               gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef,use_gpu_cuda=dtset%use_gpu_cuda)
     if (mpi_enreg%me_g0 == 1) then
       blockvectorbr(2:npw_k*my_nspinor,iblocksize)=gwavef(1,2:npw_k*my_nspinor)*sq2
       blockvectorbr(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)=gwavef(2,2:npw_k*my_nspinor)*sq2
       blockvectorbr(1,iblocksize)=gwavef(1,1)
     else
       blockvectorbr(1:npw_k*my_nspinor,iblocksize)=gwavef(1,1:npw_k*my_nspinor)*sq2
       blockvectorbr(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)=gwavef(2,1:npw_k*my_nspinor)*sq2
     end if
   end do
   ABI_DEALLOCATE(cwavef)
   ABI_DEALLOCATE(gwavef)
 else
   blockvectorbr(:,:) = blockvectorr(:,:)
 end if

!call orthonormalize(blockvectorr,blockvectorbr)
 call orthonormalize(blockvectorr,blockvectorbr,blocksize,mpi_enreg,gramrbr,vectsize)
 call dtrsm('r','u','n','n',vectsize,blocksize,one,gramrbr,blocksize,&
&              blockvectorbr,vectsize)

!compute ar
!blockvectorar=matmul(operatora,blockvectorr)
!call operatorh(blockvectorr,blockvectorar)
 ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
 ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor))
 ABI_ALLOCATE(gvnlc,(2,npw_k*my_nspinor))

 do iblocksize=1,blocksize
   iband=iblocksize
   if (mpi_enreg%me_g0 == 1) then
     cwavef(1,2:npw_k*my_nspinor)=blockvectorr(2:npw_k*my_nspinor,iblocksize)/sq2
     cwavef(2,2:npw_k*my_nspinor)=blockvectorr(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)/sq2
     cwavef(2,1)=zero
     cwavef(1,1)=blockvectorr(1,iblocksize)
   else
     cwavef(1,1:npw_k*my_nspinor)=blockvectorr(1:npw_k*my_nspinor,iblocksize)/sq2
     cwavef(2,1:npw_k*my_nspinor)=blockvectorr(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)/sq2
   end if
   tim_getghc=7 ; sij_opt=0

   if(present(vxctaulocal))then
     call getghc(-1,cwavef,cprj_dum,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
   &  kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,my_nspinor,ntypat,&
   &  nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
   else
    call getghc(-1,cwavef,cprj_dum,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
   &  kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,my_nspinor,ntypat,&
   &  nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
   end if
   if (mpi_enreg%me_g0 == 1) then
     blockvectorar(2:npw_k*my_nspinor,iblocksize)=gwavef(1,2:npw_k*my_nspinor)*sq2
     blockvectorar(npw_k*my_nspinor+1:2*npw_k*my_nspinor-1,iblocksize)=gwavef(2,2:npw_k*my_nspinor)*sq2
     blockvectorar(1,iblocksize)=gwavef(1,1)
   else
     blockvectorar(1:npw_k*my_nspinor,iblocksize)=gwavef(1,:)*sq2
     blockvectorar(npw_k*my_nspinor+1:2*npw_k*my_nspinor,iblocksize)=gwavef(2,:)*sq2
   end if
 end do

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(gwavef)
 ABI_DEALLOCATE(gvnlc)

 end subroutine start_lobpcg

 end  subroutine lobpcgIIIwf
!!***
