!{\src2tex{textfont=tt}}
!!****f* ABINIT/pw_orthon
!! NAME
!! pw_orthon
!!
!! FUNCTION
!! Normalize nvec complex vectors each of length nelem and then
!! orthogonalize by modified Gram-Schmidt.
!! Two orthogonality conditions are available:
!!  Simple orthogonality: ${<Vec_{i}|Vec_{j}>=Delta_ij}$
!!  Orthogonality with overlap S: ${<Vec_{i}|S|Vec_{j}>=Delta_ij}$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, FF, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt.
!!
!! INPUTS
!!  icg=shift to be given to the location of the data in cg(=vecnm)
!!  igsc=shift to be given to the location of the data in gsc(=ovl_vecnm)
!!  istwf_k=option parameter that describes the storage of wfs
!!  mcg=maximum size of second dimension of cg(=vecnm)
!!  mgsc=maximum size of second dimension of gsc(=ovl_vecnm)
!!  mpi_enreg=informations about MPI parallelization
!!  nelem=number of complex elements in each vector
!!  nvec=number of vectors to be orthonormalized
!!  ortalgo= option for the choice of the algorithm
!!         -1: no orthogonalization (direct return)
!!          0: old algorithm (use of buffers)
!!          1: new algorithm (use of blas)
!!  useoverlap=select the orthogonality condition
!!               0: no overlap between vectors
!!               1: vectors are overlapping
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  vecnm= input: vectors to be orthonormalized; array of nvec column
!!                vectors, each of length nelem, shifted by icg
!!                This array is complex or else real(dp) of twice length
!!         output: orthonormalized set of vectors
!!  if (useoverlap==1) only:
!!    ovl_vecnm= input: product of overlap and input vectors:
!!                      S|vecnm>, where S is the overlap operator
!!               output: updated S|vecnm> according to vecnm
!!
!! NOTES
!! Note that each vector has an arbitrary phase which is not fixed in this
!! routine.
!!
!! TODO
!! WARNING : not yet suited for nspinor=2 with $\textrm{istwf}_k/=1$
!!
!! PARENTS
!!      vtowfk,wfconv
!!
!! CHILDREN
!!      dcopy,dtrsm,orthonormalize,timab,xcomm_init,xsum_mpi,zcopy
!!      zorthonormalize,ztrsm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pw_orthon(icg,igsc,istwf_k,mcg,mgsc,mpi_enreg,nelem,nvec,&
&                    ortalgo,ovl_vecnm,useoverlap,vecnm)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pw_orthon'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_abiutil
!End of the abilint section

 implicit none
!Following iis the hand tuned interface to lapack routines
!the tuned part are written in lower case
!interface
!  SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
!   INTEGER INCX,INCY,N
!   real*8 ZX(2,*)
!   complex*16 ZY(*)
!  END SUBROUTINE
! SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!  DOUBLE COMPLEX ALPHA
!  INTEGER LDA,LDB,M,N
!  CHARACTER DIAG,SIDE,TRANSA,UPLO
!  COMPLEX*16 A( LDA,*),B( LDB,*)
! end subroutine
!end interface
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,igsc,istwf_k,mcg,mgsc,nelem,nvec,ortalgo,useoverlap
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: ovl_vecnm(2,mgsc*useoverlap),vecnm(2,mcg)

!Local variables-------------------------------
!scalars
 integer :: cgindex,gscindex,ierr,ii,ii0,ii1,ii2,ivec,ivec2,old_paral_level
 integer :: rvectsiz,spaceComm=0,vectsize
 real(dp) :: doti,dotr,sq2,sum,xnorm
! character(len=500) :: message
!arrays
 real(dp) :: buffer2(2),tsec(2)
 real(dp),allocatable :: rblockvectorbx(:,:),rblockvectorx(:,:),rgramxbx(:,:)
 complex(dpc),allocatable :: cblockvectorbx(:,:),cblockvectorx(:,:)
 complex(dpc),allocatable :: cgramxbx(:,:)
!no_abirules
!Statement functions
 cgindex(ivec)=nelem*(ivec-1)+icg+1
 gscindex(ivec)=nelem*(ivec-1)+igsc+1

! *************************************************************************
!
!DEBUG
!write(std_out,*)' pw_orthon : enter, istwf_k= ',istwf_k,'nvec',nvec,'useoverlap',useoverlap
!if(istwf_k==2)then
!Make sure imaginary part at G=0 vanishes
!do ivec=1,nvec
!if(abs(vecnm(2,1+nelem*(ivec-1)+icg))>1.0d-10)then
!write(message,'(a,a,a,a,3i5,2es16.6,a,a)')ch10,&
!&     ' pw_orthon : BUG ',&
!&     '  For istwf_k=2, observed the following element of vecnm :',ch10,&
!&     nelem,ivec,icg,vecnm(1:2,1+nelem*(ivec-1)+icg),ch10,&
!&     '  with a non-negligible imaginary part.'
!call wrtout(std_out,message,'PERS')
!call leave_new('PERS')
!end if
!end do
!end if
!ENDDEBUG

!Nothing to do if ortalgo=-1
 if(ortalgo==-1) return

 if(ortalgo==1) then
!  =======================
!  First (new) algorithm
!  =======================
!  This first algorithm seems to be more efficient especially in the parallel band-FFT mode.

   if(istwf_k==1) then
!    write(std_out,*)'zorthonormalize WFs'
     vectsize=nelem
     ABI_ALLOCATE(cgramxbx,(nvec,nvec))
     ABI_ALLOCATE(cblockvectorx,(vectsize,nvec))
     ABI_ALLOCATE(cblockvectorbx,(vectsize,nvec))
     call zcopy(nvec*vectsize,     vecnm(:, cgindex(1): cgindex(nvec+1)-1),1,cblockvectorx,1)
     if (useoverlap == 1) then
       call zcopy(nvec*vectsize,ovl_vecnm(:,gscindex(1):gscindex(nvec+1)-1),1,cblockvectorbx,1)
     else
       call zcopy(nvec*vectsize,    vecnm(:, cgindex(1): cgindex(nvec+1)-1),1,cblockvectorbx,1)
     end if
     call zorthonormalize(cblockvectorx,cblockvectorbx,nvec,mpi_enreg,cgramxbx,vectsize)
     call zcopy(nvec*vectsize,  cblockvectorx,1,    vecnm(:, cgindex(1): cgindex(nvec+1)-1),1)
     if (useoverlap == 1) then
       call ztrsm('r','u','n','n',vectsize,nvec,cone,cgramxbx,nvec,cblockvectorbx,vectsize)
       call zcopy(nvec*vectsize,cblockvectorbx,1,ovl_vecnm(:,gscindex(1):gscindex(nvec+1)-1),1)
     end if
     ABI_DEALLOCATE(cgramxbx)
     ABI_DEALLOCATE(cblockvectorx)
     ABI_DEALLOCATE(cblockvectorbx)
   else if(istwf_k==2) then
     write(std_out,*)'orthonormalize WFs'
     sq2=sqrt(2.0_dp)
     if (mpi_enreg%me_g0 == 1) then
       vectsize=2*nelem-1
     else
       vectsize=2*nelem
     end if
     ABI_ALLOCATE(rgramxbx,(nvec,nvec))
     ABI_ALLOCATE(rblockvectorx,(vectsize,nvec))
     ABI_ALLOCATE(rblockvectorbx,(vectsize,nvec))
     rvectsiz=nelem
     do ivec=1,nvec
       if (mpi_enreg%me_g0 == 1) then
         call dcopy(1         ,     vecnm(1, cgindex(ivec))                         ,1,rblockvectorx (1                  ,ivec),1)
         call dcopy(rvectsiz-1,     vecnm(1, cgindex(ivec)+1: cgindex(ivec+1)-1)*sq2,1,rblockvectorx (2:rvectsiz         ,ivec),1)
         call dcopy(rvectsiz-1,     vecnm(2, cgindex(ivec)+1: cgindex(ivec+1)-1)*sq2,1,rblockvectorx (rvectsiz+1:vectsize,ivec),1)
         if (useoverlap == 1) then
           call dcopy(1         ,ovl_vecnm(1,gscindex(ivec))                         ,1,rblockvectorbx(1                  ,ivec),1)
           call dcopy(rvectsiz-1,ovl_vecnm(1,gscindex(ivec)+1:gscindex(ivec+1)-1)*sq2,1,rblockvectorbx(2:rvectsiz         ,ivec),1)
           call dcopy(rvectsiz-1,ovl_vecnm(2,gscindex(ivec)+1:gscindex(ivec+1)-1)*sq2,1,rblockvectorbx(rvectsiz+1:vectsize,ivec),1)
         else
           call dcopy(1         ,    vecnm(1, cgindex(ivec))                         ,1,rblockvectorbx(1                  ,ivec),1)
           call dcopy(rvectsiz-1,    vecnm(1, cgindex(ivec)+1: cgindex(ivec+1)-1)*sq2,1,rblockvectorbx(2:rvectsiz         ,ivec),1)
           call dcopy(rvectsiz-1,    vecnm(2, cgindex(ivec)+1: cgindex(ivec+1)-1)*sq2,1,rblockvectorbx(rvectsiz+1:vectsize,ivec),1)
         end if
       else
         call dcopy(rvectsiz,       vecnm(1, cgindex(ivec):   cgindex(ivec+1)-1)*sq2,1,rblockvectorx (1:rvectsiz         ,ivec),1)
         call dcopy(rvectsiz,       vecnm(2, cgindex(ivec):   cgindex(ivec+1)-1)*sq2,1,rblockvectorx (rvectsiz+1:vectsize,ivec),1)
         if (useoverlap == 1) then
           call dcopy(rvectsiz,  ovl_vecnm(1,gscindex(ivec):  gscindex(ivec+1)-1)*sq2,1,rblockvectorbx(1:rvectsiz         ,ivec),1)
           call dcopy(rvectsiz,  ovl_vecnm(2,gscindex(ivec):  gscindex(ivec+1)-1)*sq2,1,rblockvectorbx(rvectsiz+1:vectsize,ivec),1)
         else
           call dcopy(rvectsiz,      vecnm(1,cgindex(ivec):    cgindex(ivec+1)-1)*sq2,1,rblockvectorbx(1:rvectsiz         ,ivec),1)
           call dcopy(rvectsiz,      vecnm(2,cgindex(ivec):    cgindex(ivec+1)-1)*sq2,1,rblockvectorbx(rvectsiz+1:vectsize,ivec),1)
         end if
       end if
     end do
     call orthonormalize(rblockvectorx,rblockvectorbx,nvec,mpi_enreg,rgramxbx,vectsize)
     do ivec=1,nvec
       if (mpi_enreg%me_g0 == 1) then
         call dcopy(1           ,rblockvectorx(1,ivec)                      ,1,vecnm(1,cgindex(ivec))                    ,1)
!        call dcopy(1           ,zero                                       ,1,vecnm(2,cgindex(ivec))                    ,1)
         vecnm(2,cgindex(ivec))=zero
         call dcopy(rvectsiz-1,rblockvectorx(2:rvectsiz         ,ivec)/sq2,1,vecnm(1,cgindex(ivec)+1:cgindex(ivec+1)-1),1)
         call dcopy(rvectsiz-1,rblockvectorx(rvectsiz+1:vectsize,ivec)/sq2,1,vecnm(2,cgindex(ivec)+1:cgindex(ivec+1)-1),1)
       else
         call dcopy(rvectsiz,  rblockvectorx(1:rvectsiz         ,ivec)/sq2,1,vecnm(1,cgindex(ivec)  :cgindex(ivec+1)-1),1)
         call dcopy(rvectsiz,  rblockvectorx(rvectsiz+1:vectsize,ivec)/sq2,1,vecnm(2,cgindex(ivec)  :cgindex(ivec+1)-1),1)
       end if

       if(useoverlap == 1) then
         call dtrsm('r','u','n','n',vectsize,nvec,one,rgramxbx,nvec,rblockvectorbx,vectsize)
         if (mpi_enreg%me_g0 == 1) then
           call dcopy(1         ,rblockvectorbx(1,ivec)                      ,1,ovl_vecnm(1,gscindex(ivec))                     ,1)
!          call dcopy(1         ,zero                                        ,1,ovl_vecnm(2,gscindex(ivec))                     ,1)
           ovl_vecnm(2,gscindex(ivec))=zero
           call dcopy(rvectsiz-1,rblockvectorbx(2:rvectsiz         ,ivec)/sq2,1,ovl_vecnm(1,gscindex(ivec)+1:gscindex(ivec+1)-1),1)
           call dcopy(rvectsiz-1,rblockvectorbx(rvectsiz+1:vectsize,ivec)/sq2,1,ovl_vecnm(2,gscindex(ivec)+1:gscindex(ivec+1)-1),1)
         else
           call dcopy(rvectsiz  ,rblockvectorbx(1:rvectsiz         ,ivec)/sq2,1,ovl_vecnm(1,gscindex(ivec)  :gscindex(ivec+1)-1),1)
           call dcopy(rvectsiz  ,rblockvectorbx(rvectsiz+1:vectsize,ivec)/sq2,1,ovl_vecnm(2,gscindex(ivec)  :gscindex(ivec+1)-1),1)
         end if
       end if
     end do
     ABI_DEALLOCATE(rgramxbx)
     ABI_DEALLOCATE(rblockvectorx)
     ABI_DEALLOCATE(rblockvectorbx)
   end if

 else !End of the first algo

!  =======================
!  Second (old) algorithm
!  =======================

   do ivec=1,nvec
!    Normalize each vecnm(n,m) in turn:
     if (useoverlap==1) then

!      Using overlap S...
       if(istwf_k/=2)then
         sum=zero;ii0=1
       else
         if (mpi_enreg%me_g0 ==1) then
           sum=half*ovl_vecnm(1,1+nelem*(ivec-1)+igsc)&
&           *vecnm(1,1+nelem*(ivec-1)+icg)
           ii0=2
         else
           sum=zero;ii0=1
         end if
       end if
!      $OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:sum) &
!      $OMP&SHARED(icg,ivec,nelem,vecnm)
       do ii=ii0+nelem*(ivec-1),nelem*ivec
         sum=sum+vecnm(1,ii+icg)*ovl_vecnm(1,ii+igsc)+vecnm(2,ii+icg)*ovl_vecnm(2,ii+igsc)
       end do
!      $OMP END PARALLEL DO
     else

!      Without overlap...
       if(istwf_k/=2)then
         sum=zero;ii0=1
       else
         if (mpi_enreg%me_g0 ==1) then
           sum=half*vecnm(1,1+nelem*(ivec-1)+icg)**2
           ii0=2
         else
           sum=zero;ii0=1
         end if
       end if
!      $OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:sum) &
!      $OMP&SHARED(icg,ivec,nelem,vecnm)
       do ii=ii0+nelem*(ivec-1)+icg,nelem*ivec+icg
         sum=sum+vecnm(1,ii)**2+vecnm(2,ii)**2
       end do
!      $OMP END PARALLEL DO
     end if

!    Init mpi_comm
!    if (mpi_enreg%paral_compil_kpt/=1) then
     old_paral_level=mpi_enreg%paral_level
     mpi_enreg%paral_level=3
     call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_3d)
     call timab(48,1,tsec)
     call xsum_mpi(sum,spaceComm ,ierr)
     call timab(48,2,tsec)
     mpi_enreg%paral_level=old_paral_level
!    end if

     if(istwf_k>=2)sum=two*sum
     xnorm = sqrt(abs(sum)) ;  sum=1.0_dp/xnorm
!    $OMP PARALLEL DO PRIVATE(ii) &
!    $OMP&SHARED(icg,ivec,nelem,sum,vecnm)
     do ii=1+nelem*(ivec-1)+icg,nelem*ivec+icg
       vecnm(1,ii)=vecnm(1,ii)*sum
       vecnm(2,ii)=vecnm(2,ii)*sum
     end do
!    $OMP END PARALLEL DO
     if (useoverlap==1) then
!      $OMP PARALLEL DO PRIVATE(ii) &
!      $OMP&SHARED(icg,ivec,nelem,sum,ovl_vecnm)
       do ii=1+nelem*(ivec-1)+igsc,nelem*ivec+igsc
         ovl_vecnm(1,ii)=ovl_vecnm(1,ii)*sum
         ovl_vecnm(2,ii)=ovl_vecnm(2,ii)*sum
       end do
!      $OMP END PARALLEL DO
     end if

!    Remove projection in all higher states.
     if (ivec<nvec) then

       if(istwf_k==1)then
!        Cannot use time-reversal symmetry

         if (useoverlap==1) then
!          ----- Using overlap -----
           do ivec2=ivec+1,nvec
!            First compute scalar product
             dotr=zero ; doti=zero
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+igsc
!            $OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:doti,dotr) &
!            $OMP&SHARED(ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               dotr=dotr+vecnm(1,ii1+ii)*ovl_vecnm(1,ii2+ii)+&
&               vecnm(2,ii1+ii)*ovl_vecnm(2,ii2+ii)
               doti=doti+vecnm(1,ii1+ii)*ovl_vecnm(2,ii2+ii)-&
&               vecnm(2,ii1+ii)*ovl_vecnm(1,ii2+ii)
             end do
!            $OMP END PARALLEL DO

!            Init mpi_comm
!            if (mpi_enreg%paral_compil_kpt/=1) then
             old_paral_level=mpi_enreg%paral_level
             mpi_enreg%paral_level=3
             call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_3d)
             buffer2(1)=doti
             buffer2(2)=dotr
             call timab(48,1,tsec)
             call xsum_mpi(buffer2,spaceComm ,ierr)
!            call xsum_mpi(doti,spaceComm ,ierr)
!            call xsum_mpi(dotr,spaceComm ,ierr)
             call timab(48,2,tsec)
             doti=buffer2(1)
             dotr=buffer2(2)
             mpi_enreg%paral_level=old_paral_level
!            end if

!            Then subtract the appropriate amount of the lower state
!            $OMP PARALLEL DO PRIVATE(ii) &
!            $OMP&SHARED(doti,dotr,ii1,ii2,nelem,vecnm)
#ifdef FC_INTEL
!            DIR$ ivdep
#endif
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)+&
&               doti*vecnm(2,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-doti*vecnm(1,ii1+ii)-&
&               dotr*vecnm(2,ii1+ii)
             end do
             ii1=nelem*(ivec-1)+igsc;ii2=nelem*(ivec2-1)+igsc
             do ii=1,nelem
               ovl_vecnm(1,ii2+ii)=ovl_vecnm(1,ii2+ii)&
&               -dotr*ovl_vecnm(1,ii1+ii)&
&               +doti*ovl_vecnm(2,ii1+ii)
               ovl_vecnm(2,ii2+ii)=ovl_vecnm(2,ii2+ii)&
               -doti*ovl_vecnm(1,ii1+ii)&
&               -dotr*ovl_vecnm(2,ii1+ii)
             end do
!            $OMP END PARALLEL DO
           end do
         else
!          ----- No overlap -----
           do ivec2=ivec+1,nvec
!            First compute scalar product
             dotr=zero ; doti=zero
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
!            $OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:doti,dotr) &
!            $OMP&SHARED(ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               dotr=dotr+vecnm(1,ii1+ii)*vecnm(1,ii2+ii)+&
&               vecnm(2,ii1+ii)*vecnm(2,ii2+ii)
               doti=doti+vecnm(1,ii1+ii)*vecnm(2,ii2+ii)-&
&               vecnm(2,ii1+ii)*vecnm(1,ii2+ii)
             end do
!            $OMP END PARALLEL DO
!            Init mpi_comm
!            if (mpi_enreg%paral_compil_kpt/=1) then
             old_paral_level=mpi_enreg%paral_level
             mpi_enreg%paral_level=3
             call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_3d)
             buffer2(1)=doti
             buffer2(2)=dotr
             call timab(48,1,tsec)
             call xsum_mpi(buffer2,spaceComm ,ierr)
!            call xsum_mpi(doti,spaceComm ,ierr)
!            call xsum_mpi(dotr,spaceComm ,ierr)
             call timab(48,2,tsec)
             doti=buffer2(1)
             dotr=buffer2(2)
             mpi_enreg%paral_level=old_paral_level
!            end if

!            Then subtract the appropriate amount of the lower state
!            $OMP PARALLEL DO PRIVATE(ii) &
!            $OMP&SHARED(doti,dotr,ii1,ii2,nelem,vecnm)
#ifdef FC_INTEL
!            DIR$ ivdep
#endif
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)+&
&               doti*vecnm(2,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-doti*vecnm(1,ii1+ii)-&
&               dotr*vecnm(2,ii1+ii)
             end do
!            $OMP END PARALLEL DO
           end do

         end if  ! Test on useoverlap

       else if(istwf_k==2)then
!        At gamma point use of time-reversal symmetry saves cpu time.

         if (useoverlap==1) then
!          ----- Using overlap -----
           do ivec2=ivec+1,nvec
!            First compute scalar product
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+igsc
             if (mpi_enreg%me_g0 ==1) then
               dotr=half*vecnm(1,ii1+1)*ovl_vecnm(1,ii2+1)
!              Avoid double counting G=0 contribution
!              Imaginary part of vecnm at G=0 should be zero, so only take real part
!              $OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) &
!              $OMP&SHARED(ii1,ii2,nelem,vecnm)
               do ii=2,nelem
                 dotr=dotr+vecnm(1,ii1+ii)*ovl_vecnm(1,ii2+ii)+&
&                 vecnm(2,ii1+ii)*ovl_vecnm(2,ii2+ii)
               end do
!              $OMP END PARALLEL DO
             else
               dotr=0._dp
!              $OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) &
!              $OMP&SHARED(ii1,ii2,nelem,vecnm)
               do ii=1,nelem
                 dotr=dotr+vecnm(1,ii1+ii)*ovl_vecnm(1,ii2+ii)+&
&                 vecnm(2,ii1+ii)*ovl_vecnm(2,ii2+ii)
               end do
!              $OMP END PARALLEL DO
             end if

             dotr=two*dotr

!            Init mpi_comm
!            if (mpi_enreg%paral_compil_kpt/=1) then
             old_paral_level=mpi_enreg%paral_level
             mpi_enreg%paral_level=3
             call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_3d)
             call timab(48,1,tsec)
             call xsum_mpi(dotr,spaceComm ,ierr)
             call timab(48,2,tsec)
             mpi_enreg%paral_level=old_paral_level
!            end if

!            Then subtract the appropriate amount of the lower state
!            $OMP PARALLEL DO PRIVATE(ii) &
!            $OMP&SHARED(dotr,ii1,ii2,nelem,vecnm)
#ifdef FC_INTEL
!            DIR$ ivdep
#endif
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-dotr*vecnm(2,ii1+ii)
             end do
             ii1=nelem*(ivec-1)+igsc;ii2=nelem*(ivec2-1)+igsc
             do ii=1,nelem
               ovl_vecnm(1,ii2+ii)=ovl_vecnm(1,ii2+ii)-dotr*ovl_vecnm(1,ii1+ii)
               ovl_vecnm(2,ii2+ii)=ovl_vecnm(2,ii2+ii)-dotr*ovl_vecnm(2,ii1+ii)
             end do
!            $OMP END PARALLEL DO
           end do
         else
!          ----- No overlap -----
           do ivec2=ivec+1,nvec
!            First compute scalar product
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
             if (mpi_enreg%me_g0 ==1) then
!              Avoid double counting G=0 contribution
!              Imaginary part of vecnm at G=0 should be zero, so only take real part
               dotr=half*vecnm(1,ii1+1)*vecnm(1,ii2+1)
!              $OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) &
!              $OMP&SHARED(ii1,ii2,nelem,vecnm)
               do ii=2,nelem
                 dotr=dotr+vecnm(1,ii1+ii)*vecnm(1,ii2+ii)+&
&                 vecnm(2,ii1+ii)*vecnm(2,ii2+ii)
               end do
!              $OMP END PARALLEL DO
             else
               dotr=0._dp
!              $OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) &
!              $OMP&SHARED(ii1,ii2,nelem,vecnm)
               do ii=1,nelem
                 dotr=dotr+vecnm(1,ii1+ii)*vecnm(1,ii2+ii)+&
&                 vecnm(2,ii1+ii)*vecnm(2,ii2+ii)
               end do
!              $OMP END PARALLEL DO
             end if
             dotr=two*dotr

!            Init mpi_comm
!            if (mpi_enreg%paral_compil_kpt/=1) then
             old_paral_level=mpi_enreg%paral_level
             mpi_enreg%paral_level=3
             call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_3d)
             call timab(48,1,tsec)
             call xsum_mpi(dotr,spaceComm ,ierr)
             call timab(48,2,tsec)
             mpi_enreg%paral_level=old_paral_level
!            end if

!            Then subtract the appropriate amount of the lower state
!            $OMP PARALLEL DO PRIVATE(ii) &
!            $OMP&SHARED(dotr,ii1,ii2,nelem,vecnm)
#ifdef FC_INTEL
!            DIR$ ivdep
#endif
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-dotr*vecnm(2,ii1+ii)
             end do
!            $OMP END PARALLEL DO
           end do
         end if  ! Test on useoverlap

       else
!        At other special points, use of time-reversal symmetry saves cpu time.

         if (useoverlap==1) then
!          ----- Using overlap -----
           do ivec2=ivec+1,nvec
!            First compute scalar product
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+igsc
!            Avoid double counting G=0 contribution
!            Imaginary part of vecnm at G=0 should be zero, so only take real part
             dotr=zero
!            $OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) &
!            $OMP&SHARED(ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               dotr=dotr+vecnm(1,ii1+ii)*ovl_vecnm(1,ii2+ii)+&
&               vecnm(2,ii1+ii)*ovl_vecnm(2,ii2+ii)
             end do
!            $OMP END PARALLEL DO
             dotr=two*dotr

!            Init mpi_comm
!            if (mpi_enreg%paral_compil_kpt/=1) then
             old_paral_level=mpi_enreg%paral_level
             mpi_enreg%paral_level=3
             call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_3d)
             call timab(48,1,tsec)
             call xsum_mpi(dotr,spaceComm ,ierr)
             call timab(48,2,tsec)
             mpi_enreg%paral_level=old_paral_level
!            end if

!            Then subtract the appropriate amount of the lower state
!            $OMP PARALLEL DO PRIVATE(ii) &
!            $OMP&SHARED(dotr,ii1,ii2,nelem,vecnm)
#ifdef FC_INTEL
!            DIR$ ivdep
#endif
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-dotr*vecnm(2,ii1+ii)
             end do
             ii1=nelem*(ivec-1)+igsc;ii2=nelem*(ivec2-1)+igsc
             do ii=1,nelem
               ovl_vecnm(1,ii2+ii)=ovl_vecnm(1,ii2+ii)-dotr*ovl_vecnm(1,ii1+ii)
               ovl_vecnm(2,ii2+ii)=ovl_vecnm(2,ii2+ii)-dotr*ovl_vecnm(2,ii1+ii)
             end do
!            $OMP END PARALLEL DO
           end do
         else
!          ----- No overlap -----
           do ivec2=ivec+1,nvec
!            First compute scalar product
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
!            Avoid double counting G=0 contribution
!            Imaginary part of vecnm at G=0 should be zero, so only take real part
             dotr=zero
!            $OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) &
!            $OMP&SHARED(ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               dotr=dotr+vecnm(1,ii1+ii)*vecnm(1,ii2+ii)+&
&               vecnm(2,ii1+ii)*vecnm(2,ii2+ii)
             end do
!            $OMP END PARALLEL DO
             dotr=two*dotr

!            Init mpi_comm
!            if (mpi_enreg%paral_compil_kpt/=1) then
             old_paral_level=mpi_enreg%paral_level
             mpi_enreg%paral_level=3
             call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_3d)
             call timab(48,1,tsec)
             call xsum_mpi(dotr,spaceComm ,ierr)
             call timab(48,2,tsec)
             mpi_enreg%paral_level=old_paral_level
!            end if

!            Then subtract the appropriate amount of the lower state
!            $OMP PARALLEL DO PRIVATE(ii) &
!            $OMP&SHARED(dotr,ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-dotr*vecnm(2,ii1+ii)
             end do
!            $OMP END PARALLEL DO
           end do
         end if

!        End use of time-reversal symmetry
       end if

     end if  ! Test on "ivec"

!    end loop over vectors (or bands) with index ivec :
   end do

 end if ! End of the second algorithm

end subroutine pw_orthon
!!***
