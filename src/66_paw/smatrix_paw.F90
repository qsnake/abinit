!{\src2tex{textfont=tt}}
!!****f* ABINIT/smatrix_paw
!! NAME
!! smatrix_paw
!!
!! FUNCTION
!! Compute the overlap matrix between the k-points k and k + dk.
!! Depending on the value of job and ddkflag, compute also its determinant,
!! its inverse and the product of its inverse with the wavefunctions at k.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT  group (MVeithen,BAmadon,FJollet,PHermet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cg(2,mcg_k) = planewave coefficients of wavefunctions at k
!! cgq(2,mcg_q) = planewave coefficients of wavefunctions at q = k + dk
!! ddkflag = 1 : compute product of the inverse overlap matrix
!!               with the wavefunction at k (job = 1 or 11)
!!           0 : do not compute the product of the inverse overlap matrix
!!               with the wavefunction at k
!! icg = shift applied to the wavefunctions in the array cg
!! icg1 = shift applied to the wavefunctions in the array cgq
!! itrs = variable that governs the use of time-reversal symmetry
!!        when converting the wavefunctions from the iBZ to the fBZ
!! job = argument of zgedi
!!        0 : update overlap matrix
!!        1 : compute inverse of the overlap matrix
!!       10 : compute determinant of the overlap matrix
!!       11 : compute determinant and inverse of the overlap matrix
!! maxbd = used in case ddkflag = 1, defines the highest band for
!!         which the ddk will be computed
!! mcg_k = second dimension of cg
!! mcg_q = second dimension of cg_q
!! mcg1_k = second dimension of cg1_k, should be equal to
!!          mpw*nsppol*nspinor*(maxbd - minbd + 1)
!! minbd = used in case ddkflag = 1, defines the lowest band for
!!         which the ddk will be computed
!! mpw = maximum dimensioned size of npw
!! nband_occ = number of (occupied) valence bands
!! npw_k1 = number of plane waves at k
!! npw_k2 = number of plane waves at k + dk
!! nspinor =number of spinorial components of the wavefunctions (on current proc)
!! pwind_k = array used to compute the overlap matrix
!! pwnsfac = phase factors for non-symmorphic translations
!! shifrbd = shift applied to the location of the WF in the cg-array
!!           after each loop over bands
!!           0 : apply no shift, this is allowed in case cg
!!               contains only the wf of one single band
!!           1 : apply a shift of npw_k1*nspinor, this is the usual option
!!               when cg contains the wf for all occupied bands
!!
!! OUTPUT
!! cg1_k(2,mcg1_k) = product of the inverse overlap matrix with the
!!                   wavefunctions at k; computed in case job = 1 or 11
!! dtm_k(2) = determinant of the overlap matrix between k and k + dk;
!!            computed in case job = 10 or 11
!! smat_inv = inverse of the overlap matrix
!!
!! SIDE EFFECTS
!! Input/Output
!! sflag_k(iband) = 1 if the elements smat_k(:,iband,:) are up to date
!!                    -> they will not be recomputed
!!                  0 the elements smat_k(:,iband,:) will be recomputed
!!      at the end of the routine, sflag_k(1:nband_occ) = 1
!!      (the whole overlap matrix is up to date)
!! smat_k = overlap matrix between k, k + dk
!!          only the lines for which sflag_k = 0 are computed
!!          smat_k(:,n,m) = < u_{n,k} | u_{m,k+dk} >
!!
!! NOTES
!! This routine is quite flexible in the way it deals with the wavefunctions:
!!  - cg (WF at k) can contain either the whole WF array (all k-points
!!    and bands), in which case the location of the WF at k is specified
!!    by icg and shiftbd = 1, or the WF of a single k-point/band, in which case
!!    shiftbd = 0 and icg = 0.
!!  - cgq (WF at k + dk) can contain either the whole WF array (all k-points
!!    and bands), in which case the location of the WF at k is specified
!!    by icg1, or the WF of a single k-point, in which case
!!    icg1 = 0. cgq must contain the WF of ALL occupied bands.
!!  - cg1_k can either be computed for all valence bands or
!!    for a group of valence bands defined by minbd and maxbd.
!!
!! PARENTS
!!
!! CHILDREN
!!      dzgedi,dzgefa,leave_new,overlap_g,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine smatrix_paw(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,job,maxbd,&
&  mcg_k,mcg_q,mcg1_k,minbd,mpw,nband_occ,npw_k1,npw_k2,nspinor,&
&  pwind_k,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_k)

 use m_profiling

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'smatrix_paw'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_28_numeric_noabirule
 use interfaces_53_spacepar
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ddkflag,icg,icg1,itrs,job,maxbd,mcg1_k,mcg_k,mcg_q
 integer,intent(in) :: minbd,mpw,nband_occ,npw_k1,npw_k2,nspinor,shiftbd
!arrays
 integer,intent(in) :: pwind_k(mpw)
 integer,intent(inout) :: sflag_k(nband_occ)
 real(dp),intent(in) :: cg(2,mcg_k),cgq(2,mcg_q),pwnsfac_k(4,mpw)
 real(dp),intent(inout) :: smat_k(2,nband_occ,nband_occ)
 real(dp),intent(out) :: cg1_k(2,mcg1_k),dtm_k(2)
 real(dp),intent(out) :: smat_inv(2,nband_occ,nband_occ)

!Local variables -------------------------
!scalars
 integer :: iband,icount,info,ipw,jband,jband1,jpw,pwmax,pwmin
 real(dp) :: doti,dotr,fac,wfi,wfr
 character(len=500) :: message
!arrays
 integer,allocatable :: ipvt(:)
 real(dp) :: det(2,2)
 real(dp),allocatable :: vect1(:,:),vect2(:,:),zgwork(:,:)

! ***********************************************************************

 DBG_ENTER("COLL")

 ABI_ALLOCATE(ipvt,(nband_occ))
 ABI_ALLOCATE(zgwork,(2,nband_occ))
 ABI_ALLOCATE(vect1,(2,0:mpw))
 ABI_ALLOCATE(vect2,(2,0:mpw))
 vect1(:,0) = zero ; vect2(:,0) = zero

!Check if the values of ddkflag and job are compatible

 if ((job /= 0).and.(job /= 1).and.(job /= 10).and.(job /= 11)) then
   write(message,'(a,a,a,a,i3,a,a)')ch10,&
&   ' smatrix_paw: BUG - ',ch10,&
&   '   job is equal to ',job,ch10,&
&   '   while only the values job = 0, 1, 10 or 11 are allowed.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

 if (ddkflag == 1) then
   if ((job == 0).or.(job == 10)) then
     write(message,'(a,a,a,a,i3,a,a)')ch10,&
&     ' smatrix_paw: BUG - ',ch10,&
&     '   job is equal to ',job,ch10,&
&     '   while ddkflag = 1. This is not allowed.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
 end if

!Check the values of sflag_k
 do iband=1,nband_occ
   if (sflag_k(iband)/=0 .and. sflag_k(iband)/=1)then
     write(message,'(6a,i4,a,i4)')ch10,&
&     ' smatrix_paw: BUG - ',ch10,&
&     '  The content of sflag_k must be 0 or 1.',ch10,&
&     '  However, for iband=',iband,', sflag_k(iband)=',sflag_k(iband)
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
 end do

!Check if shiftbd is consistent with sflag_k
 if (shiftbd == 0) then
   icount = 0
   do iband = 1, nband_occ
     if (sflag_k(iband) == 0) icount = icount + 1
   end do
   if (icount > 1) then
     write(message,'(a,a,a,a)')ch10,&
&     ' smatrix_paw: BUG - ',ch10,&
&     '   in case shiftbd = 0, only 1 element of sflag can be 0'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
 end if

!Update the lines of the overlap matrix for which sflag = 0
!MVeithen: because of sflag, it is more efficient to perform
!the loop over jband inside the loop over iband

!DEBUG
!write(std_out,*)' smatrix : smat_k(1,1,1)=',smat_k(1,1,1)
!write(std_out,*)' smatrix : sflag_k=',sflag_k
!ENDDEBUG

 do iband = 1, nband_occ

   if (sflag_k(iband) == 0) then

     pwmin = (iband-1)*npw_k1*nspinor*shiftbd
     pwmax = pwmin + npw_k1*nspinor
!    
!    old version  (*** multiply by nspinor missing??? ***)
!    vect1(:,1:npw_k1) = cg(:,icg + 1 + pwmin:icg + pwmax)
!    

!    Multiply the bra wave function by the phase factor
     if (itrs==1.or.itrs==11) then  ! take complex conjugate of bra
       do ipw = 1,npw_k1
         vect1(1,ipw) = cg(1,icg+ipw+pwmin)*pwnsfac_k(1,ipw) &
&         +cg(2,icg+ipw+pwmin)*pwnsfac_k(2,ipw)
         vect1(2,ipw) = cg(1,icg+ipw+pwmin)*pwnsfac_k(2,ipw) &
&         -cg(2,icg+ipw+pwmin)*pwnsfac_k(1,ipw)
       end do
     else
       do ipw = 1,npw_k1
         vect1(1,ipw) = cg(1,icg+ipw+pwmin)*pwnsfac_k(1,ipw) &
&         -cg(2,icg+ipw+pwmin)*pwnsfac_k(2,ipw)
         vect1(2,ipw) = cg(1,icg+ipw+pwmin)*pwnsfac_k(2,ipw) &
&         +cg(2,icg+ipw+pwmin)*pwnsfac_k(1,ipw)
       end do
     end if

     if (npw_k1 < mpw) vect1(:,npw_k1+1:mpw) = zero

     do jband = 1, nband_occ

       pwmin = (jband-1)*npw_k2*nspinor
       pwmax = pwmin + npw_k2*nspinor

       if (itrs==10.or.itrs==11) then ! take complex conjugate of ket
         do ipw = 1, npw_k2
           vect2(1,ipw) = cgq(1,icg1+ipw+pwmin)*pwnsfac_k(3,ipw) &
&           +cgq(2,icg1+ipw+pwmin)*pwnsfac_k(4,ipw)
           vect2(2,ipw) = cgq(1,icg1+ipw+pwmin)*pwnsfac_k(4,ipw) &
&           -cgq(2,icg1+ipw+pwmin)*pwnsfac_k(3,ipw)
         end do
       else
         do ipw = 1, npw_k2
           vect2(1,ipw) = cgq(1,icg1+ipw+pwmin)*pwnsfac_k(3,ipw) &
&           -cgq(2,icg1+ipw+pwmin)*pwnsfac_k(4,ipw)
           vect2(2,ipw) = cgq(1,icg1+ipw+pwmin)*pwnsfac_k(4,ipw) &
&           +cgq(2,icg1+ipw+pwmin)*pwnsfac_k(3,ipw)
         end do
       end if

       if (npw_k2 < mpw) vect2(:,npw_k2+1:mpw) = zero

!      DEBUG
!      if(iband==1 .and. jband==1 .and. itrs==0 .and. npw_k1==68 .and. npw_k2==74)then
!      write(std_out,'(a)' )' smatrix : ipw,vect1,cg,pwnsfac='
!      do ipw=1,npw_k1
!      write(std_out,'(i4,6es16.6)' )ipw,vect1(:,ipw),cg(:,icg+ipw+pwmin),pwnsfac_k(1:2,ipw)
!      end do
!      do ipw=1,npw_k2
!      write(std_out,'(i4,6es16.6)' )ipw,vect2(:,ipw),cg(:,icg1+ipw+pwmin),pwnsfac_k(3:4,ipw)
!      end do
!      end if
!      ENDDEBUG

       call overlap_g(doti,dotr,mpw,npw_k1,pwind_k,vect1,vect2)

!      smat is incremented in order to include the PAW contribution when
!      called from berryphase_new
!      if this routine is called from cgwf, no modification is expected
!      write(222,*) iband,jband,smat_k(1,iband,jband),smat_k(2,iband,jband)
       smat_k(1,iband,jband) = smat_k(1,iband,jband)+ dotr
       smat_k(2,iband,jband) = smat_k(2,iband,jband)+ doti
!      write(333,*) iband,jband,smat_k(1,iband,jband),smat_k(2,iband,jband)
!      
!      

!      DEBUG
!      if(iband==1 .and. jband==1)then
!      write(std_out,'(a,2es16.6,3i4)' )' smatrix : dotr,smat_k(1,iband,jband),mpw,npw_k1,npw_k2',dotr,smat_k(1,iband,jband),mpw,npw_k1,npw_k2
!      end if
!      ENDDEBUG

     end do   ! jband

   end if    ! sflag_k(iband) == 0

 end do   ! iband

!DEBUG
!do iband=1,nband_occ
!do jband=1,nband_occ
!write(std_out,'(a,2i4,2e20.10)') 'smat',iband,jband,smat_k(1,iband,jband),smat_k(2,iband,jband)
!end do
!end do
!write(std_out,*)' smatrix : smat_k(1,1,1)=',smat_k(1,1,1)
!ENDDEBUG

!Update sflag_k
 sflag_k(:) = 1

!Depending on the value of job, compute the determinant of the
!overlap matrix, its inverse or the product of the inverse
!overlap matrix with the WF at k.

 if (job > 0) then

   smat_inv(:,:,:) = smat_k(:,:,:)

!  DEBUG
!  write(std_out,*)' smatrix : smat_inv=',smat_inv
!  ENDDEBUG

   call dzgefa(smat_inv,nband_occ,nband_occ,ipvt,info)
   call dzgedi(smat_inv,nband_occ,nband_occ,ipvt,det,zgwork,job)

!  DEBUG
!  write(std_out,*)' smatrix : det=',det
!  ENDDEBUG

!  Compute the determinant of the overlap matrix
   dtm_k(:) = zero
   if (job > 1) then
     fac = exp(log(10._dp)*det(1,2))
     dtm_k(1) = fac*(det(1,1)*cos(log(10._dp)*det(2,2)) - &
&     det(2,1)*sin(log(10._dp)*det(2,2)))
     dtm_k(2) = fac*(det(1,1)*sin(log(10._dp)*det(2,2)) + &
&     det(2,1)*cos(log(10._dp)*det(2,2)))
   end if

!  Compute the product of the inverse overlap matrix with the WF

   if (ddkflag == 1) then

     cg1_k(:,:) = zero
     jband1 = 0

     if (itrs == 10 .or. itrs == 11) then

       do jband = minbd, maxbd
         jband1 = jband1 + 1
         do iband = 1, nband_occ
           do ipw = 1, npw_k1

             jpw = pwind_k(ipw)

             if (jpw > 0) then

               wfr = cgq(1,icg1+(iband-1)*npw_k2*nspinor+jpw)*pwnsfac_k(3,jpw)&
&               -cgq(2,icg1+(iband-1)*npw_k2*nspinor+jpw)*pwnsfac_k(4,jpw)
               wfi = cgq(1,icg1+(iband-1)*npw_k2*nspinor+jpw)*pwnsfac_k(4,jpw)&
&               +cgq(2,icg1+(iband-1)*npw_k2*nspinor+jpw)*pwnsfac_k(3,jpw)

               cg1_k(1,(jband1-1)*npw_k1*nspinor + ipw) = &
&               cg1_k(1,(jband1-1)*npw_k1*nspinor + ipw) + &
&               smat_inv(1,iband,jband)*wfr + smat_inv(2,iband,jband)*wfi

               cg1_k(2,(jband1-1)*npw_k1*nspinor + ipw) = &
&               cg1_k(2,(jband1-1)*npw_k1*nspinor + ipw) - &
&               smat_inv(1,iband,jband)*wfi + smat_inv(2,iband,jband)*wfr

             end if

           end do
         end do
       end do

     else

       do jband = minbd, maxbd
         jband1 = jband1 + 1
         do iband = 1, nband_occ
           do ipw = 1, npw_k1

             jpw = pwind_k(ipw)

             if (jpw > 0) then

               wfr = cgq(1,icg1+(iband-1)*npw_k2*nspinor+jpw)*pwnsfac_k(3,jpw)&
&               -cgq(2,icg1+(iband-1)*npw_k2*nspinor+jpw)*pwnsfac_k(4,jpw)
               wfi = cgq(1,icg1+(iband-1)*npw_k2*nspinor+jpw)*pwnsfac_k(4,jpw)&
&               +cgq(2,icg1+(iband-1)*npw_k2*nspinor+jpw)*pwnsfac_k(3,jpw)

               cg1_k(1,(jband1-1)*npw_k1*nspinor + ipw) = &
&               cg1_k(1,(jband1-1)*npw_k1*nspinor + ipw) + &
&               smat_inv(1,iband,jband)*wfr - smat_inv(2,iband,jband)*wfi

               cg1_k(2,(jband1-1)*npw_k1*nspinor + ipw) = &
&               cg1_k(2,(jband1-1)*npw_k1*nspinor + ipw) + &
&               smat_inv(1,iband,jband)*wfi + smat_inv(2,iband,jband)*wfr

             end if

           end do
         end do
       end do

     end if     ! itrs

   end if

 end if         ! job > 0

 ABI_DEALLOCATE(ipvt)
 ABI_DEALLOCATE(zgwork)
 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)

 DBG_EXIT("COLL")

end subroutine smatrix_paw
!!***
