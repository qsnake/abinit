!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawlsylm
!! NAME
!! pawlsylm
!!
!! FUNCTION
!! Compute the LS operator in the real spherical harmonics basis
!! ls_ylm(ilm1,ilm2,ispin)= <sigma, S_lm1| L.S |S_lm2, sigma_prime>
!!   ilm,1m2=(l,m1,m2) with -l<=m1<=l, -l<=m2<=l and 0<l<=lmax
!!   ispin=(sigma,sigma_prime) 1=(up,up), 2=(up,dn), 3=(dn,up), 4=(dn,dn)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!
!! OUTPUT
!!  pawang%ls_ylm(2,l_max**2*(l_max**2+1)/2,2)=LS operator in the real spherical harmonics basis
!!        ls_ylm(:,:,1)=<up, S_lm1| L.S |S_lm2, up>
!!        ls_ylm(:,:,2)=<up, S_lm1| L.S |S_lm2, down>
!!        One can deduce:
!!        <down, S_lm1| L.S |S_lm2, down>=-<up, S_lm1| L.S |S_lm2, up>
!!        <down, S_lm1| L.S |S_lm2, up>  =-Conjg[<up, S_lm1| L.S |S_lm2, down>]
!!        Also, only ilm1<=ilm2 terms are stored, because:
!!         <sigma, S_lm1| L.S |S_lm2, sigma_prime>=-<sigma_prime, S_lm1| L.S |S_lm2, sigma>
!!
!! PARENTS
!!      pawinit
!!
!! CHILDREN
!!      mat_mlms2jmj,mat_slm2ylm,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawlsylm(pawang)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawlsylm'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 type(pawang_type),intent(inout) :: pawang

!Local variables ---------------------------------------
!scalars
 integer :: ii,ilm,im,j0lm,jj,jlm,jm,klm,l_max,ll,lm0,mm,ispden
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: onem
 character(len=500) :: msg
 logical,parameter :: tso=.false. ! use true to Test Spin Orbit and
!                                   write the matrix of L.S in different basis
!arrays
 complex(dpc) :: tmp(2)
 complex(dpc),allocatable :: ls_cplx(:,:,:),slm2ylm(:,:),mat_inp_c(:,:,:)
 complex(dpc),allocatable :: mat_out_c(:,:,:)
 complex(dpc),allocatable :: mat_ls_ylm(:,:,:)
 complex(dpc),allocatable :: mat_jmj(:,:)
 character(len=9),parameter :: dspin2(2)=(/"up-up    ","up-dn    "/)
 character(len=9),parameter :: dspin6(6)=(/"dn       ","up       ","dn-dn    ","up-up    ","dn-up    ","up-dn    "/)
 character(len=9),parameter :: dspinm(6)=(/"dn       ","up       ","n        ","mx       ","my       ","mz       "/)

! *************************************************************************

 DBG_ENTER("COLL")

 if (pawang%use_ls_ylm==0) then
   msg='  ls_ylm pointer is not allocated !'
   MSG_BUG(msg)
 end if

!Initialization
 pawang%ls_ylm=zero
 l_max=pawang%l_max-1

!Nothing to do if lmax=0
 if (l_max<=0) return

!Loop on l quantum number
 do ll=1,l_max

!  Transformation matrixes: real->complex spherical harmonics
   ABI_ALLOCATE(slm2ylm,(2*ll+1,2*ll+1))
   slm2ylm=czero
   do im=1,2*ll+1
     mm=im-ll-1;jm=-mm+ll+1
     onem=dble((-1)**mm)
     if (mm> 0) then
       slm2ylm(im,im)= cmplx(onem*invsqrt2,zero,kind=dp)
       slm2ylm(jm,im)= cmplx(invsqrt2,     zero,kind=dp)
     end if
     if (mm==0) then
       slm2ylm(im,im)=cone
     end if
     if (mm< 0) then
       slm2ylm(im,im)= cmplx(zero,     invsqrt2,kind=dp)
       slm2ylm(jm,im)=-cmplx(zero,onem*invsqrt2,kind=dp)
     end if
   end do

!  Compute <sigma, Y_lm1|L.S|Y_lm2, sigma_prime> (Y_lm=complex spherical harmonics)
!  1= <up|L.S|up>  ;  2= <up|L.S|dn>
   ABI_ALLOCATE(ls_cplx,(2*ll+1,2*ll+1,2))
   ls_cplx=czero
   if(tso)  then
     ABI_ALLOCATE(mat_ls_ylm,(2*ll+1,2*ll+1,4))
     if(tso) mat_ls_ylm=czero
   end if
   if(tso)  then
     ABI_ALLOCATE(mat_jmj,(2*(2*ll+1),2*(2*ll+1)))
     if(tso) mat_jmj=czero
   end if
   do im=1,2*ll+1
     mm=im-ll-1
     ls_cplx(im,im,1)=half*mm
     if(tso) mat_ls_ylm(im,im,1)=-half*mm ! dn dn
     if(tso) mat_ls_ylm(im,im,2)=half*mm  ! up up
     if ((mm+1)<= ll) then
       ls_cplx(im,im+1,2)=half*sqrt(real((ll-mm)*(ll+mm+1),kind=dp))
       if(tso) mat_ls_ylm(im,im+1,4)=half*sqrt(real((ll-mm)*(ll+mm+1),kind=dp))  ! up dn
       if(tso) mat_ls_ylm(im+1,im,3)=half*sqrt(real((ll-mm)*(ll+mm+1),kind=dp))  ! dn up
     end if
     if ((mm-1)>=-ll) then
       ls_cplx(im-1,im,2)=half*sqrt(real((ll+mm)*(ll-mm+1),kind=dp))
       if(tso) mat_ls_ylm(im-1,im,4)=half*sqrt(real((ll+mm)*(ll-mm+1),kind=dp))  ! up dn
       if(tso) mat_ls_ylm(im,im-1,3)=half*sqrt(real((ll+mm)*(ll-mm+1),kind=dp))  ! dn up
     end if
   end do

!  test : print LS in J,M_J basis
   if(tso) then
     do ispden=1,4
       write(msg,'(3a)') ch10,"value of LS in the Ylm basis for " ,trim(dspin6(ispden+2*(4/4)))
       call wrtout(std_out,msg,'COLL')
       do im=1,ll*2+1
         write(msg,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))') (mat_ls_ylm(im,jm,ispden),jm=1,ll*2+1)
         call wrtout(std_out,msg,'COLL')
       end do
     end do
     call mat_mlms2jmj(ll,mat_ls_ylm,mat_jmj,4,1,2,3)  ! optspin=2 : dn spin are first
   end if

!  Compute <sigma, S_lm1|L.S|S_lm2, sigma_prime> (S_lm=real spherical harmonics)
!  1= <up|L.S|up>  ;  2= <up|L.S|dn>
   if(tso) then
     ABI_ALLOCATE(mat_inp_c,(2*ll+1,2*ll+1,4))
     ABI_ALLOCATE(mat_out_c,(2*ll+1,2*ll+1,4))
   end if
   lm0=ll**2
   do jm=1,2*ll+1
     jlm=lm0+jm;j0lm=jlm*(jlm-1)/2
     do im=1,jm
       ilm=lm0+im;klm=j0lm+ilm
       tmp(:)=czero
       do ii=1,2*ll+1
         do jj=1,2*ll+1
           tmp(:)=tmp(:)+ls_cplx(ii,jj,:)*CONJG(slm2ylm(ii,im))*slm2ylm(jj,jm)
         end do
       end do
       pawang%ls_ylm(1,klm,:)=REAL(tmp(:),kind=dp)
       pawang%ls_ylm(2,klm,:)=AIMAG(tmp(:))
     end do
   end do

!  Test: print LS in Slm basis
   if(tso) then
     call mat_slm2ylm(ll,mat_ls_ylm,mat_inp_c,4,2,2,3) ! from Ylm to Slm, and dn spin are first
     do ispden=1,4
       write(msg,'(3a)') ch10,"value of LS in the Slm basis for " ,trim(dspin6(ispden+2*(4/4)))
       call wrtout(std_out,msg,'COLL')
       do im=1,ll*2+1
         write(msg,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))') (mat_inp_c(im,jm,ispden),jm=1,ll*2+1)
         call wrtout(std_out,msg,'COLL')
       end do
     end do
!    change into n,m basis
     mat_ls_ylm(:,:,1)=(mat_inp_c(:,:,1)+mat_inp_c(:,:,2))
     mat_ls_ylm(:,:,2)=(mat_inp_c(:,:,3)+mat_inp_c(:,:,4))
     mat_ls_ylm(:,:,3)=-cmplx(0.d0,1.d0)*(mat_inp_c(:,:,4)-mat_inp_c(:,:,3))
     mat_ls_ylm(:,:,4)=(mat_inp_c(:,:,1)-mat_inp_c(:,:,2))
     do ispden=1,4
       write(msg,'(3a)') ch10,"value of LS in the Slm basis for " ,trim(dspinm(ispden+2*(4/4)))
       call wrtout(std_out,msg,'COLL')
       do im=1,ll*2+1
         write(msg,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))') (mat_ls_ylm(im,jm,ispden),jm=1,ll*2+1)
         call wrtout(std_out,msg,'COLL')
       end do
     end do
     ABI_DEALLOCATE(mat_inp_c)
     ABI_DEALLOCATE(mat_ls_ylm)
     ABI_DEALLOCATE(mat_jmj)
     ABI_DEALLOCATE(mat_out_c)
   end if ! tso

   ABI_DEALLOCATE(ls_cplx)
   ABI_DEALLOCATE(slm2ylm)

!  End loop on l
 end do

 DBG_EXIT("COLL")

 end subroutine pawlsylm

!!***
