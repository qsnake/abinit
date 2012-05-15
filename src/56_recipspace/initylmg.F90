!{\src2tex{textfont=tt}}
!!****f* ABINIT/initylmg
!! NAME
!! initylmg
!!
!! FUNCTION
!! Calculate the real spherical harmonics Ylm (and gradients)
!! over a set of (reciprocal space) (k+G) vectors
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive
!!              translations (b^-1)
!!  kg(3,mpw)=integer coordinates of G vectors in basis sphere
!!  kptns(3,nkpt)=k points in terms of reciprocal translations
!!  mkmem =maximum number of k points which can fit in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang=1+maximum angular momentum for nonlocal pseudopotential
!!  mpw   =maximum number of planewaves in basis sphere
!!         (large number)
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  nkpt  =number of k points
!!  npwarr(nkpt)=array holding npw for each k point
!!  nsppol=1 for unpolarized, 2 for polarized
!!  optder= 0=compute Ylm(K)
!!          1=compute Ylm(K) and dYlm/dKi
!!          2=compute Ylm(K), dYlm/dKi and d2Ylm/dKidKj
!!         -1=compute only dYlm/dKi
!!  rprimd(3,3)=dimensional primitive translations in real space
!!              (bohr)
!!  unkg=unit number for (k+G) sphere data
!!  unylm=unit number for storage of Ylm on disk
!!
!! OUTPUT
!!  if (optder>=0)
!!    ylm(mpw*mkmem,mpsang*mpsang) = real spherical harmonics
!!    for each G and k point
!!  if (optder>=1 or optder==-1)
!!    ylm_gr(mpw*mkmem,1:3,mpsang*mpsang)= gradients of real
!!    spherical harmonics wrt (G+k) in reduced coordinates
!!  if (optder>=2)
!!    ylm_gr(mpw*mkmem,4:9,mpsang*mpsang)= second gradients of
!!    real spherical harmonics wrt (G+k) in reduced coordinates
!!
!! NOTES
!! Remember the expression of complex spherical harmonics:
!! $Y_{lm}(%theta ,%phi)=sqrt{{(2l+1) over (4 %pi)}
!! {fact(l-m) over fact(l+m)} } P_l^m(cos(%theta))
!! func e^{i m %phi}$
!! Remember the expression of real spherical harmonics as
!!   linear combination of complex spherical harmonics:
!! $Yr_{lm}(%theta ,%phi)=(Re{Y_{l-m}}+(-1)^m Re{Y_{lm}})/sqrt{2}
!! $Yr_{l-m}(%theta ,%phi)=(Im{Y_{l-m}}-(-1)^m Im{Y_{lm}})/sqrt{2}
!!
!! PARENTS
!!      debug_tools,gstate,ks_ddiago,loop3dte,loper3,m_cprj_bspline,m_paw_pwij
!!      m_shirley,m_wfs,mover,nstpaw3,partial_dos_fractions,respfn,scfcv,wffile
!!
!! CHILDREN
!!      leave_new,leave_test,plm_coeff,timab,wrtout,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initylmg(gprimd,kg,kptns,mkmem,mpi_enreg,mpsang,mpw,&
  & nband,nkpt,npwarr,nsppol,optder,rprimd,unkg,unylm,ylm,ylm_gr)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initylmg'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mkmem,mpsang,mpw,nkpt,nsppol,optder
 integer,intent(in) :: unkg,unylm
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),nband(nkpt*nsppol)
 integer,intent(in) :: npwarr(nkpt)
 real(dp),intent(in) :: gprimd(3,3),kptns(3,nkpt),rprimd(3,3)
 real(dp),intent(out) :: ylm(mpw*mkmem,mpsang*mpsang)
 real(dp),intent(out) :: ylm_gr(mpw*mkmem,3+6*(optder/2),mpsang*mpsang)

!Local variables ------------------------------
!scalars
 integer :: dimgr,ia,ib,idir,ig,ii,ikg,ikpt,ilang,ilm,ipw
 integer :: jj,kk,l0,ll
 integer :: me_distrb,mm,npw_k
 real(dp),parameter :: tol=1.d-10
 real(dp) :: cphi,ctheta,fact,onem,rr,sphi,stheta,work1,work2
 real(dp) :: xx,ylmcst,ylmcst2
 real(dp) :: yy,zz
 character(len=500) :: message
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/)
 integer,parameter :: beta(6)=(/1,2,3,2,1,1/)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: dphi(3),dtheta(3),iphase(mpsang-1),kpg(3)
 real(dp) :: rphase(mpsang-1),tsec(2)
 real(dp),allocatable :: blm(:,:),ylm_disk(:,:),ylm_gr_disk(:,:,:)
 real(dp),allocatable :: ylmgr2_cart(:,:,:),ylmgr2_tmp(:,:)
 real(dp),allocatable :: ylmgr_cart(:,:)
 real(dp),allocatable :: ylmgr_red(:,:)
!no_abirules

!*****************************************************************

 if (mkmem==0.and.optder<0) then
   write(message, '(4a)' )ch10,&
&   ' initylmg : BUG -',ch10,&
&   '  mkmem=0 not compatible with optder<0 !'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   call leave_new('COLL')
 end if

!Begin executable
 call xme_init(mpi_enreg,me_distrb)

!Initialisation of spherical harmonics (and gradients)
 if (optder>=0) ylm(:,:)  =zero
 if (optder/=0) ylm_gr(:,:,:)=zero
 dimgr=3+6*(optder/2)

!In case the information are to be kept on disk,
!rewind the disk file and allocate the work space
 if (mkmem==0) then
   rewind (unit=unylm)
   ABI_ALLOCATE(ylm_disk,(mpw,mpsang*mpsang))
   ylm_disk(:,:)=zero
   if (optder/=0) then
     ABI_ALLOCATE(ylm_gr_disk,(mpw,dimgr,mpsang*mpsang))
     ylm_gr_disk(:,:,:)=zero
   end if
 end if

!Loop over k-points:
 ikg=0;if (mkmem==0) rewind (unit=unkg)
 do ikpt=1,nkpt

   if(mpi_enreg%paral_compil_kpt==1)then
     if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband(ikpt),1:nsppol)-me_distrb))/=0) then
       cycle ! Skip the rest of the k-point loop
     end if
   end if

!  Load k+G-vectors, for this k-point:
   if (mkmem==0) then
     read (unkg) npw_k
     ABI_ALLOCATE(kg_k,(3,npw_k))
     read (unkg)
     read (unkg) ((kg_k(ii,ig),ii=1,3),ig=1,npw_k)
   else
     npw_k=npwarr(ikpt)
     ABI_ALLOCATE(kg_k,(3,npw_k))
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
   end if

!  Special case for l=0
   if (mkmem==0) then
     if (optder>=0) ylm_disk   (1:npw_k,1)=1._dp/sqrt(four_pi)
     if (optder/=0) ylm_gr_disk(1:npw_k,1:dimgr,1)=zero
   else
     if (optder>=0) ylm(1+ikg:npw_k+ikg,1)=1._dp/sqrt(four_pi)
     if (optder/=0) ylm_gr(1+ikg:npw_k+ikg,1:dimgr,1)=zero
   end if
   if (mpsang>1) then

!    Loop over all k+G
     do ipw=1,npw_k

!      Load k+G
       kpg(1)=kptns(1,ikpt)+real(kg_k(1,ipw),dp)
       kpg(2)=kptns(2,ikpt)+real(kg_k(2,ipw),dp)
       kpg(3)=kptns(3,ikpt)+real(kg_k(3,ipw),dp)

!      Calculate module of k+G
       xx=gprimd(1,1)*kpg(1)+gprimd(1,2)*kpg(2)+gprimd(1,3)*kpg(3)
       yy=gprimd(2,1)*kpg(1)+gprimd(2,2)*kpg(2)+gprimd(2,3)*kpg(3)
       zz=gprimd(3,1)*kpg(1)+gprimd(3,2)*kpg(2)+gprimd(3,3)*kpg(3)
       rr=sqrt(xx**2+yy**2+zz**2)

!      Continue only for k+G<>0
       if (rr>tol) then

!        Determine theta and phi
         cphi=one
         sphi=zero
         ctheta=zz/rr
         stheta=sqrt(abs((one-ctheta)*(one+ctheta)))
         if (stheta>tol) then
           cphi=xx/(rr*stheta)
           sphi=yy/(rr*stheta)
         end if
         do mm=1,mpsang-1
           rphase(mm)=dreal(dcmplx(cphi,sphi)**mm)
           iphase(mm)=aimag(dcmplx(cphi,sphi)**mm)
         end do

!        Determine gradients of theta and phi
         if (optder/=0) then
           dtheta(1)=ctheta*cphi
           dtheta(2)=ctheta*sphi
           dtheta(3)=-stheta
           dphi(1)=-sphi
           dphi(2)=cphi
           dphi(3)=zero
         end if

!        COMPUTE Ylm(K)
!        ============================================
         if (optder>=0) then
!          Loop over angular momentum l
           do ilang=2,mpsang
             ll=ilang-1
             l0=ll**2+ll+1
             fact=1._dp/real(ll*(ll+1),dp)
             ylmcst=sqrt(real(2*ll+1,dp)/four_pi)
!            Special case m=0
             if (mkmem==0) then
               ylm_disk(ipw,l0)=ylmcst*ass_leg_pol(ll,0,ctheta)
             else
               ylm(ikg+ipw,l0)=ylmcst*ass_leg_pol(ll,0,ctheta)
             end if
!            Compute for m>0
             onem=one
             do mm=1,ll
               onem=-onem
               work1=ylmcst*sqrt(fact)*onem*ass_leg_pol(ll,mm,ctheta)*sqrt(2._dp)
               if (mkmem==0) then
                 ylm_disk(ipw,l0+mm)=work1*rphase(mm)
                 ylm_disk(ipw,l0-mm)=work1*iphase(mm)
               else
                 ylm(ikg+ipw,l0+mm)=work1*rphase(mm)
                 ylm(ikg+ipw,l0-mm)=work1*iphase(mm)
               end if
               if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
             end do ! End loop over m
           end do  ! End loop over l
         end if

!        COMPUTE dYlm/dKi
!        ============================================
         if (optder/=0) then
           ABI_ALLOCATE(ylmgr_cart,(3,2))
           ABI_ALLOCATE(ylmgr_red,(3,2))
!          Loop over angular momentum l
           do ilang=2,mpsang
             ll=ilang-1
             l0=ll**2+ll+1
             fact=1._dp/real(ll*(ll+1),dp)
             ylmcst=sqrt(real(2*ll+1,dp)/four_pi)/rr
!            === Special case m=0 ===
!            1-compute gradients in cartesian coordinates
             work1=ylmcst*plm_dtheta(ll,0,ctheta)
             ylmgr_cart(1:3,1)=work1*dtheta(1:3)
!            2-Transfer gradients into reduced coordinates
             do ii=1,3
               ylmgr_red(ii,1)=(rprimd(1,ii)*ylmgr_cart(1,1)+&
&               rprimd(2,ii)*ylmgr_cart(2,1)+&
&               rprimd(3,ii)*ylmgr_cart(3,1))
             end do
!            3-Store gradients
             if (mkmem==0) then
               ylm_gr_disk(ipw,1:3,l0)=ylmgr_red(1:3,1)
             else
               ylm_gr(ikg+ipw,1:3,l0) =ylmgr_red(1:3,1)
             end if
!            === Compute for m>0 ===
             onem=one
             do mm=1,ll
               onem=-onem
!              1-compute gradients in cartesian coordinates
               work1=ylmcst*sqrt(fact)*onem*plm_dtheta(ll,mm,ctheta)*sqrt(2._dp)
               work2=ylmcst*sqrt(fact)*onem*plm_dphi  (ll,mm,ctheta)*sqrt(2._dp)
               ylmgr_cart(1:3,1)=rphase(mm)*work1*dtheta(1:3)-iphase(mm)*work2*dphi(1:3)
               ylmgr_cart(1:3,2)=iphase(mm)*work1*dtheta(1:3)+rphase(mm)*work2*dphi(1:3)
!              2-Transfer gradients into reduced coordinates
               do kk=1,2
                 do ii=1,3
                   ylmgr_red(ii,kk)=(rprimd(1,ii)*ylmgr_cart(1,kk)+&
&                   rprimd(2,ii)*ylmgr_cart(2,kk)+&
&                   rprimd(3,ii)*ylmgr_cart(3,kk))
                 end do
               end do
!              3-Store gradients
               if (mkmem==0) then
                 ylm_gr_disk(ipw,1:3,l0+mm)=ylmgr_red(1:3,1)
                 ylm_gr_disk(ipw,1:3,l0-mm)=ylmgr_red(1:3,2)
               else
                 ylm_gr(ikg+ipw,1:3,l0+mm) =ylmgr_red(1:3,1)
                 ylm_gr(ikg+ipw,1:3,l0-mm) =ylmgr_red(1:3,2)
               end if
               if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
             end do ! End loop over m
           end do  ! End loop over l
           ABI_DEALLOCATE(ylmgr_cart)
           ABI_DEALLOCATE(ylmgr_red)
         end if

!        COMPUTE d2Ylm/dKidKj
!        ============================================
         if (optder==2) then
           ABI_ALLOCATE(ylmgr2_cart,(3,3,2))
           ABI_ALLOCATE(ylmgr2_tmp,(3,3))
           ABI_ALLOCATE(ylmgr_red,(6,2))
           ABI_ALLOCATE(blm,(5,mpsang*mpsang))
           call plm_coeff(blm,mpsang,ctheta)

!          Loop over angular momentum l
           do ilang=2,mpsang
             ll=ilang-1
             l0=ll**2+ll+1
             fact=1._dp/real(ll*(ll+1),dp)
             ylmcst=sqrt(real(2*ll+1,dp)/four_pi)/(rr**2)
!            === Special case m=0 ===
!            1-compute gradients in cartesian coordinates
             ylmgr2_cart(1,1,1)=ylmcst*(-blm(3,l0)*sphi*sphi+blm(4,l0)*cphi*cphi)
             ylmgr2_cart(2,2,1)=ylmcst*(-blm(3,l0)*cphi*cphi+blm(4,l0)*sphi*sphi)
             ylmgr2_cart(3,3,1)=ylmcst*blm(1,l0)
             ylmgr2_cart(3,1,1)=ylmcst*blm(2,l0)*cphi
             ylmgr2_cart(3,2,1)=ylmcst*blm(2,l0)*sphi
             ylmgr2_cart(2,1,1)=ylmcst*(blm(3,l0)+blm(4,l0))*sphi*cphi
             ylmgr2_cart(1,3,1)=ylmgr2_cart(3,1,1)
             ylmgr2_cart(1,2,1)=ylmgr2_cart(2,1,1)
             ylmgr2_cart(2,3,1)=ylmgr2_cart(3,2,1)
!            2-Transfer gradients into reduced coordinates
             do jj=1,3
               do ii=1,3
                 ylmgr2_tmp(ii,jj)=(rprimd(1,jj)*ylmgr2_cart(1,ii,1)+&
&                 rprimd(2,jj)*ylmgr2_cart(2,ii,1)+&
&                 rprimd(3,jj)*ylmgr2_cart(3,ii,1))
               end do
             end do
             do ii=1,6
               ia=alpha(ii);ib=beta(ii)
               ylmgr_red(ii,1)=(rprimd(1,ia)*ylmgr2_tmp(1,ib)+&
&               rprimd(2,ia)*ylmgr2_tmp(2,ib)+&
&               rprimd(3,ia)*ylmgr2_tmp(3,ib))
             end do
             if (mkmem==0) then
               ylm_gr_disk(ipw,4:9,l0)=ylmgr_red(1:6,1)
             else
               ylm_gr(ikg+ipw,4:9,l0) =ylmgr_red(1:6,1)
             end if
!            === Compute for m>0 ===
             onem=one
             do mm=1,ll
               onem=-onem;ylmcst2=ylmcst*sqrt(fact)*sqrt(two)
               ylmgr2_cart(1,1,1)=ylmcst2*((-blm(3,l0+mm)*sphi*sphi+blm(4,l0+mm)*cphi*cphi)*rphase(mm)-&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*iphase(mm))
               ylmgr2_cart(1,1,2)=ylmcst2*((-blm(3,l0+mm)*sphi*sphi+blm(4,l0+mm)*cphi*cphi)*iphase(mm)+&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*rphase(mm))
               ylmgr2_cart(2,2,1)=ylmcst2*((-blm(3,l0+mm)*cphi*cphi+blm(4,l0+mm)*sphi*sphi)*rphase(mm)+&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*iphase(mm))
               ylmgr2_cart(2,2,2)=ylmcst2*((-blm(3,l0+mm)*cphi*cphi+blm(4,l0+mm)*sphi*sphi)*iphase(mm)-&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*rphase(mm))
               ylmgr2_cart(3,3,1)=ylmcst2*blm(1,l0+mm)*rphase(mm)
               ylmgr2_cart(3,3,2)=ylmcst2*blm(1,l0+mm)*iphase(mm)
               ylmgr2_cart(3,1,1)=ylmcst2*(blm(2,l0+mm)*cphi*rphase(mm)-&
&               mm*iphase(mm)*sphi*onem*plm_dtheta(ll,mm,ctheta))
               ylmgr2_cart(3,1,2)=ylmcst2*(blm(2,l0+mm)*cphi*iphase(mm)+&
&               mm*rphase(mm)*sphi*onem*plm_dtheta(ll,mm,ctheta))
               ylmgr2_cart(3,2,1)=ylmcst2*(blm(2,l0+mm)*sphi*rphase(mm)+&
&               mm*iphase(mm)*cphi*onem*plm_dtheta(ll,mm,ctheta))
               ylmgr2_cart(3,2,2)=ylmcst2*(blm(2,l0+mm)*sphi*iphase(mm)-&
&               mm*rphase(mm)*cphi*onem*plm_dtheta(ll,mm,ctheta))
               ylmgr2_cart(2,1,1)=ylmcst2*((blm(3,l0+mm)+blm(4,l0+mm))*sphi*cphi*rphase(mm)-&
&               blm(5,l0+mm)*(sphi*sphi-cphi*cphi)*mm*iphase(mm))
               ylmgr2_cart(2,1,2)=ylmcst2*((blm(3,l0+mm)+blm(4,l0+mm))*sphi*cphi*iphase(mm)+&
&               blm(5,l0+mm)*(sphi*sphi-cphi*cphi)*mm*rphase(mm))
               ylmgr2_cart(1,3,:)=ylmgr2_cart(3,1,:)
               ylmgr2_cart(1,2,:)=ylmgr2_cart(2,1,:)
               ylmgr2_cart(2,3,:)=ylmgr2_cart(3,2,:)
!              2-Transfer gradients into reduced coordinates
               do kk=1,2
                 do jj=1,3
                   do ii=1,3
                     ylmgr2_tmp(ii,jj)=(rprimd(1,jj)*ylmgr2_cart(1,ii,kk)+&
&                     rprimd(2,jj)*ylmgr2_cart(2,ii,kk)+&
&                     rprimd(3,jj)*ylmgr2_cart(3,ii,kk))
                   end do
                 end do
                 do ii=1,6
                   ia=alpha(ii);ib=beta(ii)
                   ylmgr_red(ii,kk)=(rprimd(1,ia)*ylmgr2_tmp(1,ib)+&
&                   rprimd(2,ia)*ylmgr2_tmp(2,ib)+&
&                   rprimd(3,ia)*ylmgr2_tmp(3,ib))
                 end do
               end do
               if (mkmem==0) then
                 ylm_gr_disk(ipw,4:9,l0+mm)=ylmgr_red(1:6,1)
                 ylm_gr_disk(ipw,4:9,l0-mm)=ylmgr_red(1:6,2)
               else
                 ylm_gr(ikg+ipw,4:9,l0+mm) =ylmgr_red(1:6,1)
                 ylm_gr(ikg+ipw,4:9,l0-mm) =ylmgr_red(1:6,2)
               end if
               if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
             end do ! End loop over m
           end do  ! End loop over l
           ABI_DEALLOCATE(blm)
           ABI_DEALLOCATE(ylmgr2_cart)
           ABI_DEALLOCATE(ylmgr2_tmp)
           ABI_DEALLOCATE(ylmgr_red)
         end if

!        End condition r<>0
       end if

!      End loop over k+G
     end do

!    End condition l<>0
   end if

   ABI_DEALLOCATE(kg_k)

   ikg=ikg+npw_k

!  Eventually stores Ylm (and dYlm/dKi) on disk
   if (mkmem==0) then
     write(unylm) npw_k
     if (optder/=0) then
       write(unylm) ((ylm_disk(ipw,ilm),ipw=1,npw_k),ilm=1,mpsang*mpsang),&
&       (((ylm_gr_disk(ipw,idir,ilm),ipw=1,npw_k),idir=1,dimgr),ilm=1,mpsang*mpsang)
     else
       write(unylm) ((ylm_disk(ipw,ilm),ipw=1,npw_k),ilm=1,mpsang*mpsang)
     end if
   end if

!  End Loop over k-points
 end do

 if(mkmem==0) then
   ABI_DEALLOCATE(ylm_disk)
   if (optder/=0)  then
     ABI_DEALLOCATE(ylm_gr_disk)
   end if
 end if

 if(mpi_enreg%paral_compil_kpt==1)then
   call timab(66,1,tsec)
   call leave_test()
   call timab(66,2,tsec)
 end if

end subroutine initylmg
!!***
