!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawgylm
!! NAME
!! pawgylm
!!
!! FUNCTION
!! Compute g_l(r)*Y_lm(r) (and derivatives) on the fine (rectangular) grid
!! around one atom (g_l=radial shape function).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lm_size=number of lm components to be calculated
!!  nfgd= number of (fine grid) FFT points in the paw sphere around current atom
!!  optgr0= 1 if g_l(r)*Y_lm(r) are computed
!!  optgr1= 1 if first derivatives of g_l(r)*Y_lm(r) are computed
!!  optgr2= 1 if second derivatives of g_l(r)*Y_lm(r) are computed
!!  pawtab <type(pawtab_type)>=paw tabulated starting data for current atom
!!  rfgd(3,nfgd)= coordinates of r-r_atom on the fine rect. grid around current atom
!!  rfgd_allocated= >0 if rfgd array has been allocated (for testing purpose only)
!!
!! OUTPUT
!!  if (optgr0==1)
!!    gylm(nfgd,lm_size)= g_l(r)*Y_lm(r) around current atom
!!  if (optgr1==1)
!!    gylmgr(3,nfgd,lm_size)= derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optgr2==1)
!!    gylmgr2(6,nfgd,lm_size)= second derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!
!! PARENTS
!!      nhatgrid,paw_mknewh0,pawdij,pawfrnhat,pawgrnl,pawmknhat
!!      pawmknhat_psipsi
!!
!! CHILDREN
!!      initylmr,jbessel,sort_dp,splint
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawgylm(gylm,gylmgr,gylmgr2,lm_size,nfgd,optgr0,optgr1,optgr2,pawtab,rfgd,rfgd_allocated)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawgylm'
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lm_size,nfgd,optgr0,optgr1,optgr2
 integer,intent(in) :: rfgd_allocated
 type(pawtab_type),intent(in) :: pawtab
!arrays
 real(dp),intent(in) :: rfgd(3,nfgd)
 real(dp),intent(out) :: gylm(nfgd,optgr0*lm_size)
 real(dp),intent(out) :: gylmgr(3,nfgd,optgr1*lm_size)
 real(dp),intent(out) :: gylmgr2(6,nfgd,optgr2*lm_size)

!Local variables ------------------------------
!scalars
 integer :: argl,ic,ilm,izero,l_size,lambda,ll,normchoice,option,shape_type
 real(dp) :: arg
 real(dp) :: d2shpfunc1,d2shpfunc1_ovr2_0,d2shpfunc1_ovr2_0_2
 real(dp) :: d2shpfunc1_ovr2_0_3,d2shpfunc1_ovr2_0_4
 real(dp) :: d2shpfunc2,d2shpfunc2_ovr2_0,d2shpfunc3
 real(dp) :: dshpfunc1,dshpfunc1_ovr_0,dshpfunc1_ovr_0_2
 real(dp) :: dshpfunc2,dshpfunc2_ovr_0,dshpfunc3
 real(dp) :: jbes1,jbes2,jbesp1,jbesp2,jbespp1,jbespp2,pi_over_rshp
 real(dp) :: shapefunc1,shapefunc1_0,shapefunc2,shapefunc2_0,shapefunc3
 real(dp) :: sigma,splfact
 logical :: compute_gr0,compute_gr1,compute_gr2
 character(len=500) :: msg
!arrays
 integer,allocatable :: isort(:)
 real(dp),parameter :: ffact(1:9)=(/1._dp,3._dp,15._dp,105._dp,945._dp,10395._dp,&
&                                   135135._dp,2027025._dp,34459425._dp/)
 real(dp),parameter :: toldev=tol3
 real(dp) :: ss(3)
 real(dp),allocatable :: alpha(:,:),cc(:,:),d2gfact(:,:),d2shpfuncnum(:,:),dgfact(:,:)
 real(dp),allocatable :: dshpfuncnum(:,:),gfact(:,:)
 real(dp),allocatable :: qq(:,:),rnrm(:),rnrm_inv(:),rnrm_sort(:)
 real(dp),allocatable :: shpfuncnum(:,:),work(:),ylmr(:,:),ylmrgr(:,:,:)
!no_abirules
!Statement functions -----------------------------------
!shapefunc1 is g(x) (gaussian)
 shapefunc1(arg)=exp(-(arg/sigma)**lambda)
!shapefunc1_0 is g(x) (gaussian) for small x
 shapefunc1_0(arg)=one-(arg/sigma)**lambda+half*(arg/sigma)**(2*lambda)-sixth*(arg/sigma)**(3*lambda)
!shapefunc2 is g(x) (sinc2)
 shapefunc2(arg)=(sin(pi_over_rshp*arg)/(pi_over_rshp*arg))**2
!shapefunc2_0 is g(x) (sinc2) for small x
 shapefunc2_0(arg)=one-third*(pi_over_rshp*arg)**2+two*(pi_over_rshp*arg)**4/45._dp
!shapefunc3 is g(x) (Bessel)
 shapefunc3(jbes1,jbes2,argl)= alpha(1,1+argl)*jbes1+alpha(2,1+argl)*jbes2
!dshpfunc1(x) is g_prime(x) (gaussian)
 dshpfunc1(arg)=-lambda/sigma*(arg/sigma)**(lambda-1)*exp(-(arg/sigma)**lambda)
!dshpfunc1_ovr_0(x) is g_prime(x)/x (gaussian) for small x and lambda>2
 dshpfunc1_ovr_0(arg)=-lambda/sigma**2*((arg/sigma)**(lambda-2)-(arg/sigma)**(2*lambda-2))
!dshpfunc1_ovr_0_2(x) is g_prime(x)/x (gaussian) for small x and lambda=2
 dshpfunc1_ovr_0_2(arg)=-two/sigma**2*(one-(arg/sigma)**2+half*(arg/sigma)**4)
!dshpfunc2(x) is g_prime(x) (sinc2)
 dshpfunc2(arg)=two*pi_over_rshp*sin(pi_over_rshp*arg)/(pi_over_rshp*arg)**3&
&              *(pi_over_rshp*arg*cos(pi_over_rshp*arg)-sin(pi_over_rshp*arg))
!dshpfunc2_ovr_0(x) is g_prime(x)/x (sinc2) for small x
 dshpfunc2_ovr_0(arg)=-two_thirds*pi_over_rshp**2+eight*pi_over_rshp**4*arg**2/45._dp
!dshpfunc3(x) is g_prime(x) (Bessel)
 dshpfunc3(jbesp1,jbesp2,argl)= alpha(1,1+argl)*qq(1,1+argl)*jbesp1 &
&                              +alpha(2,1+argl)*qq(2,1+argl)*jbesp2
!d2shpfunc1(x) is g_prime_prime(x) (gaussian)
 d2shpfunc1(arg)=lambda/(sigma**2)*(lambda*(arg/sigma)**(2*lambda-2) &
&               -(lambda-1)*(arg/sigma)**(lambda-2))*exp(-(arg/sigma)**lambda)
!d2shpfunc1_ovr2_0(x) is (g_prime_prime(x)-g_prime(x)/x)/x**2 (gaussian) for small x and lambda>4
 d2shpfunc1_ovr2_0(arg)=-lambda/(sigma**4)*((lambda-2)*(arg/sigma)**(lambda-4) &
&                                          -(lambda-1)*two*(arg/sigma)**(2*lambda-4))
!d2shpfunc1_ovr2_0_2(x) is (g_prime_prime(x)-g_prime(x)/x)/x**2 (gaussian) for small x and lambda==2
 d2shpfunc1_ovr2_0_2(arg)=four/(sigma**4)*(one-(arg/sigma)**2)
!d2shpfunc1_ovr2_0_3(x) is (g_prime_prime(x)-g_prime(x)/x)/x**2 (gaussian) for small x and lambda==3
 d2shpfunc1_ovr2_0_3(arg)=-three/arg/sigma**3+12._dp*arg**2/sigma**6
!d2shpfunc1_ovr2_0_4(x) is (g_prime_prime(x)-g_prime(x)/x)/x**2 (gaussian) for small x and lambda==4
 d2shpfunc1_ovr2_0_4(arg)=-eight/(sigma**4)*(one-three*(arg/sigma)**4)
!d2shpfunc2(x) is g_prime_prime(x) (sinc2)
 d2shpfunc2(arg)=two/(pi_over_rshp**2*arg**4)* &
&               (pi_over_rshp**2*arg**2*(cos(pi_over_rshp*arg))**2 &
&               +(three-pi_over_rshp**2*arg**2)*(sin(pi_over_rshp*arg))**2 &
&               -four*pi_over_rshp*arg*cos(pi_over_rshp*arg)*sin(pi_over_rshp*arg))
!d2shpfunc2_ovr2_0(x) is (g_prime_prime(x)-g_prime(x)/x)/x**2 (sinc2) for small x
 d2shpfunc2_ovr2_0(arg)=16._dp/45._dp*pi_over_rshp**4-eight/105._dp*pi_over_rshp**6*arg**2
!d2shpfunc3(x) is g_prime_prime(x) (Bessel)
 d2shpfunc3(jbespp1,jbespp2,argl)= alpha(1,1+argl)*(qq(1,1+argl)**2)*jbespp1 &
&                                 +alpha(2,1+argl)*(qq(2,1+argl)**2)*jbespp2

! *************************************************************************

 DBG_ENTER("COLL")

 if (optgr0==0.and.optgr1==0.and.optgr2==0) return
 if (nfgd==0) return

!Compatibility test
!==========================================================
 if (rfgd_allocated==0) then
   msg='rfgd array must be allocated !'
   MSG_BUG(msg)
 end if
 if (pawtab%lcut_size>9) then
   msg='l_size>10 forbidden !'
   MSG_BUG(msg)
 end if
 if (pawtab%shape_type==1.and.pawtab%shape_lambda<2) then
   msg='Exponent lambda of gaussian shape function must be > 1 !'
   MSG_ERROR(msg)
 end if

!Initializations
!==========================================================
!Options for computation
 compute_gr0=(optgr0==1.or.optgr1==1.or.optgr2==1)
 compute_gr1=(optgr1==1.or.optgr2==1)
 compute_gr2=(optgr2==1)
 l_size=pawtab%lcut_size

!Norms of vectors around the atom
 ABI_ALLOCATE(rnrm,(nfgd))
 izero=-1
 do ic=1,nfgd
   rnrm(ic)=sqrt(rfgd(1,ic)**2+rfgd(2,ic)**2+rfgd(3,ic)**2)
   if (rnrm(ic)<=tol10) izero=ic  ! Has to be consistent with initylmr !!
!  if (rnrm(ic)<=ten*epsilon(one)) izero=ic
 end do

!Initializations
 if (optgr0==1) gylm=zero
 if (optgr1==1) gylmgr=zero
 if (optgr2==1) gylmgr2=zero

!Some definitions concerning shape function g_l(r)
 shape_type=pawtab%shape_type
 sigma=pawtab%shape_sigma;lambda=pawtab%shape_lambda
 pi_over_rshp=pi/pawtab%rshp
 if (shape_type==3) then
   ABI_ALLOCATE(alpha,(2,l_size))
   ABI_ALLOCATE(qq,(2,l_size))
   do ll=1,l_size
     alpha(1:2,ll)=pawtab%shape_alpha(1:2,ll)
     qq(1:2,ll)=pawtab%shape_q(1:2,ll)
   end do
 end if

!If needed, sort selected radii by increasing norm
 if (shape_type==-1) then
   ABI_ALLOCATE(isort,(nfgd))
   ABI_ALLOCATE(rnrm_sort,(nfgd))
   do ic=1,nfgd
     isort(ic)=ic
   end do
   rnrm_sort(1:nfgd)=rnrm(1:nfgd)
   call sort_dp(nfgd,rnrm_sort,isort,tol16)
 end if

!If shape function is "numeric", spline it onto selected radii
 if (shape_type==-1) then
   ABI_ALLOCATE(work,(nfgd))
   if (compute_gr0) then
     ABI_ALLOCATE(shpfuncnum,(nfgd,l_size))
     do ll=1,l_size
       call splint(pawtab%mesh_size,pawtab%rad_for_spline,pawtab%shapefunc(:,ll),&
&       pawtab%dshpfunc(:,ll,2),nfgd,rnrm_sort,work)
       do ic=1,nfgd
         shpfuncnum(isort(ic),ll)=work(ic)
       end do
     end do
   end if
   if(compute_gr1) then
     ABI_ALLOCATE(dshpfuncnum,(nfgd,l_size))
     do ll=1,l_size
       call splint(pawtab%mesh_size,pawtab%rad_for_spline,pawtab%dshpfunc(:,ll,1),&
&       pawtab%dshpfunc(:,ll,3),nfgd,rnrm_sort,work)
       do ic=1,nfgd
         dshpfuncnum(isort(ic),ll)=work(ic)
       end do
     end do
   end if
   if(compute_gr2) then
     ABI_ALLOCATE(d2shpfuncnum,(nfgd,l_size))
     do ll=1,l_size
       call splint(pawtab%mesh_size,pawtab%rad_for_spline,pawtab%dshpfunc(:,ll,2),&
&       pawtab%dshpfunc(:,ll,4),nfgd,rnrm_sort,work)
       do ic=1,nfgd
         d2shpfuncnum(isort(ic),ll)=work(ic)
       end do
     end do
   end if
   ABI_DEALLOCATE(work)
 end if

 if (shape_type==-1)  then
   ABI_DEALLOCATE(isort)
   ABI_DEALLOCATE(rnrm_sort)
 end if

!If needed, compute limits at r=0 of shape function and derivatives
 if (izero>0) then
   ABI_ALLOCATE(cc,(3,min(l_size,3)))
   cc=zero
   if (shape_type==-1) then
     splfact=(pawtab%rad_for_spline(4)-pawtab%rad_for_spline(1))&
&     /(pawtab%rad_for_spline(3)-pawtab%rad_for_spline(2))
   end if
   do ll=1,min(l_size,3)
     if (optgr0==1.or.optgr1==1.or.optgr2==1) then
       if (shape_type==-1) then
         ss(1:3)=pawtab%shapefunc(2:4,ll)/pawtab%rad_for_spline(2:4)**(ll-1)
         cc(1,ll)=ss(3)+(ss(1)-ss(2))*splfact
       else if (shape_type==1.or.shape_type==2) then
         cc(1,ll)=one
       else if (shape_type==3) then
         cc(1,ll)=(alpha(1,ll)*qq(1,ll)**(ll-1) &
&         +alpha(2,ll)*qq(2,ll)**(ll-1))/ffact(ll)
       end if
       cc(1,ll)=cc(1,ll)*pawtab%gnorm(ll)
     end if
     if (optgr1==1.or.optgr2==1) then
       if (shape_type==-1) then
         ss(1:3)=(ss(1:3)-cc(1,ll))/pawtab%rad_for_spline(2:4)
         cc(2,ll)=ss(3)+(ss(1)-ss(2))*splfact
       else if (shape_type==1.and.lambda==1) then
         cc(2,ll)=-one/sigma
       else
         cc(2,ll)=zero
       end if
       cc(2,ll)=cc(2,ll)*pawtab%gnorm(ll)
     end if
     if (optgr2==1) then
       if (shape_type==-1) then
         ss(1:3)=(ss(1:3)-cc(2,ll))/pawtab%rad_for_spline(2:4)
         cc(3,ll)=ss(3)+(ss(1)-ss(2))*splfact
       else if (shape_type==1) then
         if (lambda==1) cc(3,ll)=half/sigma**2
         if (lambda==2) cc(3,ll)=-one/sigma**2
         if (lambda >2) cc(3,ll)=zero
       else if (shape_type==2) then
         cc(3,ll)=-third*pi_over_rshp**2
       else if (shape_type==3) then
         cc(3,ll)=-half*(alpha(1,ll)*qq(1,ll)**(ll+1) &
&         +alpha(2,ll)*qq(2,ll)**(ll+1))/ffact(ll+1)
       end if
       cc(3,ll)=cc(3,ll)*pawtab%gnorm(ll)
     end if
   end do
 end if

!Y_lm(r) calculation
!==========================================================
 normchoice=1 ; option=max(optgr0,2*optgr1,3*optgr2)
 if(compute_gr0)  then
   ABI_ALLOCATE(ylmr,(l_size**2,nfgd))
 end if
 if(compute_gr1.and.(.not.compute_gr2))  then
   ABI_ALLOCATE(ylmrgr,(3,l_size**2,nfgd))
 end if
 if(compute_gr2)  then
   ABI_ALLOCATE(ylmrgr,(9,l_size**2,nfgd))
 end if
 if (compute_gr0.and.(.not.compute_gr1).and.(.not.compute_gr2)) then
   call initylmr(l_size,normchoice,nfgd,rnrm,option,rfgd,ylmr)
 else
   call initylmr(l_size,normchoice,nfgd,rnrm,option,rfgd,ylmr,ylmrgr)
 end if

!g_l(r) calculation (and factors for derivatives)
!==========================================================
 if (compute_gr0)  then
   ABI_ALLOCATE(gfact,(nfgd,0:l_size-1))
 end if
 if (compute_gr1)  then
   ABI_ALLOCATE(dgfact,(nfgd,0:l_size-1))
 end if
 if (compute_gr2)  then
   ABI_ALLOCATE(d2gfact,(nfgd,0:l_size-1))
 end if
 if(compute_gr1) then
   ABI_ALLOCATE(rnrm_inv,(nfgd))
   do ic=1,nfgd
     if (ic/=izero) rnrm_inv(ic)=one/rnrm(ic)
   end do
   if (izero>0) rnrm_inv(izero)=zero
 end if

!----- type -1 -----
 if (shape_type==-1) then
   if (compute_gr0) then
     do ll=0,l_size-1
       gfact(1:nfgd,ll)=shpfuncnum(1:nfgd,ll+1)
     end do
   end if
   if (compute_gr1) then
     do ll=0,l_size-1
       dgfact(1:nfgd,ll)=dshpfuncnum(1:nfgd,ll+1)*rnrm_inv(1:nfgd)
     end do
   end if
   if(compute_gr2) then
     do ll=0,l_size-1
       d2gfact(1:nfgd,ll)=(d2shpfuncnum(1:nfgd,ll+1)-dgfact(1:nfgd,ll))*rnrm_inv(1:nfgd)**2
     end do
   end if

!  ----- type 1 or 2 -----
 else if (shape_type==1.or.shape_type==2) then
   if (optgr0==1.and.optgr1==0.and.optgr2==0) then
     if (shape_type==1) then
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc1_0(arg)
         else
           gfact(ic,0)=shapefunc1(arg)
         end if
       end do
     else ! shape_type==2
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc2_0(arg)
         else
           gfact(ic,0)=shapefunc2(arg)
         end if
       end do
     end if
   else if (optgr1==1.and.optgr2==0) then
     if (shape_type==1) then
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc1_0(arg)
           if (lambda==2) then
             dgfact(ic,0)=dshpfunc1_ovr_0_2(arg)
           else ! lambda>2
             dgfact(ic,0)=dshpfunc1_ovr_0(arg)
           end if
         else
           gfact(ic,0)=shapefunc1(arg)
           dgfact(ic,0)=dshpfunc1(arg)*rnrm_inv(ic)
         end if
       end do
     else ! shape_type==2
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc2_0(arg)
           dgfact(ic,0)=dshpfunc2_ovr_0(arg)
         else
           gfact(ic,0)=shapefunc2(arg)
           dgfact(ic,0)=dshpfunc2(arg)*rnrm_inv(ic)
         end if
       end do
     end if
   else if (optgr2==1) then
     if (shape_type==1) then
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc1_0(arg)
           if (lambda==2) then
             dgfact(ic,0)=dshpfunc1_ovr_0_2(arg)
             d2gfact(ic,0)=d2shpfunc1_ovr2_0_2(arg)
           else if (lambda==3) then
             dgfact(ic,0)=dshpfunc1_ovr_0(arg)
             if (ic/=izero) then
               d2gfact(ic,0)=d2shpfunc1_ovr2_0_3(arg)
             else
               d2gfact(ic,0)=zero ! Diverging case
             end if
           else if (lambda==4) then
             dgfact(ic,0)=dshpfunc1_ovr_0(arg)
             d2gfact(ic,0)=d2shpfunc1_ovr2_0_4(arg)
           else ! lambda>4
             dgfact(ic,0)=dshpfunc1_ovr_0(arg)
             d2gfact(ic,0)=d2shpfunc1_ovr2_0(arg)
           end if
         else
           gfact(ic,0)=shapefunc1(arg)
           dgfact(ic,0)=dshpfunc1(arg)*rnrm_inv(ic)
           d2gfact(ic,0)=(d2shpfunc1(arg)-dgfact(ic,0))*rnrm_inv(ic)**2
         end if
       end do
     else ! shape_type==2
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc2_0(arg)
           dgfact(ic,0)=dshpfunc2_ovr_0(arg)
           d2gfact(ic,0)=d2shpfunc2_ovr2_0(arg)
         else
           gfact(ic,0)=shapefunc2(arg)
           dgfact(ic,0)=dshpfunc2(arg)*rnrm_inv(ic)
           d2gfact(ic,0)=(d2shpfunc2(arg)-dgfact(ic,0))*rnrm_inv(ic)**2
         end if
       end do
     end if
   end if
   if (l_size>1) then
     if (compute_gr0) then
       do ll=1,l_size-1
         do ic=1,nfgd
           gfact(ic,ll)=pawtab%gnorm(ll+1)*gfact(ic,0)*(rnrm(ic)**ll)
         end do
       end do
     end if
     if (compute_gr1) then
       do ic=1,nfgd
         dgfact(ic,1)=pawtab%gnorm(2)*(dgfact(ic,0)*rnrm(ic)+rnrm_inv(ic)*gfact(ic,0))
       end do
     end if
   end if
   if (l_size>2) then
     if (compute_gr1) then
       do ic=1,nfgd
         dgfact(ic,2)=pawtab%gnorm(3)*(dgfact(ic,0)*rnrm(ic)**2+two*gfact(ic,0))
       end do
     end if
     if (compute_gr2) then
       do ic=1,nfgd
         d2gfact(ic,2)=pawtab%gnorm(3)*(rnrm(ic)**2*d2gfact(ic,0)+four*dgfact(ic,0))
       end do
     end if
   end if
   if (l_size>3) then
     do ll=3,l_size-1
       if (compute_gr1) then
         do ic=1,nfgd
           dgfact(ic,ll)=pawtab%gnorm(ll+1)*(dgfact(ic,0)*rnrm(ic)**ll+ll*gfact(ic,0)*rnrm(ic)**(ll-2))
         end do
       end if
       if (compute_gr2) then
         do ic=1,nfgd
           d2gfact(ic,ll)=pawtab%gnorm(ll+1)*rnrm_inv(ic)**2*(rnrm(ic)**(ll+2)*d2gfact(ic,0) &
&           +two*ll*rnrm(ic)**ll*dgfact(ic,0)+ll*(ll-2)*rnrm(ic)**(ll-2)*gfact(ic,0))
         end do
       end if
     end do
   end if
   if (compute_gr0) gfact(:,0)=gfact(:,0)*pawtab%gnorm(1)
   if (compute_gr1) dgfact(:,0)=dgfact(:,0)*pawtab%gnorm(1)
   if (compute_gr2) d2gfact(:,0)=d2gfact(:,0)*pawtab%gnorm(1)

!  ----- type 3 -----
 else if (shape_type==3) then
   if (optgr0==1.and.optgr1==0.and.optgr2==0) then
     do ll=0,l_size-1
       do ic=1,nfgd
         call jbessel(jbes1,jbesp1,jbespp1,ll,0,qq(1,1+ll)*rnrm(ic))
         call jbessel(jbes2,jbesp2,jbespp2,ll,0,qq(2,1+ll)*rnrm(ic))
         gfact(ic,ll)=shapefunc3(jbes1,jbes2,ll)
       end do
     end do
   else if (optgr1==1.and.optgr2==0) then
     do ll=0,l_size-1
       do ic=1,nfgd
         call jbessel(jbes1,jbesp1,jbespp1,ll,1,qq(1,1+ll)*rnrm(ic))
         call jbessel(jbes2,jbesp2,jbespp2,ll,1,qq(2,1+ll)*rnrm(ic))
         gfact(ic,ll)=shapefunc3(jbes1,jbes2,ll)
         dgfact(ic,ll)=dshpfunc3(jbesp1,jbesp2,ll)*rnrm_inv(ic)
       end do
     end do
     if (izero>0.and.l_size>=1)  dgfact(izero,0)=-third*(alpha(1,1)*qq(1,1)+alpha(2,1)*qq(2,1))
     if (izero>0.and.l_size>=3)  dgfact(izero,2)=two/15._dp*(alpha(1,1)*qq(1,1)+alpha(2,1)*qq(2,1))
!    Note: for l=1, dgfact is diverging - d2gfact is diverging for l<4
   else if (optgr2==1) then
     do ll=0,l_size-1
       do ic=1,nfgd
         call jbessel(jbes1,jbesp1,jbespp1,ll,2,qq(1,1+ll)*rnrm(ic))
         call jbessel(jbes2,jbesp2,jbespp2,ll,2,qq(2,1+ll)*rnrm(ic))
         gfact(ic,ll)=shapefunc3(jbes1,jbes2,ll)
         dgfact(ic,ll)=dshpfunc3(jbesp1,jbesp2,ll)*rnrm_inv(ic)
         d2gfact(ic,ll)=(d2shpfunc3(jbespp1,jbespp2,ll)-dgfact(ic,ll))*rnrm_inv(ic)**2
       end do
     end do
     if (izero>0.and.l_size>=1)  dgfact(izero,0)=-third*(alpha(1,1)*qq(1,1)+alpha(2,1)*qq(2,1))
     if (izero>0.and.l_size>=3)  dgfact(izero,2)=two/15._dp*(alpha(1,1)*qq(1,1)+alpha(2,1)*qq(2,1))
!    Note: for l=1, dgfact is diverging - d2gfact is diverging for l<4
   end if
 end if

!g_l(r)*Y_lm(r) calculation
!==========================================================
 if (optgr0==1) then

   do ll=0,l_size-1
     do ilm=ll**2+1,min((ll+1)**2,lm_size)
       do ic=1,nfgd
         gylm(ic,ilm)=gfact(ic,ll)*ylmr(ilm,ic)
       end do
     end do
   end do

!  Special value at r=0  (supposing shapefunc(r)->C.r**l when r->0)
   if (izero>0) then
     gylm(izero,1:lm_size)=zero
     if (lm_size>=1) gylm(izero,1)=ylmr(1,izero)*cc(1,1)
   end if

 end if

!d/dr{g_l(r)*Y_lm(r)} calculation
!==========================================================
 if(optgr1==1) then

   do ll=0,l_size-1
     do ilm=ll**2+1,min((ll+1)**2,lm_size)
       do ic=1,nfgd
         gylmgr(1:3,ic,ilm)=gfact(ic,ll)*ylmrgr(1:3,ilm,ic)&
&         +dgfact(ic,ll)*rfgd(1:3,ic)*ylmr(ilm,ic)
       end do
     end do
   end do

!  Special values at r=0  (supposing shapefunc(r)->C.r**l when r->0)
   if (izero>0) then
     gylmgr(1:3,izero,1:lm_size)=zero
     if (lm_size>=1) then
       arg=cc(2,1)/sqrt(four_pi)
       gylmgr(1:3,izero,1)=arg
     end if
     if (lm_size>=2) then
       arg=cc(1,2)*sqrt(three/four_pi)
       gylmgr(2,izero,2)=arg
       if (lm_size>=3) gylmgr(3,izero,3)=arg
       if (lm_size>=4) gylmgr(1,izero,4)=arg
     end if
   end if

 end if

!d2/dridrj{g_l(r)*Y_lm(r)} calculation
!==========================================================
 if(optgr2==1) then

   do ll=0,l_size-1
     do ilm=ll**2+1,min((ll+1)**2,lm_size)
       do ic=1,nfgd
         gylmgr2(1,ic,ilm)=gfact(ic,ll)*ylmrgr(4,ilm,ic) &
&         +dgfact(ic,ll)*(ylmr(ilm,ic)+two*rfgd(1,ic)*ylmrgr(1,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(1,ic)*rfgd(1,ic)
         gylmgr2(2,ic,ilm)=gfact(ic,ll)*ylmrgr(5,ilm,ic) &
&         +dgfact(ic,ll)*(ylmr(ilm,ic)+two*rfgd(2,ic)*ylmrgr(2,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(2,ic)*rfgd(2,ic)
         gylmgr2(3,ic,ilm)=gfact(ic,ll)*ylmrgr(6,ilm,ic) &
&         +dgfact(ic,ll)*(ylmr(ilm,ic)+two*rfgd(3,ic)*ylmrgr(3,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(3,ic)*rfgd(3,ic)
         gylmgr2(4,ic,ilm)=gfact(ic,ll)*ylmrgr(7,ilm,ic) &
&         +dgfact(ic,ll)*(rfgd(3,ic)*ylmrgr(2,ilm,ic)+rfgd(2,ic)*ylmrgr(3,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(3,ic)*rfgd(2,ic)
         gylmgr2(5,ic,ilm)=gfact(ic,ll)*ylmrgr(8,ilm,ic) &
&         +dgfact(ic,ll)*(rfgd(3,ic)*ylmrgr(1,ilm,ic)+rfgd(1,ic)*ylmrgr(3,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(3,ic)*rfgd(1,ic)
         gylmgr2(6,ic,ilm)=gfact(ic,ll)*ylmrgr(9,ilm,ic) &
&         +dgfact(ic,ll)*(rfgd(1,ic)*ylmrgr(2,ilm,ic)+rfgd(2,ic)*ylmrgr(1,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(1,ic)*rfgd(2,ic)
       end do
     end do
   end do

!  Special values at r=0  (supposing shapefunc(r)->C.r**l when r->0)
   if (izero>0) then
     gylmgr2(1:6,izero,1:lm_size)=zero
     if (lm_size>=1) then
       arg=cc(3,1)/sqrt(pi)
       gylmgr2(1:3,izero,1)=arg
     end if
     if (lm_size>=2) then
       arg=cc(2,2)*sqrt(three/four_pi)
       gylmgr2(2,izero,2)=two*arg
       gylmgr2(4,izero,2)=    arg
       if (lm_size>=3) then
         gylmgr2(1,izero,3)=two*arg
         gylmgr2(3,izero,3)=two*arg
       end if
       if (lm_size>=4) then
         gylmgr2(5,izero,4)=arg
         gylmgr2(6,izero,4)=arg
       end if
     end if
     if (lm_size>=5) then
       arg=cc(1,3)*sqrt(15._dp/four_pi)
       gylmgr2(6,izero,5)=arg
       if (lm_size>=6) gylmgr2(4,izero,6)=arg
       if (lm_size>=7) then
         gylmgr2(1,izero,7)=   -arg/sqrt3
         gylmgr2(2,izero,7)=   -arg/sqrt3
         gylmgr2(3,izero,7)=two*arg/sqrt3
       end if
       if (lm_size>=8) gylmgr2(5,izero,8)=arg
       if (lm_size>=9) then
         gylmgr2(1,izero,9)= arg
         gylmgr2(2,izero,9)=-arg
       end if
     end if
   end if

 end if

!Memory deallocation
!==========================================================
 ABI_DEALLOCATE(rnrm)
 if (compute_gr0)  then
   ABI_DEALLOCATE(gfact)
 end if
 if (compute_gr1)  then
   ABI_DEALLOCATE(dgfact)
 end if
 if (compute_gr2)  then
   ABI_DEALLOCATE(d2gfact)
 end if
 if (compute_gr1)  then
   ABI_DEALLOCATE(rnrm_inv)
 end if
 if (shape_type==3)  then
   ABI_DEALLOCATE(alpha)
   ABI_DEALLOCATE(qq)
 end if
 if (compute_gr0)  then
   ABI_DEALLOCATE(ylmr)
 end if
 if (compute_gr1)  then
   ABI_DEALLOCATE(ylmrgr)
 end if
 if (shape_type==-1) then
   if (compute_gr0)  then
     ABI_DEALLOCATE(shpfuncnum)
   end if
   if (compute_gr1)  then
     ABI_DEALLOCATE(dshpfuncnum)
   end if
   if (compute_gr2)  then
     ABI_DEALLOCATE(d2shpfuncnum)
   end if
 end if

 DBG_EXIT("COLL")

end subroutine pawgylm
!!***
