!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkffnl
!! NAME
!! mkffnl
!!
!! FUNCTION
!! Make FFNL, nonlocal form factors, for each type of atom up to ntypat
!! and for each angular momentum.
!! When Legendre polynomials are used in the application of the
!!   nonlocal operator, FFNLs depend on (l,n) components; in this
!!   case, form factors are real and divided by |k+G|^l;
!! When spherical harmonics are used, FFNLs depend on (l,m,n)
!!   components; in this case, form factors are multiplied by Ylm(k+G).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MT, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dimekb=second dimension of ekb (see ekb)
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  ekb(dimekb,ntypat*(1-usepaw))=(Real) Kleinman-Bylander energies (hartree)
!!                                ->NORM-CONSERVING PSPS ONLY
!!  ffspl(mqgrid,2,lnmax,ntypat)=form factors and spline fit to 2nd derivative
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  ider=0=>no derivative wanted; 1=>1st derivative wanted; 2=>1st and 2nd derivatives wanted
!!  idir=ONLY WHEN YLMs ARE USED:
!!       When 1st derivative has to be computed:  (see more info below)
!!       - Determine the direction(s) of the derivatives(s)
!!       - Determine the set of coordinates (reduced or cartesians)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
!!                                                     or i=lmn (if useylm=1)
!!  kg(3,npw)=integer coordinates of planewaves in basis sphere for this k point.
!!  kpg(npw,nkpg)= (k+G) components (only if useylm=1)
!!  kpt(3)=reduced coordinates of k point
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=size of q (or |G|) grid for f(q)
!!  nkpg=second dimension of kpg_k (0 if useylm=0)
!!  npw=number of planewaves in basis sphere
!!  ntypat=number of types of atoms
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  pspso(ntypat)=spin-orbit characteristics for each atom type (1, 2, or 3)
!!  qgrid(mqgrid)=uniform grid of q values from 0 to qmax
!!  rmet(3,3)=real space metric (bohr**2)
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  ylm   (npw,mpsang*mpsang*useylm)=real spherical harmonics for each G and k point
!!  ylm_gr(npw,3,mpsang*mpsang*useylm)=gradients of real spherical harmonics wrt (k+G)
!!
!! OUTPUT
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=described below
!!
!! NOTES
!!  Uses spline fit ffspl provided by Numerical Recipes spline subroutine.
!!  Form factor $f_l(q)$ is defined by
!!   \begin{equation}
!!  \textrm{f}_l(q)=\frac{1}{dvrms} \int_0^\infty [j_l(2 \pi r q) u_l(r) dV(r) r dr]
!!   \end{equation}
!!   where u_l(r)=reference state wavefunction, dV(r)=nonlocal psp
!!   correction, j_l(arg)=spherical Bessel function for angular momentum l,
!!   and
!!   \begin{equation}
!!    \textrm{dvrms} =  \int_0^\infty [(u_l(r) dV(r))^2 dr])^{1/2}
!!   \end{equation}
!!   which is square root of mean square dV, i.e.
!!     $ (\langle (dV)^2 \rangle)^{1/2} $ .
!!   This routine is passed f_l(q) in spline form in the array ffspl and then
!!   constructs the values of $f_l(q)$ on the relevant (k+G) in array ffnl.
!!   The evaluation of the integrals defining ffspl was done in mkkbff.
!!
!!  Delivers the following (for each atom type t, or itypat):
!!   --------------------------
!!   Using Legendre polynomials in the application of nl operator:
!!     ffnl are real.
!!     ffnl(ig,1,(l,0,n),itypat) $= f_ln(k+G)/|k+G|^l $
!!     === if ider>=1
!!       ffnl(ig,2,(l,0,n),itypat) $=(fprime_ln(k+G)-l*f_ln(k+G)/|k+G|)/|k+G|^(l+1) $
!!     === if ider==2
!!       ffnl(ig,3,(l,0,n),itypat) $=(fprimeprime_ln(k+G)-(2l+1)*fprime_ln(k+G)/|k+G|
!!                                   +l(l+2)*f_ln(k+G)/|k+G|**2)/|k+G|^(l+2)
!!   --------------------------
!!   Using spherical harmonics in the application of nl operator:
!!     ffnl are real (we use REAL spherical harmonics).
!!     ffnl(ig,1,(l,m,n),itypat) = ffnl_1
!!                              $= f_ln(k+G) * Y_lm(k+G) $
!!     === if ider>=1
!!     --if (idir==0)
!!       ffnl(ig,1+i,(l,m,n),itypat) = dffnl_i = 3 reduced coord. of d(ffnl_1)/dK^cart
!!         $= fprime_ln(k+G).Y_lm(k+G).(k+G)^red_i/|k+G|+f_ln(k+G).(dY_lm/dK^cart)^red_i $
!!         for i=1..3
!!     --if (0<idir<4)
!!       ffnl(ig,2,(l,m,n),itypat)= cart. coordinate idir of d(ffnl_1)/dK^red
!!                                = Sum_(mu,nu) [ Gprim(mu,idir) Gprim(mu,nu) dffnl_nu ]
!!     --if (idir==4)
!!       ffnl(ig,1+i,(l,m,n),itypat)= 3 cart. coordinates of d(ffnl_1)/dK^red
!!                                  = Sum_(mu,nu) [ Gprim(mu,i) Gprim(mu,nu) dffnl_nu ]
!!     --if (-7<idir<0)
!!       ffnl(ig,2,(l,m,n),itypat)=1/2 [d(ffnl)/dK^cart_mu K^cart_nu + d(ffnl)/dK^cart_nu K^cart_mu]
!!                                with d(ffnl)/dK^cart_i = Sum_nu [ Gprim(nu,i) dffnl_nu ]
!!                                for |idir|->(mu,nu) (1->11,2->22,3->33,4->32,5->31,6->21)
!!     --if (idir==-7)
!!       ffnl(ig,3:8,(l,m,n),itypat)=1/2 [d(ffnl)/dK^cart_mu K^cart_nu + d(ffnl)/dK^cart_nu K^cart_mu]
!!                                with d(ffnl)/dK^cart_i = Sum_nu [ Gprim(nu,i) dffnl_nu ]
!!                                for all (mu,nu) (6 independant terms)
!!     === if ider==2
!!     --if (idir==0)
!!       ffnl(ig,2+i,(l,m,n),itypat) = d2ffnl_i = 6 reduced coord. of d2(ffnl_1)/dK^cart.dK^cart
!!        for all i=(mu,nu) (6 independant terms)
!!   --------------------------
!!
!! NOTES
!!  1) l may be 0, 1, 2, or 3 in this version.
!!
!!  2) Norm-conserving psps : only FFNL for which ekb is not zero are calculated.
!!
!!  3) Each expression above approaches a constant as $|k+G| \rightarrow 0 $.
!!     In the cases where $|k+G|$ is in the denominator, there is always a
!!     factor of $(k+G)_mu$ multiplying the ffnl term where it is actually used,
!!     so that we may replace the ffnl term by any constant when $|k+G| = 0$.
!!     Below we replace 1/0 by 1/tol10, thus creating an arbitrary constant
!!     which will later be multiplied by 0.
!!
!! PARENTS
!!      ctocprj,dyfnl3,eltfrnl3,energy,forstrnps,ks_ddiago,ladielmt,lavnl
!!      m_commutator_vkbr,m_cprj_bspline,m_shirley,m_wfs,nstpaw3,nstwf3,nstwf4
!!      prctfvw1,prctfvw2,resp3dte,rhofermi3,update_mmat,vso_realspace_nonlop
!!      vtorho,vtorho3
!!
!! CHILDREN
!!      leave_new,mkkin,splfit,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkffnl(dimekb,dimffnl,ekb,ffnl,ffspl,gmet,gprimd,ider,idir,indlmn,&
&                   kg,kpg,kpt,lmnmax,lnmax,mpsang,mqgrid,nkpg,npw,ntypat,pspso,&
&                   qgrid,rmet,usepaw,useylm,ylm,ylm_gr)

 use m_profiling

 use defs_basis
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkffnl'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimekb,dimffnl,ider,idir,lmnmax,lnmax,mpsang,mqgrid,nkpg
 integer,intent(in) :: npw,ntypat,usepaw,useylm
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),kg(3,npw),pspso(ntypat)
 real(dp),intent(in) :: ekb(dimekb,ntypat*(1-usepaw))
 real(dp),intent(in) :: ffspl(mqgrid,2,lnmax,ntypat),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: kpg(npw,nkpg),kpt(3),qgrid(mqgrid),rmet(3,3)
 real(dp),intent(in) :: ylm(npw,mpsang*mpsang*useylm)
 real(dp),intent(in) :: ylm_gr(npw,3+6*(ider/2),mpsang*mpsang*useylm)
 real(dp),intent(out) :: ffnl(npw,dimffnl,lmnmax,ntypat)

!Local variables-------------------------------
!scalars
 integer :: ider_tmp,iffnl,ig,il,ilm,ilmn,iln,iln0,itypat,mu,mua,mub,nlmn,nu
 real(dp),parameter :: renorm_factor=0.5d0/pi**2
 real(dp) :: ecut,ecutsm,effmass,kpg1,kpg2,kpg3,kpgc1,kpgc2,kpgc3,rmetab
 logical :: testnl=.false.
 character(len=500) :: message
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 real(dp) :: tsec(2)
 real(dp),allocatable :: dffnl_cart(:,:),dffnl_red(:,:),dffnl_tmp(:),kpgc(:,:)
 real(dp),allocatable :: kpgn(:,:),kpgnorm(:),kpgnorm_inv(:),wk_ffnl1(:)
 real(dp),allocatable :: wk_ffnl2(:),wk_ffnl3(:),wk_ffspl(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' mkffnl : enter'
!return
!ENDDEBUG

!Keep track of time spent in mkffnl
 call timab(16,1,tsec)

!Compatibility tests
 if (mpsang>4) then
   write(message, '(a,a,a,a,i10,a,a)' ) ch10,&
&   ' mkffnl : BUG -',ch10,&
&   '  Called with mpsang > 4, =',mpsang,ch10,&
&   '  This  subroutine will not accept lmax+1 > 4.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if
 if (idir<-7.or.idir>4) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' mkffnl : BUG -',ch10,&
&   '  Called with idir<-6 or idir>4 !'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if
 if (useylm==0) then
   iffnl=1+ider
 else
   iffnl=1
   if (ider>=1) then
     if (idir==0) iffnl=iffnl+3
     if (idir/=0) iffnl=iffnl+1
     if (idir==4) iffnl=iffnl+2
     if (idir==-7) iffnl=iffnl+5
   end if
   if (ider==2) then
     if (idir==0) iffnl=iffnl+6
   end if
 end if
 if (iffnl/=dimffnl) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' mkffnl : BUG -',ch10,&
&   '  Incompatibility between ider, idir and dimffnl !'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

!Get (k+G) and |k+G|:
 ABI_ALLOCATE(kpgnorm,(npw))
 ABI_ALLOCATE(kpgnorm_inv,(npw))
 if (useylm==1) then
   ABI_ALLOCATE(kpgc,(npw,3))
   if (ider>=1) ABI_ALLOCATE(kpgn,(npw,3))
   if (nkpg<3) then
!    $OMP PARALLEL DO PRIVATE(ig) &
!    $OMP&SHARED(kpt,kg,kpgc,kpgnorm,kpgnorm_inv,gprimd,npw,tol10)
     do ig=1,npw
       kpg1=kpt(1)+dble(kg(1,ig));kpg2=kpt(2)+dble(kg(2,ig));kpg3=kpt(3)+dble(kg(3,ig))
       kpgc1=kpg1*gprimd(1,1)+kpg2*gprimd(1,2)+kpg3*gprimd(1,3)
       kpgc2=kpg1*gprimd(2,1)+kpg2*gprimd(2,2)+kpg3*gprimd(2,3)
       kpgc3=kpg1*gprimd(3,1)+kpg2*gprimd(3,2)+kpg3*gprimd(3,3)
       kpgc(ig,1)=kpgc1;kpgc(ig,2)=kpgc2;kpgc(ig,3)=kpgc3
       kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
       if (ider>=1) then
         kpgnorm_inv(ig)=1.d0/max(kpgnorm(ig),tol10)
         kpgn(ig,1)=kpg1*kpgnorm_inv(ig)
         kpgn(ig,2)=kpg2*kpgnorm_inv(ig)
         kpgn(ig,3)=kpg3*kpgnorm_inv(ig)
       end if
     end do
!    $OMP END PARALLEL DO
   else
!    $OMP PARALLEL DO PRIVATE(ig) &
!    $OMP&SHARED(kpgc,kpgnorm,kpgnorm_inv,gprimd,npw,tol10)
     do ig=1,npw
       kpgc1=kpg(ig,1)*gprimd(1,1)+kpg(ig,2)*gprimd(1,2)+kpg(ig,3)*gprimd(1,3)
       kpgc2=kpg(ig,1)*gprimd(2,1)+kpg(ig,2)*gprimd(2,2)+kpg(ig,3)*gprimd(2,3)
       kpgc3=kpg(ig,1)*gprimd(3,1)+kpg(ig,2)*gprimd(3,2)+kpg(ig,3)*gprimd(3,3)
       kpgc(ig,1)=kpgc1;kpgc(ig,2)=kpgc2;kpgc(ig,3)=kpgc3
       kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
       if (ider>=1) then
         kpgnorm_inv(ig)=1.d0/max(kpgnorm(ig),tol10)
         kpgn(ig,1:3)=kpg(ig,1:3)*kpgnorm_inv(ig)
       end if
     end do
!    $OMP END PARALLEL DO
   end if
 else
   if (nkpg<3) then
     ecut=huge(0.0d0)*0.1d0;ecutsm=zero;effmass=one
!    Note that with ecutsm=0, the right kinetic energy is computed
     call mkkin(ecut,ecutsm,effmass,gmet,kg,kpgnorm,kpt,npw)
!    $OMP PARALLEL DO PRIVATE(ig) &
!    $OMP&SHARED(renorm_factor,npw,kpgnorm,kpgnorm_inv,tol10)
     do ig=1,npw
       kpgnorm(ig)=sqrt(renorm_factor*kpgnorm(ig))
       kpgnorm_inv(ig)=1.d0/max(kpgnorm(ig),tol10)
     end do
!    $OMP END PARALLEL DO
   else
!    $OMP PARALLEL DO PRIVATE(ig) &
!    $OMP&SHARED(kpg,kpgnorm,kpgnorm_inv,gprimd,npw,tol10)
     do ig=1,npw
       kpgc1=kpg(ig,1)*gprimd(1,1)+kpg(ig,2)*gprimd(1,2)+kpg(ig,3)*gprimd(1,3)
       kpgc2=kpg(ig,1)*gprimd(2,1)+kpg(ig,2)*gprimd(2,2)+kpg(ig,3)*gprimd(2,3)
       kpgc3=kpg(ig,1)*gprimd(3,1)+kpg(ig,2)*gprimd(3,2)+kpg(ig,3)*gprimd(3,3)
       kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
       kpgnorm_inv(ig)=1.d0/max(kpgnorm(ig),tol10)
     end do
!    $OMP END PARALLEL DO
   end if
 end if
 ABI_ALLOCATE(wk_ffnl1,(npw))
 ABI_ALLOCATE(wk_ffnl2,(npw))
 ABI_ALLOCATE(wk_ffnl3,(npw))
 ABI_ALLOCATE(wk_ffspl,(mqgrid,2))
 if (ider>=1.and.useylm==1) then
   ABI_ALLOCATE(dffnl_red,(npw,3))
   if (idir/=0)  then
     ABI_ALLOCATE(dffnl_cart,(npw,3))
   end if
   if (idir>0)   then
     ABI_ALLOCATE(dffnl_tmp,(npw))
   end if
 end if

!Loop over types of atoms
 do itypat=1,ntypat

!  Loop over (l,m,n) values
   iln0=0;nlmn=count(indlmn(3,:,itypat)>0)
   ffnl(:,:,:,itypat)=zero
   do ilmn=1,nlmn
     il=indlmn(1,ilmn,itypat)
     ilm =indlmn(4,ilmn,itypat)
     iln =indlmn(5,ilmn,itypat)
     iffnl=ilmn;if (useylm==0) iffnl=iln

!    Special case : spin-orbit calculation and no spin-orbit
!    contribution for this type of psp
!    ->no spin-orbit contribution for this type of psp
     if ((indlmn(6,ilmn,itypat)==1).or.(pspso(itypat)/=0)) then

!      Compute FFNL only if ekb>0 or paw
       if (usepaw==1) testnl=.true.
       if (usepaw==0) testnl=(abs(ekb(iln,itypat))>tol10)
       if (testnl) then

!        Store form factors (from ffspl)
!        -------------------------------
         if (iln>iln0) then
!          $OMP PARALLEL DO PRIVATE(ig) &
!          $OMP&SHARED(ffspl,iln,itypat,mqgrid,wk_ffspl)
           do ig=1,mqgrid
             wk_ffspl(ig,1)=ffspl(ig,1,iln,itypat)
             wk_ffspl(ig,2)=ffspl(ig,2,iln,itypat)
           end do
!          $OMP END PARALLEL DO
           ider_tmp=min(ider,1)
           call splfit(qgrid,wk_ffnl2,wk_ffspl,ider_tmp,kpgnorm,wk_ffnl1,mqgrid,npw)
           if(ider==2) then
             call splfit(qgrid,wk_ffnl3,wk_ffspl,ider,kpgnorm,wk_ffnl1,mqgrid,npw)
           end if
         end if

!        Store FFNL and FFNL derivatives
!        -------------------------------

!        =========================================================================
!        A-USE OF SPHER. HARMONICS IN APPLICATION OF NL OPERATOR:
!        ffnl(K,l,m,n)=fnl(K).Ylm(K)
!        --if (idir==0)
!        ffnl_prime(K,1:3,l,m,n)=3 reduced coordinates of d(ffnl)/dK^cart
!        =fnl_prime(K).Ylm(K).K^red_i/|K|+fnl(K).(dYlm/dK^cart)^red_i
!        --if (0<idir<4)
!        ffnl_prime(K,l,m,n)=cart. coordinate idir of d(ffnl)/dK^red
!        --if (idir==4)
!        ffnl_prime(K,l,m,n)=3 cart. coordinates of d(ffnl)/dK^red
!        --if (-7<=idir<0) - |idir|=(mu,nu) (1->11,2->22,3->33,4->32,5->31,6->21)
!        ffnl_prime(K,l,m,n)=1/2 [d(ffnl)/dK^cart_mu K^cart_nu + d(ffnl)/dK^cart_nu K^cart_mu]
!        ffnl_prime_prime(K,l,m,n)=6 reduced coordinates of d2(ffnl)/dK^cart.dK^cart
         if (useylm==1) then

!          $OMP PARALLEL DO PRIVATE(ig) &
!          $OMP&SHARED(ffnl,wk_ffnl1,ylm,npw,iffnl,itypat,ilm)
           do ig=1,npw
             ffnl(ig,1,iffnl,itypat)=ylm(ig,ilm)*wk_ffnl1(ig)
           end do
!          $OMP END PARALLEL DO
           if (ider>=1) then
             do mu=1,3
!              $OMP PARALLEL DO PRIVATE(ig) &
!              $OMP&SHARED(dffnl_red,wk_ffnl1,wk_ffnl2,ylm,ylm_gr,kpgn,npw,mu,ilm)
               do ig=1,npw
                 dffnl_red(ig,mu)=ylm(ig,ilm)*wk_ffnl2(ig)*kpgn(ig,mu)&
&                 +ylm_gr(ig,mu,ilm)*wk_ffnl1(ig)
               end do
!              $OMP END PARALLEL DO
             end do
             if (idir==0) then
               do mu=1,3
                 ffnl(:,1+mu,iffnl,itypat)=dffnl_red(:,mu)
               end do
             else
               dffnl_cart=zero
               do nu=1,3
                 do mu=1,3
!                  $OMP PARALLEL DO PRIVATE(ig) &
!                  $OMP&SHARED(dffnl_cart,dffnl_red,gprimd,npw,mu,nu)
                   do ig=1,npw
                     dffnl_cart(ig,mu)=dffnl_cart(ig,mu)+dffnl_red(ig,nu)*gprimd(mu,nu)
                   end do
!                  $OMP END PARALLEL DO
                 end do
               end do
               if (idir>0.and.idir<4) then
                 dffnl_tmp=zero
                 do nu=1,3
!                  $OMP PARALLEL DO PRIVATE(ig) &
!                  $OMP&SHARED(ffnl,dffnl_cart,gprimd,npw,mu,nu,iffnl,itypat,idir)
                   do ig=1,npw
                     dffnl_tmp(ig)=dffnl_tmp(ig)&
&                     +dffnl_cart(ig,nu)*gprimd(nu,idir)
                   end do
!                  $OMP END PARALLEL DO
                 end do
                 ffnl(:,2,iffnl,itypat)=dffnl_tmp(:)
               else if (idir==4) then
                 do mu=1,3
                   dffnl_tmp=zero
                   do nu=1,3
!                    $OMP PARALLEL DO PRIVATE(ig) &
!                    $OMP&SHARED(ffnl,dffnl_cart,gprimd,npw,mu,nu,iffnl,itypat,idir)
                     do ig=1,npw
                       dffnl_tmp(ig)=dffnl_tmp(ig)&
&                       +dffnl_cart(ig,nu)*gprimd(nu,mu)
                     end do
!                    $OMP END PARALLEL DO
                   end do
                   ffnl(:,1+mu,iffnl,itypat)=dffnl_tmp(:)
                 end do
               else if (idir/=-7) then
                 mu=abs(idir);mua=alpha(mu);mub=beta(mu)
!                $OMP PARALLEL DO PRIVATE(ig) &
!                $OMP&SHARED(ffnl,dffnl_cart,kpgc,npw,mua,mub,iffnl,itypat)
                 do ig=1,npw
                   ffnl(ig,2,iffnl,itypat)=0.5d0* &
&                   (dffnl_cart(ig,mua)*kpgc(ig,mub) &
&                   +dffnl_cart(ig,mub)*kpgc(ig,mua))
                 end do
!                $OMP END PARALLEL DO
               else if (idir==-7) then
                 do mu=1,6
                   mua=alpha(mu);mub=beta(mu)
!                  $OMP PARALLEL DO PRIVATE(ig) &
!                  $OMP&SHARED(ffnl,dffnl_cart,kpgc,npw,mua,mub,iffnl,itypat)
                   do ig=1,npw
                     ffnl(ig,1+mu,iffnl,itypat)=0.5d0* &
&                     (dffnl_cart(ig,mua)*kpgc(ig,mub) &
&                     +dffnl_cart(ig,mub)*kpgc(ig,mua))
                   end do
!                  $OMP END PARALLEL DO
                 end do
               end if
             end if
           end if
           if (ider==2.and.idir==0) then
             do mu=1,6
               mua=alpha(mu);mub=beta(mu)
               rmetab=rmet(mua,mub)
!              $OMP PARALLEL DO PRIVATE(ig) &
!              $OMP&SHARED(npw,mu,mua,mub,iffnl,itypat)
!              $OMP&SHARED(ffnl,kpgn,kpgnorm_inv,rmetab,ylm,ylm_gr,wf_ffnl1,wk_ffnl2,wk_ffnl3)
               do ig=1,npw
                 ffnl(ig,4+mu,iffnl,itypat)= &
&                 ylm_gr(ig,3+mu,ilm)*wk_ffnl1(ig) &
&                 + (rmetab-kpgn(ig,mua)*kpgn(ig,mub))*ylm(ig,ilm)*wk_ffnl2(ig)*kpgnorm_inv(ig) &
&                 + ylm(ig,ilm)*kpgn(ig,mua)*kpgn(ig,mub)*wk_ffnl3(ig) &
&                 + (ylm_gr(ig,mua,ilm)*kpgn(ig,mub)+ylm_gr(ig,mub,ilm)*kpgn(ig,mua))*wk_ffnl2(ig)
               end do
!              $OMP END PARALLEL DO
             end do
           end if

!          =========================================================================
!          B-USE OF LEGENDRE POLYNOMIAL IN APPLICATION OF NL OPERATOR:
!          ffnl(K,l,n)=fnl(K)/|K|^l
!          ffnl_prime(K,l,n)=(fnl_prime(K)-l*fnl(K)/|K|)/|K|^(l+1)
!          ffnl_prime_prime(K,l,n)=(fnl_prime_prime(K)-(2*l+1)*fnl_prime(K)/|K|
!          +l*(l+2)*fnl(K)/|K|^2)/|K|^(l+2)
         else if (iln>iln0) then

           if (il==0) then
!            $OMP PARALLEL DO PRIVATE(ig) &
!            $OMP&SHARED(ffnl,iffnl,itypat,npw,wk_ffnl1)
             do ig=1,npw
               ffnl(ig,1,iffnl,itypat)=wk_ffnl1(ig)
             end do
!            $OMP END PARALLEL DO
           else
!            $OMP PARALLEL DO PRIVATE(ig) &
!            $OMP&SHARED(ffnl,iffnl,il,itypat,npw,kpgnorm_inv,wk_ffnl1)
             do ig=1,npw
               ffnl(ig,1,iffnl,itypat)=wk_ffnl1(ig)*kpgnorm_inv(ig)**il
             end do
!            $OMP END PARALLEL DO
           end if
           if (ider>=1) then
!            $OMP PARALLEL DO PRIVATE(ig) &
!            $OMP&SHARED(ffnl,iffnl,il,itypat,npw,kpgnorm_inv,wk_ffnl1,wk_ffnl2,wk_ffnl3)
             do ig=1,npw
               ffnl(ig,2,iffnl,itypat)= (wk_ffnl2(ig)-&
&               dble(il)*wk_ffnl1(ig)*kpgnorm_inv(ig))&
&               *kpgnorm_inv(ig)**(il+1)
             end do
!            $OMP END PARALLEL DO
             if (ider==2) then
!              $OMP PARALLEL DO PRIVATE(ig) &
!              $OMP&SHARED(ffnl,iffnl,il,itypat,npw,kpgnorm_inv,wk_ffnl1,wk_ffnl2,wk_ffnl3)
               do ig=1,npw
                 ffnl(ig,3,iffnl,itypat)= (wk_ffnl3(ig)-&
&                 dble(2*il+1)*wk_ffnl2(ig)*kpgnorm_inv(ig)+&
&                 dble(il*(il+2))*wk_ffnl1(ig)*kpgnorm_inv(ig)**2)&
&                 *kpgnorm_inv(ig)**(il+2)
               end do
!              $OMP END PARALLEL DO
             end if
           end if

!          End if - Use of Ylm or not
         end if

!        End if - a nonlocal part exists
       end if

!      End if - special case : spin orbit calc. & no spin-orbit psp
     end if

!    End do - loop over (l,m,n) values
     if (iln>iln0) iln0=iln
   end do

!  End do - loop over atom types
 end do

!deallocate(kpgnorm,kpgnorm_inv,wk_ffnl1,wk_ffnl2,wk_ffnl3,wk_ffspl)
!if (useylm==1) deallocate(kpgc)
!if (ider>=1.and.useylm==1) then
!deallocate(kpgn,dffnl_red)
!if (idir/=0) deallocate(dffnl_cart)
!if (idir>0)  deallocate(dffnl_tmp)
!end if
 if (allocated(wk_ffspl))  then
   ABI_DEALLOCATE(wk_ffspl)
 end if
 if (allocated(kpgc))  then
   ABI_DEALLOCATE(kpgc)
 end if
 if (allocated(kpgn))  then
   ABI_DEALLOCATE(kpgn)
 end if
 if (allocated(dffnl_red))  then
   ABI_DEALLOCATE(dffnl_red)
 end if
 if (allocated(dffnl_cart))  then
   ABI_DEALLOCATE(dffnl_cart)
 end if
 if (allocated(dffnl_tmp))  then
   ABI_DEALLOCATE(dffnl_tmp)
 end if
 if (allocated(wk_ffnl3))  then
   ABI_DEALLOCATE(wk_ffnl3)
 end if
 if (allocated(wk_ffnl2))  then
   ABI_DEALLOCATE(wk_ffnl2)
 end if
 if (allocated(wk_ffnl1))  then
   ABI_DEALLOCATE(wk_ffnl1)
 end if
 if (allocated(kpgnorm_inv))  then
   ABI_DEALLOCATE(kpgnorm_inv)
 end if
 if (allocated(kpgnorm))  then
   ABI_DEALLOCATE(kpgnorm)
 end if


 call timab(16,2,tsec)

!DEBUG
!write(std_out,*)' mkffnl : exit'
!ENDDEBUG

end subroutine mkffnl
!!***
