!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxcm
!! NAME
!! pawxcm
!!
!! FUNCTION
!! PAW only
!! Start from the density or spin-density, and compute xc correlation
!! potential and energies inside a paw sphere.
!! LDA+GGA - USE A DEVELOPMENT OF THE DENSITY OVER (L,M) MOMENTS
!! Driver of XC functionals.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT, GJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!! This routine has been written from rhohxc
!!
!! INPUTS
!!  corexc(pawrad%mesh_size)=core density on radial grid
!!  exexch= choice of local exact exchange. Active if exexch=3
!!  ixc= choice of exchange-correlation scheme (see above and below)
!!  lm_size=size of density array rhor (see below)
!!  lmselect(lm_size)=select the non-zero LM-moments of input density rhor
!!  nhat(pawrad%mesh_size,lm_size,nspden)=compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nkxc=second dimension of the kxc array. If /=0, the exchange-correlation kernel must be computed
!!  nspden=number of spin-density components
!!  option=0 compute both XC energies (direct+double-counting) and potential
!!         1 compute only XC potential
!!         2 compute only XC energies (direct+double-counting)
!!         3 compute only XC energy by direct scheme
!!         4 compute only XC energy by direct scheme for spherical part of the density
!!         5 compute only XC potential for spherical part of the density
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawxcdev=order of Vxc development
!!  rhor(pawrad%mesh_size,lm_size,nspden)=electron density in real space in electrons/bohr**3
!!                                       (total in 1st half and spin-up in 2nd half if nspden=2)
!!  usecore= 1 if core density has to be used in Exc/Vxc ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in double counting energy term only
!!             2 if compensation density (nhat) has to be used in Exc/Vxc and double counting energy term
!!  xclevel= XC functional level
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!
!! OUTPUT
!!  == if option==0, 2, 3, or 4 ==
!!    enxc=returned exchange and correlation energy (hartree)
!!  == if option==0 or 2 ==
!!    enxcdc=returned exchange-cor. contribution to double-counting energy
!!  == if option==0 or 1 ==
!!    vxc(pawrad%mesh_size,lm_size,nspden)=xc potential
!!       (spin up in 1st half and spin-down in 2nd half if nspden=2)
!!  == if nkxc>0 ==
!!    kxc(pawrad%mesh_size,lm_size,nkxc)=xc kernel
!!        (see notes below for nkxc)
!!
!! NOTES
!!  Dimension of Kxc:
!!   ===== if LDA (xclevel=1) :
!!    if nspden==1: return kxc(:,1)= d2Exc/drho2
!!       that is 1/2 (d2Exc/drho_up drho_up + d2Exc/drho_up drho_dn)
!!    if nspden==1: also return kxc(:,2)= d2Exc/drho_up drho_dn
!!    if nspden>=2, return  kxc(:,1)=d2Exc/drho_up drho_up
!!                          kxc(:,2)=d2Exc/drho_up drho_dn
!!                          kxc(:,3)=d2Exc/drho_dn drho_dn
!!   ===== if GGA (xclevel=2) :
!!    Treat all cases as spin-polarized, with nkxc=23
!!    kxc(:,1)= d2Ex/drho_up drho_up
!!    kxc(:,2)= d2Ex/drho_dn drho_dn
!!    kxc(:,3)= dEx/d(abs(grad(rho_up))) / abs(grad(rho_up))
!!    kxc(:,4)= dEx/d(abs(grad(rho_dn))) / abs(grad(rho_dn))
!!    kxc(:,5)= d2Ex/d(abs(grad(rho_up))) drho_up / abs(grad(rho_up))
!!    kxc(:,6)= d2Ex/d(abs(grad(rho_dn))) drho_dn / abs(grad(rho_dn))
!!    kxc(:,7)= 1/abs(grad(rho_up)) * d/d(abs(grad(rho_up)) (dEx/d(abs(grad(rho_up))) /abs(grad(rho_up)))
!!    kxc(:,8)= 1/abs(grad(rho_dn)) * d/d(abs(grad(rho_dn)) (dEx/d(abs(grad(rho_dn))) /abs(grad(rho_dn)))
!!    kxc(:,9)= d2Ec/drho_up drho_up
!!    kxc(:,10)=d2Ec/drho_up drho_dn
!!    kxc(:,11)=d2Ec/drho_dn drho_dn
!!    kxc(:,12)=dEc/d(abs(grad(rho))) / abs(grad(rho))
!!    kxc(:,13)=d2Ec/d(abs(grad(rho))) drho_up / abs(grad(rho))
!!    kxc(:,14)=d2Ec/d(abs(grad(rho))) drho_dn / abs(grad(rho))
!!    kxc(:,15)=1/abs(grad(rho)) * d/d(abs(grad(rho)) (dEc/d(abs(grad(rho))) /abs(grad(rho)))
!!    kxc(:,16)=rho_up
!!    kxc(:,17)=rho_dn
!!    kxc(:,18)=gradx(rho_up)
!!    kxc(:,19)=gradx(rho_dn)
!!    kxc(:,20)=grady(rho_up)
!!    kxc(:,21)=grady(rho_dn)
!!    kxc(:,22)=gradz(rho_up)
!!    kxc(:,23)=gradz(rho_dn)
!!
!! PARENTS
!!      pawdenpot,psp7in
!!
!! CHILDREN
!!      mkdenpos,pawxcsph,pawxcsum,simp_gen,size_dvxc,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawxcm(corexc,enxc,enxcdc,exexch,ixc,kxc,lm_size,lmselect,nhat,nkxc,nspden,option,&
&                  pawang,pawrad,pawxcdev,rhor,usecore,usexcnhat,vxc,xclevel,xc_denpos)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_radmesh,          only : simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawxcm'
 use interfaces_18_timing
 use interfaces_56_xc, except_this_one => pawxcm
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exexch,ixc,lm_size,nkxc,nspden,option,pawxcdev,usecore
 integer,intent(in) :: usexcnhat,xclevel
 real(dp),intent(in) :: xc_denpos
 real(dp),intent(out) :: enxc,enxcdc
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size)
 real(dp),intent(in) :: corexc(pawrad%mesh_size)
 real(dp),intent(in) :: nhat(pawrad%mesh_size,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor(pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(out) :: kxc(pawrad%mesh_size,lm_size,nkxc)
 real(dp),intent(out) :: vxc(pawrad%mesh_size,lm_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ilm,ir,ir1,ir2,ispden,iwarn,jr,ndvxc,ngr2,ngrad,nrad
 integer :: nspden_updn,nspgrad,nvxcdgr,optv2,order
 integer :: nd2vxc
 real(dp),parameter :: delta=1.d-4
 real(dp) :: dvxc1,dvxc2,dvxc3,dvxc4,dvxca,dvxcb,dvxcc,dvxcd
 real(dp) :: fact,invsqfpi,invsqfpi2,m_norm,m_norm_min,sqfpi,sqfpi2
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: d1vxc(:,:),d2vxc(:,:),exc_(:),exci(:),ff(:),gg(:),m_norm_inv(:)
 real(dp),allocatable :: rho_(:,:),rho_updn(:,:,:),rhoinv(:,:),rhosph(:,:)
 real(dp),allocatable :: v0sum(:,:),v1sum(:,:),v2sum(:,:,:),vxc1(:,:),vxc2(:,:)
 real(dp),allocatable :: vxcdn1(:,:),vxcdn2(:,:),vxci(:,:)
!************************************************************************

 DBG_ENTER("COLL")

 call timab(81,1,tsec)

 if(nkxc>0) then
   msg=' Computation of Kxc not yet implemented (choose pawxcdev=0) !'
   MSG_ERROR(msg)
   kxc=zero
 end if
 if(nspden==4.and.nkxc>0) then
   msg=' Kxc for nspden=4 not implemented !'
   MSG_ERROR(msg)
 end if

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

!Arrays dimensions and constants
 iwarn=0
 order=1
 nspden_updn=min(nspden,2)
 ngrad=1;if (xclevel==2) ngrad=2 ! ngrad=1 is for LDAs or LSDs; ngrad=2 is for GGAs
 nspgrad=nspden_updn*ngrad;if(nspden_updn==2.and.ngrad==2) nspgrad=5
 nrad=pawrad%mesh_size
 sqfpi=sqrt(four_pi);sqfpi2=half*sqfpi
 invsqfpi=one/sqfpi;invsqfpi2=half*invsqfpi

!Compute sizes of arrays

!to be adjusted for the call
 nd2vxc=1
 call size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden_updn,nvxcdgr,order)

!Initializations of output arrays
 if (option/=1.and.option/=5) enxc=zero
 if (option==0.or.option==2) enxcdc=zero
 if (option/=3.and.option/=4) vxc(:,:,:)=zero

 if (xclevel==0) then ! No xc at all is applied (usually for testing)
   MSG_WARNING('  Note that no xc is applied (ixc=0).')
   return
 end if

!----------------------------------------------------------------------
!----- Build several densities
!----------------------------------------------------------------------

!rho_updn contains the effective density used for XC
!with core density and/or compensation density eventually included
!-----------------------------------------------------------------
 ABI_ALLOCATE(rho_updn,(nrad,lm_size,nspden))
 rho_updn(:,:,:)=rhor(:,:,:)
 if (usexcnhat==2) rho_updn(:,:,:)=rho_updn(:,:,:)+nhat(:,:,:)
 if (usecore==1) then
   if (nspden==1.or.nspden==4) then
     rho_updn(:,1,1)=rho_updn(:,1,1)+sqfpi*corexc(:)
   else if (nspden==2) then
     rho_updn(:,1,1)=rho_updn(:,1,1)+sqfpi*corexc(:)
     rho_updn(:,1,2)=rho_updn(:,1,2)+sqfpi2*corexc(:)
   end if
 end if

!In case of collinear magnetism, separate up and down contributions
 if (nspden==2) then
   ABI_ALLOCATE(ff,(nrad))
   do ilm=1,lm_size
     ff(:)=rho_updn(:,ilm,2)
     rho_updn(:,ilm,2)=rho_updn(:,ilm,1)-ff(:)
     rho_updn(:,ilm,1)=ff(:)
   end do
   ABI_DEALLOCATE(ff)
 end if

!rhoSPH contains the spherical part of effective density
!(including Y00 spherical harmonic)
!-----------------------------------------------------------------
 ABI_ALLOCATE(rhosph,(nrad,nspden_updn))

!Non-magnetic system: rhoSPH(;,1)=(1/2).rhoSPH_total
 if (nspden==1) then
   rhosph(:,1)=rho_updn(:,1,1)*invsqfpi2

!  Collinear magnetism: rhoSPH = (rhoSPH_up, rhoSPH_dn)
 else if (nspden==2) then
   rhosph(:,1:2)=rho_updn(:,1,1:2)*invsqfpi

!  Non-collinear magnetism: rhoSPH = (rhoSPH_up, rhoSPH_dn)
!  obtained by rotating rho_updn
 else if (nspden==4) then
   ABI_ALLOCATE(m_norm_inv,(nrad))
   m_norm_min=EPSILON(0.0_dp)**2
   do ir=1,nrad
     m_norm=sqrt(rho_updn(ir,1,2)**2+rho_updn(ir,1,3)**2+rho_updn(ir,1,4)**2)
     rhosph(ir,1)=(rho_updn(ir,1,1)+m_norm)*invsqfpi2
     rhosph(ir,2)=(rho_updn(ir,1,1)-m_norm)*invsqfpi2
     if (m_norm>m_norm_min) then
!      if (m_norm>abs(rho_updn(ir,1,1))*tol10+tol14) then
       m_norm_inv(ir)=one/m_norm
     else
       m_norm_inv(ir)=zero
     end if
   end do
 end if

!Make spherical density positive
 call mkdenpos(iwarn,nrad,nspden_updn,0,rhosph,xc_denpos)

!----------------------------------------------------------------------
!----- Compute Exc(rhoSPH) and Vxc(rhoSPH)
!----------------------------------------------------------------------

 ABI_ALLOCATE(exci,(nrad))
 ABI_ALLOCATE(vxci,(nrad,nspden_updn))
 call pawxcsph(exci,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden_updn,nspgrad,nvxcdgr,&
& order,pawrad,rhosph,vxci,xclevel)

!----------------------------------------------------------------------
!----- Compute numerical derivatives of Vxc (by finite difference scheme)
!----------------------------------------------------------------------

 if (option/=4) then
   ABI_ALLOCATE(exc_,(nrad))
   ABI_ALLOCATE(rho_,(nrad,nspden_updn))

   if (nspden_updn==2) rho_(:,2)=rhosph(:,2)

!  Compute Exc, Vxc for rho+delta_rho
   ABI_ALLOCATE(vxc1,(nrad,nspden_updn))
   rho_(:,1)=(one+delta)*rhosph(:,1)
   call pawxcsph(exc_,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden_updn,nspgrad,nvxcdgr,&
&   order,pawrad,rho_,vxc1,xclevel)

!  Compute Exc, Vxc for rho-delta_rho
   ABI_ALLOCATE(vxc2,(nrad,nspden_updn))
   rho_(:,1)=(one-delta)*rhosph(:,1)
   call pawxcsph(exc_,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden_updn,nspgrad,nvxcdgr,&
&   order,pawrad,rho_,vxc2,xclevel)

!  Additional terms for spin-polarized systems
   if (nspden_updn==2) then
     rho_(:,1)=rhosph(:,1)

!    Compute Exc, Vxc for rho+delta_rho_down
     ABI_ALLOCATE(vxcdn1,(nrad,nspden_updn))
     rho_(:,2)=(one+delta)*rhosph(:,2)
     call pawxcsph(exc_,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden_updn,nspgrad,nvxcdgr,&
&     order,pawrad,rho_,vxcdn1,xclevel)

!    Compute Exc, Vxc for rho-delta_rho_down
     ABI_ALLOCATE(vxcdn2,(nrad,nspden_updn))
     rho_(:,2)=(one-delta)*rhosph(:,2)
     call pawxcsph(exc_,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden_updn,nspgrad,nvxcdgr,&
&     order,pawrad,rho_,vxcdn2,xclevel)

   end if !nspden_updn==2
   ABI_DEALLOCATE(exc_)
   ABI_DEALLOCATE(rho_)

!  Store inverse of density finite step
   ABI_ALLOCATE(rhoinv,(nrad,nspden_updn))
   fact=one/delta;if (nspden_updn==1) fact=half*fact
   do ispden=1,nspden_updn
     do ir=1,nrad
       if (rhosph(ir,ispden)>tol14) then
         rhoinv(ir,ispden)=fact/rhosph(ir,ispden)
       else
         rhoinv(ir,ispden)=zero
       end if
     end do
   end do

!  Compute numerical first derivatives of Vxc (by finite difference scheme)
   if (option/=5) then
     ABI_ALLOCATE(d1vxc,(nrad,2*nspden_updn-1))
!    Non-magnetic system: compute dVxc/dn
     if (nspden==1) then
       d1vxc(1:nrad,1)=(vxc1(1:nrad,1)-vxc2(1:nrad,1))*half*rhoinv(1:nrad,1)
!      Collinear magnetism: compute dVxc_up/dn_up,dVxc_dn/dn_up,dVxc_dn/dn_dn
     else if (nspden==2) then
       d1vxc(1:nrad,1)=(vxc1(1:nrad,1)-vxc2(1:nrad,1))*half*rhoinv(1:nrad,1)
       d1vxc(1:nrad,2)=(vxc1(1:nrad,2)-vxc2(1:nrad,2))*half*rhoinv(1:nrad,1)
       d1vxc(1:nrad,3)=(vxcdn1(1:nrad,2)-vxcdn2(1:nrad,2))*half*rhoinv(1:nrad,2)
!      Non-collinear magnetism: compute 1/2 d(Vxc_up+Vxc_dn)/dn,1/2 d(Vxc_up-Vxc_dn)/dn
!      1/2 d(Vxc_up-Vxc_dn)/dm
     else if (nspden==4) then
       do ir=1,nrad
         fact=half*rhoinv(ir,1)
         dvxc1=(vxc1  (ir,1)-vxc2  (ir,1))*fact !dVxc_up/dn_up
         dvxc2=(vxc1  (ir,2)-vxc2  (ir,2))*fact !dVxc_dn/dn_up
         fact=half*rhoinv(ir,2)
         dvxc3=(vxcdn1(ir,2)-vxcdn2(ir,2))*fact !dVxc_dn/dn_dn
         dvxca=dvxc1+dvxc3;dvxcb=dvxc1-dvxc3;dvxcc=two*dvxc2 !Temporary terms
         d1vxc(ir,1)=quarter*(dvxca+dvxcc)  ! 1/2 d(Vxc_up+Vxc_dn)/dn
         d1vxc(ir,2)=quarter* dvxcb         ! 1/2 d(Vxc_up-Vxc_dn)/dn
         d1vxc(ir,3)=quarter*(dvxca-dvxcc)  ! 1/2 d(Vxc_up-Vxc_dn)/dm
       end do
     end if
   end if

!  Compute numerical second derivatives of Vxc (by finite difference scheme)
   if (option/=3.or.(option/=5.and.pawxcdev>=2)) then
     ABI_ALLOCATE(d2vxc,(nrad,3*nspden_updn-2))
!    Non-magnetic system: compute d2Vxc/dn2
     if (nspden==1) then
       d2vxc(1:nrad,1)=(vxc1(1:nrad,1)+vxc2(1:nrad,1)-two*vxci(1:nrad,1))*rhoinv(1:nrad,1)**2
!      Collinear magnetism: compute d2Vxc_up/dn_up2,d2Vxc_dn/dn_up2,d2Vxc_up/dn_dn2,d2Vxc_dn/dn_dn2
     else if (nspden==2) then
       d2vxc(1:nrad,1)=(vxc1(1:nrad,1)+vxc2(1:nrad,1)-two*vxci(1:nrad,1))*rhoinv(1:nrad,1)**2
       d2vxc(1:nrad,2)=(vxc1(1:nrad,2)+vxc2(1:nrad,2)-two*vxci(1:nrad,2))*rhoinv(1:nrad,1)**2
       d2vxc(1:nrad,3)=(vxcdn1(1:nrad,1)+vxcdn2(1:nrad,1)-two*vxci(1:nrad,1))*rhoinv(1:nrad,2)**2
       d2vxc(1:nrad,4)=(vxcdn1(1:nrad,2)+vxcdn2(1:nrad,2)-two*vxci(1:nrad,2))*rhoinv(1:nrad,2)**2
!      Non-collinear magnetism: compute 1/2 d2(Vxc_up+Vxc_dn)/dn2,1/2 d2(Vxc_up-Vxc_dn)/dn2
!      1/2 d2(Vxc_up+Vxc_dn)/dm2,1/2 d2(Vxc_up-Vxc_dn)/dm2
     else if (nspden==4) then
       do ir=1,nrad
         fact=rhoinv(ir,1)**2
         dvxc1=(vxc1  (ir,1)+vxc2  (ir,1)-two*vxci(ir,1))*fact !d2Vxc_up/dn_up2
         dvxc2=(vxc1  (ir,2)+vxc2  (ir,2)-two*vxci(ir,2))*fact !d2Vxc_dn/dn_up2
         fact=rhoinv(ir,2)**2
         dvxc3=(vxcdn1(ir,1)+vxcdn2(ir,1)-two*vxci(ir,1))*fact !d2Vxc_up/dn_dn2
         dvxc4=(vxcdn1(ir,2)+vxcdn2(ir,2)-two*vxci(ir,2))*fact !d2Vxc_dn/dn_dn2
         dvxca=dvxc1+dvxc4;dvxcb=dvxc1-dvxc4 !Temporary terms
         dvxcc=dvxc2+dvxc3;dvxcd=dvxc2-dvxc3 !Temporary terms
         d2vxc(ir,1)=eighth*(dvxca+three*dvxcc)  ! 1/2 d2(Vxc_up+Vxc_dn)/dn2
         d2vxc(ir,2)=eighth*(dvxcb+dvxcd)        ! 1/2 d2(Vxc_up-Vxc_dn)/dn2
         d2vxc(ir,3)=eighth*(dvxca-dvxcc)        ! 1/2 d2(Vxc_up+Vxc_dn)/dm2
         d2vxc(ir,4)=eighth*(dvxcb-three*dvxcd)  ! 1/2 d2(Vxc_up-Vxc_dn)/dm2
       end do
     end if
   end if

!  If non-collinear magnetism, store 1/2(Vxc_up+Vxc_dn) and 1/2(Vxc_up-Vxc_dn)
   if (nspden==4) then
     vxci(:,1)=half*(vxci(:,1)+vxci(:,2))
     vxci(:,2)=vxci(:,1)-vxci(:,2)
   end if

   ABI_DEALLOCATE(rhoinv)
   ABI_DEALLOCATE(vxc1)
   ABI_DEALLOCATE(vxc2)
   if (nspden_updn==2) then
     ABI_DEALLOCATE(vxcdn1)
     ABI_DEALLOCATE(vxcdn2)
   end if

 end if ! option/=4

 ABI_DEALLOCATE(rhosph)

!----------------------------------------------------------------------
!----- Compute useful sums of densities
!----------------------------------------------------------------------

 if (option/=4) then

!  Non-collinear magnetism: V0SUM=(m_0.m_L)/|m_0|
!  --------------------------------------------------
   if (nspden==4) then
     ABI_ALLOCATE(v0sum,(nrad,lm_size))
     v0sum(:,1)=zero
     do ilm=2,lm_size
       v0sum(1:nrad,ilm)=(rho_updn(1:nrad,1,2)*rho_updn(1:nrad,ilm,2) &
&       +rho_updn(1:nrad,1,3)*rho_updn(1:nrad,ilm,3) &
&       +rho_updn(1:nrad,1,4)*rho_updn(1:nrad,ilm,4))*m_norm_inv(1:nrad)
     end do
   end if

!  Non-magnetic system:
!  Compute
!  V1SUM1(r)=Sum_L{n_L(r)^2}
!  V2SUM1(r,L)=Sum_L1_L2{n_L1(r)*n_L2(r)*Gaunt_(L,L1,L2)}
!  Collinear magnetism:
!  Compute
!  V1SUM1(r)=Sum_L{n^up_L(r)^2}
!  V1SUM2(r)=Sum_L{n^up_L(r)*n^dn_L(r)}
!  V1SUM3(r)=Sum_L{n^dn_L(r)^2}
!  V2SUM1(r,L)=Sum_L1_L2{n^up_L1(r)*n^up_L2(r)*Gaunt_(L,L1,L2)}
!  V2SUM2(r,L)=Sum_L1_L2{n^up_L1(r)*n^dn_L2(r)*Gaunt_(L,L1,L2)}
!  V2SUM3(r,L)=Sum_L1_L2{n^dn_L1(r)*n^dn_L2(r)*Gaunt_(L,L1,L2)}
!  Non-collinear magnetism:
!  Compute
!  V1SUM1(r)=Sum_L{n_L(r)^2}
!  V1SUM2(r)=Sum_L{n_L(r) (m_0.m_L)}/|m_0|
!  V1SUM3(r)=Sum_L{(m_0.m_L)^2}/|m_0|^2
!  V2SUM1(r,L)=Sum_L1_L2{n_L1(r)*n_L2(r)*Gaunt_(L,L1,L2)}
!  V2SUM2(r,L)=Sum_L1_L2{n_L1(r) (m_0.m_L2)*Gaunt_(L,L1,L2)}/|m_0|
!  V2SUM3(r,L)=Sum_L1_L2{(m_0.m_L1)*(m_0.m_L2)*Gaunt_(L,L1,L2)}/|m_0|^2
   optv2=pawxcdev;if (option==5) optv2=max(pawxcdev,1)
   if (pawxcdev>=1)  then
     ABI_ALLOCATE(v1sum,(nrad,2*nspden_updn-1))
   end if
   if (optv2>=2)  then
     ABI_ALLOCATE(v2sum,(nrad,lm_size,2*nspden_updn-1))
   end if
   if (nspden/=4) then
     call pawxcsum(lmselect,lmselect,lm_size,nspden_updn,nrad,optv2,pawang,&
&     rho_updn(:,:,1),rho_updn(:,:,nspden_updn),v1sum,v2sum)
   else
     call pawxcsum(lmselect,lmselect,lm_size,nspden_updn,nrad,optv2,pawang,&
&     rho_updn(:,:,1),v0sum(:,:),v1sum,v2sum)
   end if

 end if !option

!----------------------------------------------------------------------
!----- Accumulate and store XC potential
!----------------------------------------------------------------------

 if (option/=3.and.option/=4) then

!  === First order development
!  ---------------------------
   if (pawxcdev>=1) then

!    Non-magnetic system:
     if (nspden_updn==1) then
       vxc(1:nrad,1,1)=v1sum(1:nrad,1)*d2vxc(1:nrad,1)*invsqfpi2+vxci(1:nrad,1)*sqfpi
       if (option/=5) then
         do ilm=2,lm_size
           if (lmselect(ilm)) then
             vxc(1:nrad,ilm,1)=d1vxc(1:nrad,1)*rho_updn(1:nrad,ilm,1)
           end if
         end do
       end if

!      Magnetic system:
     else if (nspden_updn==2) then
       vxc(1:nrad,1,1)=vxci(1:nrad,1)*sqfpi+invsqfpi2*(v1sum(1:nrad,1)*d2vxc(1:nrad,1) &
&       +two*v1sum(1:nrad,2)*d2vxc(1:nrad,2)+v1sum(1:nrad,3)*d2vxc(1:nrad,3))
       vxc(1:nrad,1,2)=vxci(1:nrad,2)*sqfpi+invsqfpi2*(v1sum(1:nrad,1)*d2vxc(1:nrad,2) &
&       +two*v1sum(1:nrad,2)*d2vxc(1:nrad,3)+v1sum(1:nrad,3)*d2vxc(1:nrad,4))
       if (option/=5) then
         if (nspden==2) then
           do ilm=2,lm_size
             if (lmselect(ilm)) then
               vxc(1:nrad,ilm,1)=vxc(1:nrad,ilm,1)+d1vxc(1:nrad,1)*rho_updn(1:nrad,ilm,1) &
&               +d1vxc(1:nrad,2)*rho_updn(1:nrad,ilm,2)
               vxc(1:nrad,ilm,2)=vxc(1:nrad,ilm,2)+d1vxc(1:nrad,2)*rho_updn(1:nrad,ilm,1) &
&               +d1vxc(1:nrad,3)*rho_updn(1:nrad,ilm,2)
             end if
           end do
         else if (nspden==4) then
           do ilm=2,lm_size
             if (lmselect(ilm)) then
               vxc(1:nrad,ilm,1)=vxc(1:nrad,ilm,1)+d1vxc(1:nrad,1)*rho_updn(1:nrad,ilm,1) &
&               +d1vxc(1:nrad,2)*v0sum(1:nrad,ilm)
               vxc(1:nrad,ilm,2)=vxc(1:nrad,ilm,2)+d1vxc(1:nrad,2)*rho_updn(1:nrad,ilm,1) &
&               +d1vxc(1:nrad,3)*v0sum(1:nrad,ilm)
             end if
           end do
         end if
       end if
     end if
   end if ! pawxcdev>=1

!  == 2nd order development
!  ---------------------------
   if (pawxcdev>=2.and.option/=5) then

!    Non-magnetic system:
     if (nspden_updn==1) then
       do ilm=2,lm_size
         vxc(1:nrad,ilm,1)=vxc(1:nrad,ilm,1)+half*d2vxc(1:nrad,1)*v2sum(1:nrad,ilm,1)
       end do

!      Magnetic system:
     else if (nspden_updn==2) then
       do ilm=2,lm_size
         vxc(1:nrad,ilm,1)=vxc(1:nrad,ilm,1)+d2vxc(1:nrad,2)*v2sum(1:nrad,ilm,2) &
&         +half*(d2vxc(1:nrad,1)*v2sum(1:nrad,ilm,1)+d2vxc(1:nrad,3)*v2sum(1:nrad,ilm,3))
         vxc(1:nrad,ilm,2)=vxc(1:nrad,ilm,2)+d2vxc(1:nrad,3)*v2sum(1:nrad,ilm,2) &
&         +half*(d2vxc(1:nrad,2)*v2sum(1:nrad,ilm,1)+d2vxc(1:nrad,4)*v2sum(1:nrad,ilm,3))
       end do
     end if
   end if !pawxcdev=2

!  === Pathological case: if rho(r) is negative, interpolate Vxc
!  -------------------------------------------------------------
   if (lmselect(1)) then
     do ispden=1,nspden_updn
       ir1=0;ir2=0
       do ir=1,nrad
         if (rho_updn(ir,1,ispden)<tol13) then
           if (ir1==0) ir1=ir-1
           ir2=ir+1
         else if (ir1>0) then
           if (ir1>1.or.ir2<nrad) then
             fact=(vxc(ir2,1,ispden)-vxc(ir1,1,ispden))/(pawrad%rad(ir2)-pawrad%rad(ir1))
             do jr=ir1+1,ir2-1
               vxc(jr,1,ispden)=vxc(ir1,1,ispden)+fact*(pawrad%rad(jr)-pawrad%rad(ir1))
             end do
           end if
           ir1=0;ir2=0
         end if
       end do
     end do
   end if

!  === Non-collinear magnetism: "rotate" back the XC potential
!  ------- ---------------------------------------------------
   if (nspden==4) then
     do ilm=1,lm_size
       do ir=1,nrad
         dvxca=vxc(ir,ilm,1)
         fact=vxc(ir,ilm,2)*m_norm_inv(ir)
         vxc(ir,ilm,1)=dvxca+fact*rho_updn(ir,1,4)
         vxc(ir,ilm,2)=dvxca-fact*rho_updn(ir,1,4)
         vxc(ir,ilm,3)=      fact*rho_updn(ir,1,2)
         vxc(ir,ilm,4)=     -fact*rho_updn(ir,1,3)
       end do
     end do
   end if

 end if !option/=3 and option/=4

 if (nspden==4)  then
   ABI_DEALLOCATE(m_norm_inv)
 end if

!----------------------------------------------------------------------
!----- Accumulate and store XC energies
!----------------------------------------------------------------------

!----- Calculate Exc (direct scheme) term
!----------------------------------------
 if (option/=1.and.option/=5) then
   ABI_ALLOCATE(ff,(nrad))

!  Contribution from spherical part of rho
   if (nspden==1.or.nspden==4) then
     ff(1:nrad)=rho_updn(1:nrad,1,1)*exci(1:nrad)*sqfpi
   else if (nspden==2) then
     ff(1:nrad)=(rho_updn(1:nrad,1,1)+rho_updn(1:nrad,1,2))*exci(1:nrad)*sqfpi
   end if

!  Contribution from aspherical part of rho
   if (option/=4) then

!    First order development
     if (pawxcdev>=1) then
       if (nspden_updn==1) then
         ff(1:nrad)=ff(1:nrad)+half*v1sum(1:nrad,1)*d1vxc(1:nrad,1)
       else if (nspden_updn==2) then
         ff(1:nrad)=ff(1:nrad)+v1sum(1:nrad,2)*d1vxc(1:nrad,2) &
&         +half*(v1sum(1:nrad,1)*d1vxc(1:nrad,1)+v1sum(1:nrad,3)*d1vxc(1:nrad,3))
       end if
     end if

!    Second order development
     if (pawxcdev>=2) then
       ABI_ALLOCATE(gg,(nrad))

       gg=zero
       do ilm=2,lm_size
         if (lmselect(ilm)) then
           gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,1)*rho_updn(1:nrad,ilm,1)
         end if
       end do
       ff(1:nrad)=ff(1:nrad)+sixth*gg(1:nrad)*d2vxc(1:nrad,1)

       if (nspden_updn==2) then
         gg=zero
         if (nspden==2) then
           do ilm=2,lm_size
             if (lmselect(ilm)) then
               gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,3)*rho_updn(1:nrad,ilm,2)
             end if
           end do
         else if (nspden==4) then
           do ilm=2,lm_size
             if (lmselect(ilm)) then
               gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,3)*v0sum(1:nrad,ilm)
             end if
           end do
         end if
         ff(1:nrad)=ff(1:nrad)+sixth*gg(1:nrad)*d2vxc(1:nrad,4)
         gg=zero
         do ilm=2,lm_size
           if (lmselect(ilm)) then
             gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,2)*rho_updn(1:nrad,ilm,1)
           end if
         end do
         ff(1:nrad)=ff(1:nrad)+half*gg(1:nrad)*d2vxc(1:nrad,2)
         gg=zero
         do ilm=2,lm_size
           if (lmselect(ilm)) then
             gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,3)*rho_updn(1:nrad,ilm,1)
           end if
         end do
         ff(1:nrad)=ff(1:nrad)+half*gg(1:nrad)*d2vxc(1:nrad,3)

         ABI_DEALLOCATE(gg)
       end if
     end if

   end if ! option/=4

   ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
   call simp_gen(enxc,ff,pawrad)
   ABI_DEALLOCATE(ff)
 end if ! option/=1 and option/=5

 ABI_DEALLOCATE(exci)
 ABI_DEALLOCATE(vxci)
 if (nspden==4  .and.option/=4)  then
   ABI_DEALLOCATE(v0sum)
 end if
 if (pawxcdev>=1.and.option/=4)  then
   ABI_DEALLOCATE(v1sum)
 end if
 if (pawxcdev>=2.and.option/=4.and.option/=5)  then
   ABI_DEALLOCATE(v2sum)
 end if
 if (option/=4.and.(option/=3.or.(pawxcdev>=2.and.option/=5)))  then
   ABI_DEALLOCATE(d2vxc)
 end if
 if (option/=4.and.option/=5)  then
   ABI_DEALLOCATE(d1vxc)
 end if

!----- Calculate Excdc double counting term
!------------------------------------------
 if (option==0.or.option==2) then

!  Build appropriate density
   if (usexcnhat==1) then
     if (nspden==1.or.nspden==4) then
       rho_updn(:,:,:)=rho_updn(:,:,:)+nhat(:,:,:)
     else if (nspden==2) then
       rho_updn(:,:,1)=rho_updn(:,:,1)+nhat(:,:,2)
       rho_updn(:,:,2)=rho_updn(:,:,2)+nhat(:,:,1)-nhat(:,:,2)
     end if
   end if
   if (usecore==1) then
     if (nspden==1.or.nspden==4) then
       rho_updn(:,1,1)=rho_updn(:,1,1)-sqfpi*corexc(:)
     else if (nspden==2) then
       rho_updn(:,1,1)=rho_updn(:,1,1)-sqfpi2*corexc(:)
       rho_updn(:,1,2)=rho_updn(:,1,2)-sqfpi2*corexc(:)
     end if
   end if

   ABI_ALLOCATE(ff,(nrad))
   ff(1:nrad)=zero

!  Non magnetic or collinear magnetic system:
   if (nspden/=4) then
     do ispden=1,nspden_updn
       do ilm=1,lm_size
         if (lmselect(ilm)) ff(1:nrad)=ff(1:nrad)+vxc(1:nrad,ilm,ispden)*rho_updn(1:nrad,ilm,ispden)
       end do
     end do
   else
!    Non-collinear magnetic system:
     do ilm=1,lm_size
       if (lmselect(ilm)) then
         do ir=1,nrad
           dvxca=vxc(ir,ilm,1)+vxc(ir,ilm,2);dvxcb=vxc(ir,ilm,1)-vxc(ir,ilm,2)
           ff(ir)=ff(ir)+half*(dvxca*rho_updn(ir,ilm,1)+dvxcb*rho_updn(ir,ilm,4)) &
&           +vxc(ir,ilm,3)*rho_updn(ir,ilm,2)-vxc(ir,ilm,4)*rho_updn(ir,ilm,3)
         end do
       end if
     end do
   end if

   ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
   call simp_gen(enxcdc,ff,pawrad)
   ABI_DEALLOCATE(ff)

 end if ! option

 ABI_DEALLOCATE(rho_updn)

!----- End of routine
 call timab(81,2,tsec)

 DBG_EXIT("COLL")

 end subroutine pawxcm
!!***
