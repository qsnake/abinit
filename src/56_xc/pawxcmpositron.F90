!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxcmpositron
!! NAME
!! pawxcmpositron
!!
!! FUNCTION
!! PAW only
!! Compute electron-positron correlation potential and energies inside a PAW sphere
!! LDA+GGA - USE A DEVELOPMENT OF THE DENSITY OVER (L,M) MOMENTS
!! Driver of XC functionals.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT,GJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  calctype=type of electron-positron calculation:
!!           calctype=1 : positron in electronic density
!!           calctype=2 : electrons in positronic density
!!  corexc(pawrad%mesh_size)=electron core density on radial grid
!!  ixcpositron=choice of electron-positron XC scheme
!!  lm_size=size of density array rhor (see below)
!!  lmselect   (lm_size)=select the non-zero LM-moments of input density rhor    (see below)
!!  lmselect_ep(lm_size)=select the non-zero LM-moments of input density rhor_ep (see below)
!!  nhat   (pawrad%mesh_size,lm_size,nspden)=compensation density corresponding to rhor
!!  nhat_ep(pawrad%mesh_size,lm_size,nspden)=compensation density corresponding to rhor_ep
!!  nspden=number of spin-density components
!!  option=0 compute both XC energies (direct+double-counting) and potential
!!         1 compute only XC potential
!!         2 compute only XC energies (direct+double-counting)
!!         3 compute only XC energy by direct scheme
!!         4 compute only XC energy by direct scheme for spherical part of the density
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawxcdev=order of Vxc development
!!  posdensity0_limit=True if we are in the zero positron density limit
!!  rhor(pawrad%mesh_size,lm_size,nspden)=electron (or positron) density in real space
!!                             (total in 1st half and spin-up in 2nd half if nspden=2)
!!                             Contents depends on calctype value:
!!                             calctype=1: rhor is the positronic density
!!                             calctype=2: rhor is the electronic density
!!  rhor_ep(pawrad%mesh_size,lm_size,nspden)=electron (or positron) density in real space
!!                             (total in 1st half and spin-up in 2nd half if nspden=2)
!!                             Contents depends on calctype value:
!!                             calctype=1: rhor_ep is the electronic density
!!                             calctype=2: rhor_ep is the positronic density
!!  usecore= 1 if core density has to be used in Exc/Vxc for the electronic density ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in double counting energy term only
!!             2 if compensation density (nhat) has to be used in Exc/Vxc and double counting energy term
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
!!
!! NOTES
!!
!! PARENTS
!!      pawdenpot
!!
!! CHILDREN
!!      mkdenpos,pawxcsphpositron,pawxcsum,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawxcmpositron(calctype,corexc,enxc,enxcdc,ixcpositron,lm_size,lmselect,lmselect_ep,&
 &                         nhat,nhat_ep,nspden,option,pawang,pawrad,pawxcdev,posdensity0_limit,&
 &                         rhor,rhor_ep,usecore,usexcnhat,vxc,xc_denpos)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_radmesh,          only : simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawxcmpositron'
 use interfaces_56_xc, except_this_one => pawxcmpositron
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: calctype,ixcpositron,lm_size,nspden,option,pawxcdev,usecore
 integer,intent(in) :: usexcnhat
 logical,intent(in) :: posdensity0_limit
 real(dp),intent(in) :: xc_denpos
 real(dp),intent(out) :: enxc,enxcdc
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size),lmselect_ep(lm_size)
 real(dp),intent(in) :: corexc(pawrad%mesh_size)
 real(dp),intent(in) :: nhat   (pawrad%mesh_size,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: nhat_ep(pawrad%mesh_size,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor   (pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(in) :: rhor_ep(pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(out) :: vxc(pawrad%mesh_size,lm_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ilm,ir,ir1,ir2,iwarn,iwarnp,jr,nrad
 real(dp),parameter :: delta=1.d-4
 real(dp) :: fact,invsqfpi,sqfpi
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: d1vxc(:,:),d2vxc(:,:),fxc_(:),ff(:),fxci(:),gg(:)
 real(dp),allocatable :: rho_(:),rhotot(:,:),rhotot_ep(:,:),rhoinv(:),rhoinv_ep(:)
 real(dp),allocatable :: rhosph(:),rhosph_ep(:),v1sum(:,:),v2sum(:,:,:)
 real(dp),allocatable :: vxce1(:),vxce1_ep(:),vxce2(:),vxce2_ep(:)
 real(dp),allocatable :: vxcp1(:),vxcp1_ep(:),vxcp2(:),vxcp2_ep(:)
 real(dp),allocatable :: vxcei(:),vxcpi(:)

!************************************************************************

 DBG_ENTER("COLL")

!----- Check options
 if(calctype/=1.and.calctype/=2) then
   msg='  Invalid value for calctype'
   MSG_BUG(msg)
 end if

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

!Initializations and constants
 iwarn=0;iwarnp=1
 nrad=pawrad%mesh_size
 sqfpi=sqrt(four_pi)
 invsqfpi=one/sqfpi

!Initializations of output arrays
 if (option/=1) enxc=zero
 if (option==0.or.option==2) enxcdc=zero
 if (option<3) vxc(:,:,:)=zero

 if (ixcpositron==0) then ! No xc at all is applied (usually for testing)
   msg='  Note that no xc is applied (ixc=0).'
   MSG_WARNING(msg)
   return
 end if

!----------------------------------------------------------------------
!----- Build several densities
!----------------------------------------------------------------------

!rhotot/rhotot_ep contain the effective total densities used for XC
!with core density and/or compensation density eventually included
!-----------------------------------------------------------------
!Input density
 ABI_ALLOCATE(rhotot,(nrad,lm_size))
 ABI_ALLOCATE(rhotot_ep,(nrad,lm_size))
 rhotot   (:,:)=rhor   (:,:,1)
 rhotot_ep(:,:)=rhor_ep(:,:,1)
!Eventually add compensation density
 if (usexcnhat==2) then
   rhotot   (:,:)=rhotot   (:,:)+nhat   (:,:,1)
   rhotot_ep(:,:)=rhotot_ep(:,:)+nhat_ep(:,:,1)
 end if
!Eventually add core density
 if (usecore==1) then
   if (calctype==1) rhotot_ep(:,1)=rhotot_ep(:,1)+sqfpi*corexc(:)
   if (calctype==2) rhotot   (:,1)=rhotot   (:,1)+sqfpi*corexc(:)
 end if

!rhoSPH/rhoSPH_ep contain the spherical part of effective densities
!(including Y00 spherical harmonic)
!-----------------------------------------------------------------
 ABI_ALLOCATE(rhosph,(nrad))
 ABI_ALLOCATE(rhosph_ep,(nrad))

 rhosph   (:)=rhotot   (:,1)*invsqfpi
 rhosph_ep(:)=rhotot_ep(:,1)*invsqfpi

!Make spherical densities positive
 if (calctype==1) then
   if (.not.posdensity0_limit) then
     call mkdenpos(iwarnp,nrad,1,1,rhosph,xc_denpos)
   end if
   call mkdenpos(iwarn ,nrad,1,1,rhosph_ep,xc_denpos)
 else if (calctype==2) then
   call mkdenpos(iwarn ,nrad,1,1,rhosph,xc_denpos)
   if (.not.posdensity0_limit) then
     call mkdenpos(iwarnp,nrad,1,1,rhosph_ep,xc_denpos)
   end if
 end if

!----------------------------------------------------------------------
!----- Compute Exc(rhoSPH,rhoSPH_ep) and Vxc(rhoSPH,rhoSPH_ep)
!----------------------------------------------------------------------

 ABI_ALLOCATE(fxci,(nrad))
 ABI_ALLOCATE(vxcei,(nrad))
 ABI_ALLOCATE(vxcpi,(nrad))
 call pawxcsphpositron(calctype,fxci,ixcpositron,nrad,pawrad,posdensity0_limit,rhosph,rhosph_ep,vxcei,vxcpi)

!----------------------------------------------------------------------
!----- Compute numerical derivatives of Vxc (by finite diff. scheme)
!----------------------------------------------------------------------

 if (option/=4) then

   ABI_ALLOCATE(fxc_,(nrad))
   ABI_ALLOCATE(rho_,(nrad))

!  Compute Vxc for (rho+delta_rho,rho_ep)
   ABI_ALLOCATE(vxce1,(nrad))
   ABI_ALLOCATE(vxcp1,(nrad))
   rho_(:)=(one+delta)*rhosph(:)
   call pawxcsphpositron(calctype,fxc_,ixcpositron,nrad,pawrad,posdensity0_limit,rho_,rhosph_ep,vxce1,vxcp1)

!  Compute Vxc for(rho-delta_rho,rho_ep)
   ABI_ALLOCATE(vxce2,(nrad))
   ABI_ALLOCATE(vxcp2,(nrad))
   rho_(:)=(one-delta)*rhosph(:)
   call pawxcsphpositron(calctype,fxc_,ixcpositron,nrad,pawrad,posdensity0_limit,rho_,rhosph_ep,vxce2,vxcp2)

!  Compute Vxc for (rho,rho_ep+delta_rho_ep)
   ABI_ALLOCATE(vxce1_ep,(nrad))
   ABI_ALLOCATE(vxcp1_ep,(nrad))
   rho_(:)=(one+delta)*rhosph_ep(:)
   call pawxcsphpositron(calctype,fxc_,ixcpositron,nrad,pawrad,posdensity0_limit,rhosph,rho_,vxce1_ep,vxcp1_ep)

!  Compute Vxc for (rho,rho_ep-delta_rho_ep)
   ABI_ALLOCATE(vxce2_ep,(nrad))
   ABI_ALLOCATE(vxcp2_ep,(nrad))
   rho_(:)=(one-delta)*rhosph_ep(:)
   call pawxcsphpositron(calctype,fxc_,ixcpositron,nrad,pawrad,posdensity0_limit,rhosph,rho_,vxce2_ep,vxcp2_ep)

   ABI_DEALLOCATE(fxc_)
   ABI_DEALLOCATE(rho_)

!  Store inverse of density finite step
   ABI_ALLOCATE(rhoinv,(nrad))
   ABI_ALLOCATE(rhoinv_ep,(nrad))
   fact=one/delta
   do ir=1,nrad
     if (rhosph(ir)>tol14) then
       rhoinv(ir)=fact/rhosph(ir)
     else
       rhoinv(ir)=zero
     end if
     if (rhosph_ep(ir)>tol14) then
       rhoinv_ep(ir)=fact/rhosph_ep(ir)
     else
       rhoinv_ep(ir)=zero
     end if
   end do

!  Compute numerical first derivatives of Vxc (by finite difference scheme)
   ABI_ALLOCATE(d1vxc,(nrad,3))
   if (calctype==1) then
     d1vxc(:,1)=(vxcp1   (:)-vxcp2   (:))*half*rhoinv   (:)  ! dVxc+/drho+
     d1vxc(:,2)=(vxcp1_ep(:)-vxcp2_ep(:))*half*rhoinv_ep(:)  ! dVxc+/drho-
     d1vxc(:,3)=(vxce1_ep(:)-vxce2_ep(:))*half*rhoinv_ep(:)  ! dVxc-/drho-
   else if (calctype==2) then
     d1vxc(:,1)=(vxce1   (:)-vxce2   (:))*half*rhoinv   (:)  ! dVxc-/drho-
     d1vxc(:,2)=(vxcp1   (:)-vxcp2   (:))*half*rhoinv   (:)  ! dVxc+/drho-
!    d1vxc(:,2)=(vxce1_ep(:)-vxce2_ep(:))*half*rhoinv_ep(:)  ! dVxc-/drho+
     d1vxc(:,3)=(vxcp1_ep(:)-vxcp2_ep(:))*half*rhoinv_ep(:)  ! dVxc+/drho+
   end if

!  Compute numerical second derivatives of Vxc (by finite difference scheme)
   if (option<3.or.pawxcdev>1) then
     ABI_ALLOCATE(d2vxc,(nrad,4))
     if (calctype==1) then
       d2vxc(:,1)=(vxcp1   (:)+vxcp2   (:)-two*vxcpi(:))*rhoinv   (:)**2  ! d2Vxc+/drho+_drho+
       d2vxc(:,2)=(vxce1   (:)+vxce2   (:)-two*vxcei(:))*rhoinv   (:)**2  ! d2Vxc-/drho+_drho+
       d2vxc(:,3)=(vxcp1_ep(:)+vxcp2_ep(:)-two*vxcpi(:))*rhoinv_ep(:)**2  ! d2Vxc+/drho-_drho-
       d2vxc(:,4)=(vxce1_ep(:)+vxce2_ep(:)-two*vxcei(:))*rhoinv_ep(:)**2  ! d2Vxc-/drho-_drho-
     else if (calctype==2) then
       d2vxc(:,1)=(vxce1   (:)+vxce2   (:)-two*vxcei(:))*rhoinv   (:)**2  ! d2Vxc-/drho-_drho-
       d2vxc(:,2)=(vxcp1   (:)+vxcp2   (:)-two*vxcpi(:))*rhoinv   (:)**2  ! d2Vxc+/drho-_drho-
       d2vxc(:,3)=(vxce1_ep(:)+vxce2_ep(:)-two*vxcei(:))*rhoinv_ep(:)**2  ! d2Vxc-/drho+_drho+
       d2vxc(:,4)=(vxcp1_ep(:)+vxcp2_ep(:)-two*vxcpi(:))*rhoinv_ep(:)**2  ! d2Vxc+/drho+_drho+
     end if
   end if ! option

   ABI_DEALLOCATE(rhoinv)
   ABI_DEALLOCATE(rhoinv_ep)
   ABI_DEALLOCATE(vxce1)
   ABI_DEALLOCATE(vxcp1)
   ABI_DEALLOCATE(vxce2)
   ABI_DEALLOCATE(vxcp2)
   ABI_DEALLOCATE(vxce1_ep)
   ABI_DEALLOCATE(vxcp1_ep)
   ABI_DEALLOCATE(vxce2_ep)
   ABI_DEALLOCATE(vxcp2_ep)

 end if ! option/=4

 ABI_DEALLOCATE(rhosph)
 ABI_DEALLOCATE(rhosph_ep)

!----------------------------------------------------------------------
!----- Compute useful sums of densities
!----------------------------------------------------------------------

 if (option<3.or.option/=1) then

!  Compute V1SUM1(r)=Sum_L{n^el_L(r)^2}
!  V1SUM2(r)=Sum_L{n^el_L(r)*n^pos_L(r)}
!  V1SUM3(r)=Sum_L{n^pos_L(r)^2}
!  V2SUM1(r,L)=Sum_L1_L2{n^el_L1(r)*n^el_L2(r)*Gaunt_(L,L1,L2)}
!  V2SUM2(r,L)=Sum_L1_L2{n^el_L1(r)*n^pos_L2(r)*Gaunt_(L,L1,L2)}
!  V2SUM3(r,L)=Sum_L1_L2{n^pos_L1(r)*n^pos_L2(r)*Gaunt_(L,L1,L2)}
   if (pawxcdev>=1)  then
     ABI_ALLOCATE(v1sum,(nrad,3))
   end if
   if (pawxcdev>=2)  then
     ABI_ALLOCATE(v2sum,(nrad,lm_size,3))
   end if
   call pawxcsum(lmselect,lmselect_ep,lm_size,2,nrad,pawxcdev,pawang,rhotot,rhotot_ep,v1sum,v2sum)

 end if !option

!----------------------------------------------------------------------
!----- Accumulate and store XC potential
!----------------------------------------------------------------------

 if (option<3) then

!  if (option==0.or.option==2) allocate(vxc_ep(nrad,lm_size))

!  === First order development
!  ---------------------------
   if (pawxcdev>=1) then
     if (calctype==1) vxc(:,1,1)=vxcpi(:)*sqfpi
     if (calctype==2) vxc(:,1,1)=vxcei(:)*sqfpi
     vxc(:,1,1)=vxc(:,1,1)+invsqfpi*(d2vxc(:,2)*v1sum(:,2) &
&     +half*(d2vxc(:,1)*v1sum(:,1)+d2vxc(:,3)*v1sum(:,3)))
     do ilm=2,lm_size
       if (lmselect(ilm))    vxc(:,ilm,1)=vxc(:,ilm,1)+d1vxc(:,1)*rhotot   (:,ilm)
       if (lmselect_ep(ilm)) vxc(:,ilm,1)=vxc(:,ilm,1)+d1vxc(:,2)*rhotot_ep(:,ilm)
     end do
!    if (option==0.or.option==2) then
!    if (calctype==1) vxc_ep(:,1)=vxcei(:)*sqfpi
!    if (calctype==2) vxc_ep(:,1)=vxcpi(:)*sqfpi
!    vxc_ep(:,1)=vxc_ep(:,1,1)+invsqfpi*(d2vxc(:,3)*v1sum(:,2) &
!    &             +half*(d2vxc(:,2)*v1sum(:,1)+d2vxc(:,4)*v1sum(:,3)))
!    do ilm=2,lm_size
!    if (lmselect(ilm))    vxc_ep(:,ilm)=vxc_ep(:,ilm)+d1vxc(:,2)*rhotot   (:,ilm)
!    if (lmselect_ep(ilm)) vxc_ep(:,ilm)=vxc_ep(:,ilm)+d1vxc(:,3)*rhotot_ep(:,ilm)
!    end do
!    end if
   end if ! pawxcdev>=1

!  == 2nd order development
!  ---------------------------
   if (pawxcdev>=2) then
     do ilm=2,lm_size
       vxc(:,ilm,1)=vxc(:,ilm,1)+d2vxc(:,2)*v2sum(:,ilm,2) &
&       +half*(d2vxc(:,1)*v2sum(:,ilm,1)+d2vxc(:,3)*v2sum(:,ilm,3))
     end do
!    if (option==0.or.option==2) then
!    do ilm=2,lm_size
!    vxc_ep(:,ilm)=vxc_ep(:,ilm)+d2vxc(:,3)*v2sum(:,ilm,2) &
!    &                +half*(d2vxc(:,2)*v2sum(:,ilm,1)+d2vxc(:,4)*v2sum(:,ilm,3))
!    end do
!    end if
   end if !pawxcdev=2

!  === Pathological case: if rho(r) is negative, interpolate Vxc
!  -------------------------------------------------------------
   if (lmselect(1)) then
     ir1=0;ir2=0
     do ir=1,nrad
       if (rhotot(ir,1)<tol13) then
         if (ir1==0) ir1=ir-1
         ir2=ir+1
       else if (ir1>0) then
         if (ir1>1.or.ir2<nrad) then
           fact=(vxc(ir2,1,1)-vxc(ir1,1,1))/(pawrad%rad(ir2)-pawrad%rad(ir1))
           do jr=ir1+1,ir2-1
             vxc(jr,1,1)=vxc(ir1,1,1)+fact*(pawrad%rad(jr)-pawrad%rad(ir1))
           end do
         end if
         ir1=0;ir2=0
       end if
     end do
   end if
!  if (option==0.or.option==2) then
!  if (lmselect_ep(1)) then
!  ir1=0;ir2=0
!  do ir=1,nrad
!  if (rhotot_ep(ir,1)<tol14) then
!  if (ir1==0) ir1=ir-1
!  ir2=ir+1
!  else if (ir1>0) then
!  if (ir1>1.or.ir2<nrad) then
!  fact=(vxc_ep(ir2,1)-vxc_ep(ir1,1))/(pawrad%rad(ir2)-pawrad%rad(ir1))
!  do jr=ir1+1,ir2-1
!  vxc_ep(jr,1)=vxc_ep(ir1,1)+fact*(pawrad%rad(jr)-pawrad%rad(ir1))
!  end do
!  end if
!  ir1=0;ir2=0
!  end if
!  end do
!  end if
!  end if

!  When vxc is dimensionned as polarized...
   if (nspden>=2) vxc(:,:,2)=vxc(:,:,1)
   if (nspden==4) vxc(:,:,3:4)=zero

 end if !option<3

 ABI_DEALLOCATE(vxcei)
 ABI_DEALLOCATE(vxcpi)

!----------------------------------------------------------------------
!----- Accumulate and store XC energies
!----------------------------------------------------------------------

!----- Calculate Exc (direct scheme) term
!----------------------------------------

 if (option/=1) then
   ABI_ALLOCATE(ff,(nrad))

!  Contribution from spherical part of rho
   ff(:)=fxci(:)*four_pi

!  Contribution from aspherical part of rho
   if (option/=4) then

!    First order development
     if (pawxcdev>=1) then
       ff(:)=ff(:)+v1sum(:,2)*d1vxc(:,2) &
&       +half*(v1sum(:,1)*d1vxc(:,1)+v1sum(:,3)*d1vxc(:,3))
     end if

!    Second order development
     if (pawxcdev>=2) then
       ABI_ALLOCATE(gg,(nrad))
       gg=zero
       do ilm=2,lm_size
         if (lmselect(ilm))    gg(:)=gg(:)+v2sum(:,ilm,1)*rhotot(:,ilm)
       end do
       ff(:)=ff(:)+sixth*gg(:)*d2vxc(:,1)
       gg=zero
       do ilm=2,lm_size
         if (lmselect(ilm))    gg(:)=gg(:)+v2sum(:,ilm,2)*rhotot(:,ilm)
       end do
       ff(:)=ff(:) +half*gg(:)*d2vxc(:,2)
       gg=zero
       do ilm=2,lm_size
         if (lmselect(ilm))    gg(:)=gg(:)+v2sum(:,ilm,3)*rhotot(:,ilm)
       end do
       ff(:)=ff(:) +half*gg(:)*d2vxc(:,3)
       gg=zero
       do ilm=2,lm_size
         if (lmselect_ep(ilm)) gg(:)=gg(:)+v2sum(:,ilm,3)*rhotot_ep(:,ilm)
       end do
       ff(:)=ff(:)+sixth*gg(:)*d2vxc(:,4)
       ABI_DEALLOCATE(gg)
     end if ! pawxcdev>=2

   end if ! option/=4

   ff(:)=ff(:)*pawrad%rad(:)**2
   call simp_gen(enxc,ff,pawrad)
   ABI_DEALLOCATE(ff)
 end if ! option/=1

 ABI_DEALLOCATE(fxci)
 if (pawxcdev>=1.and.(option<3.or.option/=1))  then
   ABI_DEALLOCATE(v1sum)
 end if
 if (pawxcdev>=2.and.(option<3.or.option/=1))  then
   ABI_DEALLOCATE(v2sum)
 end if
 if (option<3.or.(option/=4.and.pawxcdev>1))   then
   ABI_DEALLOCATE(d2vxc)
 end if
 if (option/=4)  then
   ABI_DEALLOCATE(d1vxc)
 end if

!----- Calculate Excdc double counting term
!------------------------------------------
 if (option==0.or.option==2) then

!  Build appropriate density
   if (usexcnhat==1) rhotot(:,:)=rhotot(:,:)+nhat(:,:,1)
   if (usecore==1.and.calctype==2) rhotot(:,1)=rhotot(:,1)-sqfpi*corexc(:)

!  Integrate with potential
   ABI_ALLOCATE(ff,(nrad))
   ff(:)=zero
   do ilm=1,lm_size
     if (lmselect(ilm)) ff(:)=ff(:)+vxc(:,ilm,1)*rhotot(:,ilm)
   end do
   ff(:)=ff(:)*pawrad%rad(:)**2
   call simp_gen(enxcdc,ff,pawrad)
   ABI_DEALLOCATE(ff)
 end if ! option

 ABI_DEALLOCATE(rhotot)
 ABI_DEALLOCATE(rhotot_ep)
!if (option==0.or.option==2) deallocate(vxc_ep)

!----- End of routine
 DBG_EXIT("COLL")
 end subroutine pawxcmpositron
!!***

