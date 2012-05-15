!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxc
!! NAME
!! pawxc
!!
!! FUNCTION
!! PAW only
!! Start from the density or spin-density, and compute xc correlation
!! potential and energies inside a paw sphere.
!! USE THE DENSITY OVER A WHOLE SPHERICAL GRID (r,theta,phi)
!! Driver of XC functionals.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!! This routine has been written from rhohxc
!!
!! INPUTS
!!  corexc(pawrad%mesh_size)=core density on radial grid
!!  ixc= choice of exchange-correlation scheme (see above and below)
!!  lm_size=size of density array rhor (see below)
!!  lmselect(lm_size)=select the non-zero LM-moments of input density rhor
!!  nhat(pawrad%mesh_size,lm_size,nspden)=compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nkxc=second dimension of the kxc array. If /=0, the exchange-correlation kernel must be computed
!!  nspden=number of spin-density components
!!  option=0  compute both XC energies (direct+double-counting) and potential
!!         1  compute only XC potential
!!         2  compute only XC energies (direct+double-counting)
!!         3  compute only XC energy by direct scheme
!!         4  compute only XC energy by direct scheme for spherical part of the density
!!         5  compute only XC potential for spherical part of the density
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
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
!!  == if option=0, 2, 3, or 4 ==
!!    enxc=returned exchange and correlation energy (hartree)
!!  == if option=0 or 2 ==
!!    enxcdc=returned exchange-cor. contribution to double-counting energy
!!  == if option=0, 1 or 5 ==
!!    vxc(pawrad%mesh_size,pawang%angl_size,nspden)=xc potential
!!       (spin up in 1st half and spin-down in 2nd half if nspden=2)
!!  == if nkxc>0 ==
!!    kxc(pawrad%mesh_size,pawang%angl_size,nkxc)=xc kernel
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
!!
!! PARENTS
!!      pawdenpot,psp7in
!!
!! CHILDREN
!!      deducer0,drivexc,mkdenpos,nderiv_gen,simp_gen,size_dvxc,timab,xcmult
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawxc(corexc,enxc,enxcdc,ixc,kxc,lm_size,lmselect,nhat,nkxc,nspden,option,&
&                pawang,pawrad,rhor,usecore,usexcnhat,vxc,xclevel,xc_denpos)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_radmesh,   only : simp_gen, nderiv_gen, deducer0
#if defined HAVE_DFT_LIBXC
 use libxc_functionals
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawxc'
 use interfaces_18_timing
 use interfaces_56_xc, except_this_one => pawxc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,lm_size,nkxc,nspden,option,usecore,usexcnhat,xclevel
 real(dp),intent(in) :: xc_denpos
 real(dp),intent(out) :: enxc,enxcdc
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size)
 real(dp),intent(in) :: corexc(pawrad%mesh_size)
 real(dp),intent(in) :: nhat(pawrad%mesh_size,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in),target :: rhor(pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(out) :: kxc(pawrad%mesh_size,pawang%angl_size,nkxc)
 real(dp),intent(out) :: vxc(pawrad%mesh_size,pawang%angl_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ilm,ipts,ir,ispden,iwarn,ndvxc,nd2vxc,ngr2,ngrad,npts,nrad
 integer :: nspden_eff,nspden_updn,nspgrad,nvxcdgr,order,ylm_size
 real(dp) :: dvdn,dvdz,enxcr,factor,m_norm_min,vxcrho
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: d2vxcar(:),dgxc(:),dnexcdn(:,:),drho(:),drhocore(:),dvxcdgr(:,:),dvxci(:,:)
 real(dp),allocatable :: exci(:),ff(:),grho2_updn(:,:),gxc(:,:,:,:),m_norm(:),rho_updn(:,:)
 real(dp),allocatable :: rhoarr(:,:),rhonow(:,:,:)
 real(dp),allocatable :: vxci(:,:)
 real(dp),allocatable,target :: rhohat(:,:,:)
 real(dp),pointer :: rho_(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(81,1,tsec)

!----------------------------------------------------------------------
!----- Check options
!----------------------------------------------------------------------
 if(nspden==4.and.nkxc>0) then
   msg='  Kxc for nspden=4 not implemented !'
   MSG_ERROR(msg)
 end if
 if(nspden==4.and.xclevel==2) then
   msg='  GGA for nspden=4 not implemented !'
   MSG_ERROR(msg)
 end if
 if(pawang%angl_size==0) then
   msg='  pawang%angl_size=0 !'
   MSG_BUG(msg)
 end if
 if(.not.associated(pawang%ylmr)) then
   msg='  pawang%ylmr must be allocated !'
   MSG_BUG(msg)
 end if
 if(xclevel==2.and.(.not.associated(pawang%ylmrgr))) then
   msg='  pawang%ylmrgr must be allocated !'
   MSG_BUG(msg)
 end if
 if(option==4.or.option==5) then
   if (pawang%angl_size/=1) then
     msg='  When option=4 or 5, pawang%angl_size must be 1 !'
     MSG_BUG(msg)
   end if
   if (pawang%ylm_size/=1) then
     msg='  When option=4 or 5, pawang%ylm_size must be 1 !'
     MSG_BUG(msg)
   end if
   if (abs(pawang%anginit(1,1)-one)>tol12.or.abs(pawang%anginit(2,1))>tol12.or. &
&   abs(pawang%anginit(3,1))>tol12) then
     msg='  When option=4 or 5, pawang%anginit must be (1 0 0) !'
     MSG_BUG(msg)
   end if
 end if

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------
 iwarn=0
 nspden_updn=min(nspden,2)
 nspden_eff=nspden_updn;if (nspden==4.and.ngrad==2) nspden_eff=4
 nrad=pawrad%mesh_size
 npts=pawang%angl_size
 ylm_size=min(lm_size,pawang%ylm_size)
 ngrad=1;if(xclevel==2)ngrad=2
 nspgrad=nspden_updn*ngrad;if(nspden_updn==2.and.ngrad==2)nspgrad=5
 if (option/=1.and.option/=5) enxc=zero
 if (option==0.or.option==2) enxcdc=zero
 if (option/=3.and.option/=4) vxc(:,:,:)=zero
 if (nkxc>0) kxc(:,:,:)=zero
 m_norm_min=EPSILON(0.0_dp)**2

 if (ixc==0) then
   MSG_WARNING('  Note that no xc is applied (ixc=0).')

 else

   ABI_ALLOCATE(rhonow,(nrad,nspden,ngrad*ngrad))
   ABI_ALLOCATE(rhoarr,(nrad,nspden))
   if (usexcnhat>0) then
     ABI_ALLOCATE(rhohat,(nrad,lm_size,nspden))
     rhohat(:,:,:)=rhor(:,:,:)+nhat(:,:,:)
   end if
   if (nspden==4)  then
     ABI_ALLOCATE(m_norm,(nrad))
   end if
   if (ngrad==2.and.usecore==1) then
     ABI_ALLOCATE(drhocore,(nrad))
     call nderiv_gen(drhocore,corexc,1,pawrad)
   end if
   if (ngrad==2) then
     if (option/=4.and.option/=5)  then
       ABI_ALLOCATE(gxc,(nrad,3,ylm_size,nspden_updn))
     end if
     if (option==4.or. option==5)  then
       ABI_ALLOCATE(gxc,(nrad,1,ylm_size,nspden_updn))
     end if
     gxc=zero
   end if

!  ----------------------------------------------------------------------
!  ----- Loop on the angular part and inits
!  ----------------------------------------------------------------------

!  Do loop on the angular part
   do ipts=1,npts

!    Copy the input density for this (theta,phi)
     rhoarr(:,:)=zero
     if (usexcnhat< 2) rho_=>rhor
     if (usexcnhat==2) rho_=>rhohat
     do ispden=1,nspden
       do ilm=1,ylm_size
         if (lmselect(ilm)) then
           rhoarr(1:nrad,ispden)=rhoarr(1:nrad,ispden)+rho_(1:nrad,ilm,ispden)*pawang%ylmr(ilm,ipts)
         end if
       end do
     end do
     if (usecore==1) then
       rhoarr(1:nrad,1)=rhoarr(1:nrad,1)+corexc(1:nrad)
       if (nspden==2) rhoarr(1:nrad,2)=rhoarr(1:nrad,2)+half*corexc(1:nrad)
     end if
     rhonow(1:nrad,1:nspden,1)=rhoarr(1:nrad,1:nspden)

!    GGA: compute gradient of density
     if (ngrad==2) then
       rhonow(:,:,2:4)=zero
       ABI_ALLOCATE(drho,(1:nrad))
       do ispden=1,nspden
         do ilm=1,ylm_size
           if (lmselect(ilm)) then
             call nderiv_gen(drho,rho_(:,ilm,ispden),1,pawrad)
             do ir=1,nrad
               rhonow(ir,ispden,2:4)=rhonow(ir,ispden,2:4) &
&               +drho(ir)*pawang%ylmr(ilm,ipts)*pawang%anginit(1:3,ipts) &
&               +rho_(ir,ilm,ispden)*pawang%ylmrgr(1:3,ilm,ipts)
             end do
           end if
         end do
       end do
       ABI_DEALLOCATE(drho)
       if (usecore==1) then
         do ir=1,nrad
           rhonow(ir,1,2:4)=rhonow(ir,1,2:4)+drhocore(ir)*pawang%anginit(1:3,ipts)
         end do
         if (nspden==2) then
           do ir=1,nrad
             rhonow(ir,2,2:4)=rhonow(ir,2,2:4)+half*drhocore(ir)*pawang%anginit(1:3,ipts)
           end do
         end if
       end if
     end if

!    The variable order indicates to which derivative of the energy
!    the computation must be done. Here, no derivative, except if Kxc is requested.
     order=1;if (nkxc>0) order=2

!    Allocation of mandatory arguments of drivexc
     ABI_ALLOCATE(exci,(nrad))
     ABI_ALLOCATE(vxci,(nrad,nspden_updn))
     ABI_ALLOCATE(rho_updn,(nrad,nspden_updn))

!    Allocation of optional arguments
     nd2vxc=1 !to be adjusted for the call to size_dvxc
     call size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden_updn,nvxcdgr,order)
     if (ndvxc/=0)  then
       ABI_ALLOCATE(dvxci,(nrad,ndvxc))
     end if
     if (nvxcdgr/=0)  then
       ABI_ALLOCATE(dvxcdgr,(nrad,nvxcdgr))
     end if
     if ((ixc==3 .or. xclevel==2) .and. order==3)  then
       ABI_ALLOCATE(d2vxcar,(nrad))
     end if
     if (ngrad==2)  then
       ABI_ALLOCATE(grho2_updn,(nrad,ngr2))
       ABI_ALLOCATE(dnexcdn,(nrad,nspgrad))
     end if

!    Storage of density (and gradient) in (up,dn) format
     if (nspden==1) then
       rho_updn(1:nrad,1)=rhonow(1:nrad,1,1)*half
       if (ngrad==2) then
         grho2_updn(1:nrad,1)=quarter*(rhonow(1:nrad,1,2)**2+rhonow(1:nrad,1,3)**2+rhonow(1:nrad,1,4)**2)
       end if
     else if (nspden==2) then
       rho_updn(1:nrad,1)=rhonow(1:nrad,2,1)
       rho_updn(1:nrad,2)=rhonow(1:nrad,1,1)-rhonow(1:nrad,2,1)
       if (ngrad==2) then
         grho2_updn(1:nrad,1)=rhonow(1:nrad,2,2)**2+rhonow(1:nrad,2,3)**2+rhonow(1:nrad,2,4)**2
         grho2_updn(1:nrad,2)=(rhonow(1:nrad,1,2)-rhonow(1:nrad,2,2))**2 +   &
&         (rhonow(1:nrad,1,3)-rhonow(1:nrad,2,3))**2 +   &
&         (rhonow(1:nrad,1,4)-rhonow(1:nrad,2,4))**2
         grho2_updn(1:nrad,3)=rhonow(1:nrad,1,2)**2+rhonow(1:nrad,1,3)**2+rhonow(1:nrad,1,4)**2
       end if
     else if (nspden==4) then
       m_norm(1:nrad)=sqrt(rhonow(1:nrad,2,1)**2+rhonow(1:nrad,3,1)**2+rhonow(1:nrad,4,1)**2)
       rho_updn(1:nrad,1)=(rhonow(1:nrad,1,1)+m_norm(1:nrad))*half
       rho_updn(1:nrad,2)=(rhonow(1:nrad,1,1)-m_norm(1:nrad))*half
     end if

!    Make the density positive everywhere (but do not care about gradients)
     call mkdenpos(iwarn,nrad,nspden_updn,0,rho_updn,xc_denpos)

#if defined HAVE_DFT_LIBXC
     if (ixc<0) then
       if (libxc_functionals_isgga()) then
         if (order**2 <= 1) then
           call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,     &
&           grho2_updn=grho2_updn,vxcgrho=dvxcdgr)
         else
           call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,     &
&           grho2_updn=grho2_updn,vxcgrho=dvxcdgr,dvxc=dvxci)
         end if
       else
         if (order**2 <= 1) then
           call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr)
         else
           call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,dvxc=dvxci)
         end if
       end if
     else
#endif
!      Cases with gradient
       if (xclevel==2)then
         if (order**2 <= 1 .or. ixc == 16 .or. ixc==17 .or. ixc==26 .or. ixc==27 ) then
           if (ixc /= 13) then
             call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
&             grho2_updn=grho2_updn,vxcgrho=dvxcdgr)
           else
             call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
&             grho2_updn=grho2_updn)
           end if
         else if (order /= 3) then
           if (ixc /= 13) then
             call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
&             dvxc=dvxci,grho2_updn=grho2_updn,vxcgrho=dvxcdgr)
           else
             call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
&             dvxc=dvxci,grho2_updn=grho2_updn)
           end if
         else if (order == 3) then
           if (ixc /= 13) then
             call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
&             dvxc=dvxci,d2vxc=d2vxcar,grho2_updn=grho2_updn,vxcgrho=dvxcdgr)
           else
             call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
&             dvxc=dvxci,d2vxc=d2vxcar,grho2_updn=grho2_updn)
           end if
         end if
!        Cases without gradient
       else
         if (order**2 <=1 .or. ixc >= 31 .and. ixc<=34) then
           call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr)
         else if (order==3 .and. (ixc==3 .or. ixc>=7 .and. ixc<=10)) then
           call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
&           dvxc=dvxci,d2vxc=d2vxcar)
         else
           call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nd2vxc,nvxcdgr,&
&           dvxc=dvxci)
         end if
       end if
#if defined HAVE_DFT_LIBXC
     end if
#endif

!    ----------------------------------------------------------------------
!    ----- Accumulate and store XC potential
!    ----------------------------------------------------------------------
     if (option/=3.and.option/=4) then
       if (nspden/=4) then
         do ispden=1,nspden
           vxc(1:nrad,ipts,ispden)=vxci(1:nrad,ispden)
         end do
       else
         do ir=1,nrad
           dvdn=(vxci(ir,1)+vxci(ir,2))*half
           dvdz=(vxci(ir,1)-vxci(ir,2))*half
           if(m_norm(ir)>m_norm_min)then
!            if(m_norm(ir)>rhonow(ir,1,1)*tol10+tol14)then
             factor=dvdz/m_norm(ir)
             vxc(ir,ipts,1)=dvdn+rhonow(ir,4,1)*factor
             vxc(ir,ipts,2)=dvdn-rhonow(ir,4,1)*factor
             vxc(ir,ipts,3)= rhonow(ir,2,1)*factor
             vxc(ir,ipts,4)=-rhonow(ir,3,1)*factor
           else
             vxc(ir,ipts,1:2)=dvdn
             vxc(ir,ipts,3:4)=zero
           end if
         end do
       end if

!      For GGAs, additional terms appear
       if(ngrad==2.and.ixc/=13)then
         dnexcdn(1:nrad,1:nspden_updn)=vxci(1:nrad,1:nspden_updn)
!        Treat explicitely spin up, spin down and total spin for spin-polarized
         do ii=1,3
           if(nspden_updn==1.and.ii>=2)exit !exit when ii=1 is finished if non-spin-polarized
           do ir=1,nrad
!            If the norm of the gradient vanishes, then the different terms vanishes
             if(grho2_updn(ir,ii)<1.0d-24) then
               dnexcdn(ir,ii+nspden_updn)=zero;cycle
             end if
!            Compute the derivative of n.e_xc wrt spin up, spin down, or total density
             if(nspden_updn==1)then
               dnexcdn(ir,ii+nspden_updn)=half*dvxcdgr(ir,1) !Definition of dvxcdgr changed in v3.3
               if (nvxcdgr==3) dnexcdn(ir,ii+nspden_updn)=dnexcdn(ir,ii+nspden_updn)+dvxcdgr(ir,3)
             else if(nspden_updn==2)then
               if (nvxcdgr==3) then
                 dnexcdn(ir,ii+nspden_updn)=dvxcdgr(ir,ii)
               else if (ii/=3) then
                 dnexcdn(ir,ii+nspden_updn)=dvxcdgr(ir,ii)
               else if (ii==3) then
                 dnexcdn(ir,ii+nspden_updn)=zero
               end if
             end if
           end do
         end do
         call xcmult(dnexcdn,nrad,ngrad,nspden_eff,nspgrad,rhonow)
         factor=one;if (nspden_updn==1) factor=half
         if (option/=4.and.option/=5) then
           factor=factor*four_pi
!          Accumulate moments of gxc
           do ispden=1,nspden_updn
             do ilm=1,ylm_size
               do ii=1,3
                 gxc(1:nrad,ii,ilm,ispden)=gxc(1:nrad,ii,ilm,ispden)+rhonow(1:nrad,ispden,1+ii) &
&                 *pawang%ylmr(ilm,ipts)*pawang%angwgth(ipts)*factor
               end do
             end do
           end do
         else
           do ispden=1,nspden_updn
             gxc(1:nrad,1,1,ispden)=factor*rhonow(1:nrad,ispden,2)
           end do
         end if
       end if

     end if !option

!    ----------------------------------------------------------------------
!    ----- Accumulate and store XC energy
!    ----------------------------------------------------------------------
     if (option/=1.and.option/=5) then
       ABI_ALLOCATE(ff,(nrad))
       ff(1:nrad)=rhoarr(1:nrad,1)*exci(1:nrad)*pawrad%rad(1:nrad)**2
       call simp_gen(enxcr,ff,pawrad)
       if (option/=4) enxc=enxc+enxcr*pawang%angwgth(ipts)
       if (option==4) enxc=enxc+enxcr
       ABI_DEALLOCATE(ff)
     end if

!    ----------------------------------------------------------------------
!    ----- If LDA, accumulate and store XC double-counting energy here
!    ----------------------------------------------------------------------
     if (ngrad==1.and.(option==0.or.option==2)) then
!      (Eventually) Modify density for this (theta,phi)
       if (usexcnhat==1) then
         rho_=>rhohat
         rhoarr(:,:)=zero
         do ispden=1,nspden
           do ilm=1,ylm_size
             if (lmselect(ilm)) then
               rhoarr(1:nrad,ispden)=rhoarr(1:nrad,ispden)+rho_(1:nrad,ilm,ispden)*pawang%ylmr(ilm,ipts)
             end if
           end do
         end do
       else if (usecore==1) then
         rhoarr(1:nrad,1)=rhoarr(1:nrad,1)-corexc(1:nrad)
         if (nspden==2) rhoarr(1:nrad,2)=rhoarr(1:nrad,2)-half*corexc(1:nrad)
       end if
!      Compute integral of Vxc*rho
       ABI_ALLOCATE(ff,(nrad))
       if (nspden/=4) then
         ff(:)=vxc(:,ipts,1)*rhoarr(:,nspden)
         if (nspden==2) ff(:)=ff(:)+vxc(:,ipts,2)*(rhoarr(:,1)-rhoarr(:,2))
       else
         ff(:)=half*(vxc(:,ipts,1)*(rhoarr(:,1)+rhoarr(:,4))+vxc(:,ipts,2)*(rhoarr(:,1)-rhoarr(:,4))) &
&         +vxc(:,ipts,3)*rhoarr(:,2)-vxc(:,ipts,4)*rhoarr(:,3)
       end if
       ff(:)=ff(:)*pawrad%rad(:)**2
       call simp_gen(vxcrho,ff,pawrad)
       enxcdc=enxcdc+vxcrho*pawang%angwgth(ipts)
       ABI_DEALLOCATE(ff)
     end if

!    ----------------------------------------------------------------------
!    ----- Accumulate and store XC kernel
!    ----------------------------------------------------------------------
     if (nkxc>0) then
!      LDA
       if(nkxc/=23)then
         if (ndvxc==15) then
           if (nkxc>=3) then
             kxc(1:nrad,ipts,1)=dvxci(1:nrad,1)+dvxci(1:nrad,9)
             kxc(1:nrad,ipts,2)=dvxci(1:nrad,10)
             kxc(1:nrad,ipts,3)=dvxci(1:nrad,2)+dvxci(1:nrad,11)
             if (nkxc>3) kxc(1:nrad,ipts,4:nkxc)=zero
           else
             kxc(1:nrad,ipts,1)=half*(dvxci(1:nrad,1)+dvxci(1:nrad,9)+dvxci(1:nrad,10))
             if (nkxc>1) kxc(1:nrad,ipts,2:nkxc)=zero
           end if
         else if (nkxc<=ndvxc) then
           kxc(1:nrad,ipts,1:nkxc)=dvxci(1:nrad,1:nkxc)
         else
           if (ndvxc>0) kxc(1:nrad,ipts,1:ndvxc)=dvxci(1:nrad,1:ndvxc)
           kxc(1:nrad,ipts,ndvxc+1:nkxc)=zero
         end if
!        GGA
       else if(nkxc==23)then
         if (ndvxc==15) then
           kxc(1:nrad,ipts,1:15)=dvxci(1:nrad,1:15)
         else
           kxc(1:nrad,ipts,1:ndvxc)=dvxci(1:nrad,1:ndvxc)
           kxc(1:nrad,ipts,ndvxc+1:15)=zero
         end if
         do ispden=1,nspden_updn
           do ii=1,4
             kxc(1:nrad,ipts,13+ispden+2*ii)=rhonow(1:nrad,ispden,ii)
           end do
         end do
       end if
     end if

!    Deallocate temporary memory space
     ABI_DEALLOCATE(exci)
     ABI_DEALLOCATE(rho_updn)
     ABI_DEALLOCATE(vxci)
     if (allocated(dvxci))  then
       ABI_DEALLOCATE(dvxci)
     end if
     if (allocated(dvxcdgr))  then
       ABI_DEALLOCATE(dvxcdgr)
     end if
     if (allocated(d2vxcar))  then
       ABI_DEALLOCATE(d2vxcar)
     end if
     if (allocated(dnexcdn))  then
       ABI_DEALLOCATE(dnexcdn)
     end if
     if (allocated(grho2_updn))  then
       ABI_DEALLOCATE(grho2_updn)
     end if

!    ----------------------------------------------------------------------
!    ----- End of the loop on npts (angular part)
!    ----------------------------------------------------------------------
   end do

!  Deallocate memory
   if (nspden==4)  then
     ABI_DEALLOCATE(m_norm)
   end if
   if (ngrad==2.and.usecore==1)  then
     ABI_DEALLOCATE(drhocore)
   end if
   ABI_DEALLOCATE(rhonow)

!  ----------------------------------------------------------------------
!  ----- If GGA, modify potential with term from density gradient
!  ----------------------------------------------------------------------
   if (ngrad==2) then
     ABI_ALLOCATE(dgxc,(nrad))
!    Need to multiply gxc by 2 in the non-polarised case
     factor=one;if (nspden==1) factor=two
     if (option/=4.and.option/=5) then
!      Compute divergence of gxc and substract it from Vxc
       do ispden=1,nspden_updn
         do ilm=1,ylm_size
           do ii=1,3
             call nderiv_gen(dgxc,gxc(:,ii,ilm,ispden),1,pawrad)
             do ipts=1,npts
               vxc(1:nrad,ipts,ispden)=vxc(1:nrad,ipts,ispden) &
&               -factor*(dgxc(1:nrad)*pawang%anginit(ii,ipts)*pawang%ylmr(ilm,ipts) &
&               +gxc(1:nrad,ii,ilm,ispden)*pawang%ylmrgr(ii,ilm,ipts))
             end do
           end do
         end do
       end do
     else ! option==4 or option==5
       do ispden=1,nspden_updn
         call nderiv_gen(dgxc,gxc(:,1,1,ispden),1,pawrad)
         vxc(2:nrad,1,ispden)=vxc(2:nrad,1,ispden) &
&         -factor*(dgxc(2:nrad)+two*gxc(2:nrad,1,1,ispden)/pawrad%rad(2:nrad))
         call deducer0(vxc(:,1,ispden),nrad,pawrad)
       end do
     end if
     ABI_DEALLOCATE(dgxc)
     ABI_DEALLOCATE(gxc)
   end if ! GGA

!  ----------------------------------------------------------------------
!  ----- If GGA, accumulate and store XC double-counting energy here
!  ----------------------------------------------------------------------
   if (ngrad==2.and.(option==0.or.option==2)) then
     ABI_ALLOCATE(ff,(nrad))
     do ipts=1,npts !  Do loop on the angular part
!      Compute density for this (theta,phi)
       rhoarr(:,:)=zero
       if (usexcnhat==0) rho_=>rhor
       if (usexcnhat/=0) rho_=>rhohat
       do ispden=1,nspden
         do ilm=1,ylm_size
           if (lmselect(ilm)) then
             rhoarr(1:nrad,ispden)=rhoarr(1:nrad,ispden)+rho_(1:nrad,ilm,ispden)*pawang%ylmr(ilm,ipts)
           end if
         end do
       end do
!      Compute integral of Vxc*rho
       if (nspden/=4) then
         ff(:)=vxc(:,ipts,1)*rhoarr(:,nspden)
         if (nspden==2) ff(:)=ff(:)+vxc(:,ipts,2)*(rhoarr(:,1)-rhoarr(:,2))
       else
         ff(:)=half*(vxc(:,ipts,1)*(rhoarr(:,1)+rhoarr(:,4))+vxc(:,ipts,2)*(rhoarr(:,1)-rhoarr(:,4))) &
&         +vxc(:,ipts,3)*rhoarr(:,2)-vxc(:,ipts,4)*rhoarr(:,3)
       end if
       ff(:)=ff(:)*pawrad%rad(:)**2
       call simp_gen(vxcrho,ff,pawrad)
       enxcdc=enxcdc+vxcrho*pawang%angwgth(ipts)
     end do ! End of the loop on npts (angular part)
     ABI_DEALLOCATE(ff)
   end if ! option

!  ----------------------------------------------------------------------
!  ----- End
!  ----------------------------------------------------------------------
!  Add the four*pi factor of the Exc and Excdc angular integration
   if (option/=1.and.option/=5) enxc=enxc*four_pi
   if (option==0.or.option==2) enxcdc=enxcdc*four_pi

!  Final memory deallocation
   nullify(rho_)
   ABI_DEALLOCATE(rhoarr)
   if (usexcnhat>0)  then
     ABI_DEALLOCATE(rhohat)
   end if

!  ------------------------------------
!  End IF a xc part has to be computed
 end if

 call timab(81,2,tsec)

 DBG_EXIT("COLL")

 end subroutine pawxc
!!***
