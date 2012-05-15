!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxcpositron
!! NAME
!! pawxcpositron
!!
!! FUNCTION
!! PAW only
!! Compute electron-positron correlation potential and energies inside a PAW sphere
!! LDA ONLY - USE THE DENSITY OVER A WHOLE SPHERICAL GRID (r,theta,phi)
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
!!  calctype=type of electronpositron calculation:
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
!!    vxc(pawrad%mesh_size,pawang%angl_size,nspden)=xc potential
!!       (spin up in 1st half and spin-down in 2nd half if nspden=2)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!
!! PARENTS
!!      pawdenpot
!!
!! CHILDREN
!!      mkdenpos,simp_gen,wrtout,xcpositron
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawxcpositron(calctype,corexc,enxc,enxcdc,ixcpositron,lm_size,lmselect,lmselect_ep,&
&                        nhat,nhat_ep,nspden,option,pawang,pawrad,posdensity0_limit,&
&                        rhor,rhor_ep,usecore,usexcnhat,vxc,xc_denpos)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors
 
 use m_radmesh,   only : simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawxcpositron'
 use interfaces_14_hidewrite
 use interfaces_56_xc, except_this_one => pawxcpositron
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: calctype,ixcpositron,lm_size,nspden,option,usecore,usexcnhat
 logical,intent(in) :: posdensity0_limit
 real(dp),intent(in) :: xc_denpos
 real(dp),intent(out) :: enxc,enxcdc
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size),lmselect_ep(lm_size)
 real(dp),intent(in) :: corexc(pawrad%mesh_size)
 real(dp),intent(in) :: nhat(pawrad%mesh_size,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: nhat_ep(pawrad%mesh_size,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor(pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(in) :: rhor_ep(pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(out) :: vxc(pawrad%mesh_size,pawang%angl_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ilm,ipts,iwarn,iwarnp,ngr,ngrad,npts,nrad,order
 real(dp) :: enxcr,vxcrho
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: ff(:),fxci(:),grho2(:),rhoarr(:),rhoarr_ep(:),rhoarrdc(:),vxci(:),vxci_ep(:),vxcigr(:)

! *************************************************************************

 DBG_ENTER("COLL")

!----- Check options
 if(ixcpositron==3.or.ixcpositron==31) then
   msg='  GGA is not implemented (use pawxcdev/=0) !'
   MSG_ERROR(msg)
 end if
 if(calctype/=1.and.calctype/=2) then
   msg='  Invalid value for calctype'
   MSG_BUG(msg)
 end if
 if(pawang%angl_size==0) then
   msg='  pawang%angl_size=0 !'
   MSG_BUG(msg)
 end if
 if(.not.associated(pawang%ylmr)) then
   msg='  pawang%ylmr must be allocated !'
   MSG_BUG(msg)
 end if

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

!Initialization and constants
 iwarn=0;iwarnp=1
 nrad=pawrad%mesh_size
 npts=pawang%angl_size
 order=1;ngr=0;ngrad=1 ! only LDA here !

!Initializations of output arrays
 if (option/=1) enxc=zero
 if (option==0.or.option==2) enxcdc=zero
 if (option<3) vxc(:,:,:)=zero

 if (ixcpositron==0) then ! No xc at all is applied (usually for testing)
   write(msg, '(a,a,a,a)' ) ch10,&
&   ' pawxcpositron : WARNING -',ch10,&
&   '  Note that no xc is applied (ixcpositron=0).'
   call wrtout(std_out,msg,'COLL')
   return
 end if

!Allocations
 ABI_ALLOCATE(fxci,(nrad))
 ABI_ALLOCATE(vxci,(nrad))
 ABI_ALLOCATE(rhoarr,(nrad))
 ABI_ALLOCATE(rhoarr_ep,(nrad))
 if (option==0.or.option==2)  then
   ABI_ALLOCATE(rhoarrdc,(nrad))
 end if

!----------------------------------------------------------------------
!----- Loop on the angular part
 do ipts=1,npts

!  ----------------------------------------------------------------------
!  ----- Build several densities
!  ----------------------------------------------------------------------

!  Eventually add compensation density to input density
   rhoarr=zero;rhoarr_ep=zero
   if (usexcnhat==2) then
     do ilm=1,lm_size
       if (lmselect(ilm)) &
&       rhoarr(:)=rhoarr(:)+(rhor(:,ilm,1)+nhat(:,ilm,1))*pawang%ylmr(ilm,ipts)
     end do
     do ilm=1,lm_size
       if (lmselect_ep(ilm)) &
&       rhoarr_ep(:)=rhoarr_ep(:)+(rhor_ep(:,ilm,1)+nhat_ep(:,ilm,1))*pawang%ylmr(ilm,ipts)
     end do
   else
     do ilm=1,lm_size
       if (lmselect(ilm)) rhoarr(:)=rhoarr(:)+rhor(:,ilm,1)*pawang%ylmr(ilm,ipts)
     end do
     do ilm=1,lm_size
       if (lmselect_ep(ilm)) rhoarr_ep(:)=rhoarr_ep(:)+rhor_ep(:,ilm,1)*pawang%ylmr(ilm,ipts)
     end do
   end if

!  Store density for use in double-counting term
   if (option==0.or.option==2) rhoarrdc(:)=rhoarr(:)

!  Eventually add core density
   if (usecore==1) then
     if (calctype==1) rhoarr_ep(:)=rhoarr_ep(:)+corexc(:)
     if (calctype==2) rhoarr   (:)=rhoarr   (:)+corexc(:)
   end if

!  Make the densities positive
   if (calctype==1) then
     if (.not.posdensity0_limit) then
       call mkdenpos(iwarnp,nrad,1,1,rhoarr,xc_denpos)
     end if
     call mkdenpos(iwarn ,nrad,1,1,rhoarr_ep,xc_denpos)
   else if (calctype==2) then
     call mkdenpos(iwarn ,nrad,1,1,rhoarr,xc_denpos)
     if (.not.posdensity0_limit) then
       call mkdenpos(iwarnp,nrad,1,1,rhoarr_ep,xc_denpos)
     end if
   end if

!  ----------------------------------------------------------------------
!  ----- Compute XC data
!  ----------------------------------------------------------------------

!  electron-positron correlation for the positron
   ABI_ALLOCATE(vxci_ep,(nrad))
   ABI_ALLOCATE(vxcigr,(ngr))
   ABI_ALLOCATE(grho2,(ngr))
   if (calctype==1) then
     call xcpositron(fxci,grho2,ixcpositron,ngr,nrad,posdensity0_limit,rhoarr_ep,rhoarr,vxci_ep,vxcigr,vxci)
   else if (calctype==2) then
     call xcpositron(fxci,grho2,ixcpositron,ngr,nrad,posdensity0_limit,rhoarr,rhoarr_ep,vxci,vxcigr,vxci_ep)
   end if
   ABI_DEALLOCATE(vxci_ep)
   ABI_DEALLOCATE(vxcigr)
   ABI_DEALLOCATE(grho2)

!  ----------------------------------------------------------------------
!  ----- Accumulate and store XC potential
!  ----------------------------------------------------------------------
   if (option<3) then
     vxc(:,ipts,1)=vxci(:)
     if (nspden>=2) vxc(:,ipts,2)=vxci(:)
     if (nspden==4) vxc(:,ipts,3:4)=zero
   end if

!  ----------------------------------------------------------------------
!  ----- Accumulate and store XC energies
!  ----------------------------------------------------------------------

!  ----- Calculate Exc term
   if (option/=1) then
     ABI_ALLOCATE(ff,(nrad))
     ff(:)=fxci(:)*pawrad%rad(:)**2
     call simp_gen(enxcr,ff,pawrad)
     ABI_DEALLOCATE(ff)
     if (option/=4) enxc=enxc+enxcr*pawang%angwgth(ipts)
     if (option==4) enxc=enxc+enxcr
   end if

!  ----- Calculate Excdc double counting term
   if (option==0.or.option==2) then
     if (usexcnhat==1) then
       do ilm=1,lm_size
         if (lmselect(ilm)) then
           rhoarrdc(:)=rhoarrdc(:)+nhat(:,ilm,1)*pawang%ylmr(ilm,ipts)
         end if
       end do
     end if
     ABI_ALLOCATE(ff,(nrad))
     ff(:)=vxci(:)*rhoarrdc(:)*pawrad%rad(:)**2
     call simp_gen(vxcrho,ff,pawrad)
     ABI_DEALLOCATE(ff)
     enxcdc=enxcdc+vxcrho*pawang%angwgth(ipts)
   end if

!  ---------------------------------------------------
!  ----- End of the loop on npts (angular part)
 end do

!Add the four*pi factor of the angular integration
 if (option/=1) enxc=enxc*four_pi
 if (option==0.or.option==2) enxcdc=enxcdc*four_pi

!Deallocations
 ABI_DEALLOCATE(fxci)
 ABI_DEALLOCATE(vxci)
 ABI_DEALLOCATE(rhoarr)
 ABI_DEALLOCATE(rhoarr_ep)
 if (option==0.or.option==2)  then
   ABI_DEALLOCATE(rhoarrdc)
 end if

 DBG_EXIT("COLL")

 end subroutine pawxcpositron
!!***

