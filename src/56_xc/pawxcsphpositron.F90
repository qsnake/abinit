!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxcsphpositron
!! NAME
!! pawxcsphpositron
!!
!! FUNCTION
!! PAW only
!! Compute electron-positron XC energy and potential for spherical densities rho_el(r) rho_pos(r)
!! Driver of XC functionals. LDA and GGA
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!! This routine has been written from rhohxc_coll
!!
!! INPUTS
!!  calctype=type of electron-positron calculation:
!!           calctype=1 : positron in electronic density
!!           calctype=2 : electrons in positronic density
!!  ixcpositron= choice of elctron-positron exchange-correlation scheme
!!  nrad= dimension of the radial mesh
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  posdensity0_limit=True if we are in the zero positron density limit
!!  rho(nrad,lm_size)=electron (or positron) density in real space
!!                    Contents depends on calctype value:
!!                    calctype=1: rho is the positronic density
!!                    calctype=2: rho is the electronic density
!!  rho_ep(nrad,lm_size)=electron (or positron) density in real space
!!                      Contents depends on calctype value:
!!                      calctype=1: rho_ep is the electronic density
!!                      calctype=2: rho_ep is the positronic density
!!
!! OUTPUT
!!  fxc(nrad)= electron-positron XC energy per unit volume
!!  vxce(nrad)= electron-positron XC potential for the electron
!!  vxcp(nrad)= electron-positron XC potential for the positron
!!
!! PARENTS
!!      pawxcmpositron
!!
!! CHILDREN
!!      deducer0,leave_new,nderiv_gen,wrtout,xcpositron
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawxcsphpositron(calctype,fxc,ixcpositron,nrad,pawrad,posdensity0_limit,rho,rho_ep,vxce,vxcp)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

 use m_radmesh,   only :  nderiv_gen, deducer0

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawxcsphpositron'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_56_xc, except_this_one => pawxcsphpositron
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: calctype,ixcpositron,nrad
 logical,intent(in) :: posdensity0_limit
 type(pawrad_type),intent(in) :: pawrad
!arrays
 real(dp),intent(in) :: rho(nrad),rho_ep(nrad)
 real(dp),intent(out) :: fxc(nrad),vxce(nrad),vxcp(nrad)

!Local variables-------------------------------
!scalars
 integer :: ngr
 character(len=500) :: message
!arrays
 real(dp),allocatable :: dff(:),rhograd(:),rhograd2(:),vxcegr(:)

! *************************************************************************

 if(nrad/=pawrad%mesh_size)then
   write(message, '(4a)' ) ch10,&
&   ' pawxcsphpositron :  BUG -',ch10,&
&   '  nrad is not equal to radial mesh size !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Need gradient of density for GGA
 ngr=0;if (ixcpositron==3.or.ixcpositron==31) ngr=nrad
 ABI_ALLOCATE(rhograd,(ngr))
 ABI_ALLOCATE(rhograd2,(ngr))
 ABI_ALLOCATE(vxcegr,(ngr))
 if (ngr==nrad) then
   if (calctype==1) then
     call nderiv_gen(rhograd,rho_ep,1,pawrad)
   else if (calctype==2) then
     call nderiv_gen(rhograd,rho,1,pawrad)
   end if
   rhograd2(:)=rhograd(:)**2
 end if

!---- Computation of Fxc and Vxc for the positron
!rho    is the positronic density
!rho_ep is the electronic density
 if (calctype==1) then
   call xcpositron(fxc,rhograd2,ixcpositron,ngr,nrad,posdensity0_limit,rho_ep,rho,vxce,vxcegr,vxcp)

!  ---- Computation of Exc and Vxc for the electron
!  rho    is the electronic density
!  rho_ep is the positronic density
 else if (calctype==2) then
   call xcpositron(fxc,rhograd2,ixcpositron,ngr,nrad,posdensity0_limit,rho,rho_ep,vxce,vxcegr,vxcp)
 end if

 ABI_DEALLOCATE(rhograd2)

!---- GGA - gradient corrections
 if (ngr==nrad) then
   ABI_ALLOCATE(dff,(nrad))
   vxcegr(1:nrad)=vxcegr(1:nrad)*rhograd(1:nrad)
   call nderiv_gen(dff,vxcegr,1,pawrad)
   vxcp(2:nrad)=vxcp(2:nrad)-(dff(2:nrad)+two*vxcegr(2:nrad)/pawrad%rad(2:nrad))
   call deducer0(vxcp,nrad,pawrad)
   ABI_DEALLOCATE(dff)
 end if

 ABI_DEALLOCATE(vxcegr)
 ABI_DEALLOCATE(rhograd)

 end subroutine pawxcsphpositron
!!***
