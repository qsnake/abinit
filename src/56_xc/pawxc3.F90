!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxc3
!! NAME
!! pawxc3
!!
!! FUNCTION
!! PAW only
!! Compute first-order change of XC potential and contribution to
!! 2nd-order change of XC energy inside a PAW sphere.
!! LDA ONLY - USE THE DENSITY OVER A WHOLE SPHERICAL GRID (r,theta,phi)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!! This routine has been written from rhohxc
!!
!! INPUTS
!!  corexc1(cplex_den*pawrad%mesh_size)=first-order change of core density on radial grid
!!  cplex_den= if 1, 1st-order densities are REAL, if 2, COMPLEX
!!  cplex_vxc= if 1, 1st-order XC potential is complex, if 2, COMPLEX
!!  kxc(pawrad%mesh_size,pawang%angl_size,nkxc)=GS xc kernel
!!  lm_size=size of density array rhor (see below)
!!  lmselect(lm_size)=select the non-zero LM-moments of input density rhor1
!!  nhat1(cplex_den*pawrad%mesh_size,lm_size,nspden)=first-order change of compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  option=0  compute both 2nd-order XC energy and 1st-order potential
!!         1  compute only 1st-order XC potential
!!         2  compute only 2nd-order XC energy, XC potential is temporary computed here
!!         3  compute only 2nd-order XC energy, XC potential is input in vxc1(:)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  rhor1(cplex_den*pawrad%mesh_size,lm_size,nspden)=first-order change of density
!!  usecore= 1 if core density has to be used in Exc/Vxc ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in d2Exc only
!!             2 if compensation density (nhat) has to be used in d2Exc and Vxc1
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  == if option=0 or 2 or 3 ==
!!    d2enxc   =returned exchange-cor. contribution to 2nd-order XC energy
!!    d2enxc_im=returned IMAGINARY PART of exchange-cor. contribution to 2nd-order XC energy
!!              (optional argument)
!!
!! SIDE EFFECTS
!!    vxc1(cplex_vxc*pawrad%mesh_size,pawang%angl_size,nspden)=1st-order XC potential
!!      Output if option==0 or 1
!!      Unused if option==2
!!      Input  if option==3
!!
!! PARENTS
!!      pawdenpot,pawnstd2e,pawxcm3
!!
!! CHILDREN
!!      simp_gen,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawxc3(corexc1,cplex_den,cplex_vxc,d2enxc,kxc,lm_size,lmselect,nhat1,nkxc,nspden,&
&                 option,pawang,pawrad,rhor1,usecore,usexcnhat,vxc1,xclevel,&
&                 d2enxc_im) ! optional

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_radmesh,          only : simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawxc3'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex_den,cplex_vxc,lm_size,nkxc,nspden,option,usecore,usexcnhat,xclevel
 real(dp),intent(out) :: d2enxc
 real(dp),intent(out),optional :: d2enxc_im
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size)
 real(dp),intent(in) :: corexc1(cplex_den*pawrad%mesh_size)
 real(dp),intent(in) :: kxc(pawrad%mesh_size,pawang%angl_size,nkxc)
 real(dp),intent(in) :: nhat1(cplex_den*pawrad%mesh_size,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor1(cplex_den*pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(inout) :: vxc1(cplex_vxc*pawrad%mesh_size,pawang%angl_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ilm,ipts,ir,ispden,jr,npts,nrad
 logical :: need_impart
 real(dp) :: rho_dn,rho_up,rhoim_dn,rhoim_up,ro11i,ro11r,ro12i,ro12r,ro21i,ro21r,ro22i,ro22r
 real(dp) :: v11i,v11r,v12i,v12r,v21i,v21r,v22i,v22r,vxcrho
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: ff(:),gg(:),rhoarr(:,:),vxc1i(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(81,1,tsec)

!----------------------------------------------------------------------
!----- Check options
!----------------------------------------------------------------------

 if(option<0.or.option>3) then
   msg='  Wrong option !'
   MSG_BUG(msg)
 end if
 if(option/=3.and.nkxc/=2*min(nspden,2)-1) then
   msg='  nkxc must be 1 or 3 !'
   MSG_BUG(msg)
 end if
 if(xclevel==2) then
   msg='  GGA is not implemented !'
   MSG_ERROR(msg)
 end if
 if(nspden==4.and.option/=3) then
   msg='  nspden=4 not implemented (for vxc) !'
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

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

 nrad=pawrad%mesh_size
 npts=pawang%angl_size
 need_impart=present(d2enxc_im)
 if (option/=1) then
   d2enxc=zero
   if (need_impart) d2enxc_im=zero
 end if
 if (option<=1) vxc1(:,:,:)=zero
 ABI_ALLOCATE(rhoarr,(cplex_den*nrad,nspden))
 ABI_ALLOCATE(vxc1i,(cplex_vxc*nrad,nspden))

!----------------------------------------------------------------------
!----- Loop on the angular part and inits
!----------------------------------------------------------------------

!Do loop on the angular part (theta,phi)
 do ipts=1,npts

!  Copy the input density for this (theta,phi)
   rhoarr(:,:)=zero
   if (usexcnhat==0.or.usexcnhat==1) then
     do ispden=1,nspden
       do ilm=1,lm_size
         if (lmselect(ilm)) rhoarr(:,ispden)=rhoarr(:,ispden) &
&         +rhor1(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
   else
     do ispden=1,nspden
       do ilm=1,lm_size
         if (lmselect(ilm)) rhoarr(:,ispden)=rhoarr(:,ispden) &
&         +(rhor1(:,ilm,ispden)+nhat1(:,ilm,ispden))*pawang%ylmr(ilm,ipts)
       end do
     end do
   end if

   if (usecore==1) then
     rhoarr(:,1)=rhoarr(:,1)+corexc1(:)
     if (nspden==2) rhoarr(:,2)=rhoarr(:,2)+half*corexc1(:)
   end if

!  
!  ----------------------------------------------------------------------
!  ----- Accumulate and store 1st-order change of XC potential
!  ----------------------------------------------------------------------

   if (option/=3) then
!    Non-spin-polarized
     if(nspden==1)then
       if (cplex_vxc==1) then
         if (cplex_den==1) then  ! cplex_vxc==1 and cplex_den==1
           vxc1i(1:nrad,1)=kxc(1:nrad,ipts,1)*rhoarr(1:nrad,1)
         else                    ! cplex_vxc==1 and cplex_den==2
           do ir=1,nrad
             vxc1i(ir,1)=kxc(ir,ipts,1)*rhoarr(2*ir-1,1)
           end do
         end if
       else
         if (cplex_den==1) then  ! cplex_vxc==2 and cplex_den==1
           do ir=1,nrad
             vxc1i(2*ir-1,1)=kxc(ir,ipts,1)*rhoarr(ir,1)
             vxc1i(2*ir  ,1)=zero
           end do
         else                    ! cplex_vxc==2 and cplex_den==2
           do ir=1,nrad
             vxc1i(2*ir-1,1)=kxc(ir,ipts,1)*rhoarr(2*ir-1,1)
             vxc1i(2*ir  ,1)=kxc(ir,ipts,1)*rhoarr(2*ir  ,1)
           end do
         end if
       end if

!      Spin-polarized
     else
       if (cplex_vxc==1) then
         if (cplex_den==1) then  ! cplex_vxc==1 and cplex_den==1
           do ir=1,nrad
             rho_up=rhoarr(ir,2);rho_dn=rhoarr(ir,1)-rho_up
             vxc1i(ir,1)=kxc(ir,ipts,1)*rho_up+kxc(ir,ipts,2)*rho_dn
             vxc1i(ir,2)=kxc(ir,ipts,2)*rho_up+kxc(ir,ipts,3)*rho_dn
           end do
         else                    ! cplex_vxc==1 and cplex_den==2
           do ir=1,nrad
             jr=2*ir-1
             rho_up=rhoarr(jr,2);rho_dn=rhoarr(jr,1)-rho_up
             vxc1i(ir,1)=kxc(ir,ipts,1)*rho_up+kxc(ir,ipts,2)*rho_dn
             vxc1i(ir,2)=kxc(ir,ipts,2)*rho_up+kxc(ir,ipts,3)*rho_dn
           end do
         end if
       else
         if (cplex_den==1) then  ! cplex_vxc==2 and cplex_den==1
           do ir=1,nrad
             jr=2*ir-1
             rho_up=rhoarr(ir,2);rho_dn=rhoarr(ir,1)-rho_up
             vxc1i(jr,1)=kxc(ir,ipts,1)*rho_up+kxc(ir,ipts,2)*rho_dn
             vxc1i(jr,2)=kxc(ir,ipts,2)*rho_up+kxc(ir,ipts,3)*rho_dn
           end do
         else                    ! cplex_vxc==2 and cplex_den==2
           do ir=1,nrad
             jr=2*ir
             rho_up  =rhoarr(jr-1,2);rho_dn  =rhoarr(jr-1,1)-rho_up
             rhoim_up=rhoarr(jr  ,2);rhoim_dn=rhoarr(jr  ,1)-rhoim_up
             vxc1i(jr-1,1)=kxc(ir,ipts,1)*rho_up  +kxc(ir,ipts,2)*rho_dn
             vxc1i(jr  ,1)=kxc(ir,ipts,1)*rhoim_up+kxc(ir,ipts,2)*rhoim_dn
             vxc1i(jr-1,2)=kxc(ir,ipts,2)*rho_up  +kxc(ir,ipts,3)*rho_dn
             vxc1i(jr  ,2)=kxc(ir,ipts,2)*rhoim_up+kxc(ir,ipts,3)*rhoim_dn
           end do
         end if
       end if
     end if

     if (option<=1) then
       vxc1(1:cplex_vxc*nrad,ipts,1:nspden)=vxc1i(1:cplex_vxc*nrad,1:nspden)
     end if

   else  ! option==3
     vxc1i(1:cplex_vxc*nrad,1:nspden)=vxc1(1:cplex_vxc*nrad,ipts,1:nspden)
   end if

!  ----------------------------------------------------------------------
!  ----- Accumulate and store 2nd-order change of XC energy
!  ----------------------------------------------------------------------
   if (option/=1) then

!    For usexnhat=1 particular case, add now compensation density
     if (usexcnhat==1) then
       do ispden=1,nspden
         do ilm=1,lm_size
           if (lmselect(ilm)) rhoarr(:,ispden)=rhoarr(:,ispden)+nhat1(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
         end do
       end do
     end if

!    ----- Calculate d2Exc=Int[Vxc^(1)^*(r).n^(1)(r).dr]
     ABI_ALLOCATE(ff,(nrad))
     if (need_impart) ABI_ALLOCATE(gg,(nrad))

!    COLLINEAR MAGNETISM
     if (nspden/=4) then
       if (cplex_vxc==1.and.cplex_den==1) then       ! cplex_vxc==1 and cplex_den==1
         ff(:)=vxc1i(:,1)*rhoarr(:,nspden)
         if (nspden==2) ff(:)=ff(:)+vxc1i(:,2)*(rhoarr(:,1)-rhoarr(:,2))
         if (need_impart) gg(:)=zero
       else if (cplex_vxc==2.and.cplex_den==2) then  ! cplex_vxc==2 and cplex_den==2
         if (.not.need_impart) then      ! Real part only
           do ir=1,nrad
             jr=2*ir;v11r=vxc1i(jr-1,1);v11i=vxc1i(jr,1)
             ro11r=rhoarr(jr-1,nspden);ro11i=rhoarr(jr,nspden)
             ff(ir)=v11r*ro11r+v11i*ro11i
           end do
           if (nspden==2) then
             do ir=1,nrad
               jr=2*ir;v22r=vxc1i(jr-1,2);v22i=vxc1i(jr,2)
               ro22r=rhoarr(jr-1,1)-rhoarr(jr-1,2)
               ro22i=rhoarr(jr  ,1)-rhoarr(jr  ,2)
               ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
             end do
           end if
         else
           do ir=1,nrad                  ! Real and imaginary parts
             jr=2*ir;v11r=vxc1i(jr-1,1);v11i=vxc1i(jr,1)
             ro11r=rhoarr(jr-1,nspden);ro11i=rhoarr(jr,nspden)
             ff(ir)=v11r*ro11r+v11i*ro11i
             gg(ir)=v11r*ro11i-v11i*ro11r
           end do
           if (nspden==2) then
             do ir=1,nrad
               jr=2*ir;v22r=vxc1i(jr-1,2);v22i=vxc1i(jr,2)
               ro22r=rhoarr(jr-1,1)-rhoarr(jr-1,2)
               ro22i=rhoarr(jr  ,1)-rhoarr(jr  ,2)
               ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
               gg(ir)=gg(ir)+v22r*ro22i-v22i*ro22r
             end do
           end if
         end if
       else                                          ! other cases for cplex_vxc and cplex_den
         v11i=zero;ro11i=zero
         do ir=1,nrad
           jr=cplex_vxc*(ir-1)+1
           v11r=vxc1i(jr,1);if (cplex_vxc==2) v11i=vxc1i(jr+1,1)
           jr=cplex_den*(ir-1)+1
           ro11r=rhoarr(jr,nspden);if (cplex_den==2) ro11i=rhoarr(jr+1,nspden)
           ff(ir)=v11r*ro11r+v11i*ro11i
           if (need_impart) gg(ir)=v11r*ro11i-v11i*ro11r
         end do
         if (nspden==2) then
           v22i=zero;ro22i=zero
           do ir=1,nrad
             jr=cplex_vxc*(ir-1)+1
             v22r=vxc1i(jr,2);if (cplex_vxc==2) v22i=vxc1i(jr+1,2)
             jr=cplex_den*(ir-1)+1
             ro22r=rhoarr(jr,1)-rhoarr(jr,2)
             if (cplex_den==2) ro22i=rhoarr(jr+1,1)-rhoarr(jr+1,2)
             ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
             gg(ir)=gg(ir)+v22r*ro22i-v22i*ro22r
           end do
         end if
       end if ! cplex_vxc and cplex_den

!      NON-COLLINEAR MAGNETISM
     else
       if (cplex_vxc==1.and.cplex_den==1) then   ! cplex_vxc==1 and cplex_den==1
         ff(:)=half*(vxc1i(:,1)*(rhoarr(:,1)+rhoarr(:,4)) &
&         +vxc1i(:,2)*(rhoarr(:,1)-rhoarr(:,4))) &
&         +vxc1i(:,3)*rhoarr(:,2) &
&         -vxc1i(:,4)*rhoarr(:,3)
         if (need_impart) gg(:)=zero
       else                                      ! other cases for cplex_vxc and cplex_den

!        V is stored as : v^11, v^22, V^12, i.V^21 (each are complex)
!        N is stored as : n, m_x, m_y, mZ          (each are complex)
         do ir=1,nrad
           jr=cplex_vxc*(ir-1)+1
           v11r= vxc1i(jr,1);v22r= vxc1i(jr,2)
           v12r= vxc1i(jr,3);v21i=-vxc1i(jr,1)
           if (cplex_vxc==2) then
             v11i= vxc1i(jr+1,1);v22i= vxc1i(jr+1,2)
             v12i= vxc1i(jr+1,3);v21r= vxc1i(jr+1,1)
           else
             v11i=zero;v22i=zero
             v12i=zero;v21i=zero
           end if
           jr=cplex_den*(ir-1)+1
           ro11r= rhoarr(jr,1)+rhoarr(jr,4)
           ro22r= rhoarr(jr,1)-rhoarr(jr,4)
           ro12r= rhoarr(jr,2);ro12i=-rhoarr(jr,3)
           ro21r= rhoarr(jr,2);ro21i= rhoarr(jr,3)
           if (cplex_den==2) then
             ro11i=rhoarr(jr+1,1)+rhoarr(jr+1,4)
             ro22i=rhoarr(jr+1,1)-rhoarr(jr+1,4)
             ro12r=ro12r+rhoarr(jr+1,3);ro12i=ro12i+rhoarr(jr+1,2)
             ro21r=ro21r-rhoarr(jr+1,3);ro21i=ro21i+rhoarr(jr+1,2)
           else
             ro11i=zero;ro22i=zero
           end if
!          Real part
           ff(ir)=half*(v11r*ro11r+v11i*ro11i+v22r*ro22r+v22i*ro22i &
&           +v12r*ro12r+v12i*ro12i+v21r*ro21r+v21i*ro21i)
!          Imaginary part
           if (need_impart) gg(ir)=half*(v11r*ro11i-v11i*ro11r+v22r*ro22i-v22i*ro22r &
&           +v12r*ro12i-v12i*ro12r+v21r*ro21i-v21i*ro21r)
         end do
       end if ! cplex_vxc and cplex_den
     end if ! nspden

     ff(:)=ff(:)*pawrad%rad(:)**2
     call simp_gen(vxcrho,ff,pawrad)
     d2enxc=d2enxc+vxcrho*pawang%angwgth(ipts)

     if (need_impart) then
       gg(:)=gg(:)*pawrad%rad(:)**2
       call simp_gen(vxcrho,gg,pawrad)
       d2enxc_im=d2enxc_im+vxcrho*pawang%angwgth(ipts)
     end if

     ABI_DEALLOCATE(ff)
     if (need_impart) ABI_DEALLOCATE(gg)

   end if

!  ----- End of the loop on npts (angular part)
 end do

!Add the four*pi factor of the angular integration
 if (option/=1) then
   d2enxc=d2enxc*four_pi
   if (need_impart) d2enxc_im=d2enxc_im*four_pi
 end if

 ABI_DEALLOCATE(rhoarr)
 ABI_DEALLOCATE(vxc1i)

 call timab(81,2,tsec)

 DBG_EXIT("COLL")

 end subroutine pawxc3
!!***
