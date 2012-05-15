!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxcm3
!! NAME
!! pawxcm3
!!
!! FUNCTION
!! PAW only
!! PAW only
!! Compute first-order change of XC potential and contribution to
!! 2nd-order change of XC energy inside a PAW sphere.
!! LDA+GGA - USE A DEVELOPMENT OF THE DENSITY OVER (L,M) MOMENTS
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
!!  kxc(pawrad%mesh_size,lm_size,nkxc)=GS xc kernel
!!  lm_size=size of density array rhor (see below)
!!  lmselect(lm_size)=select the non-zero LM-moments of input density rhor1
!!  nhat1(cplex_den*pawrad%mesh_size,lm_size,nspden)=first-order change of compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  option=0  compute both 2nd-order XC energy and 1st-order potential
!!         1  compute only 1st-order XC potential
!!         2  compute only 2nd-order XC energy, XC potential is temporary computed here
!!         3  compute only 2nd-order XC energy, XC potential is input in vxc1(:)!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawxcdev=order of Vxc development
!!  rhor1(cplex_den*pawrad%mesh_size,lm_size,nspden)=first-order change of density
!!  usecore= 1 if core density has to be used in Exc/Vxc ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in d2Exc only
!!             2 if compensation density (nhat) has to be used in d2Exc and Vxc1
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  == if option=0 or 2 or 3 ==
!!    d2enxc=returned exchange-cor. contribution to 2nd-order XC energy
!!
!! SIDE EFFECTS
!!    vxc1(cplex_vxc*pawrad%mesh_size,pawang%angl_size,nspden)=1st-order XC potential
!!      Output if option==0 or 1
!!      Unused if option==2
!!      Input  if option==3
!!
!! PARENTS
!!      pawdenpot,pawnstd2e
!!
!! CHILDREN
!!      pawxc3,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawxcm3(corexc1,cplex_den,cplex_vxc,d2enxc,kxc,lm_size,lmselect,nhat1,nkxc,nspden,&
&                   option,pawang,pawrad,pawxcdev,rhor1,usecore,usexcnhat,vxc1,xclevel,&
&                   d2enxc_im) ! optional

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawxcm3'
 use interfaces_18_timing
 use interfaces_56_xc, except_this_one => pawxcm3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex_den,cplex_vxc,lm_size,nkxc,nspden,option
 integer,intent(in) :: pawxcdev,usecore,usexcnhat,xclevel
 real(dp),intent(out) :: d2enxc
 real(dp),intent(out),optional :: d2enxc_im
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size)
 real(dp),intent(in) :: corexc1(cplex_den*pawrad%mesh_size)
 real(dp),intent(in) :: kxc(pawrad%mesh_size,lm_size,nkxc)
 real(dp),intent(in) :: nhat1(cplex_den*pawrad%mesh_size,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor1(cplex_den*pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(inout) :: vxc1(cplex_vxc*pawrad%mesh_size,lm_size,nspden)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)

!************************************************************************

 DBG_ENTER("COLL")

 call timab(81,1,tsec)

 if(pawxcdev/=0)then
   msg='  Not yet implemented (choose pawxcdev=0) !'
   MSG_ERROR(msg)
 end if

!NOTE (MT)
!lmselect and lm_size are not necessarily the samer for densities, kxc and vxc1
!This is not taken into account for the moment, but has to be programmed...

!NOTE   (LD)
!TEMPORARY solution : To call pawxc3 (LD after consulting MT)
!pawxc3 and pawxcm3 play the same role, but using two different algorithms.
!the first one computes the first-order xc potential using densities developped on (r, theta,phi);
!the second one computes this same first-order xc potential using densities developped on spherical harmonics. They do the same and the choice is imposed by the pawxcdev keyword in the input file.

 if (present(d2enxc_im)) then
   call pawxc3(corexc1,cplex_den,cplex_vxc,d2enxc,kxc,lm_size,lmselect,nhat1,nkxc,nspden,&
&   option,pawang,pawrad,rhor1,usecore,usexcnhat,vxc1,xclevel,d2enxc_im=d2enxc_im)
 else
   call pawxc3(corexc1,cplex_den,cplex_vxc,d2enxc,kxc,lm_size,lmselect,nhat1,nkxc,nspden,&
&   option,pawang,pawrad,rhor1,usecore,usexcnhat,vxc1,xclevel)
 end if

!----- End of routine
 call timab(81,2,tsec)

 DBG_EXIT("COLL")

 end subroutine pawxcm3
!!***
