!{\src2tex{textfont=tt}}
!!****f* ABINIT/predict_steepest
!! NAME
!! predict_steepest
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images.
!! Here, simple steepest descent algorithm, based on the value of the forces on the current timimage step.
!! No change of acell, rprim and vel at present.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! itimimage=number of the current time for image propagation (itimimage+1 is to be predicted here)
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example : in the NEB of string method, one expect the two end images to be fixed.
!! natom=dimension of vel_timimage and xred_timimage
!! ndynimage=number of dynamical images
!! nimage=number of images
!! ntimimage=dimension of several arrays
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! results_img(ntimimage,nimage)=datastructure that hold all the history of previous computations.
!!   results_img(:,:)%acell_timimage(3)
!!    at input, history of the values of acell for all images, up to itimimage
!!    at output, the predicted values of acell for all images, stored in acell_timimage(3,itimimage+1,nimage)
!!   results_img(:,:)%results_gs
!!    at input, history of the values of energies and forces for all images, up to itimimage
!!   results_img(:,:)%rprim_timimage(3,3)
!!    at input, history of the values of rprim for all images, up to itimimage
!!    at output, the predicted values of rprim for all images, stored in rprim_timimage(3,itimimage+1,nimage)
!!   results_img(:,:)%vel_timimage(3,natom)
!!    at input, history of the values of vel for all images, up to itimimage
!!    at output, the predicted values of vel for all images, stored in vel_timimage(3,natom,itimimage+1,nimage)
!!   results_img(:,:)%xred_timimage(3,natom)
!!    at input, history of the values of xred for all images, up to itimimage
!!    at output, the predicted values of xred for all images, stored in xred_timimage(3,natom,itimimage+1,nimage)
!!
!! PARENTS
!!      predictimg
!!
!! CHILDREN
!!      mkrdim,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine predict_steepest(fxcartfactor,itimimage,list_dynimage,natom,ndynimage,nimage,&
&                           ntimimage,results_img)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_results_img, only : results_img_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'predict_steepest'
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,natom,ndynimage,nimage,ntimimage
 real(dp),intent(in) :: fxcartfactor
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 type(results_img_type) :: results_img(nimage,ntimimage)

!Local variables-------------------------------
!scalars
 integer :: idynimage,iimage
!arrays
 real(dp) :: acell(3),rprim(3,3),rprimd(3,3)
 real(dp),allocatable :: xcart(:,:)

! *************************************************************************

 ABI_ALLOCATE(xcart,(3,natom))
 do idynimage=1,ndynimage

   iimage=list_dynimage(idynimage)

   acell(:)  =results_img(iimage,itimimage)%acell(:)
   rprim(:,:)=results_img(iimage,itimimage)%rprim(:,:)
   call mkrdim(acell,rprim,rprimd)
   call xredxcart(natom,1,rprimd,xcart,results_img(iimage,itimimage)%xred)

!  This is the core of the algorithm ... very simple ...
!  Note that one uses fcart, for which the sum of forces on all atoms vanish (this is not the case for fred).
   xcart(:,:)=xcart(:,:)+fxcartfactor*results_img(iimage,itimimage)%results_gs%fcart(:,:)

   call xredxcart(natom,-1,rprimd,xcart,results_img(iimage,itimimage+1)%xred)

   results_img(iimage,itimimage+1)%acell(:)  =results_img(iimage,itimimage)%acell(:)
   results_img(iimage,itimimage+1)%rprim(:,:)=results_img(iimage,itimimage)%rprim(:,:)
   results_img(iimage,itimimage+1)%vel(:,:)  =results_img(iimage,itimimage)%vel(:,:)

 end do  ! idynimage
 ABI_DEALLOCATE(xcart)

end subroutine predict_steepest
!!***
