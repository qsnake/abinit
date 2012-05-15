!{\src2tex{textfont=tt}}
!!****f* ABINIT/predictimg
!! NAME
!! predictimg
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! deltae=averaged energy difference used to control convergence over images
!! fxcartfactor=factor that allows to transform forces in cartesian displacements for selected algorithms.
!! iatfix(3,natom)= fix xred coordinates if 1
!! imagealgo_str=name of the algorithm (with images) used
!! imgmov=gives the algorithm to be used for prediction of new set of images
!! itimimage=number of the current time for image propagation (itimimage+1 is to be predicted here)
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example : in the NEB method, or in the string method, one expect the two end images to be fixed.
!! mpi_enreg=MPI-parallelisation information
!! natom=dimension of vel_timimage and xred_timimage
!! ndynimage=number of dynamical images
!! nimage=number of images (treated by current proc)
!! nimage_tot=total number of images
!! ntimimage=dimension of several arrays
!! pimd_param=several parameters for Path-Integral MD
!! prtvolimg=printing volume
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
!!      gstateimg
!!
!! CHILDREN
!!      predict_copy,predict_ga,predict_pimd,predict_steepest,predict_string
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine predictimg(deltae,fxcartfactor,iatfix,imagealgo_str,imgmov,itimimage,list_dynimage,&
&                     mpi_enreg,natom,ndynimage,nimage,nimage_tot,ntimimage,&
&                     pimd_param,prtvolimg,results_img)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_pimd
 use m_results_img
 use m_results_gs , only : results_gs_type
 use m_use_ga

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'predictimg'
 use interfaces_14_hidewrite
 use interfaces_45_geomoptim, except_this_one => predictimg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: imgmov,itimimage,natom,ndynimage,nimage,nimage_tot,ntimimage,prtvolimg
 character(len=60),intent(in) :: imagealgo_str
 real(dp),intent(in) :: deltae,fxcartfactor
 type(pimd_type),intent(in) :: pimd_param
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 integer,intent(in) :: iatfix(3,natom)
 type(results_img_type) :: results_img(nimage,ntimimage)

!Local variables-------------------------------
!scalars
 integer, SAVE :: idum=5
 character(len=500) :: msg
!arrays

! *************************************************************************

 write(msg,'(3a)') ch10,&
& '------------------------------------------------------------',ch10
 if (prtvolimg<2) write(msg,'(5a)') trim(msg),' ',trim(imagealgo_str),':',ch10
 if (itimimage>1) write(msg,'(2a,es11.3,2a)') trim(msg),&
& ' Average[Abs(Etotal(t)-Etotal(t-dt))]=',deltae,' Hartree',ch10
 write(msg,'(2a)') trim(msg),' Moving images of the cell...'
 call wrtout(ab_out ,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 select case(imgmov)

   case(0)

     call predict_copy(itimimage,list_dynimage,ndynimage,nimage,ntimimage,results_img)

   case(1)

     call predict_steepest(fxcartfactor,itimimage,list_dynimage,natom,ndynimage,nimage,&
&     ntimimage,results_img)

   case(2)

     call predict_string(fxcartfactor,iatfix,itimimage,list_dynimage,mpi_enreg,natom,&
&     ndynimage,nimage,nimage_tot,ntimimage,results_img)

   case(4)

     call predict_ga(itimimage,idum,natom,nimage,ntimimage,results_img)

   case(9, 13)
!    Path Integral Molecular Dynamics
     call predict_pimd(imgmov,itimimage,mpi_enreg,natom,nimage,nimage_tot,&
&     ntimimage,pimd_param,prtvolimg,results_img)

     case default

 end select

end subroutine predictimg
!!***
