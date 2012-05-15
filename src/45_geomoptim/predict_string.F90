!{\src2tex{textfont=tt}}
!!****f* ABINIT/predict_string
!! NAME
!! predict_string
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images.
!! Given the images, the changes on the geometry and others are predicted by rescaling the path
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
!! This is quite useful when ground states of the A and B states is known
!! mpi_enreg=MPI-parallelisation information
!! natom=dimension of vel_timimage and xred_timimage
!! ndynimage=number of dynamical images
!! nimage=number of images (on current proc)
!! nimage_tot=total number of images
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
!!      gather_array_img,mkrdim,spline,splint,xcast_mpi,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine predict_string(fxcartfactor,iatfix,itimimage,list_dynimage,mpi_enreg,natom,&
&                         ndynimage,nimage,nimage_tot,ntimimage,results_img)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_results_img, only : results_img_type,gather_array_img
 use m_splines
 use m_xmpi, only : xcast_mpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'predict_string'
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,natom,ndynimage,nimage,nimage_tot,ntimimage
 real(dp),intent(in) :: fxcartfactor
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in)     :: iatfix(3,natom),list_dynimage(ndynimage)
 type(results_img_type) :: results_img(nimage,ntimimage)

!Local variables-------------------------------
!scalars
 integer :: idynimage,ierr,ii,iimage,iatom
 real(dp) :: step
!arrays
 real(dp) :: acell(3),coordif(3),rprim(3,3),rprimd(3,3)
 real(dp),allocatable :: dequal(:),dimage(:)
 real(dp),allocatable :: x(:),y(:),z(:),x2(:),y2(:),z2(:)
 real(dp),allocatable :: xcart(:,:),xout(:),yout(:),zout(:)
 real(dp),allocatable,target :: xred(:,:,:)
 real(dp),pointer :: xred_all(:,:,:)

! *************************************************************************

 ABI_ALLOCATE(xred,(3,natom,nimage))

!Parallelism over image: only one process per image of the cell
 if (mpi_enreg%me_one_img==0) then

!  Compute new atomic positions in each cell
   do iimage=1,nimage
     xred(:,:,iimage)=results_img(iimage,itimimage)%xred(:,:)
   end do
   ABI_ALLOCATE(xcart,(3,natom))
   do idynimage=1,ndynimage
     iimage=list_dynimage(idynimage)
     acell(:)  =results_img(iimage,itimimage)%acell(:)
     rprim(:,:)=results_img(iimage,itimimage)%rprim(:,:)
     call mkrdim(acell,rprim,rprimd)
     call xredxcart(natom,1,rprimd,xcart,xred(:,:,iimage))
!    Note that one uses fcart, for which the sum of forces on all atoms vanish
!    (this is not the case for fred).
     xcart(:,:)=xcart(:,:)+fxcartfactor*results_img(iimage,itimimage)%results_gs%fcart(:,:)
     call xredxcart(natom,-1,rprimd,xcart,xred(:,:,iimage))
!    In case atom is fixed, we restore its previous value
     where(iatfix(:,:)==1)
     xred(:,:,iimage)=results_img(iimage,itimimage)%xred(:,:)
     end where
   end do
   ABI_DEALLOCATE(xcart)

!  Retrieve new atomic positions for all images
   if (mpi_enreg%paral_img==1) then
     ABI_ALLOCATE(xred_all,(3,natom,nimage_tot))
     call gather_array_img(xred,xred_all,mpi_enreg,allgather=.true.)
   else
     xred_all => xred
   end if

!  Reparametrize the string
!  Here the distance between images is calculated and normalized to a string length of 1.0
   ABI_ALLOCATE(dimage,(nimage_tot))
   dimage=zero
   do iimage=2,nimage_tot
     dimage(iimage)=dimage(iimage-1)
     do iatom=1,natom
       coordif=xred_all(:,iatom,iimage)-xred_all(:,iatom,iimage-1)
       dimage(iimage)=dimage(iimage)+sqrt(dot_product(coordif,coordif))
     end do
   end do
   dimage(:)=dimage(:)/dimage(nimage_tot)

!  dequal is just the parametrization on the string with equal arc-lengths
   ABI_ALLOCATE(dequal,(nimage_tot))
   dequal(1)=zero
   step=one/dble(nimage_tot-1)
   do iimage=2,nimage_tot
     dequal(iimage)=dequal(iimage-1)+step
   end do

!  New image coordinates are calculated and such that now the mesh is uniform
   ABI_ALLOCATE(x,(nimage_tot))
   ABI_ALLOCATE(y,(nimage_tot))
   ABI_ALLOCATE(z,(nimage_tot))
   ABI_ALLOCATE(x2,(nimage_tot))
   ABI_ALLOCATE(y2,(nimage_tot))
   ABI_ALLOCATE(z2,(nimage_tot))
   ABI_ALLOCATE(xout,(nimage_tot))
   ABI_ALLOCATE(yout,(nimage_tot))
   ABI_ALLOCATE(zout,(nimage_tot))
   do iatom=1,natom
     do iimage=1,nimage_tot
       x(iimage)=xred_all(1,iatom,iimage)
       y(iimage)=xred_all(2,iatom,iimage)
       z(iimage)=xred_all(3,iatom,iimage)
     end do
     call spline(dimage,x,nimage_tot,greatest_real,greatest_real,x2)
     call spline(dimage,y,nimage_tot,greatest_real,greatest_real,y2)
     call spline(dimage,z,nimage_tot,greatest_real,greatest_real,z2)
     call splint(nimage_tot,dimage,x,x2,nimage_tot,dequal,xout)
     call splint(nimage_tot,dimage,y,y2,nimage_tot,dequal,yout)
     call splint(nimage_tot,dimage,z,z2,nimage_tot,dequal,zout)
!    After a spline, new image coordinate for that particular
!    atom are generated only if they are dynamical
     do idynimage=1,ndynimage
       iimage=list_dynimage(idynimage)
       ii=mpi_enreg%index_img(iimage)
       xred(1,iatom,iimage)=xout(ii)
       xred(2,iatom,iimage)=yout(ii)
       xred(3,iatom,iimage)=zout(ii)
     end do
   end do  ! iatom

!  Free memory
   ABI_DEALLOCATE(x2)
   ABI_DEALLOCATE(y2)
   ABI_DEALLOCATE(z2)
   ABI_DEALLOCATE(x)
   ABI_DEALLOCATE(y)
   ABI_DEALLOCATE(z)
   ABI_DEALLOCATE(xout)
   ABI_DEALLOCATE(yout)
   ABI_DEALLOCATE(zout)
   ABI_DEALLOCATE(dimage)
   ABI_DEALLOCATE(dequal)
   if (mpi_enreg%paral_img==1)  then
     ABI_DEALLOCATE(xred_all)
   end if
   nullify(xred_all)

 end if ! mpi_enreg%me_one_img==0

!Store acell, rprim, xred and vel for the new iteration
 call xcast_mpi(xred,0,mpi_enreg%comm_one_img,ierr)
 do iimage=1,nimage
   results_img(iimage,itimimage+1)%xred(:,:) =xred(:,:,iimage)
   results_img(iimage,itimimage+1)%acell(:)  =results_img(iimage,itimimage)%acell(:)
   results_img(iimage,itimimage+1)%rprim(:,:)=results_img(iimage,itimimage)%rprim(:,:)
   results_img(iimage,itimimage+1)%vel(:,:)  =results_img(iimage,itimimage)%vel(:,:)
 end do
 ABI_DEALLOCATE(xred)

end subroutine predict_string
!!***
