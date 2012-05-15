!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_results_img
!! NAME
!!  m_results_img
!!
!! FUNCTION
!!  This module provides the definition of the results_img_type used
!!  to store results from GS calculations for a given image of the cell.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_results_img

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_energies
 use m_results_gs
 use m_errors
 use m_xmpi

 implicit none

 private

! public procedures.
 public :: init_results_img
 public :: destroy_results_img
 public :: nullify_results_img
 public :: copy_results_img
 public :: gather_results_img
 public :: gather_array_img
 public :: scatter_array_img
!!***

!!****t* m_results_img/results_img_type
!! NAME
!! results_img_type
!!
!! FUNCTION
!! This structured datatype contains the results of a GS calculation
!! for a given image of the cell:
!!   energy, forces, stresses,positions, velocities, cell parameter

!!
!! SOURCE

 type, public :: results_img_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar

  integer :: natom
   ! The number of atoms for this image

! Real (real(dp)) arrays

  real(dp), pointer :: acell(:)
   ! acell(3)
   ! Dimensions of the cell

  real(dp), pointer :: rprim(:,:)
   ! rprim(3,3)
   ! Primitive translations of the cell

  real(dp), pointer :: vel(:,:)
   ! vel(3,natom)
   ! Velocities of the atoms

  real(dp), pointer :: xred(:,:)
   ! xred(3,natom)
   ! Reduced coordinates of the atoms

! Other types of data
  type(results_gs_type), pointer :: results_gs
   ! Energies, forces and stresses from the GS calculation

 end type results_img_type
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_results_img/init_results_img
!! NAME
!!  init_results_img
!!
!! FUNCTION
!!  Init all scalars and pointers in an array of results_img datastructures
!!
!! INPUTS
!!  natom=number of atoms in cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_img(:)=<type(results_img_type)>=results_img datastructure array
!!
!! PARENTS
!!      gstateimg,m_results_img
!!
!! CHILDREN
!!      xallgather_mpi,xcast_mpi,xscatterv_mpi
!!
!! SOURCE

subroutine init_results_img(natom,results_img)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_results_img'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 type(results_img_type),intent(inout) :: results_img(:)
!Local variables-------------------------------
!scalars
 integer :: ii,results_img_size
!arrays

!************************************************************************

 !@results_img_type

 results_img_size=size(results_img)

 if (results_img_size>0) then

   do ii=1,results_img_size

     results_img(ii)%natom  =natom

     ABI_ALLOCATE(results_img(ii)%results_gs,)
     call init_results_gs(natom,results_img(ii)%results_gs)

     ABI_ALLOCATE(results_img(ii)%acell,(3))
     results_img(ii)%acell=zero
     ABI_ALLOCATE(results_img(ii)%rprim,(3,3))
     results_img(ii)%rprim=zero
     ABI_ALLOCATE(results_img(ii)%xred,(3,natom))
     results_img(ii)%xred =zero
     ABI_ALLOCATE(results_img(ii)%vel,(3,natom))
     results_img(ii)%vel  =zero

   end do
 end if

end subroutine init_results_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/destroy_results_img
!! NAME
!!  destroy_results_img
!!
!! FUNCTION
!!  Clean and destroy an array of results_img datastructures
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_img(:)=<type(results_img_type)>=results_img datastructure array
!!
!! PARENTS
!!      gstateimg,prtimg
!!
!! CHILDREN
!!      xallgather_mpi,xcast_mpi,xscatterv_mpi
!!
!! SOURCE

subroutine destroy_results_img(results_img)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_results_img'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(results_img_type),intent(inout) :: results_img(:)
!Local variables-------------------------------
!scalars
 integer :: ii,results_img_size

!************************************************************************

 !@results_img_type

 results_img_size=size(results_img)
 if (results_img_size>0) then

   do ii=1,results_img_size

     results_img(ii)%natom=0

     if (associated(results_img(ii)%results_gs)) then
       call destroy_results_gs(results_img(ii)%results_gs)
       ABI_DEALLOCATE(results_img(ii)%results_gs)
     end if

     if (associated(results_img(ii)%acell))  then
       ABI_DEALLOCATE(results_img(ii)%acell)
     end if
     if (associated(results_img(ii)%rprim))  then
       ABI_DEALLOCATE(results_img(ii)%rprim)
     end if
     if (associated(results_img(ii)%xred))   then
       ABI_DEALLOCATE(results_img(ii)%xred)
     end if
     if (associated(results_img(ii)%vel))    then
       ABI_DEALLOCATE(results_img(ii)%vel)
     end if
   end do

 end if

end subroutine destroy_results_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/nullify_results_img
!! NAME
!!  nullify_results_img
!!
!! FUNCTION
!!  Nullify an array of results_img datastructures
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_img(:)=<type(results_img_type)>=results_img datastructure array
!!
!! PARENTS
!!
!! CHILDREN
!!      xallgather_mpi,xcast_mpi,xscatterv_mpi
!!
!! SOURCE

subroutine nullify_results_img(results_img)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_results_img'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(results_img_type),intent(inout) :: results_img(:)
!Local variables-------------------------------
!scalars
 integer :: ii,results_img_size

!************************************************************************

 !@results_img_type

 results_img_size=size(results_img)
 if (results_img_size>0) then

   do ii=1,results_img_size
     results_img(ii)%natom=0
     nullify(results_img(ii)%acell)
     nullify(results_img(ii)%rprim)
     nullify(results_img(ii)%xred)
     nullify(results_img(ii)%vel)
     nullify(results_img(ii)%results_gs)
   end do

 end if

end subroutine nullify_results_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/copy_results_img
!! NAME
!!  copy_results_img
!!
!! FUNCTION
!!  Copy a results_img datastructure into another
!!
!! INPUTS
!!  results_img_in=<type(results_img_type)>=input results_img datastructure
!!
!! OUTPUT
!!  results_img_out=<type(results_img_type)>=output results_img datastructure
!!
!! PARENTS
!!      m_results_img
!!
!! CHILDREN
!!      xallgather_mpi,xcast_mpi,xscatterv_mpi
!!
!! SOURCE

subroutine copy_results_img(results_img_in,results_img_out)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_results_img'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(results_img_type),intent(in) :: results_img_in
 type(results_img_type),intent(out) :: results_img_out
!Local variables-------------------------------
!scalars
 integer :: natom_in,natom_out

!************************************************************************

 !@results_img_type

 natom_in =results_img_in%natom
 natom_out=results_img_out%natom

 if (natom_in>natom_out) then
   if (associated(results_img_out%acell))  then
     ABI_DEALLOCATE(results_img_out%acell)
   end if
   if (associated(results_img_out%rprim))  then
     ABI_DEALLOCATE(results_img_out%rprim)
   end if
   if (associated(results_img_out%xred))   then
     ABI_DEALLOCATE(results_img_out%xred)
   end if
   if (associated(results_img_out%vel))    then
     ABI_DEALLOCATE(results_img_out%vel)
   end if

   if (associated(results_img_in%acell))  then
     ABI_ALLOCATE(results_img_out%acell,(3))
   end if
   if (associated(results_img_in%rprim))  then
     ABI_ALLOCATE(results_img_out%rprim,(3,3))
   end if
   if (associated(results_img_in%xred))   then
     ABI_ALLOCATE(results_img_out%xred,(3,natom_in))
   end if
   if (associated(results_img_in%vel))    then
     ABI_ALLOCATE(results_img_out%vel,(3,natom_in))
   end if
 end if

 results_img_out%natom  =results_img_in%natom

 call copy_results_gs(results_img_in%results_gs,results_img_out%results_gs)

 if (associated(results_img_in%acell)) results_img_out%acell(:)=results_img_in%acell(:)
 if (associated(results_img_in%rprim)) results_img_out%rprim(:,:)=results_img_in%rprim(:,:)
 if (associated(results_img_in%xred))  results_img_out%xred(:,1:natom_in)=results_img_in%xred(:,1:natom_in)
 if (associated(results_img_in%vel))   results_img_out%vel(:,1:natom_in)=results_img_in%vel(:,1:natom_in)

end subroutine copy_results_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/gather_results_img
!! NAME
!!  gather_results_img
!!
!! FUNCTION
!!  Gather results_img datastructures using communicator over images (replicas) of the cell.
!!  Each contribution of single processor is gathered into a big array on master processor
!!
!! INPUTS
!!  allgather= --optional, default=false--  if TRUE do ALL_GATHER instead of GATHER
!!  master= --optional, default=0-- index of master proc receiving gathered data (if allgather=false)
!!  mpi_enreg=informations about MPI parallelization
!!  only_one_per_img= --optional, default=true--  if TRUE, the gather operation
!!                    is only done by one proc per image (master of the comm_one_img)
!!  results_img(:)=<type(results_img_type)>=results_img datastructure array on each proc
!!
!! SIDE EFFECTS
!!  results_img_all(:)=<type(results_img_type)>=global (gathered) results_img datastructure array
!!
!! PARENTS
!!      prtimg
!!
!! CHILDREN
!!      xallgather_mpi,xcast_mpi,xscatterv_mpi
!!
!! SOURCE

subroutine gather_results_img(mpi_enreg,results_img,results_img_all,&
&                 master,allgather,only_one_per_img) ! optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gather_results_img'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: master
 logical,optional,intent(in) :: allgather,only_one_per_img
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 type(results_img_type),intent(inout) :: results_img(:)
 type(results_img_type),intent(inout) :: results_img_all(:)

!Local variables-------------------------------
!scalars
 integer :: ibufr,ierr,iproc,jj
 integer :: master_all,master_img,master_one_img,natom,nimage,nimagetot
 integer :: rsize,rsize_img
 logical :: do_allgather,i_am_master,one_per_img,use_results_all
 character(len=500) :: msg
!arrays
 integer,allocatable :: iimg(:),nimage_all(:),rbufshft(:),rsize_img_all(:)
 real(dp),allocatable :: rbuffer(:),rbuffer_all(:)

!************************************************************************

 !@results_img_type

 one_per_img=.true.;if (present(only_one_per_img)) one_per_img=only_one_per_img
 do_allgather=.false.;if (present(allgather)) do_allgather=allgather
 master_all=0;if (present(master)) master_all=master
 i_am_master=(mpi_enreg%me==master_all)

 master_img=0;master_one_img=0
 use_results_all= &
&  (((     do_allgather).and.(     one_per_img).and.(mpi_enreg%me_one_img==master_one_img)) .or. &
&   ((     do_allgather).and.(.not.one_per_img))                                            .or. &
&   ((.not.do_allgather).and.(     one_per_img).and.(mpi_enreg%me==master_all))             .or. &
&   ((.not.do_allgather).and.(.not.one_per_img).and.(mpi_enreg%me_img==master_img)))

!Create global results_img_all datastructure
 if (use_results_all) then
   call init_results_img(results_img(1)%natom,results_img_all)
 end if

 if ((.not.one_per_img).or.(mpi_enreg%me_one_img==master_one_img)) then

!  Simple copy in case of 1 image
   if (use_results_all) then
     if (size(results_img_all,1)<=1) then
       call copy_results_img(results_img(1),results_img_all(1))
       return
     end if
   endif

!  Gather number of images treated by each proc
   nimage=size(results_img,1)
   ABI_ALLOCATE(nimage_all,(mpi_enreg%nproc_img))
   call xallgather_mpi(nimage,nimage_all,mpi_enreg%comm_img,ierr)
   nimagetot=sum(nimage_all)
   if (use_results_all) then
     if (size(results_img_all,1)/=nimagetot) then
       msg='  Wrong results_img_all size ! '
       MSG_BUG(msg)
     end if
   end if

!  Copy natom from distributed results_img to gathered one
   natom=results_img(1)%natom
   if (use_results_all) then
     do jj=1,nimagetot
       results_img_all(jj)%natom=natom
     enddo
   end if

!  Compute number of data
   rsize=29+n_energies+24*natom
   rsize_img=nimage*rsize
   ABI_ALLOCATE(rsize_img_all,(mpi_enreg%nproc_img))
   rsize_img_all(:)=rsize*nimage_all(:)
   ABI_DEALLOCATE(nimage_all)

!  Compute shifts in buffer arrays for each proc
   ABI_ALLOCATE(rbufshft,(mpi_enreg%nproc_img))
   rbufshft(1)=0
   do jj=2,mpi_enreg%nproc_img
     rbufshft(jj)=rbufshft(jj-1)+rsize_img_all(jj-1)
   end do

!  Load buffers
   ABI_ALLOCATE(rbuffer,(rsize_img))
   ibufr=0
   do jj=1,nimage
     rbuffer(ibufr+1)  =results_img(jj)%results_gs%deltae
     rbuffer(ibufr+2)  =results_img(jj)%results_gs%diffor
     rbuffer(ibufr+3)  =results_img(jj)%results_gs%entropy
     rbuffer(ibufr+4)  =results_img(jj)%results_gs%etotal
     rbuffer(ibufr+5)  =results_img(jj)%results_gs%fermie
     rbuffer(ibufr+6)  =results_img(jj)%results_gs%residm
     rbuffer(ibufr+7)  =results_img(jj)%results_gs%res2
     rbuffer(ibufr+8)  =results_img(jj)%results_gs%vxcavg
     rbuffer(ibufr+9:ibufr+11) =results_img(jj)%results_gs%pel(1:3)
     rbuffer(ibufr+12:ibufr+17)=results_img(jj)%results_gs%strten(1:6)
     rbuffer(ibufr+18:ibufr+20)=results_img(jj)%acell(1:3)
     rbuffer(ibufr+21:ibufr+23)=results_img(jj)%rprim(1:3,1)
     rbuffer(ibufr+24:ibufr+26)=results_img(jj)%rprim(1:3,2)
     rbuffer(ibufr+27:ibufr+29)=results_img(jj)%rprim(1:3,3)
     ibufr=ibufr+29
     call energies_to_array(results_img(jj)%results_gs%energies,&
&                           rbuffer(ibufr+1:ibufr+n_energies),1)
     ibufr=ibufr+n_energies
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%fcart(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%fred(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%gresid(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%grewtn(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%grxc(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%results_gs%synlgr(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%xred(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
     rbuffer(ibufr+1:ibufr+3*natom)=reshape(results_img(jj)%vel(1:3,1:natom),(/3*natom/))
     ibufr=ibufr+3*natom
   end do
   if (ibufr/=rsize_img) then
     msg='  wrong buffer size ! '
     MSG_BUG(msg)
   end if

!  Gather all data
   if (use_results_all)  then
     ABI_ALLOCATE(rbuffer_all,(rsize*nimagetot))
   end if
   if (.not.use_results_all)  then
     ABI_ALLOCATE(rbuffer_all,(0))
   end if
   if (do_allgather) then
     call xallgatherv_mpi(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                         mpi_enreg%comm_img,ierr)
   else
     call xgatherv_mpi(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                          master_img,mpi_enreg%comm_img,ierr)
   end if
   ABI_DEALLOCATE(rbuffer)
   ABI_DEALLOCATE(rsize_img_all)

!  Transfer buffers into gathered results_img_all (master proc only)
   if (use_results_all) then
     ABI_ALLOCATE(iimg,(mpi_enreg%nproc_img))
     iimg=0
     do jj=1,nimagetot
!      The following line supposes that images are sorted by increasing index
       iproc=mpi_enreg%distrb_img(jj)+1;iimg(iproc)=iimg(iproc)+1
       ibufr=rbufshft(iproc)+(iimg(iproc)-1)*rsize
       results_img_all(jj)%results_gs%deltae     =rbuffer_all(ibufr+1)
       results_img_all(jj)%results_gs%diffor     =rbuffer_all(ibufr+2)
       results_img_all(jj)%results_gs%entropy    =rbuffer_all(ibufr+3)
       results_img_all(jj)%results_gs%etotal     =rbuffer_all(ibufr+4)
       results_img_all(jj)%results_gs%fermie     =rbuffer_all(ibufr+5)
       results_img_all(jj)%results_gs%residm     =rbuffer_all(ibufr+6)
       results_img_all(jj)%results_gs%res2       =rbuffer_all(ibufr+7)
       results_img_all(jj)%results_gs%vxcavg     =rbuffer_all(ibufr+8)
       results_img_all(jj)%results_gs%pel(1:3)   =rbuffer_all(ibufr+9:ibufr+11)
       results_img_all(jj)%results_gs%strten(1:6)=rbuffer_all(ibufr+12:ibufr+17)
       results_img_all(jj)%acell(1:3)   =rbuffer_all(ibufr+18:ibufr+20)
       results_img_all(jj)%rprim(1:3,1)=rbuffer_all(ibufr+21:ibufr+23)
       results_img_all(jj)%rprim(1:3,2)=rbuffer_all(ibufr+24:ibufr+26)
       results_img_all(jj)%rprim(1:3,3)=rbuffer_all(ibufr+27:ibufr+29)
       ibufr=ibufr+29
       call energies_to_array(results_img_all(jj)%results_gs%energies,&
&                             rbuffer_all(ibufr+1:ibufr+n_energies),-1)
       ibufr=ibufr+n_energies
       results_img_all(jj)%results_gs%fcart(1:3,1:natom)= &
&             reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%fred(1:3,1:natom)= &
&             reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%gresid(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%grewtn(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%grxc(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%results_gs%synlgr(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%xred(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
       results_img_all(jj)%vel(1:3,1:natom)= &
  &            reshape(rbuffer_all(ibufr+1:ibufr+3*natom),(/3,natom/))
       ibufr=ibufr+3*natom
     end do
     ABI_DEALLOCATE(iimg)
   end if

!  Free memory
   ABI_DEALLOCATE(rbufshft)
   ABI_DEALLOCATE(rbuffer_all)

 end if

end subroutine gather_results_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/gather_array_img
!! NAME
!!  gather_array_img
!!
!! FUNCTION
!!  Gather an real 2D-array (part of a results_img datastructure) using communicator
!!  over images (replicas) of the cell.
!!  Each contribution of single processor is gathered into a big array on master processor
!!
!! INPUTS
!!  allgather= --optional, default=false--  if TRUE do ALL_GATHER instead of GATHER
!!  master= --optional, default=0-- index of master proc receiving gathered data (if allgather=false)
!!  mpi_enreg=informations about MPI parallelization
!!  only_one_per_img= --optional, default=true--  if TRUE, the gather operation
!!                    is only done by one proc per image (master of the comm_one_img)
!!  array_img(:,:,:)= (real) 2D-array distributed (has 3 dimensions; the 3rd one is nimage)
!!
!! SIDE EFFECTS
!!  array_img_all(:,:,:)= (real) global (gathered) 2D-array
!!                        (has 3 dimensions; the 3rd one is nimagetot)
!!
!! PARENTS
!!      predict_pimd,predict_string
!!
!! CHILDREN
!!      xallgather_mpi,xcast_mpi,xscatterv_mpi
!!
!! SOURCE

subroutine gather_array_img(array_img,array_img_all,mpi_enreg,&
&                           master,allgather,only_one_per_img) ! optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gather_array_img'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: master
 logical,optional,intent(in) :: allgather,only_one_per_img
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: array_img(:,:,:)
 real(dp),intent(inout) :: array_img_all(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: ibufr,ierr,iproc,jj
 integer :: master_all,master_img,master_one_img,nimage,nimagetot
 integer :: rsize,rsize_img,size1,size2
 logical :: do_allgather,i_am_master,one_per_img,use_array_all
 character(len=500) :: msg
!arrays
 integer,allocatable :: iimg(:),nimage_all(:),rbufshft(:),rsize_img_all(:)
 real(dp),allocatable :: rbuffer(:),rbuffer_all(:)

!************************************************************************

 !@results_img_type

 one_per_img=.true.;if (present(only_one_per_img)) one_per_img=only_one_per_img
 do_allgather=.false.;if (present(allgather)) do_allgather=allgather
 master_all=0;if (present(master)) master_all=master
 i_am_master=(mpi_enreg%me==master_all)

 master_img=0;master_one_img=0
 use_array_all= &
&  (((     do_allgather).and.(     one_per_img).and.(mpi_enreg%me_one_img==master_one_img)) .or. &
&   ((     do_allgather).and.(.not.one_per_img))                                            .or. &
&   ((.not.do_allgather).and.(     one_per_img).and.(mpi_enreg%me==master_all))             .or. &
&   ((.not.do_allgather).and.(.not.one_per_img).and.(mpi_enreg%me_img==master_img)))

 size1=size(array_img,1);size2=size(array_img,2)
 if (use_array_all) then
   if (size(array_img_all,1)/=size1.or.size(array_img_all,2)/=size2) then
     msg='  Wrong array_img_all size (1) ! '
     MSG_BUG(msg)
   end if
 end if

 if ((.not.one_per_img).or.(mpi_enreg%me_one_img==master_one_img)) then

!  Simple copy in case of 1 image
   if (use_array_all) then
     if (size(array_img_all,3)<=1) then
       array_img_all(:,:,1)=array_img(:,:,1)
       return
     end if
   endif

!  Gather number of images treated by each proc
   nimage=size(array_img,3)
   ABI_ALLOCATE(nimage_all,(mpi_enreg%nproc_img))
   call xallgather_mpi(nimage,nimage_all,mpi_enreg%comm_img,ierr)
   nimagetot=sum(nimage_all)
   if (use_array_all) then
     if (size(array_img_all,3)/=nimagetot) then
       msg='  Wrong array_img_all size (2) ! '
       MSG_BUG(msg)
     endif
   end if

!  Compute number of data
   rsize=size1*size2;rsize_img=nimage*rsize
   ABI_ALLOCATE(rsize_img_all,(mpi_enreg%nproc_img))
   rsize_img_all(:)=rsize*nimage_all(:)
   ABI_DEALLOCATE(nimage_all)

!  Compute shifts in buffer arrays for each proc
   ABI_ALLOCATE(rbufshft,(mpi_enreg%nproc_img))
   rbufshft(1)=0
   do jj=2,mpi_enreg%nproc_img
     rbufshft(jj)=rbufshft(jj-1)+rsize_img_all(jj-1)
   end do

!  Load buffers
   ABI_ALLOCATE(rbuffer,(rsize_img))
   ibufr=0
   do jj=1,nimage
     rbuffer(ibufr+1:ibufr+rsize)=reshape(array_img(1:size1,1:size2,jj),(/rsize/))
     ibufr=ibufr+rsize
   end do

!  Gather all data
   if (use_array_all)  then
     ABI_ALLOCATE(rbuffer_all,(rsize*nimagetot))
   end if
   if (.not.use_array_all)  then
     ABI_ALLOCATE(rbuffer_all,(0))
   end if
   if (do_allgather) then
     call xallgatherv_mpi(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                         mpi_enreg%comm_img,ierr)
   else
     call xgatherv_mpi(rbuffer,rsize_img,rbuffer_all,rsize_img_all,rbufshft,&
&                      master_img,mpi_enreg%comm_img,ierr)
   end if
   ABI_DEALLOCATE(rbuffer)
   ABI_DEALLOCATE(rsize_img_all)

!  Transfer buffers into gathered array_img_all (master proc only)
   if (use_array_all) then
     ABI_ALLOCATE(iimg,(mpi_enreg%nproc_img))
     iimg=0
     do jj=1,nimagetot
!      The following line supposes that images are sorted by increasing index
       iproc=mpi_enreg%distrb_img(jj)+1;iimg(iproc)=iimg(iproc)+1
       ibufr=rbufshft(iproc)+(iimg(iproc)-1)*rsize
       array_img_all(1:size1,1:size2,jj)=reshape(rbuffer_all(ibufr+1:ibufr+rsize),(/size1,size2/))
     end do
     ABI_DEALLOCATE(iimg)
   end if

!  Free memory
   ABI_DEALLOCATE(rbufshft)
   ABI_DEALLOCATE(rbuffer_all)

 end if

end subroutine gather_array_img
!!***

!----------------------------------------------------------------------

!!****f* m_results_img/scatter_array_img
!! NAME
!!  scatter_array_img
!!
!! FUNCTION
!!  Scatter an real 2D-array (part of a results_img datastructure) using communicator
!!  over images (replicas) of the cell.
!!  A big array on master processor is scattered into each contribution
!!  of single processor is gathered
!!
!! INPUTS
!!  master= --optional, default=0-- index of master proc sending data
!!  mpi_enreg=informations about MPI parallelization
!!  only_one_per_img= --optional, default=true--  if TRUE, the scatter operation
!!                    is only done by one proc per image (master of the comm_one_img)
!!  array_img_all(:,:,:)= (real) global 2D-array (has 3 dimensions; the 3rd one is nimagetot)
!!
!! SIDE EFFECTS
!!  array_img(:,:,:)= (real) distributed 2D-array
!!                    (has 3 dimensions; the 3rd one is nimage)
!!
!! PARENTS
!!      predict_pimd
!!
!! CHILDREN
!!      xallgather_mpi,xcast_mpi,xscatterv_mpi
!!
!! SOURCE

subroutine scatter_array_img(array_img,array_img_all,mpi_enreg,&
&                            master,only_one_per_img) ! optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scatter_array_img'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: master
 logical,optional,intent(in) :: only_one_per_img
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: array_img(:,:,:)
 real(dp),intent(in) :: array_img_all(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: ibufr,ierr,iproc,jj



 integer :: master_all,master_img,master_one_img,nimage,nimagetot
 integer :: rsize,rsize_img,size1,size2
 logical :: i_am_master,one_per_img,use_array_all
 character(len=500) :: msg
!arrays
 integer,allocatable :: iimg(:),nimage_all(:),rbufshft(:),rsize_img_all(:)
 real(dp),allocatable :: rbuffer(:),rbuffer_all(:)

!************************************************************************

 !@results_img_type

 one_per_img=.true.;if (present(only_one_per_img)) one_per_img=only_one_per_img
 master_all=0;if (present(master)) master_all=master
 i_am_master=(mpi_enreg%me==master_all)

 use_array_all=i_am_master
 master_img=0;master_one_img=0

 size1=size(array_img,1);size2=size(array_img,2)
 if (use_array_all) then
   if (size(array_img_all,1)/=size1.or.size(array_img_all,2)/=size2) then
     msg='  Wrong array_img_all size (1) ! '
     MSG_BUG(msg)
   end if
 end if

 if ((.not.one_per_img).or.(mpi_enreg%me_one_img==master_one_img)) then

!  Compute (by gather operation) total number of images
   nimage=size(array_img,3)
   ABI_ALLOCATE(nimage_all,(mpi_enreg%nproc_img))
   call xallgather_mpi(nimage,nimage_all,mpi_enreg%comm_img,ierr)
   nimagetot=sum(nimage_all)
   if (use_array_all) then
     if (size(array_img_all,3)/=nimagetot) then
       msg='  Wrong array_img_all size (2) ! '
       MSG_BUG(msg)
     endif
   end if

!  Simple copy in case of 1 image
   if (nimagetot<=1) then
     if (use_array_all) array_img(:,:,1)=array_img_all(:,:,1)

   else

!    Compute number of data
     rsize=size1*size2;rsize_img=nimage*rsize
     ABI_ALLOCATE(rsize_img_all,(mpi_enreg%nproc_img))
     rsize_img_all(:)=rsize*nimage_all(:)

!    Compute shifts in buffer arrays for each proc
     ABI_ALLOCATE(rbufshft,(mpi_enreg%nproc_img))
     rbufshft(1)=0
     do jj=2,mpi_enreg%nproc_img
       rbufshft(jj)=rbufshft(jj-1)+rsize_img_all(jj-1)
     end do

!    Load buffer
     if (use_array_all)  then
       ABI_ALLOCATE(rbuffer_all,(rsize*nimagetot))
     end if
     if (.not.use_array_all)  then
       ABI_ALLOCATE(rbuffer_all,(0))
     end if
     if (use_array_all) then
       ABI_ALLOCATE(iimg,(mpi_enreg%nproc_img))
       iimg=0
       do jj=1,nimagetot
!        The following line supposes that images are sorted by increasing index
         iproc=mpi_enreg%distrb_img(jj)+1;iimg(iproc)=iimg(iproc)+1
         ibufr=rbufshft(iproc)+(iimg(iproc)-1)*rsize
         rbuffer_all(ibufr+1:ibufr+rsize)=reshape(array_img_all(1:size1,1:size2,jj),(/rsize/))
       end do
       ABI_DEALLOCATE(iimg)
      end if

!    Scatter all data
     ABI_ALLOCATE(rbuffer,(rsize_img))
     call xscatterv_mpi(rbuffer_all,rsize_img_all,rbufshft,rbuffer,rsize_img,&
&                       master_img,mpi_enreg%comm_img,ierr)
     ABI_DEALLOCATE(rbuffer_all)
     ABI_DEALLOCATE(rbufshft)
     ABI_DEALLOCATE(rsize_img_all)

!    Transfered distributed buffers into array_img (master proc only)
     ibufr=0
     do jj=1,nimage
       array_img(1:size1,1:size2,jj)=reshape(rbuffer(ibufr+1:ibufr+rsize),(/size1,size2/))
       ibufr=ibufr+rsize
     end do
     ABI_DEALLOCATE(rbuffer)

   end if ! nimagetot<=1
   ABI_DEALLOCATE(nimage_all)

!  Now, is requested dispatch data inside each image
   if (.not.one_per_img) then
     call xcast_mpi(array_img,0,mpi_enreg%comm_one_img,ierr)
   end if

 end if

end subroutine scatter_array_img
!!***

END MODULE m_results_img
!!***
