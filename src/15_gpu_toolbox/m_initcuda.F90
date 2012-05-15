!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_initcuda
!! NAME
!! m_initcuda
!!
!! FUNCTION
!!  Module containing all variables concerning GPU device
!!  and the functions needed to extract them
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2012 ABINIT group (MMancini, MT, FDahm)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  Is an experimental development
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#if defined HAVE_GPU_CUDA
#include "cuda_common.h"
#endif

#include "abi_common.h"

module m_initcuda

 use m_profiling

 use defs_basis

 implicit none

#if defined HAVE_GPU_CUDA
 integer,parameter,public :: cudap=kind(CUDA_KIND)
#endif

!Structures
!!***

!!****t* m_initcuda/devGPU_type
!! NAME
!! devGPU_type
!!
!! FUNCTION
!! This structured datatype used to contains GPU properties
!!
!!
!! SOURCE
 type,public :: devGPU_type
  integer :: ndevice  !--number of available devices
  real,pointer  :: maxmemdev(:)  !--max global memory on any device
 end type devGPU_type
!!***

 private

 private ::            &
   prt_device_info       ! To print information about GPU

 public ::             &
   InitGPU,            & ! Initialise GPU
   Get_Mem_Dev,        & ! To obtain the max memory availeble on GPU device
   Get_ndevice,        & ! Number of devices of Capability>1.2
   CleanGPU,           & ! Clean devGPU_type variables
   setdevice_cuda,     & ! Set device, print info, ...
   unsetdevice_cuda      ! Unset device

CONTAINS !===========================================================
!!***


!!****f* m_initcuda/prt_device_info
!! NAME
!! prt_device_info
!!
!! FUNCTION
!! Print information about GPU device
!!
!! PARENTS
!!      m_initcuda
!!
!! CHILDREN
!!      unset_dev
!!
!! SOURCE

 subroutine prt_device_info(device)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prt_device_info'
 use interfaces_14_hidewrite
!End of the abilint section

  implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: device
!Local variables ------------------------------
!scalars
 integer :: gflop,constmem,sharemem
 integer :: ii,regist,lenname
 real 	  :: globalmem,clockRate
 character(20)  :: name
 character(20)  :: formatdev
 character(500) :: msg
!arrays
 integer :: vers(0:1)
! *********************************************************************
#if defined HAVE_GPU_CUDA
 write(msg,'(a,80a)')' ',('_',ii=1,80); call wrtout(std_out,msg,'PERS')
 write(msg,'(a25,a25,a31,a)')'________________________',' Graphic Card Properties ','_______________________________' ,ch10
 call wrtout(std_out,msg,'PERS')

 call get_dev_info(device,name,lenname,vers,globalmem,clockRate,gflop,constmem,sharemem,regist)

 write(formatdev,'(a12,i4,a)'),'(a23,i4,a3,a',lenname,')'
 write (msg,formatdev)&
       & '  Device             ',device,' : ',name(1:lenname)
 call wrtout(std_out,msg,'PERS')
 write (msg,'(a39,2(i1,a),2(a35,f7.1,2a),3(a35,i7,2a),a35,i7,a)')&
       & ' Revision number:                   ',vers(0),'.',vers(1),ch10, &
       & ' Total amount of global memory: ',globalmem,' Mbytes',ch10, &
       & ' Clock rate:                    ',clockRate,' GHz',ch10, &
       & ' Max GFLOP:                     ',gflop,' GFP',ch10, &
       & ' Total  constant memory:        ',constmem,' bytes',ch10, &
       & ' Shared memory per block:       ',sharemem,' bytes',ch10, &
       & ' Number of registers per block: ',regist,ch10
 call wrtout(std_out,msg,'PERS')
 if(device == -1)then
   write(msg,'(a)')' no cuda-GPU devices found'; call wrtout(std_out,msg,'PERS')
 end if
 write(msg,'(a,80a)')' ',('_',ii=1,80); call wrtout(std_out,msg,'PERS')
#endif
 end subroutine prt_device_info
!!***


!!****f* m_initcuda/InitGPU
!! NAME
!! InitGPU
!!
!! FUNCTION
!! Print information about GPU device
!!
!! PARENTS
!!      m_hidecudarec,m_initcuda
!!
!! CHILDREN
!!      unset_dev
!!
!! SOURCE

 subroutine InitGPU(gpuinfo,device)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'InitGPU'
#if defined HAVE_GPU_CUDA
#endif
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)              :: device
 type(devGPU_type),intent(inout) :: gpuinfo
!Local variables ------------------------------
!scalars
  real :: locmax
! *********************************************************************
 gpuinfo%ndevice = 0
 nullify(gpuinfo%maxmemdev)
#if defined HAVE_GPU_CUDA
!--Initialization
 if(device>-1)then
   !--Get the number of device for this proc
   gpuinfo%ndevice = 1
   ABI_ALLOCATE(gpuinfo%maxmemdev,(0:1))
   call get_GPU_max_mem(device,locmax)
   gpuinfo%maxmemdev(0:1) = locmax
   call  prt_device_info(device)
 endif
#endif
 end subroutine InitGPU
!!***


!****f* m_initcuda/Get_ndevice
!! NAME
!! Get_ndevice
!!
!! FUNCTION
!! Give the number of device with capability>=1.2
!!
!! PARENTS
!!      invars0,m_gpu_detect
!!
!! CHILDREN
!!      unset_dev
!!
!! SOURCE

 subroutine Get_ndevice(ndevice)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Get_ndevice'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ndevice
! *********************************************************************
#if defined HAVE_GPU_CUDA
!--Get the number of device for this proc
 call c_get_ndevice(ndevice)
#endif
 end subroutine Get_ndevice
!!***



!!****f* m_initcuda/Get_Mem_Dev
!! NAME
!! Get_Mem_Dev
!!
!! FUNCTION
!! Get the max memory availeble on device
!!
!! INPUTS
!! device  device number
!!
!! OUTPUT
!! max_mem_dev
!!
!! PARENTS
!!
!! CHILDREN
!!      unset_dev
!!
!! SOURCE

subroutine Get_Mem_Dev(device,max_mem_dev)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Get_Mem_Dev'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: device
 real,intent(out) :: max_mem_dev
!Local variables ------------------------------
! *********************************************************************
#if defined HAVE_GPU_CUDA
 call get_GPU_max_mem(device,max_mem_dev)
#endif
end subroutine Get_Mem_Dev
!!***


!!****f* m_initcuda/CleanGPU
!! NAME
!! CleanGPU
!!
!! FUNCTION
!! Print information about GPU device
!!
!! PARENTS
!!      m_hidecudarec,m_initcuda
!!
!! CHILDREN
!!      unset_dev
!!
!! SOURCE

 subroutine CleanGPU(gpuinfo)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'CleanGPU'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(devGPU_type),intent(inout) :: gpuinfo
! *********************************************************************
#if defined HAVE_GPU_CUDA
 if (associated(gpuinfo%maxmemdev))  then
   ABI_DEALLOCATE(gpuinfo%maxmemdev)
 end if
#endif

 end subroutine CleanGPU
!!***

!!****f* m_initcuda/setdevice_cuda
!! NAME
!! setdevice_cuda
!!
!! FUNCTION
!! Detect and activate a GPU device from current CPU core
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      unset_dev
!!
!! SOURCE

 subroutine setdevice_cuda(use_gpu_cuda)

 use m_xmpi, only: xmpi_world,xcomm_rank,xmpi_abort

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setdevice_cuda'
 use interfaces_14_hidewrite
#if defined HAVE_GPU_CUDA
#endif
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: use_gpu_cuda
!Local variables ------------------------------
!scalars
 integer :: device,nb_devices
 logical :: testopen
 character(len=500) :: msg
 type(devGPU_type) :: gpuinfo
! *********************************************************************

 if (use_gpu_cuda==0) return

#if defined HAVE_GPU_CUDA
 device=-1
 call c_get_ndevice(nb_devices)
 if(nb_devices>0) then
   device = MOD(xcomm_rank(xmpi_world),nb_devices)
   call set_dev(device)
   call check_context(nb_devices,msg)
   if(nb_devices==1) then !allocation succeed
     write(msg, '(4a,i1,2a)' ) ch10,&
&     ' setdevice_cuda : COMMENT -',ch10,&
&     '  GPU ',device,' has been properly initialized, continuing...',ch10
     call wrtout(std_out,msg,'PERS')
   else !gpu allocation failed we print error message returned and exit
     device=-1
     call wrtout(std_out,msg,'COLL')
     call xmpi_abort()
     inquire(std_out,OPENED=testopen)
     if (testopen) close(std_out)
#if defined FC_NAG
     call exit(-1)
#elif defined HAVE_FC_EXIT
     call exit(1)
#else
      stop 1
#endif
   end if
   call InitGPU(gpuinfo,device)
   call CleanGPU(gpuinfo)
 else
   use_gpu_cuda=0
 end if
#endif
 end subroutine setdevice_cuda
!!***


!!****f* m_initcuda/unsetdevice_cuda
!! NAME
!! setdevice_cuda
!!
!! FUNCTION
!! Deactivate a GPU device from current CPU core
!!
!! PARENTS
!!      m_hidecudarec
!!
!! CHILDREN
!!
!! SOURCE

 subroutine unsetdevice_cuda(use_gpu_cuda)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unsetdevice_cuda'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: use_gpu_cuda
!Local variables ------------------------------
!scalars
 character(len=500) :: msg
! *********************************************************************

 if (use_gpu_cuda==0) return

#if defined HAVE_GPU_CUDA
 call unset_dev()
#endif
 end subroutine unsetdevice_cuda
!!***

end module m_initcuda
!!***
