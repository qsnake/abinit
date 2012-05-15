!{\src2tex{textfont=tt}}
!!****f* ABINIT/prttagm_images
!!
!! NAME
!! prttagm_images
!!
!! FUNCTION
!! Extension to prttagm to include the printing of
!! images information, in those cases the same variable
!! is printed several times for each dataset 
!!
!! Cases where images information are relevant includes
!! xcart, xred, acell, fcart.
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      outvar1,outvars
!!
!! CHILDREN
!!      appdig,write_var_netcdf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prttagm_images(dprarr_images,iout,jdtset_,&
& marr,narrm,ncid,ndtset_alloc,token,&
& mxnimage,nimage,ndtset,prtimg,strimg)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prttagm_images'
 use interfaces_32_util
 use interfaces_57_iovars, except_this_one => prttagm_images
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,marr,ndtset_alloc,ncid
 integer,intent(in) :: mxnimage,ndtset
 character(len=*),intent(in) :: token
!arrays
 integer,intent(in) :: prtimg(mxnimage,0:ndtset_alloc)
 integer,intent(in) :: jdtset_(0:ndtset_alloc)
 integer,intent(in) :: nimage(0:ndtset_alloc)
 character(len=8),intent(in) :: strimg(mxnimage)
 integer,intent(in) :: narrm(0:ndtset_alloc)
 real(dp),intent(in) :: dprarr_images(marr,mxnimage,0:ndtset_alloc)

!Local variables-------------------------------
 integer :: idtset,iimage,jdtset
 integer :: intarr_images(marr,mxnimage,0:ndtset_alloc)
 character(len=4) :: appen
 character(len=16) :: keywd
 character(len=*), parameter :: format01160 ="(1x,a16,1x,(t22,3es18.10)) "
 character(len=*), parameter :: format01160a="(1x,a16,a,1x,(t22,3es18.10)) "

! *************************************************************************

 do idtset=1,ndtset_alloc

   if (narrm(idtset)>0)then
     do iimage=1,nimage(idtset)
       keywd=token//trim(strimg(iimage))

       if ((prtimg(iimage,idtset)==1).or.(ncid<0)) then
         if(ndtset>0)then

           jdtset=jdtset_(idtset)
           call appdig(jdtset,'',appen)
           if (prtimg(iimage,idtset)==1) write(iout,format01160a)trim(keywd),appen,dprarr_images(1:narrm(idtset),iimage,idtset)
           
#if defined HAVE_TRIO_ETSF_IO
           call write_var_netcdf(intarr_images(1:narrm(idtset),iimage,idtset),&
&           dprarr_images(1:narrm(idtset),iimage,idtset),&
&           marr,narrm(idtset),ncid,'DPR',trim(keywd)//appen)
#endif
         else
           if (prtimg(iimage,idtset)==1) write(iout,format01160)trim(keywd),dprarr_images(1:narrm(idtset),iimage,idtset)

#if defined HAVE_TRIO_ETSF_IO
           call write_var_netcdf(intarr_images(1:narrm(idtset),iimage,idtset),&
&           dprarr_images(1:narrm(idtset),iimage,idtset),&
&           marr,narrm(idtset),abs(ncid),'DPR',trim(keywd))
#endif
           
         end if
       end if
     end do
   end if
 end do

end subroutine prttagm_images
!!***
