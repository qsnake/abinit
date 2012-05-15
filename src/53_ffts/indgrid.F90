!{\src2tex{textfont=tt}}
!!****f* ABINIT/indgrid
!!
!! NAME
!! indgrid
!!
!! FUNCTION
!! Calculate the correspondance between the coarse grid and
!! the fine grid for PAW calculations.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! nfftc=total number of FFt grid=n1*n2*n3 for the coarse grid
!! nfftf=total number of FFt grid=n1*n2*n3 for the fine grid
!! ngfftc(18)=contain all needed information about 3D FFT, for the coarse grid,
!!        see ~abinit/doc/input_variables/vargs.htm#ngfft
!! ngfftf(18)=contain all needed information about 3D FFT, for the fine grid,
!!        see ~abinit/doc/input_variables/vargs.htm#ngfft
!!
!! OUTPUT
!! coatofin(nfftc)= index of the points of the coarse grid on the fine grid
!! fintocoa(nfftf)=index of the points of the fine grid on the
!!   coarse grid (=0 if the point of the fine grid does not belong to
!!   the coarse grid).
!!
!! PARENTS
!!      fourier_interpol,m_paw_toolbox,m_rec
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine indgrid(coatofin,fintocoa,nfftc,nfftf,ngfftc,ngfftf)

 use m_profiling

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'indgrid'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftc,nfftf
!arrays
 integer,intent(in) :: ngfftc(18),ngfftf(18)
 integer,intent(out) :: coatofin(nfftc),fintocoa(nfftf)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,if1,if2,if3,ii,ing,n1c,n1f,n2c,n2f,n3c,n3f,narg1,narg2
!arrays
 integer :: id(3)
 integer,allocatable :: gc(:,:),gf(:,:)
 character(len=500) :: msg

! *************************************************************************
!

 DBG_ENTER("COLL")

 n1c=ngfftc(1);n2c=ngfftc(2);n3c=ngfftc(3)
 n1f=ngfftf(1);n2f=ngfftf(2);n3f=ngfftf(3)

 ABI_ALLOCATE(gc,(3,max(n1c,n2c,n3c)))
 do ii=1,3
   id(ii)=ngfftc(ii)/2+2
   do ing=1,ngfftc(ii)
     gc(ii,ing)=ing-(ing/id(ii))*ngfftc(ii)-1
   end do
 end do

 ABI_ALLOCATE(gf,(3,max(n1f,n2f,n3f)))
 do ii=1,3
   id(ii)=ngfftf(ii)/2+2
   do ing=1,ngfftf(ii)
     gf(ii,ing)=ing-(ing/id(ii))*ngfftf(ii)-1
   end do
 end do

 coatofin=0;fintocoa=0
 do i1=1,n1c
   do if1=1,n1f
     if(gc(1,i1)==gf(1,if1)) then
       do i2=1,n2c
         do if2=1,n2f
           if(gc(2,i2)==gf(2,if2)) then
             do i3=1,n3c
               do if3=1,n3f
                 if(gc(3,i3)==gf(3,if3)) then
                   narg1=i1+n1c*(i2-1+n2c*(i3-1))
                   narg2=if1+n1f*(if2-1+n2f*(if3-1))
                   coatofin(narg1)=narg2
                   fintocoa(narg2)=narg1
                   exit
                 end if
               end do
             end do
             exit ! To avoid N_fine * N_coarse scaling
           end if
         end do
       end do
       exit ! To avoid N_fine * N_coarse scaling
     end if
   end do
 end do

!Check coatofin to make sure there are no zeros!
 do ii=1,ubound(coatofin,1)
   if (coatofin(ii)==0) then
     msg = 'A zero was found in coatofin. Check that the fine FFT mesh is finer in each dimension than the coarse FFT mesh.'
     MSG_ERROR(msg)
   end if
 end do

 ABI_DEALLOCATE(gf)
 ABI_DEALLOCATE(gc)

 DBG_EXIT("COLL")

end subroutine indgrid
!!***
