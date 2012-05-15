!{\src2tex{textfont=tt}}
!!****f* ABINIT/indfftrisc
!!
!! NAME
!! indfftrisc
!!
!! FUNCTION
!! Take the data for sphere boundary and list of planewave in sphere (kg_k), manipulate them
!! for convenient use in fourwf, and output them in indpw_k
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gbound(2*mgfft+4)=sphere boundary data
!! kg_k(3,npw_k)=reduced planewave coordinates
!! mgfft=maximum size of 1D FFTs
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! npw_k=number of G vectors in basis at this k point
!!
!! OUTPUT
!! indpw_k(4,npw_k)=array which gives fft box index for given basis sphere
!!   in a representation that is directly usable by sg_fftrisc.f
!! ngb=number of FFTs along z
!!
!! PARENTS
!!      sg_fftrisc,sg_fftrisc_2
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine indfftrisc(gbound,indpw_k,kg_k,mgfft,ngb,ngfft,npw_k)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'indfftrisc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,npw_k
 integer,intent(out) :: ngb
!arrays
 integer,intent(in) :: gbound(2*mgfft+4),kg_k(3,npw_k),ngfft(18)
 integer,intent(out) :: indpw_k(4,npw_k)

!Local variables-------------------------------
!scalars
 integer :: g1,g2,i1,i2,i3,igb,index,ipw,n1,n2,n3
!arrays
 integer,allocatable :: index2d(:,:)
!character(len=500) :: message

! *************************************************************************
!
!DEBUG
!write(std_out,*)' indfftrisc : enter, debug '
!ENDDEBUG

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!First, generate a 2d index for each column of data
 ABI_ALLOCATE(index2d,(n1,n2))
 index2d(:,:)=0
 index=1
 igb=3
 do g2=0,gbound(2)
   do g1=0,gbound(igb+1)
     index2d(g1+1,g2+1)=index
     index=index+1
   end do
   if(gbound(igb)<=-1)then
     do g1=gbound(igb)+n1,n1-1
       index2d(g1+1,g2+1)=index
       index=index+1
     end do
   end if
   igb=igb+2
 end do
 if(gbound(1)<=-1)then
   do g2=gbound(1)+n2,n2-1
     do g1=0,gbound(igb+1)
       index2d(g1+1,g2+1)=index
       index=index+1
     end do
     if(gbound(igb)<=-1)then
       do g1=gbound(igb)+n1,n1-1
         index2d(g1+1,g2+1)=index
         index=index+1
       end do
     end if
     igb=igb+2
   end do
 end if

 ngb=index-1

!DEBUG
!write(std_out,*)' index2d '
!do i2=1,n2
!do i1=1,n1
!if(index2d(i1,i2)/=0)write(std_out,*)i1,i2,index2d(i1,i2)
!end do
!end do
!ENDDEBUG

!The 2d index has been generated
!Now, contract indpw_k(1,ipw) and indpw_k(2,ipw) into indpw_k(4,ipw)
!indpw_k(1,ipw) and indpw_k(2,ipw) are used to hold inverse of index2d,
!and for them, the second index does not fill 1:npw . It is only
!the number of z-transform FFTs.

!$OMP PARALLEL DO PRIVATE(ipw,i1,i2,i3) SHARED(n1,n2,n3,index2d,indpw_k,kg_k,npw_k)
 do ipw=1,npw_k
   i1=kg_k(1,ipw); if(i1<0)i1=i1+n1 ; i1=i1+1
   i2=kg_k(2,ipw); if(i2<0)i2=i2+n2 ; i2=i2+1
   i3=kg_k(3,ipw); if(i3<0)i3=i3+n3 ; i3=i3+1
   indpw_k(4,ipw)=index2d(i1,i2)
   indpw_k(3,ipw)=i3
 end do
!$OMP END PARALLEL DO

!DEBUG
!write(std_out,*)' indfftrisc : indpw,indpw_k(1:4) for each pw '
!write(std_out,*)' n1,n2,n3 ',n1,n2,n3
!do ipw=1,npw_k
!write(std_out,'(6i6)' )ipw,indpw(ipw),indpw_k(1:4,ipw)
!end do
!stop
!ENDDEBUG

 do i1=1,n1
   do i2=1,n2
     index=index2d(i1,i2)
     if(index/=0)then
       indpw_k(1,index)=i1
       indpw_k(2,index)=i2
     end if
   end do
 end do

!DEBUG
!write(std_out,*)' indfftrisc : indpw,indpw_k(1:4) for each pw '
!write(std_out,*)' n1,n2,n3 ',n1,n2,n3
!do ipw=1,npw_k
!write(std_out,'(6i6)' )ipw,indpw(ipw),indpw_k(1:4,ipw)
!end do
!stop
!ENDDEBUG

 ABI_DEALLOCATE(index2d)

end subroutine indfftrisc
!!***
