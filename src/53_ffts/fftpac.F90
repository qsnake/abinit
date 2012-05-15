!{\src2tex{textfont=tt}}
!!****f* ABINIT/fftpac
!! NAME
!! fftpac
!!
!! FUNCTION
!! Allow for data copying to modify the stride (dimensioning) of a three-
!! dimensional array, for more efficient three dimensional fft.
!! Note that arrays aa and bb may be the same array (start at the same address).
!! The array aa also incorporate a spin variable.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, MF, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ispden=actual spin-density of interest
!!  nspden=number of spin-density components
!!  n1,n2,n3=actual data dimensions, dimensions of complex array a
!!  nd1,nd2,nd3=array dimensions of (larger) array b
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  option= see description of side effects
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  aa & bb arrays are treated as input or output depending on option: 
!!  option=1  aa(n1*n2*n3,ispden) <-- bb(nd1,nd2,nd3) real case
!!  option=2  aa(n1*n2*n3,ispden) --> bb(nd1,nd2,nd3) real case
!!  option=10 aa(n1*n2*n3,ispden) <-- bb(nd1,nd2,nd3) complex case like option 1 real part
!!  option=11 aa(n1*n2*n3,ispden) <-- bb(nd1,nd2,nd3) complex case like option 1 imag part
!!
!! PARENTS
!!      dens_in_sph,energy,ks_ddiago,ladielmt,lavnl,mkrho,mkrho3,nstpaw3
!!      prctfvw1,prctfvw2,resp3dte,rhofermi3,suscep_dyn,suscep_kxc_dyn
!!      suscep_stat,vtorho,vtorho3
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fftpac(ispden,nspden,n1,n2,n3,nd1,nd2,nd3,ngfft,aa,bb,option)

 use m_profiling

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftpac'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ispden,n1,n2,n3,nd1,nd2,nd3,nspden,option
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: aa(n1*n2*n3/ngfft(10),nspden),bb(nd1,nd2,nd3)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,index,me_fft,nproc_fft
 character(len=500) :: message

! *************************************************************************

 me_fft=ngfft(11) ; nproc_fft=ngfft(10)


 if(option==1.or.option==2) then
   if (nd1<n1.or.nd2<n2.or.nd3<n3) then
     write(message, '(4a,3i8,2a,3i8,a)' ) ch10,&
&     ' fftpac: BUG -',ch10,&
&     '  Each of nd1,nd2,nd3=',nd1,nd2,nd3,ch10,&
&     '  must be >=      n1, n2, n3 =',n1,n2,n3,'.'
     MSG_BUG(message)
   end if
 else
   if (2*nd1<n1.or.nd2<n2.or.nd3<n3) then
     write(message, '(4a,3i8,2a,3i8,a)' ) ch10,&
&     ' fftpac: BUG -',ch10,&
&     '  Each of 2*nd1,nd2,nd3=',2*nd1,nd2,nd3,ch10,&
&     '  must be >=      n1, n2, n3 =',n1,n2,n3,'.'
     MSG_BUG(message)
   end if
 end if

 if (option==1) then
   do i3=1,n3
     if (((i3-1)/(n3/nproc_fft))==me_fft) then
       do i2=1,n2
         do i1=1,n1
           aa(i1+n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1)),ispden)=bb(i1,i2,i3)
         end do
       end do
     end if
   end do

 else if (option==2) then
!  Here we avoid corrupting the data in a while writing to b in the
!  case in which a and b are same array.
!  Also: replace "trash" data with 0 s to avoid floating point
!  exceptions when this data is actually manipulated in fft.

   do i3=nd3,n3+1,-1
     do i2=nd2,1,-1
       do i1=nd1,1,-1
         bb(i1,i2,i3)=0.d0
       end do
     end do
   end do
   do i3=n3,1,-1
     if (((i3-1)/(n3/nproc_fft))==me_fft) then
       do i2=nd2,n2+1,-1
         do i1=nd1,1,-1
           bb(i1,i2,i3)=0.d0
         end do
       end do
       do i2=n2,1,-1
         do i1=nd1,n1+1,-1
           bb(i1,i2,i3)=0.d0
         end do
         do i1=n1,1,-1
           bb(i1,i2,i3)=aa(i1+n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1)),ispden)
         end do
       end do
     end if
   end do
!  MF
 else if (option==10 .or. option==11) then
   index=1
   if(option==11) index=2
   do i3=1,n3
     do i2=1,n2
       do i1=1,n1/2
         aa(index,ispden)=bb(i1,i2,i3)
         index=index+2
       end do
     end do
   end do
!  MF
 else
   write(message, '(a,a,a,a,i12,a)' ) ch10,&
&   ' fftpac: BUG -',ch10,&
&   '  Bad option =',option,'.'
   MSG_BUG(message)
 end if

end subroutine fftpac
!!***
