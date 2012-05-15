!{\src2tex{textfont=tt}}
!!****f* ABINIT/bigbx9
!!
!! NAME
!! bigbx9
!!
!! FUNCTION
!! Generation of a Big Box containing all the R points in the
!! cartesian real space needed to Fourier Transforms the dynamical
!! matrix into its corresponding interatomic force.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (JCC,XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! brav= Bravais Lattice (1=S.C.;2=F.C.C.;3=BCC;4=Hex.)
!! choice= if 0, simply count nrpt ; if 1, checks that the input mrpt
!!  is the same as nrpt, and generate rpt(3,mrpt)
!! mrpt=dimension of rpt
!! ngqpt(3)= Numbers used to generate the q points to sample the
!!  Brillouin zone using an homogeneous grid
!! nqshft= number of q-points in the repeated cell for
!!  the Brillouin zone sampling
!!  When nqshft is not 1, but 2 or 4 (only other allowed
!!  values), the limits for the big box have to be extended by a factor of 2.
!! rprim(3,3)= Normalized coordinates in real space  !!! IS THIS CORRECT?
!!
!! OUTPUT
!! nprt= Number of R points in the Big Box
!! rpt(3,mrpt)= Canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3)!!)
!!  (output only if choice=1)
!!
!! PARENTS
!!      mkifc9
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine bigbx9(brav,choice,mrpt,ngqpt,nqshft,nrpt,rprim,rpt)

 use m_profiling

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bigbx9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav,choice,mrpt,nqshft
 integer,intent(out) :: nrpt
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: rprim(3,3)
 real(dp),intent(out) :: rpt(3,mrpt)

!Local variables -------------------------
!In some cases, the atoms coordinates are not packed in the
! [0,1]^3 cube. Then, the parameter "buffer" might be increased,
!to search relevant pairs of atoms in bigger boxes than usual.
!scalars
 integer,parameter :: buffer=1
 integer :: irpt,lim1,lim2,lim3,lqshft,r1,r2,r3
 character(len=500) :: msg

! *********************************************************************

 lqshft=1
 if(nqshft/=1)lqshft=2

!Simple Cubic Lattice

 if (brav==1) then
   lim1=((ngqpt(1)/2)+1)*lqshft+buffer
   lim2=((ngqpt(2)/2)+1)*lqshft+buffer
   lim3=((ngqpt(3)/2)+1)*lqshft+buffer
   nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)
   if(choice/=0)then
     if (nrpt/=mrpt) then
       write(msg,'(2(a,i0))')' nrpt=',nrpt,' is not equal to mrpt= ',mrpt
       MSG_BUG(msg)
     end if
     irpt=0
     do r1=-lim1,lim1
       do r2=-lim2,lim2
         do r3=-lim3,lim3
           irpt=irpt+1
           rpt(1,irpt)=r1*rprim(1,1)+r2*rprim(1,2)+r3*rprim(1,3)
           rpt(2,irpt)=r1*rprim(2,1)+r2*rprim(2,2)+r3*rprim(2,3)
           rpt(3,irpt)=r1*rprim(3,1)+r2*rprim(3,2)+r3*rprim(3,3)
         end do
       end do
     end do
   end if

!  Face Centered Cubic Lattice

 else if (brav==2) then
   lim1=((ngqpt(1)+3)/4)*lqshft+buffer
   lim2=((ngqpt(2)+3)/4)*lqshft+buffer
   lim3=((ngqpt(3)+3)/4)*lqshft+buffer
   nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)*4
   if(choice/=0)then
     if (nrpt/=mrpt) then
       write(msg,'(2(a,i0))')' nrpt=',nrpt,' is not equal to mrpt= ',mrpt
       MSG_BUG(msg)
     end if
     irpt=0
     do r1=-lim1,lim1
       do r2=-lim2,lim2
         do r3=-lim3,lim3
           irpt=irpt+4
           rpt(1,irpt-3)=r1
           rpt(2,irpt-3)=r2
           rpt(3,irpt-3)=r3
           rpt(1,irpt-2)=r1
           rpt(2,irpt-2)=r2+0.5
           rpt(3,irpt-2)=r3+0.5
           rpt(1,irpt-1)=r1+0.5
           rpt(2,irpt-1)=r2
           rpt(3,irpt-1)=r3+0.5
           rpt(1,irpt)=r1+0.5
           rpt(2,irpt)=r2+0.5
           rpt(3,irpt)=r3
         end do
       end do
     end do
   end if

!  Body Centered Cubic Lattice

 else if (brav==3) then
   lim1=((ngqpt(1)+3)/4)*lqshft+buffer
   lim2=((ngqpt(2)+3)/4)*lqshft+buffer
   lim3=((ngqpt(3)+3)/4)*lqshft+buffer
   nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)*2
   if(choice/=0)then
     if(nrpt/=mrpt) then
       write(msg,'(2(a,i0))')' nrpt= ',nrpt,' is not equal to mrpt= ',mrpt
       MSG_BUG(msg)
     end if
     irpt=0
     do r1=-lim1,lim1
       do r2=-lim2,lim2
         do r3=-lim3,lim3
           irpt=irpt+2
           rpt(1,irpt-1)=r1
           rpt(2,irpt-1)=r2
           rpt(3,irpt-1)=r3
           rpt(1,irpt)=r1+0.5
           rpt(2,irpt)=r2+0.5
           rpt(3,irpt)=r3+0.5
         end do
       end do
     end do
   end if

!  Hexagonal Lattice

 else if (brav==4) then
   lim1=(ngqpt(1)+1)*lqshft+buffer
   lim2=(ngqpt(2)+1)*lqshft+buffer
   lim3=((ngqpt(3)/2)+1)*lqshft+buffer
   nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)
   if(choice/=0)then
     if(nrpt/=mrpt)then
       write(msg,'(2(a,i0))')' nrpt=',nrpt,' is not equal to mrpt=',mrpt
       MSG_BUG(msg)
     end if
     irpt=0
     do r1=-lim1,lim1
       do r2=-lim2,lim2
         do r3=-lim3,lim3
           irpt=irpt+1
           rpt(1,irpt)=r1*rprim(1,1)+r2*rprim(1,2)+r3*rprim(1,3)
           rpt(2,irpt)=r1*rprim(2,1)+r2*rprim(2,2)+r3*rprim(2,3)
           rpt(3,irpt)=r1*rprim(3,1)+r2*rprim(3,2)+r3*rprim(3,3)
         end do
       end do
     end do
   end if

 else
   write(msg,'(a,i0,a)')' The value of brav= ',brav,' is not allowed (should be 1, 2 or 4).'
   MSG_BUG(msg)
 end if

end subroutine bigbx9
!!***
