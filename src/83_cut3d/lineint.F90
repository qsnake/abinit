!{\src2tex{textfont=tt}}
!!****f* ABINIT/lineint
!! NAME
!! lineint
!!
!! FUNCTION
!! Computes the values along a line
!! defined by two points
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (GMR,RC,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gridtt(nr1,nr2,nr3)=Total density
!! gridux(nr1,nr2,nr3)=spin-Up density, or magnetization density in X direction
!! griddy(nr1,nr2,nr3)=spin-Down density, or magnetization density in Y direction
!! gridmz(nr1,nr2,nr3)=spin-polarization density or magnetization density in Z direction
!! nr1=grid size along x
!! nr2=grid size along y
!! nr3=grid size along z
!! nspden=number of spin-density components
!! rprimd(3,3)=orientation of the unit cell in 3D
!!
!! OUTPUT
!!  only writing
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!      interpol3d,normalize,reduce
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


 subroutine lineint(gridtt,gridux,griddy,gridmz,nr1,nr2,nr3,nspden,rprimd)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lineint'
 use interfaces_32_util
 use interfaces_83_cut3d, except_this_one => lineint
!End of the abilint section

 implicit none

!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: nr1,nr2,nr3,nspden
!arrays
 real(dp),intent(in) :: griddy(nr1,nr2,nr3)
 real(dp),intent(in) :: gridmz(nr1,nr2,nr3),gridtt(nr1,nr2,nr3)
 real(dp),intent(in) :: gridux(nr1,nr2,nr3),rprimd(3,3)

!Local variables--------------------------------------------------------
!scalars
 integer :: inpopt,inpopt2,k2,nresol,okline
 real(dp) :: denvaldy,denvalmz,denvaltt,denvalux,dx,dy,dz,length
 character(len=fnlen) :: filnam
!arrays
 real(dp) :: cent(3),r1(3),r2(3),rcart(3),rr(3),x1(3),x2(3)

! *********************************************************************

 okline=0
 do while (okline==0)
   write(std_out,*) ' Type 1) for a line between two cartesian-defined points'
   write(std_out,*) '   or 2) for a line between two crystallographic-defined points '
   write(std_out,*) '   or 3) for a line defined by its direction in cartesion coordinates'
   write(std_out,*) '   or 4) for a line defined by its direction in crystallographic coordinates'
   read(*,*) inpopt
   write(std_out,*) ' You typed ',inpopt,ch10
   if (inpopt==1 .or. inpopt ==2 .or. inpopt==3 .or. inpopt==4) okline=1
 end do

!In the case of a line defined by its two extreme points
 if (inpopt==1) then
   write(std_out,*) ' Type the first point coordinates (Bohrs):'
   write(std_out,*) '    -> X-dir   Y-dir   Z-dir:'
   read(*,*) x1
   write(std_out,'(a,3es16.6,a)') ' You typed ',x1,ch10
   call reduce(r1,x1,rprimd)

   write(std_out,*) ' Type the second point coordinates (Bohrs):'
   write(std_out,*) '    -> X-dir   Y-dir   Z-dir:'
   read(*,*) x2
   write(std_out,'(a,3es16.6,a)') ' You typed ',x2,ch10
   call reduce(r2,x2,rprimd)
 end if

 if (inpopt==2) then
   write(std_out,*) ' Type the first point coordinates (fractional):'
   write(std_out,*) '    -> X-dir   Y-dir   Z-dir:'
   read(*,*) r1
   write(std_out,'(a,3es16.6,a)') ' You typed ',r1,ch10

   write(std_out,*) ' Type the second point coordinates (fractional):'
   write(std_out,*) '    -> X-dir   Y-dir   Z-dir:'
   read(*,*) r2
   write(std_out,'(a,3es16.6,a)') ' You typed ',r2,ch10
 end if

 if(inpopt==3 .or. inpopt==4 )then

   write(std_out,*) 'Please enter now the line direction:'
   write(std_out,*) '    -> X-dir   Y-dir   Z-dir:'
   read(*,*) x2
   write(std_out,'(a,3es16.6,a)') 'The line direction is:',x2(1),x2(2),x2(3),ch10

   if (inpopt == 4) then
     rcart=matmul(x2,rprimd)
     x2(:)=rcart(:)
     write(std_out,'(a,3es16.6,a)') 'Expressed in cartesian coordinates: ',x2(1),x2(2),x2(3),ch10
   end if

   call normalize(x2)

   write(std_out,*) 'Enter now the central point of line:'
   write(std_out,*) 'Type 1) for cartesian coordinates'
   write(std_out,*) '  or 2) for crystallographic coordinates'
   read(*,*) inpopt2
   if (inpopt2==1 .or. inpopt2==2) then
     write(std_out,*) 'Type the point coordinates:'
     write(std_out,*) '    -> X-Coord   Y-Coord   Z-Coord:'
     read(*,*) cent
     write(std_out,'(a,3es16.6,a)') 'Central point coordinates:', cent(1),cent(2),cent(3),ch10
     if (inpopt2==2) then
       rcart=matmul(cent,rprimd)
       cent(:)=rcart(:)
       write(std_out,'(a,3es16.6,a)') 'Expressed in cartesian coordinates:',cent(1),cent(2),cent(3),ch10
     end if
     write(std_out,*) 'Enter line length (in cartesian coordinates, in Bohr):'
     read(*,*) length

!    Compute the extremal points in cartesian coordinates
     x1(:)=cent(:)-length*x2(:)*half
     x2(:)=cent(:)+length*x2(:)*half

!    Transfer to crystallographic coordinates
     call reduce(r1,x1,rprimd)
     call reduce(r2,x2,rprimd)

   end if

 end if ! inpopt

 write(std_out,*)
 write(std_out,'(a,3es16.6)' ) ' Crystallographic coordinates of the first point  :',r1
 write(std_out,'(a,3es16.6)' ) ' Crystallographic coordinates of the second point :',r2
 write(std_out,*)

 write(std_out,*) '  Enter line resolution:   (integer, number of points on the line)'
 read(*,*) nresol
 write(std_out,*) ' You typed',nresol,ch10

!At this moment the code knows everything about the geometric input, the data and
!the line direction. It will further calculate the values along this line using
!an interpolation

 write(std_out,*) ch10,'  Enter the name of an output file:'
 read(*,*) filnam
 write(std_out,*) '  The name of your file is : ',trim(filnam),ch10

 open(unit=31,file=trim(filnam),status='unknown')

 dx=(r2(1)-r1(1))/nresol
 dy=(r2(2)-r1(2))/nresol
 dz=(r2(3)-r1(3))/nresol

!DEBUG
!write(std_out,*)' nspden=',nspden
!ENDDEBUG

 if(nspden==1)then
   write(std_out,*)' Index of point   value '
 else if (nspden==2)then
   write(std_out,*)' Index of point   non-spin-polarized   spin up       spin down     difference '
 else if (nspden==4)then
   write(std_out,*)' Index of point   non-spin-polarized      x              y              z '
 end if

 do k2=0,nresol

   rr(1)=r1(1)+k2*dx
   rr(2)=r1(2)+k2*dy
   rr(3)=r1(3)+k2*dz

   rr(1)=mod(mod(rr(1),1._dp)+1._dp,1._dp)
   rr(2)=mod(mod(rr(2),1._dp)+1._dp,1._dp)
   rr(3)=mod(mod(rr(3),1._dp)+1._dp,1._dp)

   call interpol3d(rr,nr1,nr2,nr3,denvaltt,gridtt)
   if(nspden==1)then
     write(31, '(i13,es22.12)' ) k2,denvaltt
     write(std_out,'(i13,es22.12)' ) k2,denvaltt
   else if(nspden==2 .or. nspden==4)then
     call interpol3d(rr,nr1,nr2,nr3,denvalux,gridux)
     call interpol3d(rr,nr1,nr2,nr3,denvaldy,griddy)
     call interpol3d(rr,nr1,nr2,nr3,denvalmz,gridmz)
     write(31, '(i13,4(es22.12))' ) k2,denvaltt,denvalux,denvaldy,denvalmz
     write(std_out,'(i13,4es22.12)' ) k2,denvaltt,denvalux,denvaldy,denvalmz
   end if
 end do

 close(31)

end subroutine lineint
!!***
