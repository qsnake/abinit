!{\src2tex{textfont=tt}}
!!****f* ABINIT/shellin
!! NAME
!! shellin
!!
!! FUNCTION
!! TO BE DESCRIBED 090903 : not very explicit ...
!! Compute the shells q+b for each q point
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell = lattice parameters
!! nqpt = number of q points
!! qpoint = coordinates of q point
!! rprim = orientation of the lattice parameters
!!
!! OUTPUT
!! qneigh = matrix whose record index the neighbouring 6 k+b points for each k point
!!
!! NOTES
!! Now, works only for primitive cartesian cells.
!!
!! PARENTS
!!      lwf
!!
!! CHILDREN
!!      matr3inv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine shellin(acell,nqpt,qneigh,qpoint,rprim)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shellin'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqpt
!arrays
 integer,intent(out) :: qneigh(nqpt,6)
 real(dp),intent(in) :: acell(3),qpoint(nqpt,3),rprim(3,3)

!Local variables-------------------------------
!scalars
 integer :: contx,conty,contz,icorr,ii,iqpt,jj,jqpt
 real(dp) :: distmininit,distminx,distminy,distminz,redx,redy,redz,tmpcont
 real(dp) :: tmpvar,xi,xj,xx,yi,yj,yy,zi,zj,zz
!arrays
 real(dp) :: distx(3*nqpt,2),disty(3*nqpt,2),distz(3*nqpt,2),gprimd(3,3)
 real(dp) :: rprimd(3,3)

!******************************************************************
!BEGIN EXECUTABLE SECTION

!DEBUG
!write(std_out,*) ' shellin: enter'
!write(std_out,*) ' acell=',acell
!write(std_out,*) ' rprim=',rprim
!write(std_out,*) ' nqpt=',nqpt
!ENDDEBUG

!initialize the maximum distance
 do ii=1,3
   do jj=1,3
     rprimd(ii,jj)=rprim(ii,jj)*acell(jj)
   end do
 end do
 xx=acell(1)*rprimd(1,1)+acell(2)*rprimd(2,1)+acell(3)*rprimd(3,1)
 yy=acell(1)*rprimd(1,2)+acell(2)*rprimd(2,2)+acell(3)*rprimd(3,2)
 zz=acell(1)*rprimd(1,3)+acell(2)*rprimd(2,3)+acell(3)*rprimd(3,3)
 distmininit=xx*xx+yy*yy+zz*zz
 call matr3inv(rprimd,gprimd)


!DEBUG
!write(std_out,*) 'distmininit=',distmininit
!write(std_out,*) ' gprimd, rec space=',gprimd
!ENDDEBUG

!compute the overlaps
 do iqpt=1,nqpt          ! loop over all q points (=q)
!  DEBUG
!  write(std_out,'(a,3f10.5)') ' shellin : current qpoint :',qpoint(iqpt,:)
!  ENDDEBUG
   contx=0
   conty=0
   contz=0
   distx=zero
   disty=zero
   distz=zero
   distminx=distmininit
   distminy=distmininit
   distminz=distmininit
   xi=qpoint(iqpt,1)*gprimd(1,1)+qpoint(iqpt,2)*gprimd(2,1)+qpoint(iqpt,3)*gprimd(3,1)
   yi=qpoint(iqpt,1)*gprimd(1,2)+qpoint(iqpt,2)*gprimd(2,2)+qpoint(iqpt,3)*gprimd(3,2)
   zi=qpoint(iqpt,1)*gprimd(1,3)+qpoint(iqpt,2)*gprimd(2,3)+qpoint(iqpt,3)*gprimd(3,3)

!  DEBUG
!  write(std_out,*) ' xi,yi,zi=',xi,yi,zi
!  ENDDEBUG

   do jqpt=1,nqpt         ! loop over all q points (=q+b)
!    DEBUG
!    write(std_out,'(a,i4,3f10.5)') ' comparing with :',jqpt,qpoint(jqpt,:)
!    ENDDEBUG

     do icorr=1,3          ! loop over the three possible phase corrections
!      looking along z red coord
       if (qpoint(iqpt,1)==qpoint(jqpt,1) .and. qpoint(iqpt,2)==qpoint(jqpt,2)) then
!        DEBUG
!        write(std_out,*) 'along z'
!        ENDDEBUG
         contz=contz+1
         redz=qpoint(jqpt,3)-two+icorr
         xj=qpoint(jqpt,1)*gprimd(1,1)+qpoint(jqpt,2)*gprimd(2,1)+redz*gprimd(3,1)
         yj=qpoint(jqpt,1)*gprimd(1,2)+qpoint(jqpt,2)*gprimd(2,2)+redz*gprimd(3,2)
         zj=qpoint(jqpt,1)*gprimd(1,3)+qpoint(jqpt,2)*gprimd(2,3)+redz*gprimd(3,3)
         distz(contz,1)=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))
         distz(contz,2)=jqpt
!        DEBUG
!        write(std_out,*) 'distz=',distz((icorr-1)*nqpt+jqpt,:)
!        ENDDEBUG
       end if

!      looking along x red coord
       if (qpoint(iqpt,2)==qpoint(jqpt,2) .and. qpoint(iqpt,3)==qpoint(jqpt,3)) then
!        DEBUG
!        write(std_out,*) 'along x'
!        ENDDEBUG
         contx=contx+1
         redx=qpoint(jqpt,1)-two+icorr
         xj=redx*gprimd(1,1)+qpoint(jqpt,2)*gprimd(2,1)+qpoint(jqpt,3)*gprimd(3,1)
         yj=redx*gprimd(1,2)+qpoint(jqpt,2)*gprimd(2,2)+qpoint(jqpt,3)*gprimd(3,2)
         zj=redx*gprimd(1,3)+qpoint(jqpt,2)*gprimd(2,3)+qpoint(jqpt,3)*gprimd(3,3)
         distx(contx,1)=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))
         distx(contx,2)=jqpt
!        DEBUG
!        write(std_out,*) 'distx=',distx((icorr-1)*nqpt+jqpt,:)
!        ENDDEBUG
       end if

!      looking along y red coord
       if (qpoint(iqpt,3)==qpoint(jqpt,3) .and. qpoint(iqpt,1)==qpoint(jqpt,1)) then
!        DEBUG
!        write(std_out,*) 'along y'
!        ENDDEBUG
         conty=conty+1
         redy=qpoint(jqpt,2)-two+icorr
         xj=qpoint(jqpt,1)*gprimd(1,1)+redy*gprimd(2,1)+qpoint(jqpt,3)*gprimd(3,1)
         yj=qpoint(jqpt,1)*gprimd(1,2)+redy*gprimd(2,2)+qpoint(jqpt,3)*gprimd(3,2)
         zj=qpoint(jqpt,1)*gprimd(1,3)+redy*gprimd(2,3)+qpoint(jqpt,3)*gprimd(3,3)
         disty(conty,1)=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))
         disty(conty,2)=jqpt
!        DEBUG
!        write(std_out,*) 'disty=',disty((icorr-1)*nqpt+jqpt,:)
!        ENDDEBUG
       end if

     end do         ! three possible phase corrections
   end do          ! q points (=q+b)

!  ordering distances
!  along x,y,z
   do ii=1,contx
     do jj=ii+1,contx
       if (distx(ii,1)>distx(jj,1)) then
         tmpvar=distx(jj,1)
         distx(jj,1)=distx(ii,1)
         distx(ii,1)=tmpvar
         tmpcont=distx(jj,2)
         distx(jj,2)=distx(ii,2)
         distx(ii,2)=tmpcont
       end if
     end do
   end do
   do ii=1,conty
     do jj=ii+1,conty
       if (disty(ii,1)>disty(jj,1)) then
         tmpvar=disty(jj,1)
         disty(jj,1)=disty(ii,1)
         disty(ii,1)=tmpvar
         tmpcont=disty(jj,2)
         disty(jj,2)=disty(ii,2)
         disty(ii,2)=tmpcont
       end if
     end do
   end do
   do ii=1,contz
     do jj=ii+1,contz
       if (distz(ii,1)>distz(jj,1)) then
         tmpvar=distz(jj,1)
         distz(jj,1)=distz(ii,1)
         distz(ii,1)=tmpvar
         tmpcont=distz(jj,2)
         distz(jj,2)=distz(ii,2)
         distz(ii,2)=tmpcont
       end if
     end do
   end do


   qneigh(iqpt,1)=distx(2,2)
   qneigh(iqpt,2)=distx(3,2)
   qneigh(iqpt,3)=disty(2,2)
   qneigh(iqpt,4)=disty(3,2)
   qneigh(iqpt,5)=distz(2,2)
   qneigh(iqpt,6)=distz(3,2)

!  DEBUG
!  !write(std_out,'(a,3f10.5)') ' shellin : current qpoint :',qpoint(iqpt,:)
!  write(std_out,'(a,a,a,i4,a,3f8.5,a,a,i4,a,3f8.5,a,a,i4,a,3f8.5,a,a,i4,a,3f8.5,a,a,i4,a,3f8.5,a,a,i4,a,3f8.5,a)')&
!  & ' Neighbouring q-points :',ch10,&
!  &'along the x axis: no.',qneigh(iqpt,1),' at ',qpoint(qneigh(iqpt,1),:),ch10,&
!  &'along the x axis: no.',qneigh(iqpt,2),' at ',qpoint(qneigh(iqpt,2),:),ch10,&
!  &'along the y axis: no.',qneigh(iqpt,3),' at ',qpoint(qneigh(iqpt,3),:),ch10,&
!  &'along the y axis: no.',qneigh(iqpt,4),' at ',qpoint(qneigh(iqpt,4),:),ch10,&
!  &'along the z axis: no.',qneigh(iqpt,5),' at ',qpoint(qneigh(iqpt,5),:),ch10,&
!  &'along the z axis: no.',qneigh(iqpt,6),' at ',qpoint(qneigh(iqpt,6),:),ch10
!  ENDDEBUG

 end do           ! q points (=q)

end subroutine shellin
!!***
