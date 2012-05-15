!{\src2tex{textfont=tt}}
!!****f* ABINIT/bldlwf
!! NAME
!! bldlwf
!!
!! FUNCTION
!! Builds the lattice Wannier functions:
!! WARNING : does not work ... wannvect is not yet initialized ...
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! grdsize= size of the grid of q points = limit of the shells of the LWF
!! natom= number of atoms per unit cell
!! nqpt= number of q points in the whole BZ
!! qpoint(nqpt,3)= coordinates of the q points
!! rcenter(3)= center of the lattice wannier functions
!! nwnn=no of wannier functions to be generated
!!
!! OUTPUT
!! wannvect(nqpt,nwnn,natom,3,2)= eigenvector of the Wannier states
!!
!! PARENTS
!!      lwf
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine bldlwf(grdsize,iout,natom,nqpt,nwnn,qpoint,rcenter,wannvect)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bldlwf'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,natom,nqpt,nwnn
!arrays
 integer,intent(in) :: grdsize(3)
 real(dp) :: wannvect(nqpt,nwnn,natom,3,2)
 real(dp),intent(in) :: qpoint(nqpt,3),rcenter(3)

!Local variables-------------------------------
!scalars
 integer :: aa,iqpt,iwann,shx,shy,shz,tao
 real(dp) :: phi
!arrays
 real(dp) :: LWF(nwnn,2*grdsize(1)+1,2*grdsize(2)+1,2*grdsize(3)+1,natom,3,2)

!******************************************************************
!BEGIN EXECUTABLE SECTION
!DEBUG
!write(std_out,*) 'bldlwf : enter'
!write(std_out,*) ' nwnn=',nwnn
!write(std_out,*) ' grdsize=',grdsize
!write(std_out,*) ' nqpt=',nqpt
!write(std_out,*) ' natom=',natom
!do iqpt=1,nqpt
!write(std_out,'(a,i4,3f8.5)') 'qpt',iqpt,qpoint(iqpt,:)
!end do
!ENDDEBUG

 LWF(:,:,:,:,:,:,:)=zero
 do iwann=1,nwnn                    !LWF
   do shx=1,2*grdsize(1)+1             !LWF shell along x
     do shy=1,2*grdsize(2)+1            !LWF shell along y
       do shz=1,2*grdsize(3)+1           !LWF shell along z

         do iqpt=1,nqpt                 !loop over q points
           phi=2*pi* ( qpoint(iqpt,1)*(rcenter(1)+dble(shx-grdsize(1)-1)) + &
&           qpoint(iqpt,2)*(rcenter(2)+dble(shy-grdsize(2)-1))  + qpoint(iqpt,3)*(rcenter(3)+dble(shz-grdsize(3)-1)) )

!          write(std_out,'(a,i4,3f8.5,f12.8)') 'qpt',iqpt,qpoint(iqpt,:),phi

           do tao=1,natom                !loop over atoms
             do aa=1,3                    !loop over cartesian directions
               LWF(iwann,shx,shy,shz,tao,aa,1)=LWF(iwann,shx,shy,shz,tao,aa,1)+&
&               (cos(phi)*wannvect(iqpt,iwann,tao,aa,1)+sin(phi)*wannvect(iqpt,iwann,tao,aa,2))/nqpt
               LWF(iwann,shx,shy,shz,tao,aa,2)=LWF(iwann,shx,shy,shz,tao,aa,2)+&
&               (cos(phi)*wannvect(iqpt,iwann,tao,aa,2)-sin(phi)*wannvect(iqpt,iwann,tao,aa,1))/nqpt
             end do
           end do                        !atoms
         end do                         ! qpoints

!        DEBUG
!        do tao=1,natom
!        write(std_out,'(a,5i4,3f12.8)') 'ATOM - SHELL X Y Z - LWF - REAL',tao,shx-grdsize(1)-1,shy-grdsize(2)-1,shz-grdsize(3)-1,iwann,&
!        LWF(iwann,shx,shy,shz,tao,1,1),LWF(iwann,shx,shy,shz,tao,2,1),&
!        LWF(iwann,shx,shy,shz,tao,3,1)
!        write(std_out,'(a,5i4,3f12.8)') 'ATOM - SHELL X Y Z - LWF - IMAG',tao,shx-grdsize(1)-1,shy-grdsize(2)-1,shz-grdsize(3)-1,iwann,&
!        LWF(iwann,shx,shy,shz,tao,1,2),LWF(iwann,shx,shy,shz,tao,2,2),&
!        LWF(iwann,shx,shy,shz,tao,3,2)
!        write(std_out,'(a,5i4,3f12.8)') 'ATOM - SHELL X Y Z - LWF - ABS ',tao,shx-grdsize(1)-1,shy-grdsize(2)-1,shz-grdsize(3)-1,iwann,&
!        sqrt(LWF(iwann,shx,shy,shz,tao,1,1)*LWF(iwann,shx,shy,shz,tao,1,1)+LWF(iwann,shx,shy,shz,tao,1,2)*LWF(iwann,shx,shy,shz,tao,1,2)),&
!        sqrt(LWF(iwann,shx,shy,shz,tao,2,1)*LWF(iwann,shx,shy,shz,tao,2,1)+LWF(iwann,shx,shy,shz,tao,2,2)*LWF(iwann,shx,shy,shz,tao,2,2)),&
!        sqrt(LWF(iwann,shx,shy,shz,tao,3,1)*LWF(iwann,shx,shy,shz,tao,3,1)+LWF(iwann,shx,shy,shz,tao,3,2)*LWF(iwann,shx,shy,shz,tao,3,2))
!        end do
!        ENDDEBUG

       end do                           !LWF z
     end do                            !LWF y
   end do                             !LWF x
 end do                              ! Wannier functions


!write in the output files
 do iwann=1,nwnn
   do tao=1,natom
     do shx=1,2*grdsize(1)+1
       do shy=1,2*grdsize(2)+1
         do shz=1,2*grdsize(3)+1
!          write(std_out,'(a,5i4,3f12.8)') 'ATOM - SHELL X Y Z - LWF - REAL',tao,shx-grdsize(1)-1,shy-grdsize(2)-1,shz-grdsize(3)-1,iwann,&
!          LWF(iwann,shx,shy,shz,tao,1,1),LWF(iwann,shx,shy,shz,tao,2,1),LWF(iwann,shx,shy,shz,tao,3,1)
           write(iout,'(a,5i4,3f12.8)')&
&           'ATOM - SHELL X Y Z - LWF - REAL',tao,shx-grdsize(1)-1,shy-grdsize(2)-1,shz-grdsize(3)-1,iwann,&
&           LWF(iwann,shx,shy,shz,tao,1,1),LWF(iwann,shx,shy,shz,tao,2,1),LWF(iwann,shx,shy,shz,tao,3,1)

         end do
       end do
     end do
   end do

   do tao=1,natom
     do shx=1,2*grdsize(1)+1
       do shy=1,2*grdsize(2)+1
         do shz=1,2*grdsize(3)+1
!          write(std_out,'(a,5i4,3f12.8)') 'ATOM - SHELL X Y Z - LWF - IMAG',tao,shx-grdsize(1)-1,shy-grdsize(2)-1,shz-grdsize(3)-1,iwann,&
!          LWF(iwann,shx,shy,shz,tao,1,2),LWF(iwann,shx,shy,shz,tao,2,2),LWF(iwann,shx,shy,shz,tao,3,2)
           write(iout,'(a,5i4,3f12.8)')&
&           'ATOM - SHELL X Y Z - LWF - IMAG',tao,shx-grdsize(1)-1,shy-grdsize(2)-1,shz-grdsize(3)-1,iwann,&
&           LWF(iwann,shx,shy,shz,tao,1,2),LWF(iwann,shx,shy,shz,tao,2,2),LWF(iwann,shx,shy,shz,tao,3,2)

         end do
       end do
     end do
   end do

   do tao=1,natom
     do shx=1,2*grdsize(1)+1
       do shy=1,2*grdsize(2)+1
         do shz=1,2*grdsize(3)+1
!          write(std_out,'(a,5i4,3f12.8)') 'ATOM - SHELL X Y Z - LWF - ABS ',tao,shx-grdsize(1)-1,shy-grdsize(2)-1,shz-grdsize(3)-1,iwann,&
!          sqrt(LWF(iwann,shx,shy,shz,tao,1,1)*LWF(iwann,shx,shy,shz,tao,1,1)+LWF(iwann,shx,shy,shz,tao,1,2)*LWF(iwann,shx,shy,shz,tao,1,2)),&
!          sqrt(LWF(iwann,shx,shy,shz,tao,2,1)*LWF(iwann,shx,shy,shz,tao,2,1)+LWF(iwann,shx,shy,shz,tao,2,2)*LWF(iwann,shx,shy,shz,tao,2,2)),&
!          sqrt(LWF(iwann,shx,shy,shz,tao,3,1)*LWF(iwann,shx,shy,shz,tao,3,1)+LWF(iwann,shx,shy,shz,tao,3,2)*LWF(iwann,shx,shy,shz,tao,3,2))
           write(iout,'(a,5i4,3f12.8)') 'ATOM - SHELL X Y Z - LWF - ABS ',&
&           tao,shx-grdsize(1)-1,shy-grdsize(2)-1,shz-grdsize(3)-1,iwann,&
&           sqrt(LWF(iwann,shx,shy,shz,tao,1,1)**2+LWF(iwann,shx,shy,shz,tao,1,2)**2),&
&           sqrt(LWF(iwann,shx,shy,shz,tao,2,1)**2+LWF(iwann,shx,shy,shz,tao,2,2)**2),&
&           sqrt(LWF(iwann,shx,shy,shz,tao,3,1)**2+LWF(iwann,shx,shy,shz,tao,3,2)**2)

         end do
       end do
     end do
   end do

 end do

!DEBUG
!write(std_out,*) 'bldlwf : exit'
!ENDDEBUG

end subroutine bldlwf
!!***
