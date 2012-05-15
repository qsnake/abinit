!{\src2tex{textfont=tt}}
!!****f* ABINIT/wanvec
!! NAME
!! wanvec
!!
!! FUNCTION
!! Builds the eigenvectors corresponding to the Wannier (interpolated) bands
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! atmass(natom)= square root of atomic masses
!! g_subsp(nqpt,3*natom)= index of the bands from the global window
!! iwf= identifier for the Wannier file
!! lambda(nqpt,nwnn,maxqsize,2)= interpolation coefficients of the Wannier bands
!! maxqsize= maximum size of the global window
!! natom= number of atoms per unit cell
!! nqpt= number of q points in the whole BZ
!! qpoint(nqpt,3)= coordinates of the q points
!! qsize(nqpt,3)= size of the windows
!! nwnn=no of wannier functions to be generated
!!
!! OUTPUT
!! wannvect(nqpt,nwnn,natom,3,2)= eigenvector of the Wannier states
!!
!! SIDE EFFECTS
!! eigvect(nqpt,3*natom,natom,3,2)=
!!
!! PARENTS
!!      lwf
!!
!! CHILDREN
!!      rwwan
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wanvec(atmass,eigvect,g_subsp,iwf,lambda,maxqsize,&
& natom,nqpt,nwnn,qpoint,qsize,wannvect)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wanvec'
 use interfaces_77_lwf, except_this_one => wanvec
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iwf,maxqsize,natom,nqpt,nwnn
!arrays
 integer,intent(in) :: g_subsp(nqpt,3*natom),qsize(nqpt,3)
 real(dp),intent(in) :: atmass(natom),lambda(nqpt,nwnn,maxqsize,2)
 real(dp),intent(inout) :: eigvect(nqpt,3*natom,natom,3,2),qpoint(nqpt,3)
 real(dp),intent(out) :: wannvect(nqpt,nwnn,natom,3,2)

!Local variables-------------------------------
!scalars
 integer :: aa,iband,iqpt,iwann,tao

!******************************************************************
!BEGIN EXECUTABLE SECTION

!DEBUG
!write(std_out,*) ' wanvec : enter'
!ENDDEBUG

 wannvect(:,:,:,:,:)=zero
 do iqpt=1,nqpt                        ! loop over Q points
!  write(std_out,'(a,i4,a,3f8.5)') ' Q point no.',iqpt,' at ',qpoint(iqpt,:)
!  BEAUTIFICATION NOTE : iout has been removed from the input variables
!  write(iout,'(a,i4,a,3f8.5)') ' Q point no.',iqpt,' at ',qpoint(iqpt,:)
   do iwann=1,nwnn                         ! loop over Wannier states
!    write(std_out,'(a,i4)') 'Wannier state no.',iwann
!    write(iout,'(a,i4)') 'Wannier state no.',iwann
     do iband=1,qsize(iqpt,1)            ! loop over global bands

!      write(std_out,'(a,3f12.7)') 'lambdas in the end',lambda(iqpt,iwann,iband,1),lambda(iqpt,iwann,iband,2),&
!      lambda(iqpt,iwann,iband,1)*lambda(iqpt,iwann,iband,1)+lambda(iqpt,iwann,iband,2)*lambda(iqpt,iwann,iband,2)
!      write(iout,'(a,3f12.7,f10.5)') 'lambdas in the end',lambda(iqpt,iwann,iband,1),lambda(iqpt,iwann,iband,2),&
!      lambda(iqpt,iwann,iband,1)*lambda(iqpt,iwann,iband,1)+lambda(iqpt,iwann,iband,2)*lambda(iqpt,iwann,iband,2),eigval(iqpt,g_subsp(iqpt,iband))

       do tao=1,natom
!        write(std_out,'(a,i5,3f14.9)') ' corresponding displacem',g_subsp(iqpt,iband),eigvect(iqpt,g_subsp(iqpt,iband),tao,1,1),&
!        eigvect(iqpt,g_subsp(iqpt,iband),tao,2,1),eigvect(iqpt,g_subsp(iqpt,iband),tao,3,1)
!        write(std_out,'(a,i5,3f14.9)') ' corresponding displacem',g_subsp(iqpt,iband),eigvect(iqpt,g_subsp(iqpt,iband),tao,1,2),&
!        eigvect(iqpt,g_subsp(iqpt,iband),tao,2,2),eigvect(iqpt,g_subsp(iqpt,iband),tao,3,2)
         eigvect(iqpt,g_subsp(iqpt,iband),tao,1,1)=eigvect(iqpt,g_subsp(iqpt,iband),tao,1,1)/atmass(tao)
         eigvect(iqpt,g_subsp(iqpt,iband),tao,2,1)=eigvect(iqpt,g_subsp(iqpt,iband),tao,2,1)/atmass(tao)
         eigvect(iqpt,g_subsp(iqpt,iband),tao,3,1)=eigvect(iqpt,g_subsp(iqpt,iband),tao,3,1)/atmass(tao)
         eigvect(iqpt,g_subsp(iqpt,iband),tao,1,2)=eigvect(iqpt,g_subsp(iqpt,iband),tao,1,2)/atmass(tao)
         eigvect(iqpt,g_subsp(iqpt,iband),tao,2,2)=eigvect(iqpt,g_subsp(iqpt,iband),tao,2,2)/atmass(tao)
         eigvect(iqpt,g_subsp(iqpt,iband),tao,3,2)=eigvect(iqpt,g_subsp(iqpt,iband),tao,3,2)/atmass(tao)
!        write(std_out,'(a,i5,3f14.9)') ' corresponding eigvector',g_subsp(iqpt,iband),eigvect(iqpt,g_subsp(iqpt,iband),tao,1,1),&
!        eigvect(iqpt,g_subsp(iqpt,iband),tao,2,1),eigvect(iqpt,g_subsp(iqpt,iband),tao,3,1)
!        write(std_out,'(a,i5,3f14.9)') ' corresponding eigvector',g_subsp(iqpt,iband),eigvect(iqpt,g_subsp(iqpt,iband),tao,1,2),&
!        eigvect(iqpt,g_subsp(iqpt,iband),tao,2,2),eigvect(iqpt,g_subsp(iqpt,iband),tao,3,2)
       end do

       do tao=1,natom                     ! loop over atoms
         do aa=1,3                         ! loop over directions
           wannvect(iqpt,iwann,tao,aa,1)=wannvect(iqpt,iwann,tao,aa,1)+lambda(iqpt,iwann,iband,1)*&
&           eigvect(iqpt,g_subsp(iqpt,iband),tao,aa,1)-lambda(iqpt,iwann,iband,2)*&
&           eigvect(iqpt,g_subsp(iqpt,iband),tao,aa,2)
           wannvect(iqpt,iwann,tao,aa,2)=wannvect(iqpt,iwann,tao,aa,2)+lambda(iqpt,iwann,iband,1)*&
&           eigvect(iqpt,g_subsp(iqpt,iband),tao,aa,2)+lambda(iqpt,iwann,iband,2)*&
&           eigvect(iqpt,g_subsp(iqpt,iband),tao,aa,1)
         end do                             !directions
       end do                              !atoms
     end do                               !global bands
   end do                                !Wannier states
 end do                                 !q points


!write the eigenvectors of the Wannier states in a separate file
 call rwwan(1,iwf,natom,nqpt,nwnn,qpoint,wannvect)



!ENDDEBUG
!write(std_out,*) ' wanvec : exit'
!ENDDEBUG

end subroutine wanvec
!!***
