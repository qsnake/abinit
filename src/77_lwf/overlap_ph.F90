!{\src2tex{textfont=tt}}
!!****f* ABINIT/overlap_ph
!! NAME
!! overlap_ph
!!
!! FUNCTION
!! Calculate the ovelaps between the phonons in the global window
!!  at each q point and its neighbours:
!!  <u_m,k|u_n,k+b>
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  eigvect = full eigenvectors matrix
!!  g_subsp = full global subspace (index of the bands)
!!  maxqsize = maxim number of bands in the global window
!!  natom = number of atoms
!!  nqpt = number of q points
!!  qneigh = index of the 6 neighbouring k+b q point
!!  qsize = size of the global subspace at each q point
!!
!! OUTPUT
!!  mmnkb=the M=<m,k|n,k+b> matrix
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


subroutine overlap_ph(eigvect,g_subsp,maxqsize,mmnkb,natom,nqpt,qneigh,qsize)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'overlap_ph'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: maxqsize,natom,nqpt
!arrays
 integer,intent(in) :: g_subsp(nqpt,3*natom),qneigh(nqpt,6),qsize(nqpt,3)
 real(dp),intent(in) :: eigvect(nqpt,3*natom,natom,3,2)
 real(dp),intent(out) :: mmnkb(nqpt,6,maxqsize,maxqsize,2)

!Local variables-------------------------------
!scalars
 integer :: aa,iqn,iqpt,mband,nband,tao
 real(dp) :: normqm,normqn

!******************************************************************
!BEGIN EXECUTABLE SECTION

!DEBUG
!write(std_out,*)
!write(std_out,*)
!write(std_out,*)
!write(std_out,*) 'overlap_ph : starting'
!write(std_out,*) ' nqpt=',nqpt
!write(std_out,*) ' maxqsize=',maxqsize
!ENDDEBUG

 do iqpt=1,nqpt            ! loop over all q-points
!  DEBUG
!  write(std_out,*) ' overlap_ph: no of bands G Z F',qsize(iqpt,:)
!  ENDDEBUG
   do iqn=1,6               ! loop over all q+b neighbours
     do mband=1,qsize(iqpt,1)       ! loop over the m global bands at Q
!      ! COMMENT
!      ! The overlap matrix is constructed on all the global bands, while in zmnbld.f only the free bands part is used

       do nband=1,qsize(qneigh(iqpt,iqn),1)  !loop over the n global bands, of the neighbours, Q+b
         mmnkb(iqpt,iqn,mband,nband,1)=zero
         mmnkb(iqpt,iqn,mband,nband,2)=zero
         normqn=zero
         normqm=zero

         do tao=1,natom        ! loop over atoms
           do aa=1,3            ! loop over directions
             mmnkb(iqpt,iqn,mband,nband,1)=mmnkb(iqpt,iqn,mband,nband,1)+&
&             eigvect(iqpt,g_subsp(iqpt,mband),tao,aa,1)*&
&             eigvect(qneigh(iqpt,iqn),g_subsp(qneigh(iqpt,iqn),nband),tao,aa,1)+&
&             eigvect(iqpt,g_subsp(iqpt,mband),tao,aa,2)*&
&             eigvect(qneigh(iqpt,iqn),g_subsp(qneigh(iqpt,iqn),nband),tao,aa,2)

             mmnkb(iqpt,iqn,mband,nband,2)=mmnkb(iqpt,iqn,mband,nband,2)+&
&             eigvect(iqpt,g_subsp(iqpt,mband),tao,aa,1)*&
&             eigvect(qneigh(iqpt,iqn),g_subsp(qneigh(iqpt,iqn),nband),tao,aa,2)-&
&             eigvect(iqpt,g_subsp(iqpt,mband),tao,aa,2)*&
&             eigvect(qneigh(iqpt,iqn),g_subsp(qneigh(iqpt,iqn),nband),tao,aa,1)

             normqn=normqn+&
&             eigvect(qneigh(iqpt,iqn),g_subsp(qneigh(iqpt,iqn),nband),tao,aa,1)*&
&             eigvect(qneigh(iqpt,iqn),g_subsp(qneigh(iqpt,iqn),nband),tao,aa,1)+&
&             eigvect(qneigh(iqpt,iqn),g_subsp(qneigh(iqpt,iqn),nband),tao,aa,2)*&
&             eigvect(qneigh(iqpt,iqn),g_subsp(qneigh(iqpt,iqn),nband),tao,aa,2)

             normqm=normqm+&
&             eigvect(iqpt,g_subsp(iqpt,mband),tao,aa,1)*eigvect(iqpt,g_subsp(iqpt,mband),tao,aa,1)+&
&             eigvect(iqpt,g_subsp(iqpt,mband),tao,aa,2)*eigvect(iqpt,g_subsp(iqpt,mband),tao,aa,2)

           end do                ! directions
         end do                 ! atoms

         mmnkb(iqpt,iqn,mband,nband,1)=mmnkb(iqpt,iqn,mband,nband,1)/dble(sqrt(normqn*normqm))
         mmnkb(iqpt,iqn,mband,nband,2)=mmnkb(iqpt,iqn,mband,nband,2)/dble(sqrt(normqn*normqm))

       end do                  ! n_bands at q+b
     end do                   ! m_bands at q

   end do                    ! b neighbours
 end do                     ! qpoints

!DEBUG
!do iqpt=1,nqpt            ! loop over all q-points
!do iqn=1,6               ! loop over all q+b neighbours
!do mband=1,qsize(iqpt,1)       ! loop over the m bands
!do nband=1,qsize(qneigh(iqpt,iqn),1)  !loop over the n bands, of the neighbours
!write(std_out,'(a,a,i4,a,i4,a,i4,a,i4,a,2i3,3f20.16)') ' mmnkb','band no',g_subsp(iqpt,mband),'at',iqpt,&
!'band no',g_subsp(qneigh(iqpt,iqn),nband),'at',qneigh(iqpt,iqn),' m,n ',mband,nband,&
!mmnkb(iqpt,iqn,mband,nband,1),mmnkb(iqpt,iqn,mband,nband,2),&
!mmnkb(iqpt,iqn,mband,nband,1)*mmnkb(iqpt,iqn,mband,nband,1)+mmnkb(iqpt,iqn,mband,nband,2)*mmnkb(iqpt,iqn,mband,nband,2)
!end do
!end do
!end do
!end do
!write(std_out,*) 'overlap_ph : end'
!ENDDEBUG

end subroutine overlap_ph

!!***
