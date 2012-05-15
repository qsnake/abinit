!{\src2tex{textfont=tt}}
!!****f* ABINIT/zmnbld
!! NAME
!! zmnbld
!!
!! FUNCTION
!! TO BE DESCRIBED 090903 : not very explicit ...
!! Builds the Z_mn(q) matrices:
!! Z_mn(q) = Sum_w ( Sum_b ( < u_m(q) | lambda_w(q+b) > < lambda_w(q+b) | u_n(q) >  ) )
!!         = Sum_w ( Sum_b ( Sum_ii ( < u_m(q) | u_ii(q+b) > | lambda_w(q+b) > < lambda_w(q+b) | < u_ii(q+b) | | u_n(q) >  ) ) )
!!         = Sum_w ( Sum_b ( Sum_ii (   M_m,ii,q,b           | lambda_w(q+b) > < lambda_w(q+b) | M_ii,n,b,q                ) ) )
!!         = Sum_w ( Sum_b ( Sum_ii (   M_m,ii,q,b           | lambda_w(q+b) > < lambda_w(q+b) | (M_n,ii,q,b)*             ) ) )
!! bl1 = M_m,ii,q,b | lambda_w(q+b) >    = (mr + i mi)(lr + i li) = (lr.mr - li.mi) + i (lr.mi + li.mr)
!! bl2 = < lambda_w(q+b) | (M_n,ii,q,b)* = (lr - i li)(mr - i mi) = (lr.mr - li.mi) - i (lr.mi + li.mr)
!! Z_mn(q) = bl1.bl2
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! f_subsp = indexes of the bands in the free window
!! lambda(nqpt,nwnn,maxqsize,2)= output of secinit in the first time, then from the
!!                               previous diagonalization of the Zmn matrix
!! maxqsize= max number of free bands
!! mmnkb(nqpt,6,maxqsize,maxqsize,2) = <u_m(q)|u_n(q+b)>
!! natom= number of atoms per unit cell
!! nqpt= number of q points in the whole BZ
!! nwnn=no of wannier functions to be generated
!! nwnz(nqpt)= number of Wannier states in the frozen window
!! qneigh(nqpt,6)= index of the neighbouring b points for each q point
!! qsize(nqpt,3)= number of the bands in the energy window : G - Z - F
!!
!! OUTPUT
!! Zmn((nqpt,nwnn,maxqsize,maxqsize,2)=Z matrix
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


subroutine zmnbld(f_subsp,lambda,maxqsize,mmnkb,natom,nqpt,nwnn,nwnz,qneigh,qsize,Zmn)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'zmnbld'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: maxqsize,natom,nqpt,nwnn
!arrays
 integer,intent(in) :: f_subsp(nqpt,3*natom),nwnz(nqpt)
 integer,intent(in) :: qneigh(nqpt,6),qsize(nqpt,3)
 real(dp),intent(in) :: lambda(nqpt,nwnn,maxqsize,2)
 real(dp),intent(in) :: mmnkb(nqpt,6,maxqsize,maxqsize,2)
 real(dp),intent(out) :: Zmn(nqpt,maxqsize,maxqsize,2)

!Local variables-------------------------------
!scalars
 integer :: ii,iqpt,ishell,iwann,mband,mf,nband,nf
 real(dp) :: nshell
!arrays
 real(dp) :: bl1(2),bl2(2)

!******************************************************************
!BEGIN EXECUTABLE SECTION

!DEBUG
!write(std_out,*) ' zmnbld : enter'
!write(std_out,*) ' zmnbld : nqpt',nqpt
!write(std_out,*) ' zmnbld : maxqsize',maxqsize
!ENDDEBUG

 Zmn=zero

!number of neighbouring q+b vectors
!NOTE: in the current version, the no. of b's is always 6        (')

 nshell=6.0

 do iqpt=1,nqpt                 ! loop over Q points

!  DEBUG
!  write(std_out,'(a,2i4,a)') ' zmnbld : current Q point',iqpt,nwnz(iqpt)
!  write(std_out,'(a,i8)') ' zmnbld : current Q point',iqpt
!  do ishell=1,6
!  write(std_out,'(a,2i6)') ' b vector:',ishell,qneigh(iqpt,ishell)
!  end do
!  write(std_out,*) ' zmnbld : with qsize G Z F ',qsize(iqpt,:)
!  do ii=1,nwnn
!  write(std_out,*) ' zmnbld : current Wannier state:', ii
!  do mband=1,qsize(iqpt,1)
!  write(std_out,'(a,2f12.7)') ' lambda(iqpt,nwnn)',lambda(iqpt,ii,mband,:)
!  end do
!  end do
!  ENDDEBUG

   if (nwnz(iqpt)<nwnn) then
     mf=0

     do mband=1,qsize(iqpt,1)         ! loop over the global (m)bands at Q
!      write(std_out,*) ' mband =',mband,f_subsp(iqpt,mband)
       if (f_subsp(iqpt,mband)>0) then    !select free bands
         mf=mf+1
         nf=0
!        do nband=1,qsize(iqpt,1)        ! loop over the global (n)bands at Q
         do nband=1,mband                ! accelerate the computation of Z matrix
!          write(std_out,*) ' nband =',nband,f_subsp(iqpt,nband)
           if (f_subsp(iqpt,nband)>0) then  !select free bands
             nf=nf+1
             do iwann=1,nwnn                ! loop over Wannier states
!              write(std_out,*) ' iwann =',iwann
               do ishell=1,6
!                write(std_out,'(a,i3,a,i3,3i4)') ' ishell=',ishell,'q point:',qneigh(iqpt,ishell),qsize(qneigh(iqpt,ishell),:)
                 do ii=1,qsize(qneigh(iqpt,ishell),1)   ! loop over j(global) bands at Q+b
!                  bl1 = M_m,ii,q,b | lambda_w(q+b) >    = (mr + i mi)(lr + i li) = (lr.mr - li.mi) + i (lr.mi + li.mr)
!                  bl2 = < lambda_w(q+b) | (M_n,ii,q,b)* = (lr - i li)(mr - i mi) = (lr.mr - li.mi) - i (lr.mi + li.mr)
                   bl1=zero
                   bl2=zero
!                  BEAUTIFICATION NOTE : g_subsp has been removed from the input variables
!                  write(std_out,'(a,3i4,2f10.6,4i4,2f10.6,4i4)') 'mmnkbs',mband,nband,ii,&
!                  &mmnkb(iqpt,ishell,mband,ii,1),mmnkb(iqpt,ishell,mband,ii,2),&
!                  &iqpt,ishell,f_subsp(iqpt,mband),g_subsp(qneigh(iqpt,ishell),ii),&
!                  &mmnkb(iqpt,ishell,nband,ii,1),mmnkb(iqpt,ishell,nband,ii,2),&
!                  &iqpt,ishell,f_subsp(iqpt,nband),g_subsp(qneigh(iqpt,ishell),ii)
!                  write(std_out,'(a,2f10.6,a)') 'lambdas',lambda(qneigh(iqpt,ishell),iwann,ii,1),lambda(qneigh(iqpt,ishell),iwann,ii,2)
!                  write(std_out,'(a,i3,i4,i3,i4)') 'm,fsub_m and n,fsub_n',mband,f_subsp(iqpt,mband),nband,f_subsp(iqpt,nband)
!                  write(std_out,'(a,5i4)') 'ii,ishell,iqpt,qnei,gsubsp',ii,ishell,iqpt,qneigh(iqpt,ishell),g_subsp(qneigh(iqpt,ishell),ii)

                   bl1(1)=bl1(1)+lambda(qneigh(iqpt,ishell),iwann,ii,1)*mmnkb(iqpt,ishell,mband,ii,1)-&
&                   lambda(qneigh(iqpt,ishell),iwann,ii,2)*mmnkb(iqpt,ishell,mband,ii,2)
                   bl1(2)=bl1(2)+lambda(qneigh(iqpt,ishell),iwann,ii,1)*mmnkb(iqpt,ishell,mband,ii,2)+&
&                   lambda(qneigh(iqpt,ishell),iwann,ii,2)*mmnkb(iqpt,ishell,mband,ii,1)

                   bl2(1)=bl2(1)+lambda(qneigh(iqpt,ishell),iwann,ii,1)*mmnkb(iqpt,ishell,nband,ii,1)-&
&                   lambda(qneigh(iqpt,ishell),iwann,ii,2)*mmnkb(iqpt,ishell,nband,ii,2)
                   bl2(2)=bl2(2)-lambda(qneigh(iqpt,ishell),iwann,ii,1)*mmnkb(iqpt,ishell,nband,ii,2)-&
&                   lambda(qneigh(iqpt,ishell),iwann,ii,2)*mmnkb(iqpt,ishell,nband,ii,1)

!                  write(std_out,'(a,4f10.6)') 'after, bl1,bl2',bl1(1),bl1(2),bl2(1),bl2(2)
                   Zmn(iqpt,mf,nf,1)=Zmn(iqpt,mf,nf,1)+(1/nshell)*(bl1(1)*bl2(1)-bl1(2)*bl2(2))
                   Zmn(iqpt,mf,nf,2)=Zmn(iqpt,mf,nf,2)+(1/nshell)*(bl1(2)*bl2(1)+bl1(1)*bl2(2))

                 end do                                  !j bands at Q+b
               end do                         !(Q+)b points
             end do                          !Wannier states
             Zmn(iqpt,nband,mband,1)=Zmn(iqpt,mband,nband,1)
             Zmn(iqpt,nband,mband,2)=-Zmn(iqpt,mband,nband,2)
!            write(std_out,'(a,2i5,a,i3,a,2i5,3f12.8)') 'Zmn',mf,nf,' qpt',iqpt,' bands',mband,nband,Zmn(iqpt,mf,nf,1),Zmn(iqpt,mf,nf,2),&
!            &Zmn(iqpt,mf,nf,1)*Zmn(iqpt,mf,nf,1)+Zmn(iqpt,mf,nf,2)*Zmn(iqpt,mf,nf,2)
           end if
         end do                           !(n)bands at Q
       end if
     end do                            !(m)bands at Q

   end if                             !end if nwnz<nwnn

 end do                             !Q points

!DEBUG
!write(std_out,*) 'zmnbld : exit'
!ENDDEBUG

end subroutine zmnbld
!!***
