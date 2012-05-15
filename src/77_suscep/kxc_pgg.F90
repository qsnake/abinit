!{\src2tex{textfont=tt}}
!!****f* ABINIT/kxc_pgg
!! NAME
!! kxc_pgg
!!
!! FUNCTION
!! Compute the PGG-exchange kernel in reciprocal space
!! (Phys. Rev. Lett. 76, 1212 (1996)).
!!
!! COPYRIGHT
!! Copyright (C) 2001-2012 ABINIT group (MF, XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gmet=reciprocal space metrix (bohr**-2)
!!  npw=number of plane waves
!!  rcut_coulomb=real space cutoff radius for Coulomb interaction (bohr)
!!  susmat(2,npw,npw)=density weighted squared density matrix in reciprocal space
!!  ucvol=unit cell volume (bohr**3)
!!
!! OUTPUT
!!  khxcg(2,npwdiel,nspden,npwdiel,nspden)=PGG-exhange kernel in G space, at
!!       full interaction strength
!!
!! NOTES
!! The density weighted squared density matrix (actually the reduced=summed-over-spin
!! density matrix) is convolved with the spherically cutoff Coulomb interaction.
!!
!! WARNINGS
!! a - 'rcut_coulomb' should be chosen consistently with cutoffs elsewhere,
!!     for instance dieltcel8.f
!! b - applicable for spin-polarized case as well, through input 'susmat',
!!     but this has not been checked
!!
!! TODO
!! If simply the squared density matrix is input through 'susmat' the
!! exchange energy is recovered as the zero-G component of the resulting 'khxcg'
!! (then not the kernel of course). This could help to check convergence
!! with respect to 'npw'. See +ex_pgg comment.
!!
!! PARENTS
!!      xcacfd
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine kxc_pgg(gmet,kg,khxcg,npw,rcut_coulomb,susmat,ucvol)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kxc_pgg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw
 real(dp),intent(in) :: rcut_coulomb,ucvol
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gmet(3,3),susmat(2,npw,npw)
 real(dp),intent(out) :: khxcg(2,npw,npw)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,ig,ii,ikg11,ikg12,ikg13,ikg21,ikg22,ikg23,ipw1,ipw2
 integer :: j1,j2,j3,jg,jj
 real(dp),parameter :: diffgsq=1.d-2
 real(dp) :: gred1,gred2,gred3,gsquar,tpisq
!arrays
 integer :: kgmax(3)
 integer,allocatable :: index_g_inv(:,:,:),jgarr(:)
 real(dp),allocatable :: gsq(:),sumg(:),vcoul(:)

! *************************************************************************

!DEBUG
!write(std_out,*) '%kxc_pgg: enter'
!write(std_out,*) 'npw=',npw
!ENDDEBUG

!tpisq is (2 Pi) **2:
 tpisq=(two_pi)**2

 kgmax(:)=0
 do ipw1=1,npw
   do jj=1,3
     kgmax(jj)=max( kg(jj,ipw1), kgmax(jj) )
   end do
 end do

!DEBUG
!write(std_out,*) 'kgmax:',kgmax(1:3)
!ENDDEBUG

!Perform allocations
 ABI_ALLOCATE(index_g_inv,(-2*kgmax(1):2*kgmax(1),-2*kgmax(2):2*kgmax(2),-2*kgmax(3):2*kgmax(3)))
 ABI_ALLOCATE(jgarr,(npw))
 ABI_ALLOCATE(gsq,(npw))
 ABI_ALLOCATE(sumg,(2))
 ABI_ALLOCATE(vcoul,(npw))

!DEBUG
!write(std_out,*) '%kxc_pg: creating plane wave index and coulomb potential'
!ENDDEBUG
 index_g_inv(:,:,:)=0
 do ipw1=1,npw
   index_g_inv(kg(1,ipw1),kg(2,ipw1),kg(3,ipw1))=ipw1

!  DEBUG
!  write(std_out,'(i5,2x,3i3,2x,i4)') ipw1,kg(1,ipw1),kg(2,ipw1),kg(3,ipw1)
!  ENDDEBUG

   gred1=dble(kg(1,ipw1))
   gred2=dble(kg(2,ipw1))
   gred3=dble(kg(3,ipw1))
   gsquar=tpisq*(gmet(1,1)*gred1**2+gmet(2,2)*gred2**2+gmet(3,3)*gred3**2 &
&   +2.0_dp*( (gmet(1,2)*gred2+gmet(1,3)*gred3)* gred1 +      &
&   gmet(2,3)*gred2*gred3) )
!  Distinguish G=0 from other elements
   if(gsquar > 1.0d-12)then
     vcoul(ipw1)=four_pi/gsquar*(1._dp-cos(sqrt(gsquar)*rcut_coulomb))
   else
     vcoul(ipw1)=four_pi*0.5_dp*rcut_coulomb**2
   end if

 end do

!DEBUG
!write(std_out,*) '%kxc_pg: starting convolution integral'
!ENDDEBUG

!loop over G1,G2 components of the density matrix
 do ipw2=1,npw
   ikg21=kg(1,ipw2)
   ikg22=kg(2,ipw2)
   ikg23=kg(3,ipw2)

   do ii=1,npw
     j1=ikg21-kg(1,ii)
     j2=ikg22-kg(2,ii)
     j3=ikg23-kg(3,ii)
     jgarr(ii)=index_g_inv(j1,j2,j3)
   end do

   do ipw1=1,ipw2
     ikg11=kg(1,ipw1)
     ikg12=kg(2,ipw1)
     ikg13=kg(3,ipw1)

!    do the convolution integral
     sumg(:)=0._dp
     do ii=1,npw

       if( jgarr(ii) /= 0 ) then

         i1=ikg11-kg(1,ii)
         i2=ikg12-kg(2,ii)
         i3=ikg13-kg(3,ii)

!        j1=ikg21-kg(1,ii)
!        j2=ikg22-kg(2,ii)
!        j3=ikg23-kg(3,ii)

         ig=index_g_inv(i1,i2,i3)
!        jg=index_g_inv(j1,j2,j3)

         if( ig /= 0 ) then

           jg=jgarr(ii)

!          DEBUG
!          write(std_out,'(i5,2x,3i3,1x,3i3,2x,2i4)') ii,i1,i2,i3,&
!          &                                             kg(1,jg),kg(2,jg),kg(3,jg),&
!          &                                             ig,jg
!          ENDDEBUG

           sumg(1)=sumg(1)+susmat(1,ig,jg)*vcoul(ii)
           sumg(2)=sumg(2)+susmat(2,ig,jg)*vcoul(ii)

         end if

       end if

     end do
     khxcg(:,ipw1,ipw2)=-sumg(:)*ucvol

!    DEBUG
!    if(ipw1==ipw2) write(std_out,'(2i4,2(1x,es14.6))') ipw1,ipw2,khxcg(1,ipw1,ipw1),vcoul(ipw1)
!    write(std_out,'(2i4,3(1x,es14.6))') ipw1,ipw2,khxcg(1:2,ipw1,ipw2),vcoul(ipw1)
!    ENDDEBUG

   end do
 end do

!DEBUG
!verify hermiticity, note: ipw1 loop above must end at npw
!write(std_out,*) '%kxc_pgg: check hermiticity of pgg kernel'
!do ipw2=1,npw,max(2,npw/10)
!do ipw1=ipw2,npw,max(2,npw/10)
!write(std_out,'(2i4,2(1x,es14.6))') ipw1,ipw2,&
!&   khxcg(1,ipw1,ipw2)-khxcg(1,ipw2,ipw1),&
!&   khxcg(2,ipw1,ipw2)+khxcg(2,ipw2,ipw1)
!end do
!end do
!ENDDEBUG

!Impose hermiticity
 write(std_out,*) '%kxc_pg: imposing hermiticity'
 do ipw2=1,npw
   do ipw1=ipw2+1,npw
     khxcg(1,ipw1,ipw2)= khxcg(1,ipw2,ipw1)
     khxcg(2,ipw1,ipw2)=-khxcg(2,ipw2,ipw1)
   end do
 end do

!DEBUG
!write(std_out,'(a10,2(1x,es20.12))') '+ex_pgg? ', 0.5_dp*khxcg(1,1,1)/ucvol
!ENDDEBUG

 ABI_DEALLOCATE(index_g_inv)
 ABI_DEALLOCATE(jgarr)
 ABI_DEALLOCATE(gsq)
 ABI_DEALLOCATE(sumg)
 ABI_DEALLOCATE(vcoul)

!DEBUG
!write(std_out,*) '%kxc_pgg: done'
!ENDDEBUG

end subroutine kxc_pgg
!!***
