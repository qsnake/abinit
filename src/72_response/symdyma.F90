!{\src2tex{textfont=tt}}
!!****f* ABINIT/symdyma
!!
!! NAME
!! symdyma
!!
!! FUNCTION
!! Symmetrize the dynamical matrices
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! indsym(4,msym*natom)=indirect indexing array : for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!! natom=number of atoms in unit cell
!! nsym=number of space group symmetries
!! qptn(3)=normalized phonon wavevector
!! rprimd(3,3)=dimensional primitive translations (bohr)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! Input/Output
!! dmati(2*3*natom*3*natom)=dynamical matrices relative to the q
!!  points of the B.Z. sampling
!!
!! NOTES
!! the procedure of the symmetrization of the dynamical matrix follows the
!! equations in: Hendrikse et al., Computer Phys. Comm. 86, 297 (1995)
!!
!! TODO
!! A full description of the equations should be included
!!
!! PARENTS
!!      phfrq3,relaxpol,scphon_new_frequencies
!!
!! CHILDREN
!!      mati3inv,matr3inv,symq3
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symdyma(dmati,indsym,natom,nsym,qptn,rprimd,symrel)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symdyma'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrel(3,3,nsym)
 real(dp),intent(in) :: qptn(3),rprimd(3,3)
 real(dp),intent(inout) :: dmati(2*3*natom*3*natom)

!Local variables -------------------------
!scalars
 integer :: i1,i2,iat,idir,ii,index,isign,isym,itirev,jat,jdir,jj,kk,ll
! integer :: bigloop ! commented out below
 integer :: niat,njat,timrev
 real(dp) :: arg1,arg2,dmint,im,re,sumi,sumr
!arrays
 integer :: indij(natom,natom),symq(4,2,nsym),symrec(3,3,nsym)
 real(dp) :: TqR(3,3),TqS_(3,3),dynmat(2,3,natom,3,natom)
 real(dp) :: dynmatint(2*nsym,2,3,natom,3,natom),gprimd(3,3)
 real(dp) :: symcart(3,3,nsym)

! *********************************************************************

!0) initializations

!DEBUG
!write(std_out,*) ' enter : symdyma DEBUG'
!write(std_out,*) ' natom =',natom
!write(std_out,*) ' nsym  =',nsym
!write(std_out,'(a,3f9.5)') ' qptn =',qptn
!do isym=1,nsym
!write(std_out,'(a,i4,2i6)') 'sym.no.',isym,symq(4,1,isym),symq(4,2,isym)
!end do
!ENDDEBUG

 call matr3inv(rprimd,gprimd)
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

 TqR=zero
 TqS_=zero
 dynmat=zero

!Note: dynmat is used as work space here
 i1=0
 do iat=1,natom
   do idir=1,3
     i1=i1+1
     i2=0
     do jat=1,natom
       do jdir=1,3
         i2=i2+1
         index=i1+3*natom*(i2-1)
         dynmat(1,idir,iat,jdir,jat)=dmati(2*index-1)
         dynmat(2,idir,iat,jdir,jat)=dmati(2*index  )
       end do
     end do
   end do
 end do

!Transform symrel to cartesian coordinates (RC coding)
!do isym=1,nsym
!symcart(:,:,isym)=matmul(rprimd,matmul(dble(symrel(:,:,isym)),gprimd))
!end do

!Coding from symdm9
 do isym=1,nsym
   do jj=1,3
     symcart(:,jj,isym)=zero
     do kk=1,3
       do ll=1,3
         symcart(:,jj,isym)=symcart(:,jj,isym)+rprimd(:,kk)*gprimd(jj,ll)*symrel(kk,ll,isym)
       end do
     end do
   end do
 end do


!Get the symq of the CURRENT Q POINT
!mjv: set prtvol=0 for production runs.
 call symq3(nsym,qptn,symq,symrec,timrev,prtvol=0)

!DEBUG
!do isym=1,nsym
!write(std_out,'(a,i4,2i6)') 'sym.no.',isym,symq(4,1,isym),symq(4,2,isym)
!write(std_out,'(a,3e16.6)') 'cartesian symmetry',symcart(1,:,isym)
!write(std_out,'(a,3e16.6)') '                  ',symcart(2,:,isym)
!write(std_out,'(a,3e16.6)') '                  ',symcart(3,:,isym)
!end do
!ENDDEBUG

!do bigloop=1,3


 indij(:,:)=0
 dynmatint=zero

 do isym=1,nsym                                 ! loop over all the symmetries
!  write(std_out,*) 'current symmetry',isym
   do itirev=1,2                                 ! loop over the time-reversal symmetry
     isign=3-2*itirev                             ! to take into accont the time-reversal
!    write(std_out,*) 'timereversal',isign
     if (symq(4,itirev,isym)==1) then             ! isym belongs to the wave vector point group
!      write(std_out,*) 'isym belongs to the wave vector point group'
       do iat=1,natom                              ! loop over the atoms
         do jat=1,natom                            ! loop over the atoms
           niat=indsym(4,isym,iat)                   ! niat={R|t}iat
           njat=indsym(4,isym,jat)                   ! njat={R|t}jat
           indij(niat,njat)=indij(niat,njat)+1
!          write(std_out,'(a,5i5)') 'current status:',iat,jat,niat,njat,indij(niat,njat)
!          phase calculation, arg1 and arg2 because of two-atom derivative
           arg1=two_pi*( qptn(1)*indsym(1,isym,iat)+&
&           qptn(2)*indsym(2,isym,iat)+&
&           qptn(3)*indsym(3,isym,iat) )
           arg2=two_pi*( qptn(1)*indsym(1,isym,jat)+&
&           qptn(2)*indsym(2,isym,jat)+&
&           qptn(3)*indsym(3,isym,jat) )
           re=cos(arg1)*cos(arg2)+sin(arg1)*sin(arg2)
           im=isign*(cos(arg2)*sin(arg1)-cos(arg1)*sin(arg2))

           do idir=1,3                               ! loop over displacements
             do jdir=1,3                              ! loop over displacements
!              we pick the (iat,jat) (3x3) block of the dyn.mat.
               sumr=zero
               sumi=zero
               do ii=1,3
                 do jj=1,3
                   sumr=sumr+symcart(idir,ii,isym)*dynmat(1,ii,niat,jj,njat)*symcart(jdir,jj,isym)
                   sumi=sumi+symcart(idir,ii,isym)*dynmat(2,ii,niat,jj,njat)*symcart(jdir,jj,isym)
                 end do
               end do
               sumi=isign*sumi

               dynmatint(nsym*(itirev-1)+isym,1,idir,iat,jdir,jat)=re*sumr-im*sumi
               dynmatint(nsym*(itirev-1)+isym,2,idir,iat,jdir,jat)=re*sumi+im*sumr

             end do
           end do
         end do
       end do                        ! end treatment of the (iat,jat) (3x3) block of dynmat
     end if                                        ! symmetry check
   end do                                         ! time-reversal
 end do                                          ! symmetries



!4) make the average, get the final symmetric dynamical matrix
 do iat=1,natom
   do jat=1,natom
!    write(std_out,*) 'symmetrizing iat,jat:',iat,jat
!    write(std_out,*) ' indij ',dble(indij(iat,jat))
!    write(std_out,'(a,3f20.16)') 'dynmat r',dynmat(1,1,iat,:,jat)
!    write(std_out,'(a,3f20.16)') 'dynmat r',dynmat(1,2,iat,:,jat)
!    write(std_out,'(a,3f20.16)') 'dynmat r',dynmat(1,3,iat,:,jat)
!    write(std_out,'(a,3f20.16)') 'dynmat i',dynmat(2,1,iat,:,jat)
!    write(std_out,'(a,3f20.16)') 'dynmat i',dynmat(2,2,iat,:,jat)
!    write(std_out,'(a,3f20.16)') 'dynmat i',dynmat(2,3,iat,:,jat)
     do idir=1,3
       do jdir=1,3
         dmint=zero
         do isym=1,2*nsym
           dmint=dmint+dynmatint(isym,1,idir,iat,jdir,jat)
         end do
         dynmat(1,idir,iat,jdir,jat)=dmint/dble(indij(iat,jat))
         dmint=zero
         do isym=1,2*nsym
           dmint=dmint+dynmatint(isym,2,idir,iat,jdir,jat)
         end do
         dynmat(2,idir,iat,jdir,jat)=dmint/dble(indij(iat,jat))
       end do
     end do
!    write(std_out,'(a,3f20.16)') 'dynmat r',dynmat(1,1,iat,:,jat)
!    write(std_out,'(a,3f20.16)') 'dynmat r',dynmat(1,2,iat,:,jat)
!    write(std_out,'(a,3f20.16)') 'dynmat r',dynmat(1,3,iat,:,jat)
!    write(std_out,'(a,3f20.16)') 'dynmat i',dynmat(2,1,iat,:,jat)
!    write(std_out,'(a,3f20.16)') 'dynmat i',dynmat(2,2,iat,:,jat)
!    write(std_out,'(a,3f20.16)') 'dynmat i',dynmat(2,3,iat,:,jat)
   end do
 end do

!end do

!do iat=1,natom
!do jat=1,natom
!write(std_out,*) 'symmetrizing iat,jat:',iat,jat
!write(std_out,'(a,3f20.16)') 'dynmat r',dynmat(1,1,iat,:,jat)
!write(std_out,'(a,3f20.16)') 'dynmat r',dynmat(1,2,iat,:,jat)
!write(std_out,'(a,3f20.16)') 'dynmat r',dynmat(1,3,iat,:,jat)
!write(std_out,'(a,3f20.16)') 'dynmat i',dynmat(2,1,iat,:,jat)
!write(std_out,'(a,3f20.16)') 'dynmat i',dynmat(2,2,iat,:,jat)
!write(std_out,'(a,3f20.16)') 'dynmat i',dynmat(2,3,iat,:,jat)
!end do
!end do

 i1=0
 do iat=1,natom
   do idir=1,3
     i1=i1+1
     i2=0
     do jat=1,natom
       do jdir=1,3
         i2=i2+1
         index=i1+3*natom*(i2-1)
         dmati(2*index-1)=dynmat(1,idir,iat,jdir,jat)
         dmati(2*index  )=dynmat(2,idir,iat,jdir,jat)
       end do
     end do
   end do
 end do

end subroutine symdyma
!!***
