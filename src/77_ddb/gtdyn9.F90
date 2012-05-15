!{\src2tex{textfont=tt}}
!!****f* ABINIT/gtdyn9
!!
!! NAME
!! gtdyn9
!!
!! FUNCTION
!! Generates a dynamical matrix from interatomic force
!! constants and long-range electrostatic interactions.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales by which rprim is to be multiplied
!! atmfrc(2,3,natom,3,natom,nrpt)
!!  = Interatomic Forces in real space
!!  (imaginary part only for debugging)
!! dielt(3,3) = dielectric tensor
!! dipdip= if 0, no dipole-dipole interaction was subtracted in atmfrc
!!  if 1, atmfrc has been build without dipole-dipole part
!! dyewq0(3,3,natom)= Ewald part of the dynamical matrix, at q=0.
!! gmet(3,3)= metric tensor in reciprocal space.
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! mpert =maximum number of ipert
!! natom= Number of atoms in the unit cell
!! nrpt= Number of R points in the Big Box
!! qphnrm= Normalisation coefficient for qpt
!! qpt(3)= Reduced coordinates of the q vectors in reciprocal space
!! rmet(3,3)= Metric tensor in real space.
!! rprim(3,3)= dimensionless primitive translations in real space
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3)!!)
!! trans(3,natom)= Atomic translations : xred = rcan + trans
!! ucvol= unit cell volume
!! wghatm(natom,natom,nrpt)
!!  = Weights associated to a pair of atoms and to a R vector
!! xred(3,natom)= relative coords of atoms in unit cell (dimensionless)
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement
!!
!! OUTPUT
!! d2cart(2,3,mpert,3,mpert)=
!!  dynamical matrix obtained for the wavevector qpt
!!  (normalized using qphnrm)
!!
!! PARENTS
!!      anaddb,inpphon,interpolate_phfrq,m_phdos,mkifc9,mkphbs,refineblk,thm9
!!
!! CHILDREN
!!      dymfz9,ewald9,ftifc_r2q,nanal9,q0dy3_apply,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine gtdyn9(acell,atmfrc,dielt,dipdip,&
& dyewq0,d2cart,gmet,gprim,mpert,natom,&
& nrpt,qphnrm,qpt,rmet,rprim,rpt,&
& trans,ucvol,wghatm,xred,zeff)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gtdyn9'
 use interfaces_14_hidewrite
 use interfaces_72_response
 use interfaces_77_ddb, except_this_one => gtdyn9
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: dipdip,mpert,natom,nrpt
 real(dp),intent(in) :: qphnrm,ucvol
!arrays
 real(dp),intent(in) :: acell(3),dielt(3,3),gmet(3,3),gprim(3,3),qpt(3)
 real(dp),intent(in) :: rmet(3,3),rprim(3,3),rpt(3,nrpt)
 real(dp),intent(in) :: trans(3,natom),wghatm(natom,natom,nrpt),xred(3,natom)
 real(dp),intent(in) :: zeff(3,3,natom)
 real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
 real(dp),intent(in) :: dyewq0(3,3,natom)
 real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer, save :: enough=50
 integer :: i1,i2,ib,iqpt,nqpt,nsize,option,plus
 integer :: sumg0
 character(len=500) :: message
!arrays
 real(dp) :: qphon(3)
 real(dp),allocatable :: dq(:,:,:,:,:),dyew(:,:,:,:,:)

! *********************************************************************

 if(enough==1)then
   write(message,'(a)' )' gtdyn9 : suppress messages '
   call wrtout(std_out,message,'COLL')
 end if
 if(enough/=0)enough=enough-1

!call timein(tcpui,twalli)

 ABI_ALLOCATE(dq,(2,3,natom,3,natom))

!Get the normalized wavevector
!write(std_out,*)' gtdyn9 : wavevector entirely blocked on zero'
 if(abs(qphnrm)<1.0d-7)then
   qphon(1:3)=0.0_dp
 else
   qphon(1:3)=qpt(1:3)/qphnrm
 end if

!Generate the analytical part from the interatomic forces
 nqpt=1

!Tests
 if(enough/=0)then
   write(message, '(a,3es15.5)' )&
&   '  gtdyn9 : enter ftiaf9 with q =',qphon(1:3)
   call wrtout(std_out,message,'COLL')
 end if

 call ftifc_r2q (atmfrc,dq,gprim,natom,nqpt,nrpt,rpt,qphon,wghatm)

!The analytical dynamical matrix dq has been generated
!in the normalized canonical coordinate system. Now, the
!phase is modified, in order to recover the usual (xred)
!coordinate of atoms.

 option=2
 nqpt=1
 call dymfz9(dq,natom,nqpt,gprim,option,qphon,trans)

 if(dipdip==1)then
!  Add the non-analytical part
   sumg0=0

!  Compute dyew(2,3,natom,3,natom)= Ewald part of the dynamical matrix,
!  second energy derivative wrt xred(3,natom) in Hartrees
!  (Denoted A-bar in the notes)
   ABI_ALLOCATE(dyew,(2,3,natom,3,natom))
   call ewald9(acell,dielt,dyew,gmet,gprim,natom,&
&   qphon,rmet,rprim,sumg0,ucvol,xred,zeff)
   call q0dy3_apply(natom,dyewq0,dyew)
   plus=1
   iqpt=1
   call nanal9(dyew,dq,iqpt,natom,nqpt,plus)
   ABI_DEALLOCATE(dyew)
 end if

!Copy the dynamical matrix in the proper location

!First zero all the elements
 nsize=2*(3*mpert)**2
 d2cart(:,:,:,:,:)=0.0_dp

!Copy the elements from dq to d2cart
 d2cart(:,:,1:natom,:,1:natom)=dq(:,:,1:natom,:,1:natom)

!In case we have the gamma point,
 if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-14)then

!  Copy the effective charge and dielectric constant
!  in the final array
   do i1=1,3
     do i2=1,3
       d2cart(1,i1,natom+2,i2,natom+2)=dielt(i1,i2)
       do ib=1,natom
         d2cart(1,i1,natom+2,i2,ib)=zeff(i1,i2,ib)
         d2cart(1,i2,ib,i1,natom+2)=zeff(i1,i2,ib)
       end do
     end do
   end do

   if(enough/=0)then
     write(message, '(a)' )' gtdyn9 : finished '
     call wrtout(std_out,message,'COLL')
   end if

 end if

 ABI_DEALLOCATE(dq)

end subroutine gtdyn9
!!***
