!{\src2tex{textfont=tt}}
!!****f* ABINIT/phfrq3
!!
!! NAME
!! phfrq3
!!
!! FUNCTION
!! Get the phonon frequencies and eigenvectors (as well as the corresponding displacements)
!! If q is at Gamma, the non-analytical behaviour can be included.
!! Then, the effective dielectric tensor, the effective charges
!! and oscillator strengths for the limiting direction are also returned
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  amu(ntypat)=mass of the atoms (atomic mass unit)
!!  matrix (diagonal in the atoms)
!!  d2cart(2,3,mpert,3,mpert)=
!!   dynamical matrix, effective charges, dielectric tensor,....
!!   all in cartesian coordinates
!!  indsym(4,msym*natom)=indirect indexing array : for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!!  mpert =maximum number of ipert
!!  msym=maximum number of symmetries
!!  natom=number of atoms in unit cell
!!  nsym=number of space group symmetries
!!  ntypat=number of atom types
!!  qphnrm=(described above)
!!  qphon(3)= to be divided by qphnrm, give the phonon wavevector;
!!     if qphnrm==0.0_dp, then the wavevector is zero (Gamma point)
!!     and qphon gives the direction of
!!     the induced electric field; in the latter case, if qphon is
!!     zero, no non-analytical contribution is included.
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!  symdynmat=if 1, (re)symmetrize the dynamical matrix, except if Gamma wavevector with electric field added.
!!  symrel(3,3,nsym)=matrices of the group symmetries (real space)
!!  typat(natom)=integer label of each type of atom (1,2,...)
!!  ucvol=unit cell volume
!!
!! OUTPUT
!!  displ(2*3*natom*3*natom)= at the end, contains the displacements of atoms in cartesian coordinates.
!!    The first index means either the real or the imaginary part,
!!    The second index runs on the direction and the atoms displaced
!!    The third index runs on the modes.
!!  eigval(3*natom)=contains the eigenvalues of the dynamical matrix
!!  eigvec(2*3*natom*3*natom)= at the end, contains the eigenvectors of the dynamical matrix.
!!  phfrq(3*natom)=phonon frequencies (square root of the dynamical matrix eigenvalues,
!!    except if these are negative, and in this case, give minus the square root of the absolute value
!!    of the matrix eigenvalues). Hartree units.
!!
!! NOTES
!! one makes the dynamical matrix hermitian...
!! in case of q=Gamma, only the real part is used.
!!
!! PARENTS
!!      anaddb,inpphon,interpolate_phfrq,m_phdos,mkifc9,mkphbs,respfn,thm9
!!      thmeig
!!
!! CHILDREN
!!      fxphas,initmpi_seq,mkherm,symdyma,wrtout,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
& mpert,msym,natom,nsym,ntypat,phfrq,qphnrm,qphon,rprimd,&
& symdynmat,symrel,typat,ucvol)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phfrq3'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_66_wfs
 use interfaces_72_response, except_this_one => phfrq3
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,msym,natom,nsym,ntypat,symdynmat
 real(dp),intent(in) :: qphnrm,ucvol
!arrays
 integer,intent(in) :: indsym(4,msym*natom),symrel(3,3,nsym),typat(natom)
 real(dp),intent(in) :: amu(ntypat),d2cart(2,3,mpert,3,mpert),rprimd(3,3)
 real(dp),intent(inout) :: qphon(3)
 real(dp),intent(out) :: displ(2*3*natom*3*natom),eigval(3*natom)
 real(dp),intent(out) :: eigvec(2*3*natom*3*natom),phfrq(3*natom)

!Local variables -------------------------
!scalars
 integer :: analyt,enough=50,i1,i2,idir1,idir2,ier,ii,imode,index,ipert1,ipert2
 real(dp),parameter :: break_symm=1.0d-12
 real(dp) :: epsq,fac,norm,qphon2
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg_dum
!arrays
 real(dp) :: nearidentity(3,3),qptn(3)
 real(dp),allocatable :: dum(:,:),matrx(:,:),zeff(:,:),zhpev1(:,:),zhpev2(:)

! *********************************************************************

!DEBUG
!write(message, '(a)' )' phfrq3 : enter '
!call wrtout(std_out,message,'COLL')
!write(std_out,*)' qphon =',qphon
!ENDDEBUG

 if(enough==1)then
   write(message, '(a)' )' phfrq3 : suppress messages '
   call wrtout(std_out,message,'COLL')
 end if
 if(enough/=0)enough=enough-1

!Prepare the diagonalisation : analytical part.
!Note: displ is used as work space here
 i1=0
 do ipert1=1,natom
   do idir1=1,3
     i1=i1+1
     i2=0
     do ipert2=1,natom
       do idir2=1,3
         i2=i2+1
         index=i1+3*natom*(i2-1)
         displ(2*index-1)=d2cart(1,idir1,ipert1,idir2,ipert2)
         displ(2*index  )=d2cart(2,idir1,ipert1,idir2,ipert2)
       end do
     end do
   end do
 end do

!Determine the analyticity of the matrix.
 analyt=1
 if(abs(qphnrm)<tol8) analyt=0
 if(abs(qphon(1))<tol8.and.abs(qphon(2))<tol8.and.abs(qphon(3))<tol8) analyt=2

!In case of q=Gamma, only the real part is used
 if(analyt==0 .or. analyt==2)then
   do i1=1,3*natom
     do i2=1,3*natom
       index=i1+3*natom*(i2-1)
       displ(2*index)=0.0_dp
     end do
   end do
 end if

!In the case the non-analyticity is required :
 if(analyt==0)then

!  Normalize the limiting direction
   qphon2=qphon(1)**2+qphon(2)**2+qphon(3)**2
   qphon(:)=qphon(:)/sqrt(qphon2)

!  Get the dielectric constant for the limiting direction
   epsq=0.0_dp
   do idir1=1,3
     do idir2=1,3
       epsq=epsq+qphon(idir1)*qphon(idir2)*&
&       d2cart(1,idir1,natom+2,idir2,natom+2)
     end do
   end do

   ABI_ALLOCATE(zeff,(3,natom))

!  Get the effective charges for the limiting direction
   do idir1=1,3
     do ipert1=1,natom
       zeff(idir1,ipert1)=0.0_dp
       do idir2=1,3
         zeff(idir1,ipert1)=zeff(idir1,ipert1)+qphon(idir2)*&
&         d2cart(1,idir1,ipert1,idir2,natom+2)
       end do
     end do
   end do

!  Get the non-analytical part of the dynamical matrix,
!  and suppress its imaginary part.
   i1=0
   do ipert1=1,natom
     do idir1=1,3
       i1=i1+1
       i2=0
       do ipert2=1,natom
         do idir2=1,3
           i2=i2+1
           index=i1+3*natom*(i2-1)
           displ(2*index-1)=displ(2*index-1)+four_pi/ucvol*&
&           zeff(idir1,ipert1)*zeff(idir2,ipert2)/epsq
           displ(2*index  )=0.0_dp
         end do
       end do
     end do
   end do

   ABI_DEALLOCATE(zeff)

!  End of the non-analyticity treatment :
 end if

!This slight breaking of the symmetry allows the
!results to be more portable between machines
 nearidentity(:,:)=one
 nearidentity(1,1)=one+break_symm
 nearidentity(3,3)=one-break_symm

!Include the masses in the dynamical matrix
 do ipert1=1,natom
   do ipert2=1,natom
     fac=1.0_dp/sqrt(amu(typat(ipert1))*amu(typat(ipert2)))/amu_emass
     do idir1=1,3
       do idir2=1,3
         i1=idir1+(ipert1-1)*3
         i2=idir2+(ipert2-1)*3
         index=i1+3*natom*(i2-1)
         displ(2*index-1)=displ(2*index-1)*fac*nearidentity(idir1,idir2)
         displ(2*index  )=displ(2*index  )*fac*nearidentity(idir1,idir2)
!        This is to break slightly the translation invariance, and make
!        the automatic tests more portable
         if(ipert1==ipert2 .and. idir1==idir2)then
           displ(2*index-1)=displ(2*index-1)+break_symm*natom/amu_emass/idir1*0.01_dp
         end if
       end do
     end do
   end do
 end do

!Make the dynamical matrix hermitian
 call mkherm(displ,3*natom)

!DEBUG
!write(message, '(a)' )' phfrq3 : before making the dyn mat symmetric '
!call wrtout(std_out,message,'COLL')
!ENDDEBUG

!***********************************************************************
!Diagonalize the dynamical matrix

!Symmetrize the dynamical matrix
!FIXME: swap the next 2 lines and update test files to include symmetrization for Gamma point too (except in non-analytic case)
!if (symdynmat==1 .and. analyt > 0) then
 if (symdynmat==1 .and. analyt == 1) then
   qptn(:)=qphon(:)
   if (analyt==1) then
     qptn(:)=qphon(:)/qphnrm
   end if
   call symdyma(displ,indsym,natom,nsym,qptn,rprimd,symrel)
 end if

 ier=0
 ii=1
 ABI_ALLOCATE(matrx,(2,(3*natom*(3*natom+1))/2))
 do i2=1,3*natom
   do i1=1,i2
     matrx(1,ii)=displ(1+2*(i1-1)+2*(i2-1)*3*natom)
     matrx(2,ii)=displ(2+2*(i1-1)+2*(i2-1)*3*natom)
     ii=ii+1
   end do
 end do
 ABI_ALLOCATE(zhpev1,(2,2*3*natom-1))
 ABI_ALLOCATE(zhpev2,(3*3*natom-2))
 call ZHPEV ('V','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,zhpev2,ier)
 ABI_DEALLOCATE(matrx)
 ABI_DEALLOCATE(zhpev1)
 ABI_DEALLOCATE(zhpev2)

!***********************************************************************

!Get the phonon frequencies (negative by convention, if
!the eigenvalue of the dynamical matrix is negative)
 do imode=1,3*natom
   if(eigval(imode)>=1.0d-16)then
     phfrq(imode)=sqrt(eigval(imode))
   else if(eigval(imode)>=-1.0d-16)then
     phfrq(imode)=zero
   else
     phfrq(imode)=-sqrt(-eigval(imode))
   end if
 end do

!Fix the phase of the eigenvectors
 ABI_ALLOCATE(dum,(2,0))
 call initmpi_seq(mpi_enreg_dum)
 call fxphas(eigvec,dum,0,0,1,3*natom*3*natom,0,mpi_enreg_dum,3*natom,3*natom,0)
!call fxphas(eigvec,dum,0,0,1,2*3*natom*3*natom,0,mpi_enreg_dum,3*natom,3*natom,0)
 ABI_DEALLOCATE(dum)

!Normalise the eigenvectors
 do imode=1,3*natom
   norm=0.0_dp
   do idir1=1,3
     do ipert1=1,natom
       i1=idir1+(ipert1-1)*3
       index=i1+3*natom*(imode-1)
       norm=norm+eigvec(2*index-1)**2+eigvec(2*index)**2
     end do
   end do
   norm=sqrt(norm)
   do idir1=1,3
     do ipert1=1,natom
       i1=idir1+(ipert1-1)*3
       index=i1+3*natom*(imode-1)
       eigvec(2*index-1)=eigvec(2*index-1)/norm
       eigvec(2*index)=eigvec(2*index)/norm
     end do
   end do
 end do

!Get the phonon displacements
 do imode=1,3*natom
   do idir1=1,3
     do ipert1=1,natom
       i1=idir1+(ipert1-1)*3
       index=i1+3*natom*(imode-1)
       displ(2*index-1)=eigvec(2*index-1)&
&       /  sqrt(amu(typat(ipert1))*amu_emass)
       displ(2*index  )=eigvec(2*index  )&
&       /  sqrt(amu(typat(ipert1))*amu_emass)
     end do
   end do
 end do

end subroutine phfrq3
!!***
