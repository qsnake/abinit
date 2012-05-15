!{\src2tex{textfont=tt}}
!!****f* ABINIT/corrmetalwf1
!!
!! NAME
!! corrmetalwf1
!!
!! FUNCTION
!! Response function calculation only:
!! Correct 1st-order wave-function, taking into account "metallic" occupations.
!! 1st-order WF orthogonal to C_n,k+q, restore the "active space" content of the first-order WF.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcgq)=planewave coefficients of wavefunctions at k+q
!!  cprjq(natom,mcprjq)= wave functions at k+q projected with non-local projectors
!!  cwavef(2,npw1*nspinor)= 1st-order wave-function before correction
!!  cwaveprj(natom,nspinor)= 1st-order wave-function before correction
!!                           projected on NL projectors (PAW)
!!  eig1(2*nband**2)=first-order eigenvalues (hartree)
!!  fermie1=derivative of fermi energy wrt (strain) perturbation
!!  ghc(2,npw1*nspinor)=<G|H0-eig0_k.I|C1 band,k> (NCPP) or <G|H0-eig0_k.S0|C1 band,k> (PAW)
!!                      (C1 before correction)
!!  iband=index of current band
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  icgq=shift to be applied on the location of data in the array cgq
!!  istwf_k=option parameter that describes the storage of wfs
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands
!!  npw1=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  occ(nband)=occupation number for each band for each k.
!!  rocceig(nband,nband)= (occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n)),
!!    if this ratio has been attributed to the band n (second argument), zero otherwise
!!  timcount=index used to accumulate timing (0 from vtowfk3, 1 from nstwf3)
!!  usepaw=flag for PAW
!!
!! OUTPUT
!!  cwave1(2,npw1*nspinor)= 1st-order wave-function after correction
!!  cwaveprj1(natom,nspinor)= 1st-order wave-function after correction
!!                            projected on NL projectors (PAW)
!!  edocc(nband)=correction to 2nd-order total energy coming from changes of occupations
!!  wf_corrected=flag put to 1 if input cwave1 is effectively different from output cwavef
!!
!! NOTES
!!  Was part of vtowfk3 before.
!!
!! PARENTS
!!      vtowfk3
!!
!! CHILDREN
!!      cprj_copy,cprj_zaxpby,dotprod_g,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine corrmetalwf1(cgq,cprjq,cwavef,cwave1,cwaveprj,cwaveprj1,edocc,eig1,fermie1,ghc,iband, &
&          ibgq,icgq,istwf_k,mcgq,mcprjq,mpi_enreg,natom,nband,npw1,nspinor,occ,rocceig,timcount,&
&          usepaw,wf_corrected)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'corrmetalwf1'
 use interfaces_18_timing
 use interfaces_44_abitypes_defs
 use interfaces_53_spacepar
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,ibgq,icgq,istwf_k,mcgq,mcprjq,natom,nband,npw1,nspinor,timcount,usepaw
 integer,intent(out) :: wf_corrected
 real(dp),intent(in) :: fermie1
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: cgq(2,mcgq),cwavef(2,npw1*nspinor)
 real(dp),intent(in) :: eig1(2*nband**2),ghc(2,npw1*nspinor),occ(nband),rocceig(nband,nband)
 real(dp),intent(out) :: cwave1(2,npw1*nspinor),edocc(nband)
 type(cprj_type),intent(in) :: cprjq(natom,mcprjq),cwaveprj(natom,nspinor*usepaw)
 type(cprj_type),intent(out) :: cwaveprj1(natom,nspinor*usepaw)

!Local variables-------------------------------
!scalars
 integer :: ibandkq,index_cgq,index_cprjq,index_eig1,ii
 real(dp) :: facti,factr,invocc
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwcorr(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 call timab(214+timcount,1,tsec)

!At this stage, the 1st order function cwavef is orthogonal to cgq (unlike when it is input to cgwf3).
!Here, restore the "active space" content of the 1st-order wavefunction, to give cwave1 .

!First copy input WF into output WF
 wf_corrected=0
!$OMP PARALLEL DO PRIVATE(ii) SHARED(cwave1,cwavef,npw1,nspinor)
 do ii=1,npw1*nspinor
   cwave1(1,ii)=cwavef(1,ii)
   cwave1(2,ii)=cwavef(2,ii)
 end do
!$OMP END PARALLEL DO
 if (usepaw==1) then
   call cprj_copy(cwaveprj,cwaveprj1)
 end if

!Correct WF only for occupied states
 if (abs(occ(iband)) > tol8) then
   invocc=one/occ(iband)

   edocc(iband)=zero

!  Loop over WF at k+q subspace
   do ibandkq=1,nband

!    Select bands with variable occupation
     if (abs(rocceig(ibandkq,iband))>tol8) then

       wf_corrected=1

       index_eig1=2*ibandkq-1+(iband-1)*2*nband
       index_cgq=npw1*nspinor*(ibandkq-1)+icgq

       if(ibandkq==iband) then
         factr=rocceig(ibandkq,iband)*invocc*(eig1(index_eig1)-fermie1)
       else
         factr=rocceig(ibandkq,iband)*invocc*eig1(index_eig1)
       end if
       facti= rocceig(ibandkq,iband)*invocc*eig1(index_eig1+1)

!      Apply correction to 1st-order WF
!      $OMP PARALLEL DO PRIVATE(ii) SHARED(cgq,cwave1,cwavef,cwcorr,facti,factr,index_cgq,npw1,nspinor)
       do ii=1,npw1*nspinor
         cwave1(1,ii)=cwave1(1,ii)+(factr*cgq(1,ii+index_cgq)-facti*cgq(2,ii+index_cgq))
         cwave1(2,ii)=cwave1(2,ii)+(facti*cgq(1,ii+index_cgq)+factr*cgq(2,ii+index_cgq))
       end do
!      $OMP END PARALLEL DO

!      In the PAW case, also apply correction to projected WF
       if (usepaw==1) then
         index_cprjq=nspinor*(ibandkq-1)+ibgq
         call cprj_zaxpby((/factr,facti/),(/one,zero/),cprjq(:,index_cprjq+1:index_cprjq+nspinor),cwaveprj1)
       end if

!      The factor of two is needed because we compute the 2DTE, and not E(2)
       edocc(iband)=edocc(iband)-two*(factr*eig1(index_eig1)+facti*eig1(index_eig1+1))

     end if ! Variable occupations
   end do ! Loop over k+q subspace
 end if ! occupied states

!In the PAW case, compute <Psi^(1)_ortho|H-Eig0_k.S|Psi^(1)_parallel> contribution to 2DTE
 if (usepaw==1.and.wf_corrected==1) then
   ABI_ALLOCATE(cwcorr,(2,npw1*nspinor))
   cwcorr=cwave1-cwavef
   call dotprod_g(factr,facti,istwf_k,mpi_enreg,npw1*nspinor,1,cwcorr,ghc)
   edocc(iband)=edocc(iband)+four*factr
   ABI_DEALLOCATE(cwcorr)
 end if

 call timab(213+timcount,2,tsec)

 DBG_EXIT("COLL")

end subroutine corrmetalwf1
!!***
