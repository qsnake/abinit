!{\src2tex{textfont=tt}}
!!****f* ABINIT/getdc1
!!
!! NAME
!! getdc1
!!
!! FUNCTION
!! Compute |delta_C^(1)> from one wave function C - PAW ONLY
!! Compute <G|delta_C^(1)> and eventually <P_i| delta_C^(1)> (P_i= non-local projector)
!! delta_C^(1) is the variation of wavefunction only due to variation of overlap operator S.
!! delta_C^(1)=-1/2.Sum_j [ <C_j|S^(1)|C>.C_j
!!         see PRB 78, 035105 (2008), Eq. (42)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cgq(2,mcgq)=wavefunction coefficients for ALL bands at k+Q
!!  cprjq(natom,mcprjq)= wave functions at k+q projected with non-local projectors: cprjq=<P_i|Cnk+q>
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  icgq=shift to be applied on the location of data in the array cgq
!!  istwfk=option parameter that describes the storage of wfs
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mpi_enreg=informations about MPI parallelization
!!  natom= number of atoms in cell
!!  nband=number of bands
!!  npw1=number of planewaves in basis sphere at k+Q
!!  nspinor=number of spinorial components of the wavefunctions
!!  opt_cprj=flag governing the computation of <P_i|delta_C^(1)> (P_i= non-local projector)
!!  ortalg=governs the choice of the algorithm for orthogonalisation.
!!  s1cwave0(2,npw1*nspinor)=<G|S^(1)|C> where S^(1) is the first-order overlap operator
!!
!! OUTPUT
!!  dcwavef(2,npw1*nspinor)=change of wavefunction due to change of overlap PROJECTED ON PLANE-WAVES:
!!         dcwavef is delta_C(1)=-1/2.Sum_{j}[<C0_k+q_j|S(1)|C0_k_i>.|C0_k+q_j>]
!!  === if optcprj=1 ===
!!  dcwaveprj(natom,nspinor*optcprj)=change of wavefunction due to change of overlap PROJECTED ON NL-PROJECTORS:
!!
!! PARENTS
!!      cgwf3,nstpaw3
!!
!! CHILDREN
!!      cprj_axpby,cprj_lincom,projbd
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getdc1(cgq,cprjq,dcwavef,dcwaveprj,ibgq,icgq,istwfk,mcgq,mcprjq,&
&                 mpi_enreg,natom,nband,npw1,nspinor,optcprj,ortalg,s1cwave0)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getdc1'
 use interfaces_44_abitypes_defs
 use interfaces_66_wfs, except_this_one => getdc1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ibgq,icgq,istwfk,mcgq,mcprjq,natom,nband,npw1,nspinor,optcprj,ortalg
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) ::cgq(2,mcgq),s1cwave0(2,npw1*nspinor)
 real(dp),intent(out) :: dcwavef(2,npw1*nspinor)
 type(cprj_type),intent(in) :: cprjq(natom,mcprjq)
 type(cprj_type),intent(out) :: dcwaveprj(natom,nspinor*optcprj)

!Local variables-------------------------------
!scalars
 integer, parameter :: iprint=0,tim_projbd=0
 integer :: ipw
 real(dp),parameter :: scal=-half
!arrays
 real(dp), allocatable :: dummy(:,:),scprod(:,:)
 type(cprj_type),allocatable :: tmpcprj(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!$OMP PARALLEL DO PRIVATE(ipw) &
!$OMP SHARED(dcwavef,s1cwave0,npw1,nspinor)
 do ipw=1,npw1*nspinor
   dcwavef(1:2,ipw)=s1cwave0(1:2,ipw)
 end do
!$OMP END PARALLEL DO
 ABI_ALLOCATE(dummy,(0,0))
 ABI_ALLOCATE(scprod,(2,nband))

!=== 1- COMPUTE: <G|S^(1)|C_k> - Sum_j [<C_k+q,j|S^(1)|C_k>.<G|C_k+q,j>]
!!               using the projb routine
!Note the subtlety: projbd is called with useoverlap=0 and s1cwave0
!in order to get Sum[<cgq|s1|c>|cgq>]=Sum[<cgq|gs1>|cgq>]
 call projbd(cgq,dcwavef,-1,icgq,0,istwfk,mcgq,mpi_enreg,0,nband,npw1,nspinor,&
& ortalg,iprint,dummy,scprod,0,tim_projbd,0)

!=== 2- COMPUTE: <G|delta_C^(1)> = -1/2.Sum_j [<C_k+q,j|S^(1)|C_k>.<G|C_k+q,j>]
!by substraction
!$OMP PARALLEL DO PRIVATE(ipw) &
!$OMP SHARED(dcwavef,s1cwave0,npw1,nspinor)
 do ipw=1,npw1*nspinor
   dcwavef(1:2,ipw)=scal*(s1cwave0(1:2,ipw)-dcwavef(1:2,ipw))
 end do
!$OMP END PARALLEL DO

!=== 3- COMPUTE: <P_i|delta_C^(1)> = -1/2.Sum_j [<C_k+q,j|S^(1)|C_k>.<P_i|C_k+q,j>]
 if (optcprj==1.and.mcprjq>0) then
   ABI_ALLOCATE(tmpcprj,(natom,nspinor))
   call cprj_lincom(scprod,cprjq(:,ibgq+1:ibgq+nspinor*nband),dcwaveprj,nband)
   call cprj_axpby(zero,scal,tmpcprj,dcwaveprj)
   ABI_DEALLOCATE(tmpcprj)
 end if

 ABI_DEALLOCATE(dummy)
 ABI_DEALLOCATE(scprod)

 DBG_EXIT("COLL")

end subroutine getdc1
!!***
