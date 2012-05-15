!{\src2tex{textfont=tt}}
!!****f* ABINIT/mag_loc_k
!! NAME
!! mag_loc_k
!!
!! FUNCTION
!! This routine computes the on-site local magnetization.
!!
!! COPYRIGHT
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=inverse index table for atoms (see scfcv.f)
!!  cg(2,mcg)=input wavefunctions
!!  cprj(natom,nband_k)=cprj for all bands at current k point, ordered as in input file
!!  dimffnl = 2nd dimension of ffnl (1 + number of derivatives)
!!  ffnl(npw_k,dimffnl,lmnmax,ntypat) = nonlocal form factors
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  icg=shift to be applied on the location of data in the array cg
!!  ikpt = index of present kpoint, used for storing result in detefield
!!  indlmn(6,lmnmax,ntypat) = array giving l,m,n,lm,ln,s for i = lmn
!!  istwfk_k = option parameter governing storage of wavefunctions
!!  kg_k(3,npw_k) = integer coords of planewave in basis sphere
!!  kpg_k(npw_k,nkpg) = (k+G) components and related data
!!  kpt(3) = k point in terms of recip translations
!!  lmnmax = max number of (l,m,n) components over all atom types
!!  matblk = dimension of array ph3d
!!  mcg=second dimension of the cg array
!!  mgfft = max size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+max(spin*angular momentum) for nonlocal pseudopotentials
!!  natom = number of atoms in cell
!!  nattyp(ntypat) = number of atoms of each type
!!  nband_k = number of bands at this k point
!!  ngfft(18) = data about 3D FFT
!!  nkpg = second dimension size of kpg
!!  nloalg(5) = data concerning choice of nonlocal operator algorithm
!!  nspinor = number of spinors
!!  npw_k = number of plane waves for this k point
!!  ntypat = number of types of atoms
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  phkxred(2,natom) = phase factors exp(2\pi i kpt.xred)
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  ph3d(2,npw_k,matblk) = 3D structure factors for each atom and planewave
!!  ucvol = unit cell volume
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtefield <type(efield_type)> = variables related to Berry phase. The contribution to the local
!!   magnetization at the given k point is stored in dtefield%mag_local_k(idir,ikpt)
!!
!! NOTES
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      cprj_alloc,cprj_free,cprj_get,nonlop
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mag_loc_k(atindx,atindx1,cg,cprj,dimffnl,dtefield,ffnl,gmet,gprimd,icg,&
&                       ikpt,indlmn,istwfk_k,kg_k,kpg_k,kpt,lmnmax,matblk,mcg,mgfft,mpi_enreg,&
&                       mpsang,mpssoang,natom,nattyp,nband_k,ngfft,nkpg,nloalg,npw_k,&
&                       nspinor,ntypat,pawtab,phkxred,ph1d,ph3d,ucvol)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_efield

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mag_loc_k'
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,dimffnl,ikpt,istwfk_k,lmnmax,matblk,mcg,mgfft
 integer,intent(in) :: mpsang,mpssoang,natom,nband_k,nkpg,npw_k,nspinor
 integer,intent(in) :: ntypat
 real(dp),intent(in) :: ucvol
 type(efield_type),intent(inout) :: dtefield
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),indlmn(6,lmnmax,ntypat),kg_k(3,npw_k)
 integer,intent(in) :: ngfft(18),nattyp(ntypat),nloalg(5)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kpg_k(npw_k,nkpg),kpt(3)
 real(dp),intent(in) :: phkxred(2,natom),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(inout) :: ph3d(2,npw_k,matblk)
 type(cprj_type),intent(in) :: cprj(natom,nband_k)
 type(pawtab_type),intent(in)  :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,dimenl1,dimenl2,iat,iatom,iband,idir,itypat
 integer :: klmn,nnlout,nspinortot,only_SO,paw_opt,signs
 integer :: tim_nonlop,useylm
 real(dp) :: lambda

!arrays
 integer,allocatable :: dimlmn_srt(:)
 real(dp) :: cwavef(2,npw_k),enlout(1),sij(1,1),svectout(2,1)
 real(dp),allocatable :: Lij(:,:,:)
 type(cprj_type),allocatable :: cprj_band_srt(:,:)
! ************************************************************************

!======================================================================
!compute magnetization due to on-site local velocity at current k point
!======================================================================

 choice = 1 ! apply non-local potential Lij only
 ABI_ALLOCATE(cprj_band_srt,(natom,1))
 ABI_ALLOCATE(dimlmn_srt,(natom))
 iatom = 0
 do itypat = 1, ntypat
   do iat = 1, nattyp(itypat)
     iatom = iatom + 1
     dimlmn_srt(iatom) = dtefield%lmn_size(itypat)
   end do
 end do
 call cprj_alloc(cprj_band_srt,0,dimlmn_srt)

 cpopt = 2 ! use cprj in memory
 dimenl2 = natom
 lambda = zero
 nnlout = 1
 only_SO = 0
 paw_opt = 1 ! use Lij only in nonlop
 signs = 1 ! get contracted output <u|Lij|u>
 tim_nonlop = 0 ! no timing call, maybe should change this
 useylm = 1
 nspinortot=min(2,(1+mpi_enreg%paral_spin)*nspinor)

 dimenl1 = 2*lmnmax*(lmnmax+1)/2
 ABI_ALLOCATE(Lij,(dimenl1,dimenl2,1))
 Lij(:,:,:) = zero


 dtefield%mag_local_k(:,ikpt) = zero
 do idir = 1, 3

   iatom = 0
   do itypat = 1, ntypat
     do iat = 1, nattyp(itypat)
       iatom = iatom + 1
       do klmn = 1, pawtab(itypat)%lmn2_size
         Lij(2*klmn-1,iatom,1) = dtefield%Lij(1,klmn,itypat,idir)
         Lij(2*klmn,iatom,1) = dtefield%Lij(2,klmn,itypat,idir)
       end do ! end loop over lmn2_size
     end do ! end loop over iat
   end do ! end loop over ntypat

   do iband = 1, nband_k
     cwavef(:,:)=cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
     call cprj_get(atindx,cprj_band_srt,cprj,natom,iband,0,ikpt,1,1,&
&     nband_k,1,mpi_enreg,natom,1,nband_k,1,1,0)
     call nonlop(atindx1,choice,cpopt,cprj_band_srt,dimenl1,dimenl2,dimffnl,dimffnl,Lij,enlout,&
&     ffnl,ffnl,gmet,gprimd,idir,indlmn,istwfk_k,kg_k,kg_k,kpg_k,kpg_k,kpt,kpt,&
&     lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,nattyp,ngfft,&
&     nkpg,nkpg,nloalg,nnlout,npw_k,npw_k,nspinor,nspinortot,ntypat,only_SO,paw_opt,&
&     phkxred,phkxred,ph1d,ph3d,ph3d,signs,sij,svectout,tim_nonlop,ucvol,useylm,&
&     cwavef,cwavef)
     dtefield%mag_local_k(idir,ikpt) = dtefield%mag_local_k(idir,ikpt) + enlout(1)

   end do ! end loop over bands

 end do ! end loop over idir

 call cprj_free(cprj_band_srt)
 ABI_DEALLOCATE(Lij)
 ABI_DEALLOCATE(cprj_band_srt)
 ABI_DEALLOCATE(dimlmn_srt)

end subroutine mag_loc_k
!!***
