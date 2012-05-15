!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmkrho
!! NAME
!! pawmkrho
!!
!! FUNCTION
!! PAW only:
!! Build total pseudo (compensated) density (\tild_rho + \hat_rho)
!! Build compensation charge density (\hat_rho)
!! Build occupation matrix (packed storage)
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!         1 for GS calculations
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  ipert=index of perturbation if pawrhoij is a pertubed rhoij
!!        no meaning for ground-state calculations (should be 0)
!!  idir=direction of atomic displacement (in case of atomic displ. perturb.)
!!  lmnmax=maximum number of PAW radial wavefunctions
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nspden=number of spin-density components
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  paral_kgb=option for (kpt,g vectors,bands) parallelism
!!  pawang <type(pawang_type)>=angular mesh discretization and related data
!!  pawang_sym <type(pawang_type)>=angular data used for symmetrization
!!                                 optional parameter only needed for RF calculations
!!  pawfgr <type(paw_fgr_type)>=fine rectangular grid parameters
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrhoij0(natom) <type(pawrhoij_type)>= GS paw rhoij occupancies and related data (used only if ipert>0)
!!                                          optional parameter only needed for RF calculations
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  qphon(3)=wavevector of the phonon (RF only)
!!  rhopsg(2,pawfgr%nfftc)= pseudo density given on the coarse grid in reciprocal space
!!  rhopsr(pawfgr%nfftc,nspden)= pseudo density given on the coarse grid in real space
!!  rprimd(3,3)=real space primitive translations.
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!                   reciprocal space primitive translations
!!  typat(natom)=type for each atom
!!  ucvol=volume of the unit cell
!!  xred(3,natom)= reduced atomic coordinates
!!
!! OUTPUT
!!  compch_fft=compensation charge inside spheres integrated over fine fft grid
!!  pawnhat(pawfgr%nfft,nspden)=compensation charge density on fine rectangular grid (optional argument)
!!  rhog(2,pawfgr%nfft)= compensated pseudo density given on the fine grid in reciprocal space
!!                       This output is optional
!!  rhor(pawfgr%nfft,nspden)= compensated pseudo density given on the fine grid in real space
!!
!! SIDE EFFECTS
!!  pawrhoij(natom)= PAW occupancies
!!                   At input : in unpacked storage (pawrhoij()%rhoij_)
!!                   At output: in   packed storage (pawrhoij()%rhoijp)
!!
!! PARENTS
!!      nstpaw3,vtorho,vtorho3
!!
!! CHILDREN
!!      fourdp,pawmknhat,rhoij_destroy_unpacked,symrhoij,timab,transgrid
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawmkrho(compch_fft,cplex,gprimd,idir,indlmn,indsym,ipert,lmnmax,mpi_enreg,&
&          natom,nspden,nsym,ntypat,paral_kgb,pawang,pawfgr,pawfgrtab,pawprtvol,pawrhoij,&
&          pawtab,qphon,rhopsg,rhopsr,rhor,rprimd,symafm,symrec,typat,ucvol,xred,&
&          pawang_sym,pawnhat,pawrhoij0,rhog) ! optional arguments

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawmkrho'
 use interfaces_18_timing
 use interfaces_44_abitypes_defs
 use interfaces_53_ffts
 use interfaces_66_paw, except_this_one => pawmkrho
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,lmnmax,natom,nspden,nsym,ntypat,paral_kgb,pawprtvol
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: compch_fft
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(pawang_type),intent(in),optional :: pawang_sym
 type(pawfgr_type),intent(in) :: pawfgr
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),indsym(4,nsym,natom)
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: gprimd(3,3),qphon(3),rprimd(3,3),xred(3,natom)
 real(dp),intent(out),target,optional :: pawnhat(cplex*pawfgr%nfft,nspden)
 real(dp),intent(out) :: rhor(cplex*pawfgr%nfft,nspden)
 real(dp),intent(out),optional :: rhog(2,pawfgr%nfft)
 real(dp),intent(inout) :: rhopsg(2,pawfgr%nfftc),rhopsr(cplex*pawfgr%nfftc,nspden)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawrhoij_type),intent(inout),target :: pawrhoij(natom)
 type(pawrhoij_type),intent(in),target,optional :: pawrhoij0(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: choice,ider,izero,option
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: rhodum(:,:,:)
 real(dp),pointer :: pawnhat_ptr(:,:)
 type(pawrhoij_type),pointer :: pawrhoij0_ptr(:)

! ***********************************************************************

 DBG_ENTER("COLL")
 call timab(556,1,tsec)

!Compatibility tests
 if (pawrhoij(1)%use_rhoij_==0) then
   msg='  rhoij_ field must be allocated in pawrhoij !'
   MSG_BUG(msg)
 end if
 if (ipert>0.and.(.not.present(pawrhoij0))) then
   msg='  pawrhoij0 must be present when ipert>0 !'
   MSG_BUG(msg)
 end if

!Symetrize PAW occupation matrix and store it in packed storage
 call timab(557,1,tsec)
 option=1;choice=1
 if (present(pawang_sym)) then
   call symrhoij(choice,gprimd,indlmn,indsym,ipert,lmnmax,natom,natom,nsym,ntypat,option,&
&   pawang_sym,pawprtvol,pawrhoij,rprimd,symafm,symrec,typat)
 else
   call symrhoij(choice,gprimd,indlmn,indsym,ipert,lmnmax,natom,natom,nsym,ntypat,option,&
&   pawang,pawprtvol,pawrhoij,rprimd,symafm,symrec,typat)
 end if
 call rhoij_destroy_unpacked(pawrhoij)
 call timab(557,2,tsec)

!Compute compensation charge density
 ider=0;izero=0
 if (present(pawnhat)) then
   pawnhat_ptr => pawnhat
 else
   ABI_ALLOCATE(pawnhat_ptr,(pawfgr%nfft,nspden))
 end if
 if (present(pawrhoij0)) then
   pawrhoij0_ptr => pawrhoij0
 else
   pawrhoij0_ptr => pawrhoij
 end if
 call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,mpi_enreg,natom,natom,&
& pawfgr%nfft,pawfgr%ngfft,ider,nspden,ntypat,paral_kgb,pawang,pawfgrtab,&
& rhodum,pawnhat_ptr,pawrhoij,pawrhoij0_ptr,pawtab,qphon,rprimd,ucvol,xred)
!Transfer pseudo density from coarse grid to fine grid
 call transgrid(cplex,mpi_enreg,nspden,+1,1,0,paral_kgb,pawfgr,rhopsg,rhodum,rhopsr,rhor)

!Add pseudo density and compensation charge density (on fine grid)
 rhor(:,:)=rhor(:,:)+pawnhat_ptr(:,:)

 if (.not.present(pawnhat))  then
   ABI_DEALLOCATE(pawnhat_ptr)
 end if
 nullify(pawnhat_ptr)

!Compute compensated pseudo density in reciprocal space
 if (present(rhog)) then
   call fourdp(cplex,rhog,rhor(:,1),-1,mpi_enreg,pawfgr%nfft,pawfgr%ngfft,paral_kgb,0)
 end if

 call timab(556,2,tsec)
 DBG_EXIT("COLL")

end subroutine pawmkrho
!!***
