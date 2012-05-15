!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfcv_init2
!! NAME
!! scfcv_init2
!!
!! FUNCTION
!! WARNING : Temporary wrapper to scfcv
!! Self-consistent-field convergence.
!! Conducts set of passes or overall iterations of preconditioned
!! conjugate gradient algorithm to converge wavefunctions to
!! ground state and optionally to compute forces and energy.
!! This routine is called to compute forces for given atomic
!! positions or else to do non-SCF band structures.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG, GMR, AR, MKV, MT, FJ, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mcg)=updated wavefunctions; if mkmem>=nkpt, these are kept in a disk file.
!!  cpus= cpu time limit in seconds
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtpawuj(ndtpawuj)= data used for the automatic determination of U
!!     (relevant only for PAW+U) calculations (see initberry.f)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs for the "coarse" grid (see NOTES below)
!!   | mkmem =number of k points which can fit in memory; set to 0
!!   |        if use disk
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in cell.
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  electronpositron <type(electronpositron_type)>=quantities for the
!!      electron-positron annihilation
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  iapp=indicates the eventual suffix to be appended to the generic
!!       output root
!!       if 0 : no suffix to be appended (called directly from gstate)
!!       if positive : append "_TIM//iapp" (called from move or brdmin)
!!       if -1 : append "_TIM0" (called from brdmin)
!!       if -2, -3, -4, -5: append "_TIMA", ... ,"_TIMD", (called from move)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialization of the gstate run is not yet
!!     finished
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible
!!     zone data
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  ndtpawuj=size of dtpawuj
!!  npwarr(nkpt)=number of planewaves in basis at this k point!
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!     for the "fine" grid (see NOTES below)
!!  occ(mband*nkpt*nsppol)=occupation number for each band (often 2)
!!     at each k point
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and
!!        related data
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic
!!     occupancies
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic
!!     translation phases
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points
!!     and spins
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!     forces and its components, the stress tensor) of a ground-state
!!     computation (should be made a pure output quantity)
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in el./bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  scf_history <type(scf_history_type)>=arrays obtained from previous
!!     SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic
!!     energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for
!!     each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real
!!     spherical harmonics
!!  wffnew,wffnow=struct info for wf disk files.
!!  wvl <type(wvl_data)>=all wavelets data.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)= at input, previous reduced dimensionless atomic
!!     coordinates at output, current xred is transferred to xred_old
!!
!! OUTPUT
!!  ab_scfcv <type (ab_scfcvargs)> = Datatype for most of the
!!     arguments in scfcv
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! NOTES
!!
!! Arguments not in init:
!!
!!        electronpositron,paw_dmft,rprimd,xred,xred_old,wffnew,wffnow
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scfcv_init2(scfcv_args,&
&  ab_scfcv_in,ab_scfcv_inout,&
&  dtset,paw_dmft,wffnew,wffnow)

 use m_profiling
! &  dtset,electronpositron,paw_dmft,wffnew,wffnow)

 use defs_basis
 use defs_datatypes
 use m_wffile
 use defs_abitypes
 use defs_scftypes
! use defs_wvltypes
! use defs_parameters
! use defs_rectypes
 use defs_scfcvargs
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
! use m_rec
! use m_io_tools, only : flush_unit
 use m_paw_dmft, only: paw_dmft_type

#if defined HAVE_TRIO_ETSF_IO
! use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scfcv_init2'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ab_scfcvargs),intent(out) :: scfcv_args
 type(ab_scfcv_args_in),intent(in),target :: ab_scfcv_in
 type(ab_scfcv_args_inout),intent(in),target :: ab_scfcv_inout
 type(dataset_type),intent(in),target :: dtset
! type(electronpositron_type),intent(in),target :: electronpositron
 type(paw_dmft_type), intent(in),target :: paw_dmft
 type(wffile_type),intent(in),target :: wffnew,wffnow
!arrays

!Local variables -------------------------
!scalars

! ****************************************************************

 scfcv_args%ab_scfcv_in=>ab_scfcv_in
 scfcv_args%ab_scfcv_inout=>ab_scfcv_inout
 scfcv_args%dtset=>dtset
!scfcv_args%electronpositron=>electronpositron
 scfcv_args%paw_dmft=>paw_dmft
 scfcv_args%wffnew=>wffnew
 scfcv_args%wffnow=>wffnow

end subroutine scfcv_init2
!!***
