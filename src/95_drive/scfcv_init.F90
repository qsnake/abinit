!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfcv_init
!! NAME
!! scfcv_init
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
!!   | nspinor=number of spinorial components of the wavefunctions
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
!!  ab_scfcv <type (ab_scfcv_args)> = Datatype for most of the
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


subroutine scfcv_init(ab_scfcv_in,ab_scfcv_inout,atindx,atindx1,cg,cpus,&
&  dtefield,dtfil,dtpawuj,dtset,ecore,eigen,hdr,iapp,&
&  indsym,initialized,irrzon,kg,mcg,mpi_enreg,nattyp,ndtpawuj,&
&  nfftf,npwarr,occ,pawang,pawfgr,pawrad,pawrhoij,&
&  pawtab,phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,&
&  resid,results_gs,scf_history,fatvshift,&
&  symrec,taug,taur,wvl,ylm,ylmgr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_wffile
 use defs_abitypes
 use defs_scftypes
 use defs_wvltypes
 use defs_parameters
 use defs_rectypes
 use m_scf_history, only : scf_history_type
 use m_results_gs , only : results_gs_type
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 use m_rec
 use m_io_tools, only : flush_unit
 use m_paw_dmft, only: paw_dmft_type
 use m_efield
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scfcv_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),target :: iapp,mcg,ndtpawuj,pwind_alloc
 integer,intent(in),target :: initialized,nfftf
 real(dp),intent(in),target :: cpus,ecore
 real(dp),intent(in),target :: fatvshift
 type(ab_scfcv_args_in),intent(out) :: ab_scfcv_in
 type(ab_scfcv_args_inout),intent(out) :: ab_scfcv_inout
 type(MPI_type),intent(in),target :: mpi_enreg
 type(datafiles_type),intent(in),target :: dtfil
 type(dataset_type),intent(in),target :: dtset
 type(efield_type),intent(in),target :: dtefield
! type(electronpositron_type),pointer :: electronpositron
 type(hdr_type),intent(in),target :: hdr
 type(pawang_type),intent(in),target :: pawang
 type(pawfgr_type),intent(in),target :: pawfgr
 type(pseudopotential_type),intent(in),target :: psps
 type(recursion_type),intent(in),target :: rec_set
 type(results_gs_type),intent(in),target :: results_gs
 type(scf_history_type),intent(in),target :: scf_history
! type(wffile_type),intent(in),target :: wffnew,wffnow
 type(wvl_data),intent(in),target :: wvl
!arrays
 integer,intent(in),target :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(in),target :: indsym(4,dtset%nsym,dtset%natom)
!no_abirules
 integer, intent(in),target :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer, intent(in),target :: kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in),target :: nattyp(psps%ntypat),npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
 integer, intent(in),target :: symrec(3,3,dtset%nsym)
 real(dp), intent(in),target :: cg(2,mcg)
 real(dp), intent(in),target :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in),target :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in),target :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), intent(in),target :: pwnsfac(2,pwind_alloc)
! real(dp), intent(in),target :: rprimd(3,3)
! real(dp), pointer :: rhog(:,:),rhor(:,:)
 real(dp), pointer :: taug(:,:),taur(:,:)
 real(dp), intent(in),target :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
! real(dp), intent(in),target :: xred(3,dtset%natom),xred_old(3,dtset%natom)
 real(dp), intent(in),target :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in),target :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(macro_uj_type),intent(in),target :: dtpawuj(0:ndtpawuj)
 type(pawrhoij_type), intent(in),target :: pawrhoij(mpi_enreg%natom*psps%usepaw)
 type(pawrad_type), intent(in),target :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type), intent(in),target :: pawtab(psps%ntypat*psps%usepaw)
! type(paw_dmft_type), intent(in),target :: paw_dmft

!Local variables -------------------------
!scalars
 logical :: DEBUG=.FALSE.

! ****************************************************************

 if (DEBUG) then
   write(std_out,*) 'INTENT(IN) ARGUMENTS ON SCFCV'
   write(std_out,*) 'atindx=',ab_scfcv_in%atindx
   write(std_out,*) 'atindx1=',ab_scfcv_in%atindx1
   write(std_out,*) 'cpus=',ab_scfcv_in%cpus
   write(std_out,*) 'ecore=',ab_scfcv_in%ecore
   write(std_out,*) 'fatvshift=',ab_scfcv_in%fatvshift
   write(std_out,*) 'iapp=',ab_scfcv_in%iapp
   write(std_out,*) 'indsym=',ab_scfcv_in%indsym
   write(std_out,*) 'kg=',ab_scfcv_in%kg
   write(std_out,*) 'nattyp=',ab_scfcv_in%nattyp
   write(std_out,*) 'ndtpawuj=',ab_scfcv_in%ndtpawuj
   write(std_out,*) 'npwarr=',ab_scfcv_in%npwarr
   write(std_out,*) 'phnons=',ab_scfcv_in%phnons
   write(std_out,*) 'pwind=',ab_scfcv_in%pwind
   write(std_out,*) 'pwind_alloc=',ab_scfcv_in%pwind_alloc
   write(std_out,*) 'pwnsfac=',ab_scfcv_in%pwnsfac
   write(std_out,*) 'ylm=',ab_scfcv_in%ylm
   write(std_out,*) 'ylmgr=',ab_scfcv_in%ylmgr
!  write(std_out,*) 'pawang=',ab_scfcv_in%pawang
!  write(std_out,*) 'pawrad=',ab_scfcv_in%pawrad
!  write(std_out,*) 'pawtab=',ab_scfcv_in%pawtab
!  write(std_out,*) 'psps=',ab_scfcv_in%psps
 end if

 ab_scfcv_in%atindx=>atindx
 ab_scfcv_in%atindx1=>atindx1
 ab_scfcv_in%cpus=>cpus
 ab_scfcv_in%ecore=>ecore
 ab_scfcv_in%fatvshift=>fatvshift
 ab_scfcv_in%iapp=>iapp
 ab_scfcv_in%indsym=>indsym
 ab_scfcv_in%kg=>kg
 ab_scfcv_in%mcg=>mcg
 ab_scfcv_in%nattyp=>nattyp
 ab_scfcv_in%ndtpawuj=>ndtpawuj
 ab_scfcv_in%npwarr=>npwarr
 ab_scfcv_in%pawang=>pawang
 ab_scfcv_in%pawrad=>pawrad
 ab_scfcv_in%pawtab=>pawtab
 ab_scfcv_in%phnons=>phnons
 ab_scfcv_in%psps=>psps
 ab_scfcv_in%pwind=>pwind
 ab_scfcv_in%pwind_alloc=>pwind_alloc
 ab_scfcv_in%pwnsfac=>pwnsfac
 ab_scfcv_in%ylm=>ylm
 ab_scfcv_in%ylmgr=>ylmgr

 ab_scfcv_inout%cg=>cg
 ab_scfcv_inout%dtefield=>dtefield
 ab_scfcv_inout%dtfil=>dtfil
 ab_scfcv_inout%dtpawuj=>dtpawuj
 ab_scfcv_inout%eigen=>eigen
 ab_scfcv_inout%hdr=>hdr
 ab_scfcv_inout%initialized=>initialized
 ab_scfcv_inout%irrzon=>irrzon
 ab_scfcv_inout%mpi_enreg=>mpi_enreg
 ab_scfcv_inout%nfftf=>nfftf
 ab_scfcv_inout%occ=>occ
 ab_scfcv_inout%pawfgr=>pawfgr
 ab_scfcv_inout%pawrhoij=>pawrhoij
 ab_scfcv_inout%rec_set=>rec_set
 ab_scfcv_inout%resid=>resid
 ab_scfcv_inout%results_gs=>results_gs
 ab_scfcv_inout%scf_history=>scf_history
 ab_scfcv_inout%symrec=>symrec
 ab_scfcv_inout%taug=>taug
 ab_scfcv_inout%taur=>taur
 ab_scfcv_inout%wvl=>wvl

!ab_scfcv_inout%dtset=>dtset
!ab_scfcv_inout%electronpositron=>electronpositron
!ab_scfcv_inout%paw_dmft=>paw_dmft
!ab_scfcv_inout%rhog=>rhog
!ab_scfcv_inout%rhor=>rhor
!ab_scfcv_inout%rprimd=>rprimd
!ab_scfcv_inout%wffnew=>wffnew
!ab_scfcv_inout%wffnow=>wffnow
!ab_scfcv_inout%xred=>xred
!ab_scfcv_inout%xred_old=>xred_old

end subroutine scfcv_init
!!***
