!{\src2tex{textfont=tt}}
!!****f* ABINIT/mover
!! NAME
!! mover
!!
!! FUNCTION
!! Move ion or change acell acording to forces and stresses
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amass(natom)=mass of each atom, in unit of electronic mass (=amu*1822...)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
!!   |  angular momentum for nonlocal pseudopotential
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in unit cell
!!   |  except on first call (hartree/bohr); updated on output
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points.
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=informations about MPI parallelization
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  npwarr(nkpt)=number of planewaves in basis and boundary at this k point.
!!  nattyp(ntypat)= # atoms of each type.
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!
!! OUTPUT
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points.
!!
!! SIDE EFFECTS
!! Rest of i/o is related to lda
!!  acell(3)=length scales of primitive translations (bohr)
!!  cg(2,mcg)=array for planewave coefficients of wavefunctions.
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  initialized= if 0 the initialisation of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  occ(mband*nkpt*nsppol=occupation number for each band (usually 2) at each k point.
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in electrons/bohr**3.
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  wffnew,wffnow= struct info for wf disk files.
!!  vel(3,natom)=old value of velocity; updated on output
!!  xred(3,natom)=reduced dimensionless atomic coordinates; updated on output
!!  xred_old(3,natom)=work space for old xred
!!
!! NOTES
!! This subroutine uses the arguments natom, xred, vel, fcart, amass,
!! vis, and dtion (the last two contained in dtset) to make
!! molecular dynamics updates.  The rest of the lengthy
!! argument list supports the underlying lda computation
!! of forces, returned from subroutine scfcv
!!
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      ab_movehistory_bcast,ab_movehistory_end,ab_movehistory_nullify
!!      ab_movetype_end,ab_movetype_nullify,fcart2fred,fconv,fred2fcart
!!      hist2var,hist_compare,initylmg,metric,mkradim,prec_simple,pred_bfgs
!!      pred_delocint,pred_diisrelax,pred_isokinetic,pred_isothermal
!!      pred_langevin,pred_moldyn,pred_nose,pred_simple,pred_srkna14
!!      pred_steepdesc,pred_verlet,prtxfase,read_md_hist,scfcv_new2,symzat
!!      var2hist,vel2hist,write_md_hist,wrt_moldyn_netcdf,wrtout,xcomm_world
!!      xfh_update,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mover(scfcv_args,&
  & ab_xfh,acell,amass,dtfil,&
  & electronpositron,&
  & rhog,rhor,&
  & rprimd,&
  & vel,&
  & xred,xred_old)

 use m_profiling

use defs_basis
use defs_datatypes
use defs_abitypes
!use defs_scftypes
!use defs_rectypes
!use defs_wvltypes
!use m_paw_dmft, only: paw_dmft_type
use m_electronpositron, only : electronpositron_type
!use m_wffile
use defs_scfcvargs
use defs_mover

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mover'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
 use interfaces_45_geomoptim
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
 use interfaces_61_ionetcdf
 use interfaces_67_common
 use interfaces_95_drive, except_this_one => mover
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
type(ab_scfcvargs),intent(inout) :: scfcv_args
!type(MPI_type),intent(inout) :: mpi_enreg
type(datafiles_type),intent(inout),target :: dtfil
!type(dataset_type),intent(inout),target :: dtset
type(electronpositron_type),pointer :: electronpositron
!type(paw_dmft_type) :: paw_dmft
!type(ab_scfcv_args_in),intent(inout) :: ab_scfcv_in
!type(ab_scfcv_args_inout),intent(inout) :: ab_scfcv_inout
!type(wffile_type),intent(inout) :: wffnew,wffnow
type(ab_xfh_type),intent(inout) :: ab_xfh
!arrays
!no_abirules
real(dp),intent(inout) :: acell(3)
real(dp), intent(in),target :: amass(scfcv_args%dtset%natom)
real(dp), pointer :: rhog(:,:),rhor(:,:)
real(dp), intent(inout) :: xred(3,scfcv_args%dtset%natom),xred_old(3,scfcv_args%dtset%natom)
real(dp), intent(inout) :: vel(3,scfcv_args%dtset%natom),rprimd(3,3)

!Local variables-------------------------------
!scalars
type(ab_movehistory) :: hist,hist_prev
type(ab_movetype) :: ab_mover
type(ab_forstr) :: forstr_precon ! Preconditioned forces and stress
type(mttk_type) :: mttk_vars
integer,parameter :: level=102
integer :: itime,icycle,iexit,ncycle,kk,jj,master,me,nloop,ntime,option,paral_compil,spaceworld
character(len=500) :: message
character(len=60) :: method
character(len=10) :: type4xml
character(len=8) :: crit4xml
character(len=8) :: stat4xml
character(len=35) :: fmt
character(len=fnlen) :: filename
real(dp) :: ucvol,favg
logical :: isFconv ! If the convergence is needed
logical :: DEBUG=.FALSE.
logical :: change
integer :: minIndex,ii,similar
real(dp) :: minE,tolerance

!arrays
real(dp) :: xred_tmp(3,scfcv_args%dtset%natom)
real(dp) :: xcart(3,scfcv_args%dtset%natom)

real(dp) :: fcart2(3,scfcv_args%dtset%natom)
real(dp) :: fred2(3,scfcv_args%dtset%natom)
real(dp) :: favg2(3)

real(dp) :: fred_corrected(3,scfcv_args%dtset%natom)
real(dp) :: rprim(3,3)
real(dp) :: gprimd(3,3)
real(dp) :: gmet(3,3)
real(dp) :: rmet(3,3)

! ***************************************************************

!Table of contents
!(=>) Refers to an important call (scfcv,pred_*)
!
!01. Initialization of indexes and allocations of arrays
!02. Particularities of each predictor
!03. Set the number of iterations ntime
!04. Try to read history of previous calculations
!05. Allocate the hist structure
!06. First output before any itime or icycle
!07. Compute xcart and fill the history of the first SCFCV
!08. Loop for itime (From 1 to ntime)
!09. Loop for icycle (From 1 to ncycles)
!10. Output for each icycle (and itime)
!11. Symmetrize atomic coordinates over space group elements
!12. => Call to SCFCV routine and fill history with forces
!13. Write the history into the _HIST file
!14. Output after SCFCV
!15. => Test Convergence of forces and stresses
!16. => Precondition forces, stress and energy
!17. => Call to each predictor
!18. Use the history  to extract the new values
!19. End loop icycle
!20. End loop itime
!21. Set the final values of xcart and xred
!22. XML Output at the end
!23. Deallocate hist and ab_mover datatypes
!
!Conditions:
!
!This routine (mover.F90) should produce output, no mathematical
!expresions should appear.
!
!History:
!
!History of several variables are under the type "ab_movehistory"
!It contains:
!* mxhist                  : Maximum size of history
!* ihist                   : index of history
!* histA(3,mxhist)         : Acell
!* histE(mxhist)           : Energy
!* histEk(mxhist)          : Ionic Kinetic Energy
!* histT(mxhist)           : Time (Or iteration number for GO)
!* histR(3,3,mxhist)       : Rprimd
!* histS(6,mxhist)         : Strten
!* histV(3,natom,mxhist)   : Velocity
!* histXF(3,natom,4,mxhist): Xcart,Xred,Fcart,Fred
!
!IONMOV values:
!
!1.  Molecular dynamics without viscosity (vis=0)
!1.  Molecular dynamics with viscosity (vis/=0)
!2.  Broyden-Fletcher-Goldfard-Shanno method (forces)
!3.  Broyden-Fletcher-Goldfard-Shanno method (forces,Tot energy)
!4.  Conjugate gradient of potential and ionic degrees of freedom
!5.  Simple relaxation of ionic positions
!6.  Verlet algorithm for molecular dynamics
!7.  Verlet algorithm blocking every atom where dot(vel,force)<0
!8.  Verlet algorithm with a nose-hoover thermostat
!9.  Langevin molecular dynamics
!10. BFGS with delocalized internal coordinates
!11. Conjugate gradient algorithm
!12. Isokinetic ensemble molecular dynamics
!13. Isothermal/isenthalpic ensemble molecular dynamics
!14. Symplectic algorithm Runge-Kutta-Nyström SRKNa14
!20. Ionic positions relaxation using DIIS
!21. Steepest descent algorithm
!30. Self consistent phonon structure using a supercell

!write(std_out,*) 'mover 01'
!###########################################################
!### 01. Initialization of indexes and allocations of arrays

 call ab_movetype_nullify(ab_mover)

!Copy the information from the Dataset (dtset)
!to the ab_mover structure
!Delay of Permutation (Used by pred_langevin only)
 ab_mover%delayperm=>scfcv_args%dtset%delayperm
!DIIS memory (Use by pred_diisrelax only)
 ab_mover%diismemory=>scfcv_args%dtset%diismemory
!Geometry Optimization Precondition option
 ab_mover%goprecon=>scfcv_args%dtset%goprecon
!include a JELLium SLAB in the cell
 ab_mover%jellslab=>scfcv_args%dtset%jellslab
!Number of ATOMs
 ab_mover%natom=>scfcv_args%dtset%natom
!Number of CONstraint EQuations
 ab_mover%nconeq=>scfcv_args%dtset%nconeq
!Use by pred_isothermal only
 ab_mover%nnos=>scfcv_args%dtset%nnos
!Number of SYMmetry operations
 ab_mover%nsym=>scfcv_args%dtset%nsym
!Number of Type of atoms
 ab_mover%ntypat=>scfcv_args%dtset%ntypat
!OPTimize the CELL shape and dimensions
 ab_mover%optcell=>scfcv_args%dtset%optcell
!RESTART Xcart and Fred
 ab_mover%restartxf=>scfcv_args%dtset%restartxf
!Sign of Permutation (Used by pred_langevin only)
 ab_mover%signperm=>scfcv_args%dtset%signperm
!Ion movements
 ab_mover%ionmov=>scfcv_args%dtset%ionmov

!Use by pred_isothermal only
 ab_mover%bmass=>scfcv_args%dtset%bmass
!Delta Time for IONs
 ab_mover%dtion=>scfcv_args%dtset%dtion
!Used by pred_langevin only
 ab_mover%friction=>scfcv_args%dtset%friction
!Used by pred_langevin only
 ab_mover%mdwall=>scfcv_args%dtset%mdwall
!NOT DOCUMENTED
 ab_mover%noseinert=>scfcv_args%dtset%noseinert
!STRess PRECONditioner
 ab_mover%strprecon=>scfcv_args%dtset%strprecon
!VIScosity
 ab_mover%vis=>scfcv_args%dtset%vis

!Indices of AToms that are FIXed
 ab_mover%iatfix=>scfcv_args%dtset%iatfix
!SYMmetry in REaL space
 ab_mover%symrel=>scfcv_args%dtset%symrel
!TYPe of ATom
 ab_mover%typat=>scfcv_args%dtset%typat
!PRTint QTom LIST
 ab_mover%prtatlist=>scfcv_args%dtset%prtatlist

!Mass of each atom (NOT IN DTSET)
 ab_mover%amass=>amass
!Geometry Optimization Preconditioner PaRaMeters
 ab_mover%goprecprm=>scfcv_args%dtset%goprecprm
!Molecular Dynamics Temperatures
 ab_mover%mdtemp=>scfcv_args%dtset%mdtemp
!STRess TARGET
 ab_mover%strtarget=>scfcv_args%dtset%strtarget
!Use by pred_isothermal only
 ab_mover%qmass=>scfcv_args%dtset%qmass
!Z number of each NUCLeus
 ab_mover%znucl=>scfcv_args%dtset%znucl

!Filename for Hessian matrix (NOT IN DTSET)
 ab_mover%fnameabi_hes=>dtfil%fnameabi_hes
!Filename for _HIST file
 ab_mover%filnam_ds=>dtfil%filnam_ds

 call ab_movehistory_nullify(hist)
 call ab_movehistory_nullify(hist_prev)

!Init MPI data
 master=0;paral_compil=scfcv_args%ab_scfcv_inout%mpi_enreg%paral_compil
 call xcomm_world(scfcv_args%ab_scfcv_inout%mpi_enreg,spaceworld,myrank=me)

!!DEBUG
!call ab_movetype_print(ab_mover,ab_out)
!!DEBUG

!write(std_out,*) 'mover 02'
!###########################################################
!### 02. Particularities of each predictor

!Default values first
!--------------------

!acell and rprimd are never change except if optcell/=0
 if (ab_mover%optcell/=0)then
   hist%isARused=.TRUE.
 else
   hist%isARused=.FALSE.
 end if

!Velocities are never change except for ionmov=1,6,7,8
 hist%isVused=.FALSE.

!In general convergence is needed
 isFconv=.TRUE.

!ncycle is 1 by default except for ionmov=1,9,14
 ncycle=1

!This is the initialization for ionmov==1
!-----------------------------------------
 if(ab_mover%ionmov==1) then
   ncycle=4 ! Number of internal cycles for first itime
   isFconv=.FALSE.     ! Convergence is not used for MD
   hist%isVused=.TRUE. ! Velocities are used
!  TEMPORARLY optcell is not allow
   hist%isARused=.FALSE.
!  Values use in XML Output
   type4xml='moldyn'
   crit4xml='none'
!  Name of method
   if (abs(ab_mover%vis)<=1.d-8) then
     write(method,'(a)')&
     'Molecular dynamics without viscosity (vis=0)'
   else
     write(method,'(a,1p,e12.5,a)')&
     'Molecular dynamics with viscosity (vis=',&
&     ab_mover%vis,')'
   end if
!  This is the initialization for ionmov==2,3
!  -------------------------------------------
 else if (ab_mover%ionmov==2.or.ab_mover%ionmov==3)then
!  Values use in XML Output
   type4xml='bfgs'
   crit4xml='tolmxf'
!  Name of method
   if (ab_mover%ionmov==2) then
     write(method,'(a)')&
     'Broyden-Fletcher-Goldfard-Shanno method (forces)'
   else
     write(method,'(a)')&
     'Broyden-Fletcher-Goldfard-Shanno method (forces,Tot energy)'
   end if
!  This is the initialization for ionmov==4,5
!  -------------------------------------------
 else if (ab_mover%ionmov==4.or.ab_mover%ionmov==5)then
!  Values use in XML Output
   type4xml='simple'
   crit4xml='tolmxf'
!  Name of method
   if (ab_mover%ionmov==4) then
     write(method,'(a)')&
     'Conjugate gradient of potential and ionic degrees of freedom'
   else
     write(method,'(a)')&
     'Simple relaxation of ionic positions'
   end if
!  This is the initialization for ionmov==6
!  ------------------------------------------
 else if (ab_mover%ionmov==6)then
   isFconv=.FALSE.     ! Convergence is not used for MD
!  TEMPORARLY optcell is not allow
   hist%isARused=.FALSE.
   hist%isVused=.TRUE. ! Velocities are used
!  Values use in XML Output
   type4xml='verlet'
   crit4xml='tolmxf'
!  Name of method
   write(method,'(a)')&
   'Verlet algorithm for molecular dynamics'
!  This is the initialization for ionmov==7
!  ------------------------------------------
 else if (ab_mover%ionmov==7)then
!  TEMPORARLY optcell is not allow
   hist%isARused=.FALSE.
   hist%isVused=.TRUE. ! Velocities are used
!  Values use in XML Output
   type4xml='verlet'
   crit4xml='tolmxf'
!  Name of method
   write(method,'(a)')&
   'Verlet algorithm blocking every atom where dot(vel,force)<0'
!  This is the initialization for ionmov==8
!  ------------------------------------------
 else if (ab_mover%ionmov==8)then
   hist%isVused=.TRUE.
!  TEMPORARLY optcell is not allow
   hist%isARused=.FALSE.
!  Values use in XML Output
   type4xml='nose'
   crit4xml='tolmxf'
!  Name of method
   write(method,'(a)')&
   'Verlet algorithm with a nose-hoover thermostat'
!  This is the initialization for ionmov==9
!  ------------------------------------------
 else if (ab_mover%ionmov==9)then
!  TEMPORARLY optcell is not allow
   hist%isARused=.FALSE.
   hist%isVused=.TRUE.  ! Velocities are used
   ncycle=3
!  Values use in XML Output
   type4xml='langevin'
   crit4xml='tolmxf'
!  Name of method
   write(method,'(a)')&
   'Langevin molecular dynamics'
!  This is the initialization for ionmov==10
!  -------------------------------------------
 else if (ab_mover%ionmov==10)then
!  TEMPORARLY optcell is not allow
   hist%isARused=.FALSE.
!  Values use in XML Output
   type4xml='delocint'
   crit4xml='tolmxf'
!  Name of method
   write(method,'(a)')&
   'BFGS with delocalized internal coordinates'
!  This is the initialization for ionmov==11
!  -------------------------------------------
 else if (ab_mover%ionmov==11)then
!  TEMPORARLY optcell is not allow
   hist%isARused=.FALSE.
!  Values use in XML Output
   type4xml='cg'
   crit4xml='tolmxf'
!  Name of method
   write(method,'(a)')&
   'Conjugate gradient algorithm'
!  This is the initialization for ionmov==12
!  -------------------------------------------
 else if (ab_mover%ionmov==12)then
!  TEMPORARLY optcell is not allow
   hist%isARused=.FALSE.
   hist%isVused=.TRUE.  ! Velocities are used
!  Values use in XML Output
   isFconv=.FALSE.      ! Convergence is not used for MD
   type4xml='isokin'
   crit4xml='tolmxf'
!  Name of method
   write(method,'(a)')&
   'Isokinetic ensemble molecular dynamics'
!  This is the initialization for ionmov==13
!  -------------------------------------------
 else if (ab_mover%ionmov==13)then
!  optcell is allow
   hist%isARused=.TRUE. ! RPRIMD and ACELL may change
   hist%isVused=.TRUE.  ! Velocities are used
   isFconv=.FALSE.      ! Convergence is not used for MD
!  Values use in XML Output
   type4xml='isother'
   crit4xml='tolmxf'
   ABI_ALLOCATE(mttk_vars%glogs,(ab_mover%nnos))
   ABI_ALLOCATE(mttk_vars%vlogs,(ab_mover%nnos))
   ABI_ALLOCATE(mttk_vars%xlogs,(ab_mover%nnos))
!  Name of method
   write(method,'(a)')&
   'Isothermal/isenthalpic ensemble molecular dynamics'
!  This is the initialization for ionmov==14
!  -------------------------------------------
 else if (ab_mover%ionmov==14)then
   ncycle=16
   isFconv=.FALSE.     ! Convergence is not used for MD
   hist%isVused=.TRUE. ! Velocities are used
!  TEMPORARLY optcell is not allow
   hist%isARused=.FALSE.
!  Values use in XML Output
   type4xml='srkna14'
   crit4xml='tolmxf'
!  Name of method
   write(method,'(a)')&
   'Symplectic algorithm Runge-Kutta-Nyström SRKNa14'
!  This is the initialization for ionmov==20
!  -------------------------------------------
 else if (ab_mover%ionmov==20)then
!  TEMPORARLY optcell is not allow
   hist%isARused=.FALSE.
!  Values use in XML Output
   type4xml='diisrelax'
   crit4xml='tolmxf'
!  Name of method
   write(method,'(a)')&
   'Ionic positions relaxation using DIIS'
!  This is the initialization for ionmov==21
!  -------------------------------------------
 else if (ab_mover%ionmov==21)then
   hist%isARused=.TRUE.
!  Values use in XML Output
   type4xml='steepdesc'
   crit4xml='tolmxf'
!  Name of method
   write(method,'(a)')&
   'Steepest descend algorithm'
!  This is the initialization for ionmov==30
!  -------------------------------------------
 else if (ab_mover%ionmov==30)then
!  TEMPORARLY optcell is not allow
   hist%isARused=.FALSE.
!  Values use in XML Output
   type4xml='scphon'
   crit4xml='tolmxf'
!  Name of method
   write(method,'(a)')&
   'Self consistent phonon structure using a supercell'
 end if

!write(std_out,*) 'mover 03'
!###########################################################
!### 03. Set the number of iterations ntime
!###     By default ntime==1 but if the user enter a lower
!###     value mover will execute at least one iteration

 if (scfcv_args%dtset%ntime<1)then
   ntime=1
 else
   ntime=scfcv_args%dtset%ntime
 end if

!write(std_out,*) 'mover 04'
!###########################################################
!### 04. Try to read history of previous calculations
!###     It requires access to the NetCDF library

#if defined HAVE_TRIO_NETCDF
 filename=trim(ab_mover%filnam_ds(4))//'_HIST'

 if (ab_mover%restartxf<=0)then
!  Read history from file (and broadcast if MPI)
   if (me==0) then
     call read_md_hist(filename,hist_prev)
   end if
   if (paral_compil==1) then
     call ab_movehistory_bcast(hist_prev,master,spaceworld)
   end if

!  If restartxf specifies to reconstruct the history
   if (hist_prev%mxhist>0.and.ab_mover%restartxf==-1)then
     ntime=ntime+hist_prev%mxhist
   end if

!  If restartxf specifies to start from the lowest energy
   if (hist_prev%mxhist>0.and.ab_mover%restartxf==-2)then
     minE=hist_prev%histE(1)
     minIndex=1
     do ii=1,hist_prev%mxhist
       write(std_out,*) 'Iteration:',ii,' Total Energy:',hist_prev%histE(ii)
       if (minE>hist_prev%histE(ii))then
         minE=hist_prev%histE(ii)
         minIndex=ii
       end if
     end do
     write(std_out,*) 'The lowest energy occurs at iteration:',minIndex,'etotal=',minE
     acell(:)   =hist_prev%histA(:,minIndex)
     xcart(:,:) =hist_prev%histXF(:,:,1,minIndex)
     xred(:,:)  =hist_prev%histXF(:,:,2,minIndex)
     rprimd(:,:)=hist_prev%histR(:,:,minIndex)
   end if

 end if !if (ab_mover%restartxf<=0)

#endif

!write(std_out,*) 'mover 05'
!###########################################################
!### 05. Allocate the hist structure

 hist%mxhist=ncycle*ntime
 if (scfcv_args%dtset%nctime>0) hist%mxhist=hist%mxhist+1

!Initialize indexes
 hist%ihist=1
 iexit=0

!DEBUG
!write(std_out,*) 'MXHIST:',hist%mxhist
!write(std_out,*) 'NTIME :',ntime
!write(std_out,*) 'NCYCLE:',ncycle
!DEBUG

!Allocate all the histories
 ABI_ALLOCATE(hist%histA,(3,hist%mxhist))
 ABI_ALLOCATE(hist%histE,(hist%mxhist))
 ABI_ALLOCATE(hist%histEk,(hist%mxhist))
 ABI_ALLOCATE(hist%histT,(hist%mxhist))
 ABI_ALLOCATE(hist%histR,(3,3,hist%mxhist))
 ABI_ALLOCATE(hist%histS,(6,hist%mxhist))
 ABI_ALLOCATE(hist%histV,(3,ab_mover%natom,hist%mxhist))
 ABI_ALLOCATE(hist%histXF,(3,ab_mover%natom,4,hist%mxhist))

 ABI_ALLOCATE(forstr_precon%fcart,(3,ab_mover%natom))
 ABI_ALLOCATE(forstr_precon%fred,(3,ab_mover%natom))

!write(std_out,*) 'mover 06'
!###########################################################
!### 06. First output before any itime or icycle

 write(message,'(a,a,i2,a,a,a,80a)')&
& ch10,'=== [ionmov=',ab_mover%ionmov,'] ',method,&
& ch10,('=',kk=1,80)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

!Format for printing on each cycle
 write(fmt,'(a6,i2,a4,i2,a4,i2,a4,i2,a9)')&
& '(a,a,i',int(log10(real(ntime))+1),&
& ',a,i',int(log10(real(ntime))+1),&
& ',a,i',int(log10(real(ncycle))+1),&
& ',a,i',int(log10(real(ncycle))+1),&
& ',a,a,80a)'
!write(std_out,*) fmt
!write(std_out,*) HUGE(1)
!write(std_out,*) int(log10(real(HUGE(1)))+1)

!write(std_out,*) 'mover 07'
!###########################################################
!### 07. Compute xcart and fill the history of the first SCFCV

!Compute xcart from xred, and rprimd
 call xredxcart(ab_mover%natom,1,rprimd,xcart,xred)

!Compute rprim from rprimd and acell
 do kk=1,3
   do jj=1,3
     rprim(jj,kk)=rprimd(jj,kk)/acell(kk)
   end do
 end do

!Fill history with the values of xred,xcart,acell and rprimd
 call var2hist(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,DEBUG)

!Fill velocities and ionic kinetic energy
 call vel2hist(ab_mover%amass,hist,ab_mover%natom,vel)
 hist%histT(hist%ihist)=0.0

!write(std_out,*) 'mover 08'
!###########################################################
!### 08. Loop for itime (From 1 to ntime)
 do itime=1,ntime
!  call status(1000*itime,dtfil%filstat,iexit,level,'Loop itime')

!  write(std_out,*) 'mover 09'
!  ###########################################################
!  ### 09. Loop for icycle (From 1 to ncycles)
   do icycle=1,ncycle
!    call status(1000*itime+icycle,dtfil%filstat,iexit,level,'Loop icycle')

!    write(std_out,*) 'mover 10'
!    ###########################################################
!    ### 10. Output for each icycle (and itime)

     write(message,fmt)&
&     ch10,'--- Iteration: (',itime,'/',ntime,&
&     ') Internal Cycle: (',icycle,'/',ncycle,')',ch10,&
&     ('-',kk=1,80)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

     call prtxfase(ab_mover,hist,std_out,mover_BEFORE)

!    write(std_out,*) 'mover 11'
!    ###########################################################
!    ### 11. Symmetrize atomic coordinates over space group elements

     change=.FALSE.
     xred_tmp(:,:)=xred(:,:)
     call symzat(scfcv_args%ab_scfcv_in%indsym,&
&     ab_mover%natom,&
&     scfcv_args%dtset%nsym,&
&     scfcv_args%dtset%symrel,&
&     scfcv_args%dtset%tnons,&
&     xred)
     do kk=1,ab_mover%natom
       do jj=1,3
         if (xred(jj,kk)/=xred_tmp(jj,kk)) change=.TRUE.
       end do
     end do

     if (change)then
       call xredxcart(ab_mover%natom,1,rprimd,xcart,xred)
       hist%histXF(:,:,1,hist%ihist)=xcart(:,:)
       hist%histXF(:,:,2,hist%ihist)=xred(:,:)
       write(std_out,*) 'WARNING: ATOMIC COORDINATES WERE SYMMETRIZED'
       write(std_out,*) 'DIFFERENCES:'

       do kk=1,ab_mover%natom
         write(std_out,*) xred(:,kk)-xred_tmp(:,kk)
       end do

       xred_tmp(:,:)=xred(:,:)

     end if

!    write(std_out,*) 'mover 12'
!    ###########################################################
!    ### 12. => Call to SCFCV routine and fill history with forces
     write(message,'(a,3a,33a,44a)')&
&     ch10,('-',kk=1,3),&
&     'SELF-CONSISTENT-FIELD CONVERGENCE',('-',kk=1,44)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

     if(hist_prev%mxhist>0.and.ab_mover%restartxf==-1.and.hist_prev%ihist<=hist_prev%mxhist)then

       call hist_compare(hist_prev,hist,ab_mover%natom,similar,tolerance)
       hist_prev%ihist=hist_prev%ihist+1

     else

       scfcv_args%ab_scfcv_in%ndtpawuj=0
       scfcv_args%ab_scfcv_in%iapp=itime
       if(icycle>1.and.icycle/=ncycle) scfcv_args%ab_scfcv_in%iapp=-1
       if(itime==1 .and. icycle/=ncycle ) scfcv_args%ab_scfcv_in%iapp=-icycle-1
       if (ab_mover%ionmov==14.and.(icycle<ncycle)) scfcv_args%ab_scfcv_in%iapp=-1

       call scfcv_new2(scfcv_args,&
&       electronpositron,&
&       rhog,rhor,&
&       rprimd,&
&       xred,xred_old)

!      ANOMALOUS SITUATION
!      This is the only case where rprimd could change inside scfcv
!      It generates an weird condition, we start with a certain
!      value for rprimd before scfcv and after we finish with
!      a different value.
!      Notice that normally scfcv should not change rprimd
!      And even worse if optcell==0
!      The solution here is to recompute acell and rprim
!      and store those values in the present record
!      even if initially those values were not exactly
!      the values entering in scfcv
!      
!      One test case with these condition is bigdft/t10
       if (scfcv_args%dtset%usewvl == 1) then
         call mkradim(acell,rprim,rprimd)
         hist%histA(:,hist%ihist)=acell(:)
         hist%histR(:,:,hist%ihist)=rprimd(:,:)
       end if

!      ANOMALOUS SITUATIONS
!      * In ionmov 4 & 5 xred could change inside SCFCV
!      So we need to take the values from the output
!      
!      * Inside scfcv.F90 there is a call to symzat.F90
!      for the first SCF cycle symzat could change xred
!      so we need always take xred convert to xcart and
!      store in the history

       if (ab_mover%ionmov<10)then

         change=.FALSE.

         do kk=1,ab_mover%natom
           do jj=1,3
             if (xred(jj,kk)/=xred_tmp(jj,kk)) change=.TRUE.
           end do
         end do

         if (change)then
           call xredxcart(ab_mover%natom,1,rprimd,xcart,xred)
           hist%histXF(:,:,1,hist%ihist)=xcart(:,:)
           hist%histXF(:,:,2,hist%ihist)=xred(:,:)
           write(std_out,*) 'WARNING: ATOMIC COORDINATES WERE SYMMETRIZED AFTER SCFCV'
           write(std_out,*) 'DIFFERENCES:'

           do kk=1,ab_mover%natom
             write(std_out,*) xred(:,kk)-xred_tmp(:,kk)
           end do

         end if !if (change)

!        call xredxcart(ab_mover%natom,1,rprimd,xcart,xred)
!        hist%histXF(:,:,1,hist%ihist)=xcart(:,:)
!        hist%histXF(:,:,2,hist%ihist)=xred(:,:)
       end if

!      Fill velocities and ionic kinetic energy
       call vel2hist(ab_mover%amass,hist,ab_mover%natom,vel)
       hist%histE(hist%ihist)       =scfcv_args%ab_scfcv_inout%results_gs%etotal
       hist%histXF(:,:,3,hist%ihist)=scfcv_args%ab_scfcv_inout%results_gs%fcart(:,:)
       hist%histXF(:,:,4,hist%ihist)=scfcv_args%ab_scfcv_inout%results_gs%fred(:,:)
       hist%histS(:,hist%ihist)     =scfcv_args%ab_scfcv_inout%results_gs%strten(:)
       hist%histT(hist%ihist)       =itime

!      !######################################################################
!      ! Test of convertion fcart and fred

       call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

       call fred2fcart(favg2,fcart2,scfcv_args%ab_scfcv_inout%results_gs%fred,gprimd,ab_mover%jellslab,ab_mover%natom)
       call fcart2fred(scfcv_args%ab_scfcv_inout%results_gs%fcart,fred2,rprimd,ab_mover%natom)

!      write (std_out,*) 'FCART'
!      do kk=1,ab_mover%natom
!      write (std_out,*) scfcv_args%ab_scfcv_inout%results_gs%fcart(:,kk)
!      end do
!      write (std_out,*) 'FCART converted from fred'
!      do kk=1,ab_mover%natom
!      write (std_out,*) fcart2(:,kk)
!      end do
!      write (std_out,*) 'FRED'
!      do kk=1,ab_mover%natom
!      write (std_out,*) scfcv_args%ab_scfcv_inout%results_gs%fred(:,kk)
!      end do
!      write (std_out,*) 'FRED converted from fcart'
!      do kk=1,ab_mover%natom
!      write (std_out,*) fred2(:,kk)
!      end do

!      !###################################################################### 

     end if ! if(hist_prev%mxhist>0.and.ab_mover%restartxf==-1.and.hist_prev%ihist<=hist_prev%mxhist)then

     if(ab_xfh%nxfh==0.or.itime/=1) then

       call mkradim(acell,rprim,rprimd)
!      Get rid of mean force on whole unit cell, but only if no
!      generalized constraints are in effect
!      hist%histXF(:,:,4,hist%ihist) are reduced forces
       if(ab_mover%nconeq==0)then
         do kk=1,3
           favg=sum(hist%histXF(kk,:,4,hist%ihist))/dble(ab_mover%natom)
           fred_corrected(kk,:)=hist%histXF(kk,:,4,hist%ihist)-favg
           if(ab_mover%jellslab/=0.and.kk==3)&
&           fred_corrected(kk,:)=hist%histXF(kk,:,4,hist%ihist)
         end do
       else
         fred_corrected(:,:)=hist%histXF(:,:,4,hist%ihist)
       end if

       if (ncycle<10.and.ab_mover%restartxf>=0) call xfh_update(ab_xfh,acell,fred_corrected,ab_mover%natom,rprim,&
&       hist%histS(:,hist%ihist),xred)

     end if

!    write(std_out,*) 'mover 13'
!    ###########################################################
!    ### 13. Write the history into the _HIST file
!    ###

#if defined HAVE_TRIO_NETCDF
!    write(std_out,*) 'ihist @ mover',hist%ihist
!    write(std_out,*) 'mxhist @ mover',hist%mxhist
     if (me==0) then
       call write_md_hist(filename,hist,icycle,itime,ab_mover%natom)
     end if
#endif

!    write(std_out,*) 'mover 14'
!    ###########################################################
!    ### 14. Output after SCFCV
     write(message,'(a,3a,a,72a)')&
&     ch10,('-',kk=1,3),&
&     'OUTPUT',('-',kk=1,71)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

     call prtxfase(ab_mover,hist,ab_out,mover_AFTER)
     call prtxfase(ab_mover,hist,std_out,mover_AFTER)

!    !    DEBUG (XRA AFTER SCFCV)
!    if(DEBUG)then
!    write (std_out,*) '---XRA AFTER SCFCV---'
!    write (std_out,*) 'XCART'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xcart(:,kk)
!    end do
!    write (std_out,*) 'XRED'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xred(:,kk)
!    end do
!    if (ab_mover%ionmov==1)then
!    write (std_out,*) 'VEL'
!    do kk=1,ab_mover%natom
!    write (std_out,*) hist%histV(:,kk,hist%ihist)
!    end do
!    end if
!    write(std_out,*) 'RPRIMD'
!    do kk=1,3
!    write(std_out,*) rprimd(:,kk)
!    end do
!    write(std_out,*) 'ACELL'
!    write(std_out,*) acell(:)
!    end if

!    write(std_out,*) 'mover 15'
!    ###########################################################
!    ### 15. => Test Convergence of forces and stresses

     if (itime==ntime.and.icycle==ncycle)then
       iexit=1
       stat4xml="Failed"
     else
       stat4xml="Succeded"
     end if

!    Only if convergence is needed
     if(isFconv)then
       if ((ab_mover%ionmov/=4.and.ab_mover%ionmov/=5).or.mod(itime,2)==1)then
         call fconv(hist%histXF(:,:,3,hist%ihist),&
&         scfcv_args%dtset%iatfix, &
&         iexit, itime,&
&         ab_mover%natom,&
&         ntime,&
&         ab_mover%optcell,&
&         scfcv_args%dtset%strfact,&
&         scfcv_args%dtset%strtarget,&
&         hist%histS(:,hist%ihist),&
&         scfcv_args%dtset%tolmxf)
       end if
     end if

     if (itime==ntime.and.icycle==ncycle) iexit=1

!    write(std_out,*) 'mover 16'
!    ###########################################################
!    ### 16. => Precondition forces, stress and energy

     write(std_out,*) 'Geometry Optimization Precondition:',ab_mover%goprecon
     if (ab_mover%goprecon>0)then
       call prec_simple(ab_mover,forstr_precon,hist,icycle,itime,0)
     end if

!    write(std_out,*) 'mover 16'
!    ###########################################################
!    ### 17. => Call to each predictor

!    MT->GAF: dirty trick to predict vel(t)
!    do a double loop: 1- compute vel, 2- exit
     nloop=1
     if (scfcv_args%dtset%nctime>0.and.iexit==1) then
       iexit=0;nloop=2
     end if
     do ii=1,nloop
       if (ii==2) iexit=1

       if(ab_mover%ionmov==1) then
         call pred_moldyn(ab_mover,hist,icycle,itime,ncycle,ntime,DEBUG,iexit)
       else if (ab_mover%ionmov==2.or.ab_mover%ionmov==3)then
         call pred_bfgs(ab_mover,ab_xfh,forstr_precon,hist,ab_mover%ionmov,itime,DEBUG,iexit)
       else if (ab_mover%ionmov==4.or.ab_mover%ionmov==5)then
         call pred_simple(ab_mover,hist,iexit)
       else if (ab_mover%ionmov==6.or.ab_mover%ionmov==7)then
         call pred_verlet(ab_mover,hist,ab_mover%ionmov,itime,ntime,DEBUG,iexit)
       else if (ab_mover%ionmov==8)then
         call pred_nose(ab_mover,hist,itime,ntime,DEBUG,iexit)
       else if (ab_mover%ionmov==9)then
         call pred_langevin(ab_mover,hist,icycle,itime,ncycle,ntime,DEBUG,iexit)
       else if (ab_mover%ionmov==10.or.ab_mover%ionmov==11)then
         call pred_delocint(ab_mover,ab_xfh,forstr_precon,hist,ab_mover%ionmov,itime,DEBUG,iexit)
       else if (ab_mover%ionmov==12)then
         call pred_isokinetic(ab_mover,hist,itime,ntime,DEBUG,iexit)
       else if (ab_mover%ionmov==13)then
         call pred_isothermal(ab_mover,hist,itime,mttk_vars,ntime,DEBUG,iexit)
       else if (ab_mover%ionmov==14)then
         call pred_srkna14(ab_mover,hist,icycle,ncycle,DEBUG,iexit)
       else if (ab_mover%ionmov==20)then
         call pred_diisrelax(ab_mover,hist,itime,ntime,DEBUG,iexit)
       else if (ab_mover%ionmov==21)then
         call pred_steepdesc(ab_mover,forstr_precon,hist,itime,DEBUG,iexit)
       end if

     end do

!    Write MOLDYN netcdf and POSABIN files (done every dtset%nctime time step)
     if (scfcv_args%dtset%nctime>0) then
       jj=itime;if (hist_prev%mxhist>0.and.ab_mover%restartxf==-1) jj=jj-hist_prev%mxhist
       if (jj>0) then
         option=3
         call wrt_moldyn_netcdf(amass,scfcv_args%dtset,jj,option,dtfil%fnameabo_moldyn,&
&         scfcv_args%ab_scfcv_inout%mpi_enreg,scfcv_args%ab_scfcv_inout%results_gs,&
&         hist%histR(:,:,hist%ihist-1),dtfil%unpos,hist%histV(:,:,hist%ihist),&
&         hist%histXF(:,:,1,hist%ihist-1),hist%histXF(:,:,2,hist%ihist-1))
       end if
       if (iexit==1) hist%ihist=hist%ihist-1
     end if

     if(iexit/=0) exit

!    write(std_out,*) 'mover 17'
!    ###########################################################
!    ### 18. Use the history  to extract the new values
!    ###     acell, rprimd, xcart and xred

!    !    DEBUG (XRA BEFORE USE OF HISTORY)
!    if(DEBUG)then
!    write (std_out,*) '---XRA BEFORE USE OF HISTORY---',hist%ihist
!    write (std_out,*) 'XCART'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xcart(:,kk)
!    end do
!    write (std_out,*) 'XRED'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xred(:,kk)
!    end do
!    if (ab_mover%ionmov==1)then
!    write (std_out,*) 'VEL'
!    do kk=1,ab_mover%natom
!    write (std_out,*) hist%histV(:,kk,hist%ihist)
!    end do
!    end if
!    write(std_out,*) 'RPRIMD'
!    do kk=1,3
!    write(std_out,*) rprimd(:,kk)
!    end do
!    write(std_out,*) 'ACELL'
!    write(std_out,*) acell(:)
!    end if

     call hist2var(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,DEBUG)

     if(ab_mover%optcell/=0)then

       call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!      If metric has changed since the initialization,
!      update the Ylm's
       if (scfcv_args%ab_scfcv_in%psps%useylm==1)then
         option=0;
         if (scfcv_args%dtset%iscf>0) option=1
         call initylmg(gprimd,&
&         scfcv_args%ab_scfcv_in%kg,&
&         scfcv_args%dtset%kptns,&
&         scfcv_args%dtset%mkmem,&
&         scfcv_args%ab_scfcv_inout%mpi_enreg,&
&         scfcv_args%ab_scfcv_in%psps%mpsang,&
&         scfcv_args%dtset%mpw,&
&         scfcv_args%dtset%nband,&
&         scfcv_args%dtset%nkpt,&
&         scfcv_args%ab_scfcv_in%npwarr,&
&         scfcv_args%dtset%nsppol,&
&         option,rprimd,dtfil%unkg,&
&         dtfil%unylm,&
&         scfcv_args%ab_scfcv_in%ylm,&
&         scfcv_args%ab_scfcv_in%ylmgr)
       end if

     end if

     vel(:,:)=hist%histV(:,:,hist%ihist)

!    !    DEBUG (XRA AFTER USE OF HISTORY)
!    if(DEBUG)then
!    write (std_out,*) '---XRA AFTER USE OF HISTORY---',hist%ihist
!    write (std_out,*) 'XCART'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xcart(:,kk)
!    end do
!    write (std_out,*) 'XRED'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xred(:,kk)
!    end do
!    if (ab_mover%ionmov==1)then
!    write (std_out,*) 'VEL'
!    do kk=1,ab_mover%natom
!    write (std_out,*) hist%histV(:,kk,hist%ihist)
!    end do
!    end if
!    write(std_out,*) 'RPRIMD'
!    do kk=1,3
!    write(std_out,*) rprimd(:,kk)
!    end do
!    write(std_out,*) 'ACELL'
!    write(std_out,*) acell(:)
!    end if

!    This is needed for some compilers such as
!    pathscale, g95, xlf that do not exit
!    from a loop if you change the upper limit
!    inside
     if (icycle>=ncycle .and. scfcv_args%ab_scfcv_inout%mpi_enreg%me == 0) then
!      write(std_out,*) 'exit ICYCLE',icycle
!      write(std_out,*) 'exit NCYCLE',ncycle
       exit
     end if

!    write(std_out,*) 'mover 18'
!    ###########################################################
!    ### 19. End loop icycle
   end do ! do icycle=1,ncycle

   if(iexit/=0)exit
!  write(std_out,*) 'mover 19'
!  ###########################################################
!  ### 20. End loop itime
 end do ! do itime=1,ntime

!write(std_out,*) 'mover 20'
!###########################################################
!### 21. Set the final values of xcart and xred with the last
!###     computed values (not the last predicted)

 call hist2var(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,DEBUG)

 vel(:,:)=hist%histV(:,:,hist%ihist)

 if(DEBUG)then
   if (ab_mover%ionmov==1)then
     write (std_out,*) 'vel'
     do kk=1,ab_mover%natom
       write (std_out,*) hist%histV(:,kk,hist%ihist)
     end do
   end if
 end if

!write(std_out,*) 'mover 21'
!###########################################################
!### 22. XML Output at the end

!XML output of the status
 if (scfcv_args%ab_scfcv_inout%mpi_enreg%me == 0 .and. scfcv_args%dtset%prtxml == 1) then
   write(ab_xml_out, "(3a)") '    <geometryMinimisation type="',&
&   trim(type4xml),'">'
   write(ab_xml_out, "(5a)") '      <status cvState="',&
&   trim(stat4xml) ,&
&   '" stop-criterion="',trim(crit4xml),'" />'
   write(ab_xml_out, "(3a)") '    </geometryMinimisation>'
 end if

!write(std_out,*) 'mover 22'
!###########################################################
!### 23. Deallocate hist and ab_mover datatypes

 if (ab_mover%goprecon>0)then
   call prec_simple(ab_mover,forstr_precon,hist,1,1,1)
 end if

 if (ab_mover%ionmov==13)then
   ABI_DEALLOCATE(mttk_vars%glogs)
   ABI_DEALLOCATE(mttk_vars%vlogs)
   ABI_DEALLOCATE(mttk_vars%xlogs)
 end if

 call ab_movehistory_end(hist)
 call ab_movehistory_end(hist_prev)
 call ab_movetype_end(ab_mover)

 ABI_DEALLOCATE(forstr_precon%fcart)
 ABI_DEALLOCATE(forstr_precon%fred)

end subroutine mover
!!***
