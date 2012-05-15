!{\src2tex{textfont=tt}}
!!****f* ABINIT/gstate
!! NAME
!! gstate
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations by CG minimization.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, JYR, MKV, MT, FJ, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  codvsn=code version
!!  cpui=initial CPU time
!!
!! OUTPUT
!!  npwtot(nkpt) = total number of plane waves at each k point
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!
!! SIDE EFFECTS
!!  acell(3)=unit cell length scales (bohr)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | mband =maximum number of bands (IN)
!!   | mgfft =maximum single fft dimension (IN)
!!   | mkmem =maximum number of k points which can fit in core memory (IN)
!!   | mpw   =maximum number of planewaves in basis sphere (large number) (IN)
!!   | natom =number of atoms in unit cell (IN)
!!   | nfft  =(effective) number of FFT grid points (for this processor) (IN)
!!   | nkpt  =number of k points (IN)
!!   | nspden=number of spin-density components (IN)
!!   | nsppol=number of channels for spin-polarization (1 or 2) (IN)
!!   | nsym  =number of symmetry elements in space group
!!  iexit= exit flag
!!  initialized= 0 for the first GS calculation (not initialized), else 1
!!  mpi_enreg=MPI-parallelisation information (some already initialized,
!!   some others to be initialized here)
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in gstate, a significant part of
!!   psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,
!!     ntypat,n1xccc,usepaw,useylm, and the arrays dimensioned to npsp
!!   All the remaining components of psps are to be initialized in the call
!!   to pspini .
!!   The next time the code enters gstate, psps might be identical to the
!!   one of the previous dtset, in which case, no reinitialisation is scheduled
!!   in pspini.f .
!!  rprim(3,3)=dimensionless real space primitive translations
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  vel(3,natom)=value of velocity
!!  xred(3,natom) = reduced atomic coordinates
!!
!! NOTES
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
!! In case of wavelets:
!! --------------------
!!    - Only the usual FFT grid (defined by wvl_crmult) is used.
!!      It is defined by nfft, ngfft, mgfft, ... This is strictly not
!!      an FFT grid since its dimensions are not suited for FFTs. They are
!!      defined by wvl_setngfft().
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! TODO
!! Not yet possible to use restartxf in parallel when localrdwf==0
!!
!! PARENTS
!!      gstateimg
!!
!! CHILDREN
!!      alloc_hamilt_gpu,bstruct_clean,bstruct_init,chkexi,cleanrec,clnmpi_atom
!!      clnmpi_bandfft,clnmpi_fft,clnmpi_gs,clnup1,clnup2,create_ddb_blk
!!      dealloc_hamilt_gpu,destroy_ddb_blk,destroy_efield
!!      destroy_electronpositron,destroy_sc_dmft,energies_init,fixsym,fourdp
!!      getph,handle_ncerr,hdr_clean,hdr_init,hdr_update,init_electronpositron
!!      init_pawfgr,init_sc_dmft,init_scf_history,initberry,initmpi_atom
!!      initmpi_fft,initmpi_gs,initrec,initrhoij,initro,initylmg,inwffil,ioarr
!!      ioddb8_out,jellium,kpgio,leave_new,mkradim,mkrho,mover,newocc
!!      nullify_ddb_blk,nullify_efield,outqmc,outwf,outxfhist,pawinit
!!      pawpuxinit,pawuj_drive,prep_kpgio,print_sc_dmft,prtene,psddb8
!!      psolver_kernel,pspini,reportgap,rhoij_copy,rhoij_free,scfcv_init
!!      scfcv_init2,scfcv_new2,scphon,setsym,setsymrhoij,setup1,setup2,status
!!      timab,transgrid,wffclose,wffdelete,wffopen,wffreadskiprec,write_blok8
!!      wrtout,wvl_descr_atoms_set,wvl_descr_free,wvl_descr_psp_set,wvl_mkrho
!!      wvl_projectors_free,wvl_projectors_set,wvl_setboxgeometry,wvl_setngfft
!!      wvl_timing,wvl_wfs_free,wvl_wfs_lr_copy,wvl_wfs_set,xcomm_world
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine gstate(acell,codvsn,cpui,dtfil,dtset,iexit,initialized,&
&                 mpi_enreg,npwtot,occ,pawang,pawrad,pawtab,&
&                 psps,results_gs,rprim,scf_history,vel,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_scftypes
 use defs_wvltypes
 use m_energies, only : energies_type, energies_init
 use m_results_gs , only : results_gs_type
 use m_scf_history
 use defs_parameters
 use defs_rectypes
 use defs_mover
 use m_xmpi
 use m_errors
 use m_paw_toolbox
 use m_wffile
 use m_rec
 use m_efield
 use defs_scfcvargs

#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

 use m_paw_dmft,         only : init_sc_dmft,destroy_sc_dmft,print_sc_dmft,paw_dmft_type
 use m_electronpositron, only : electronpositron_type,init_electronpositron,destroy_electronpositron, &
&                               electronpositron_calctype
 use m_header,           only : hdr_init, hdr_clean
 use m_ebands,           only : bstruct_init, bstruct_clean, ReportGap
 use m_ddb_blk

#if defined HAVE_DFT_BIGDFT
 use BigDFT_API,         only : wvl_timing => timing
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gstate'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_43_wvl_wrappers
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
#if defined HAVE_GPU_CUDA
 use interfaces_52_manage_cuda
#endif
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_57_iovars
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_62_poisson
 use interfaces_62_wvl_wfs
 use interfaces_65_psp
 use interfaces_66_paw
 use interfaces_66_wfs
 use interfaces_67_common
 use interfaces_72_response
 use interfaces_79_seqpar_mpi
 use interfaces_95_drive, except_this_one => gstate
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: iexit,initialized
 real(dp),intent(in) :: cpui
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(inout) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
 type(results_gs_type),intent(inout) :: results_gs
 type(scf_history_type),intent(inout) :: scf_history
!arrays
 integer,intent(out) :: npwtot(dtset%nkpt)
 real(dp),intent(inout) :: acell(3),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: rprim(3,3),vel(3,dtset%natom),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!Define file format for different type of files. Presently,
!only one file format is supported for each type of files, but this might
!change soon ...
!2   for wavefunction file, new format (version 2.0 and after)    (fform)   NOT USED
!52  for density rho(r)       (fformr)
!102 for potential V(r) file. (fformv)  NOT USED
!scalars
 integer,parameter :: formeig=0,level=101,response=0
 integer,save :: nsym_old=-1
 integer :: accessfil,ask_accurate,bantot,choice,fformr=52,fullinit
 integer :: gscase,iapp,iatom,idir,ierr,ii,indx,jj,kk,ios,ir,itypat,comm
 integer :: ixfh,izero,master,mcg,me,mgfftf,mpert,msize,mu,my_nspinor,nblok
 integer :: ncerr,ncid_hdr,nfftf,nfftot,nproc,nspden_rhoij
 integer :: openexit,option,optorth,psp_gencond
 integer :: pwind_alloc,rdwr,rdwrpaw,spaceworld,tim_mkrho
 integer :: vrsddb,ndtpawuj=0
 real(dp) :: cpus,diecut_eff,ecore,ecut_eff,ecutdg_eff,etot,fermie
 real(dp) :: gsqcut_eff,gsqcutc_eff,residm,tolwfr,ucvol
 logical :: read_wf_or_den,has_to_init
 character(len=500) :: message
 character(len=fnlen) :: ddbnm,dscrpt,filnam
!character(len=60) :: method
 real(dp) :: fatvshift
 type(bandstructure_type) :: bstruct,tmp_Bst
 type(efield_type) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(hdr_type) :: hdr
 type(macro_uj_type) :: dtpawuj(0) ! ndtpawuj=0
 type(paw_dmft_type) :: paw_dmft
 type(pawfgr_type) :: pawfgr
 type(recursion_type) ::rec_set
 type(wffile_type) :: wff1,wffnew,wffnow
 type(wvl_data) :: wvl
 type(ab_scfcv_args_in) :: ab_scfcv_in
 type(ab_scfcv_args_inout) :: ab_scfcv_inout
 type(ab_xfh_type) :: ab_xfh
 type(ab_scfcvargs) :: scfcv_args

 type(ddb_blk_type), pointer :: ddb_blk

!arrays
 integer,save :: paw_gencond(6)=(/-1,-1,-1,-1,-1,-1/)
 integer :: ngfft(18),ngfftf(18)
 integer,allocatable :: atindx(:),atindx1(:),indsym(:,:,:)
 integer,allocatable :: irrzon(:,:,:),kg(:,:),nattyp(:),symrec(:,:,:)
 integer,allocatable,target :: npwarr(:)
 integer,pointer :: npwarr_(:),pwind(:,:,:)
!real(dp) :: xcart(3,dtset%natom)
 real(dp) :: gmet(3,3),gprimd(3,3)
 real(dp) :: rmet(3,3),rprimd(3,3),rprimd_orig(3,3),tsec(2)
 real(dp),allocatable :: amass(:),cg(:,:),doccde(:)
 real(dp),allocatable :: eigen(:),ph1df(:,:),phnons(:,:,:),resid(:),rhowfg(:,:)
 real(dp),allocatable :: rhowfr(:,:),spinat_dum(:,:),start(:,:),work(:)
 real(dp),allocatable :: ylm(:,:),ylmgr(:,:,:)
 real(dp),pointer :: kernel_dummy(:),pwnsfac(:,:),rhog(:,:),rhor(:,:)
 real(dp),pointer :: taug(:,:),taur(:,:),xred_old(:,:)
 type(pawrhoij_type),pointer :: pawrhoij(:)

! ***********************************************************************

 DBG_ENTER("COLL")

 call timab(32,1,tsec)
 call timab(33,3,tsec)

!###########################################################
!### 01. Initializations XML, MPI, WVL, etc

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Init MPI data
 master =0
 call xcomm_world(mpi_enreg,spaceworld,myrank=me,mysize=nproc)

!Set up MPI informations from the dataset
 if (dtset%usewvl == 0) then
   call initmpi_gs(dtset,mpi_enreg)
   call initmpi_fft(dtset,mpi_enreg)
   call initmpi_atom(dtset,mpi_enreg)
 else
!  Some checks.
   if (dtset%npsp /= dtset%ntypat) then
     write(message, '(a,a,a,a,I0,a,I0,a,a,a)' ) ch10,&
&     ' wvl_wfs_set:  consistency checks failed,', ch10, &
&     '  dtset%npsp (', dtset%npsp, ') /= dtset%ntypat (', dtset%ntypat, ').', ch10, &
&     '  No alchemy pseudo are allowed with wavelets.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
!  Some default values, to be set-up elsewhere.
   wvl%descr%h(:)                 = dtset%wvl_hgrid
#if defined HAVE_DFT_BIGDFT
   wvl%descr%orthpar%directDiag   = .true.
   wvl%descr%orthpar%norbpInguess = 5
   wvl%descr%orthpar%bsLow        = 300
   wvl%descr%orthpar%bsUp         = 800
   wvl%descr%orthpar%methOrtho    = 0
   wvl%descr%orthpar%iguessTol    = 1.d-4
#endif
!  We set the atom-related internal wvl variables.
   call wvl_descr_atoms_set(acell, dtset%icoulomb, dtset%natom, &
&   dtset%ntypat, dtset%typat, wvl%descr)
 end if

 if (me == 0 .and. dtset%prtxml == 1) then
!  gstate() will handle a dataset, so we output the dataSet markup.
   write(ab_xml_out, "(A)") '  <dataSet>'
!  We output the variables of the dataset given in argument.
!  call outvarsXML()
 end if

!Define FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 call init_pawfgr(dtset,pawfgr,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)

!If dtset%accesswff == 2 set all array outputs to netcdf format
 accessfil = 0
 if (dtset%accesswff == IO_MODE_MPI) then
   accessfil = 4
 end if
 if (dtset%accesswff == IO_MODE_NETCDF) then
   accessfil = 1
 end if
 if (dtset%accesswff == IO_MODE_ETSF) then
   accessfil = 3
 end if

!Structured debugging if prtvol==-level
 if(dtset%prtvol==-level)then
   write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&   ' gstate : enter , debug mode '
   call wrtout(std_out,message,'COLL')
 end if

!###########################################################
!### 02. Calls setup1, kpgio, initylmg

 call status(0,dtfil%filstat,iexit,level,'call setup1   ')

 ecore=zero
 results_gs%grewtn(:,:)=zero
 call energies_init(results_gs%energies)
 results_gs%pel(1:3)   =zero

!Set up for iterations
 ABI_ALLOCATE(amass,(dtset%natom))
 call setup1(acell,amass,bantot,dtset,&
& ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,&
& dtset%natom,ngfftf,ngfft,dtset%nkpt,dtset%nsppol,&
& response,rmet,rprim,rprimd,ucvol,psps%usepaw)

 ABI_ALLOCATE(npwarr,(dtset%nkpt))
 if (dtset%usewvl == 0 .and. dtset%tfkinfunc /= 2) then
   call status(0,dtfil%filstat,iexit,level,'call kpgio    ')
!  Set up the basis sphere of planewaves
   ABI_ALLOCATE(kg,(3,dtset%mpw*dtset%mkmem))
   if (mpi_enreg%mode_para=='b') then
     call prep_kpgio(dtset%accesswff,ecut_eff,dtset%exchn2n3d,gmet,&
&     dtset%istwfk,kg,dtset%kptns,dtfil%fnametmp_kg,dtset%mgfft,&
&     dtset%mkmem,'PERS',mpi_enreg,dtset%mpw,dtset%nband,dtset%nkpt,&
&     npwarr,npwtot,dtset%nsppol,dtfil%unkg)
   else
     call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg,dtfil%fnametmp_kg, &
&     dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,&
&     dtset%mpw,npwarr,npwtot,dtset%nsppol,dtfil%unkg)
   end if
 else
   npwarr(:) = 0
   npwtot(:) = 0
 end if

!Set up the Ylm for each k point
 if ( dtset%tfkinfunc /= 2) then
   ABI_ALLOCATE(ylm,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
   ABI_ALLOCATE(ylmgr,(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm))
   if (psps%useylm==1) then
     if(dtset%mkmem==0) open(dtfil%unylm,file=dtfil%fnametmp_ylm,form='unformatted',status='unknown')
     call status(0,dtfil%filstat,iexit,level,'call initylmg ')
     option=0
     if (dtset%prtstm==0.and.dtset%iscf>0.and.dtset%positron/=1) option=1
     if (abs(dtset%berryopt)==5) option=1
     call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,&
&     psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&     npwarr,dtset%nsppol,option,rprimd,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
   end if
 end if

!SCF history management (allocate it at first call)
 has_to_init=(initialized==0.or.scf_history%history_size<0)
 if (initialized==0) then
!  This call has to be done before any use of SCF history
   call init_scf_history(dtset,mpi_enreg,scf_history)
 end if

 call timab(33,2,tsec)
 call timab(701,3,tsec)

!###########################################################
!### 03. Calls pspini

!Open and read pseudopotential files
 call status(0,dtfil%filstat,iexit,level,'call pspini   ')

 call pspini(dtset,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,level,pawrad,pawtab,psps,rprimd)

 call timab(701,2,tsec)
 call timab(33,3,tsec)

!In case of isolated computations, ecore must set to zero
!because its contribution is counted in the ewald energy
!as the ion-ion interaction.
 if (dtset%icoulomb == 1) ecore = zero

!WVL - Now that psp data are available, we compute rprimd, acell...
!from the atomic positions.
 if (dtset%usewvl == 1) then
   call wvl_descr_psp_set(dtset%nsppol, psps, wvl%descr)
   call wvl_setBoxGeometry(me, dtset%prtvol, psps%gth_params%radii_cf, rprimd, xred, &
&   wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult)
   call mkradim(acell,rprim,rprimd)
!  rprim(:, :)    = reshape((/ &
!  &   real(1., dp), real(0., dp), real(0., dp), &
!  &   real(0., dp), real(1., dp), real(0., dp), &
!  &   real(0., dp), real(0., dp), real(1., dp) /), (/ 3, 3 /))
   call wvl_setngfft(dtset%ixc, dtset%mgfft, mpi_enreg, dtset%natom, dtset%nfft, &
&   dtset%ngfft, dtset%nsppol, psps, rprimd, wvl%descr, dtset%wvl_crmult, &
&   dtset%wvl_frmult, xred)
   nfftf          = dtset%nfft
   mgfftf         = dtset%mgfft
   ngfftf(:)      = dtset%ngfft(:)
 end if

!Initialize band structure datatype
 ABI_ALLOCATE(doccde,(bantot))
 ABI_ALLOCATE(eigen,(bantot))
 doccde(:)=zero ; eigen(:)=zero
 if (dtset%paral_kgb/=0) then     !  We decide to store total npw in bstruct,
   ABI_ALLOCATE(npwarr_,(dtset%nkpt))
   npwarr_(:)=npwarr(:)
   call xsum_mpi(npwarr_,mpi_enreg%commcart,ierr)
 else
   npwarr_ => npwarr
 end if
 call bstruct_init(bantot,bstruct,dtset%nelect,doccde,eigen,dtset%istwfk,dtset%kptns,&
& dtset%nband,dtset%nkpt,npwarr_,dtset%nsppol,dtset%nspinor,dtset%tphysel,&
& dtset%tsmear,dtset%occopt,occ,dtset%wtk)
 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(eigen)
 if (dtset%paral_kgb/=0)  then
   ABI_DEALLOCATE(npwarr_)
 end if
 nullify(npwarr_)

!Initialize PAW atomic occupancies
 if (scf_history%history_size>=0) then
   pawrhoij => scf_history%pawrhoij_last
 else
   ABI_ALLOCATE(pawrhoij,(mpi_enreg%natom*psps%usepaw))
 end if
 if (psps%usepaw==1.and.has_to_init) then
   nspden_rhoij=dtset%nspden;if (dtset%pawspnorb>0.and.dtset%nspinor==2) nspden_rhoij=4
   call initrhoij(dtset%pawcpxocc,psps%indlmn,dtset%lexexch,psps%lmnmax,&
&   dtset%lpawu,mpi_enreg,dtset%natom,mpi_enreg%natom,nspden_rhoij,&
&   dtset%nspinor,dtset%nsppol,dtset%ntypat,pawrhoij,pawtab,dtset%spinat,dtset%typat)
 end if

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps,wvl%descr)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 etot=hdr%etot ; fermie=hdr%fermie ; residm=hdr%residm
 call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
& residm,rprimd,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)

!Clean band structure datatype (should use it more in the future !)
 call bstruct_clean(bstruct)

!###########################################################
!### 04. Calls inwffil

 call status(0,dtfil%filstat,iexit,level,'call inwffil  ')

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
 mcg=dtset%mpw*my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
 ABI_ALLOCATE(cg,(2,mcg))
 ABI_ALLOCATE(eigen,(dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_ALLOCATE(resid,(dtset%mband*dtset%nkpt*dtset%nsppol))
 eigen(:)=zero ; resid(:)=zero
!mpi_enreg%paralbd=0 ; ask_accurate=0
 ask_accurate=0
!WVL - Branching, allocating wavefunctions as wavelets.
 if (dtset%usewvl == 1) then
   call xcomm_world(mpi_enreg,comm,myrank=me,mysize=nproc)
   call wvl_wfs_lr_copy(wvl%wfs, wvl%descr)
!  Create access arrays for wavefunctions and allocate wvl%wfs%psi (other arrays
!  are left unallocated).
   call wvl_wfs_set(dtset%fixmom, dtset%kpt, me, dtset%natom, sum(dtset%nband), &
&   dtset%nkpt, nproc, dtset%nspinor, dtset%nsppol, dtset%nwfshist, dtset%occ_orig, &
&   psps, rprimd, wvl%wfs, dtset%wtk, wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult, &
&   xred)
!  We transfer wavelets informations to the hdr structure.
#if defined HAVE_DFT_BIGDFT
   hdr%nwvlarr(1) = wvl%wfs%Glr%wfd%nvctr_c
   hdr%nwvlarr(2) = 7 * wvl%wfs%Glr%wfd%nvctr_f
#endif
!  Create access arrays for projectors and allocate them.
!  Compute projectors from each atom.
   call wvl_projectors_set(me, dtset%natom, wvl%projectors, psps, rprimd, &
&   wvl%wfs, wvl%descr, dtset%wvl_frmult, xred)
 end if

 read_wf_or_den=(dtset%iscf<=0.or.dtfil%ireadwf/=0.or.(dtfil%ireadden/=0.and.dtset%positron<=0))
 read_wf_or_den=(read_wf_or_den.and.has_to_init)

!RECURSION -  initialization
 if(has_to_init .and. dtset%userec==1) then
   call InitRec(dtset,mpi_enreg,rec_set,rmet,maxval(psps%indlmn(3,:,:)))
 end if

!Initialize wavefunctions.
 if(dtset%tfkinfunc /=2) then
   wff1%unwff=dtfil%unwff1
   optorth=1   !if (psps%usepaw==1) optorth=0
   if(psps%usepaw==1 .and. dtfil%ireadwf==1)optorth=0
   call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen,&
&   dtset%exchn2n3d,formeig,gmet,hdr,dtfil%ireadwf,dtset%istwfk,kg,&
&   dtset%kptns,dtset%localrdwf,dtset%mband,mcg,dtset%mkmem,mpi_enreg,&
&   dtset%mpw,dtset%nband,ngfft,dtset%nkpt,npwarr,&
&   dtset%nsppol,dtset%nsym,occ,optorth,rprimd,dtset%symafm,&
&   dtset%symrel,dtset%tnons,dtfil%unkg,wff1,wffnow,dtfil%unwff1,&
&   dtfil%unwft1,dtfil%fnamewffk,dtfil%fnametmp_wf1,wvl)

 end if

 if (psps%usepaw==1.and.dtfil%ireadwf==1)then
   call rhoij_copy(hdr%pawrhoij,pawrhoij,mpi_enreg=mpi_enreg)
!  Has to update header again (because pawrhoij has changed)
!  MT 2007-10-22: Why ?
   call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
&   residm,rprimd,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)
 end if

!DEBUG
!write(std_out,*)' gstate : stop after inwffil for test memory leak '
!call hdr_clean(hdr)
!return
!ENDDEBUG

!###########################################################
!### 05. Operations related to restartxf (Old version)

!Initialize xf history (should be put in inwffil)
 ab_xfh%nxfh=0
 if(dtset%restartxf>=1 .and. dtfil%ireadwf==1)then

!  Should exchange the data about history in parallel localrdwf==0
   if(mpi_enreg%paral_compil_kpt==1 .and. dtset%localrdwf==0)then
     write(message, '(a,a,a,a,a,a)' )ch10,&
&     ' gstate : BUG -',ch10,&
&     '  It is not yet possible to use non-zero restartxf,',ch10,&
&     '  in parallel, when localrdwf=0. Sorry for this ...'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  DEBUG
!  write(std_out,*)'gstate before outxfhist'
!  END DEBUG

   ABI_ALLOCATE(ab_xfh%xfhist,(3,dtset%natom+4,2,0))
   call outxfhist(ab_xfh,dtset%natom,2,wff1,ios)
   ABI_DEALLOCATE(ab_xfh%xfhist)

!  DEBUG
!  write(std_out,*)'gstate after outxfhist'
!  END DEBUG

   if(ios>0)then
     write(message, '(a,a,a,a,a,a)' )ch10,&
&     ' gstate : BUG -',ch10,&
&     '  An error occurred reading the input wavefunction file,',ch10,&
&     '  with restartxf=1.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   else if(ios==0)then
     write(message, '(a,a,i4,a)' )ch10,&
&     ' gstate : reading',ab_xfh%nxfh,' (x,f) history pairs from input wf file.'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
   end if
!  WARNING : should check that restartxf is not negative
!  WARNING : should check that restartxf /= only when dtfil%ireadwf is activated
 end if

!Allocate the xf history array : takes into account the existing
!pairs, minus those that will be discarded, then those that will
!be computed, governed by dtset%ntime, and some additional pairs
!(needed when it will be possible to use xfhist for move.f)
 ab_xfh%mxfh=(ab_xfh%nxfh-dtset%restartxf+1)+dtset%ntime+5
 ABI_ALLOCATE(ab_xfh%xfhist,(3,dtset%natom+4,2,ab_xfh%mxfh))
!WARNING : should check that the number of atoms in the wf file and natom are the same

!Initialize the xf history array
 if(ab_xfh%nxfh>=dtset%restartxf .and. ab_xfh%nxfh>0)then
!  Eventually skip some of the previous history
   if(dtset%restartxf>=2)then
     do ixfh=1,dtset%restartxf-1
       call WffReadSkipRec(ios,1,wff1)
     end do
   end if

!  Read and store the relevant history
   ab_xfh%nxfhr=ab_xfh%nxfh-dtset%restartxf+1
   call outxfhist(ab_xfh,dtset%natom,3,wff1,ios)
 end if

!Close wff1, if it was ever opened (in inwffil)
 if (dtfil%ireadwf==1) then
   call WffClose(wff1,ierr)
 end if

!Initialize second wavefunction file if needed
 if(dtset%mkmem==0 .and. dtset%nstep/=0) then
   write(message, '(a,i4,a,a)' )&
&   ' gstate about to open unit',dtfil%unwft2,' for file=',trim(dtfil%fnametmp_wf2)
   call wrtout(std_out,message,'PERS')

#if defined HAVE_TRIO_NETCDF
   if(dtset%accesswff==IO_MODE_NETCDF) then
!    Create empty netCDF file
     ncerr = nf90_create(path=trim(dtfil%fnametmp_wf2), cmode=NF90_CLOBBER, ncid=ncid_hdr)
     call handle_ncerr(ncerr," create netcdf wavefunction file")
     ncerr = nf90_close(ncid_hdr)
     call handle_ncerr(ncerr," close netcdf wavefunction file")
   else if(dtset%accesswff==IO_MODE_ETSF) then
     write (std_out,*) "FIXME: ETSF I/O support in gstate"
   end if
#endif
   call WffOpen(dtset%accesswff,spaceworld,dtfil%fnametmp_wf2,ierr,wffnew,master,me,dtfil%unwft2)
 end if

 call status(0,dtfil%filstat,iexit,level,'call setup2   ')

!###########################################################
!### 06. Calls setup2

!DEBUG
!write(std_out,*)'gstate before setup2'
!END DEBUG

!Further setup
 ABI_ALLOCATE(start,(3,dtset%natom))
 call setup2(dtset,npwtot,start,wvl%wfs,xred)

!Allocation of previous atomic positions
 if (scf_history%history_size>=0) then
   xred_old => scf_history%xred_last
 else
   ABI_ALLOCATE(xred_old,(3,dtset%natom))
 end if
 if (has_to_init) xred_old=xred

!###########################################################
!### 07. Symmetry operations when nsym>1

!Do symmetry stuff only for nsym>1
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 ABI_ALLOCATE(irrzon,(nfftot**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(phnons,(2,nfftot**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(indsym,(4,dtset%nsym,dtset%natom))
 ABI_ALLOCATE(symrec,(3,3,dtset%nsym))
 irrzon(:,:,:)=zero
 phnons(:,:,:)=zero
 indsym(:,:,:)=zero
 symrec(:,:,:)=zero

 if (dtset%nsym>1) then

   call status(0,dtfil%filstat,iexit,level,'call setsym   ')

   call setsym(indsym,irrzon,dtset%iscf,dtset%natom,&
&   nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
&   phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

!  Make sure dtset%iatfix does not break symmetry
   call status(0,dtfil%filstat,iexit,level,'call fixsym   ')

   call fixsym(dtset%iatfix,indsym,dtset%natom,dtset%nsym)

 else

!  The symrec array is used by initberry even in case nsym = 1
   symrec(:,:,1) = 0
   symrec(1,1,1) = 1 ; symrec(2,2,1) = 1 ; symrec(3,3,1) = 1

 end if

!Timing for initialisation period
 call timab(33,2,tsec)
 call timab(34,3,tsec)

!###########################################################
!### 08. Compute new occupation numbers

!Compute new occupation numbers, in case wavefunctions and eigenenergies
!were read from disk, occupation scheme is metallic (this excludes iscf=-1),
!and occupation numbers are required by iscf
 if( dtfil%ireadwf==1 .and. &
& (dtset%occopt>=3.and.dtset%occopt<=8) .and. &
& (dtset%iscf>0 .or. dtset%iscf==-3) .and. dtset%positron/=1 ) then

   call status(0,dtfil%filstat,iexit,level,'call newocc   ')
   ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
!  Warning : ideally, results_gs%entropy should not be set up here XG 20011007
!  Warning : ideally, results_gs%e_fermie should not be set up here XG 20011007
!  Do not take into account the possible STM bias
   call newocc(doccde,eigen,results_gs%energies%entropy,&
&   results_gs%energies%e_fermie,&
&   dtset%fixmom,dtset%mband,dtset%nband,&
&   dtset%nelect,dtset%nkpt,dtset%nspinor,dtset%nsppol,occ,&
&   dtset%occopt,dtset%prtvol,zero,dtset%tphysel,dtset%tsmear,dtset%wtk)
   ABI_DEALLOCATE(doccde)

 else
!  Warning : ideally, results_gs%entropy should not be set up here XG 20011007
   results_gs%energies%entropy=zero
 end if

!###########################################################
!### 09. Generate an index table of atoms

!Generate an index table of atoms, in order for them to be used
!type after type.
 ABI_ALLOCATE(atindx,(dtset%natom))
 ABI_ALLOCATE(atindx1,(dtset%natom))
 ABI_ALLOCATE(nattyp,(psps%ntypat))
 indx=1
 do itypat=1,psps%ntypat
   nattyp(itypat)=0
   do iatom=1,dtset%natom
     if(dtset%typat(iatom)==itypat)then
       atindx(iatom)=indx
       atindx1(indx)=iatom
       indx=indx+1
       nattyp(itypat)=nattyp(itypat)+1
     end if
   end do
 end do

!Compute structure factor phases for current atomic pos:
 if ((.not.read_wf_or_den).or.(scf_history%history_size>0.and.has_to_init)) then
   ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*dtset%natom))
   call status(0,dtfil%filstat,iexit,level,'call getph    ')
   call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
 end if

!Here allocation of GPU for vtorho calculations
#if defined HAVE_GPU_CUDA
 if (dtset%use_gpu_cuda==1) then
   call alloc_hamilt_gpu(atindx1,dtset,gprimd,mpi_enreg,nattyp,2,psps,dtset%use_gpu_cuda)
 end if
#endif

!###########################################################
!### 10. PAW related operations

!Initialize paw_dmft, even if neither dmft not paw are used
!write(std_out,*) "dtset%usedmft",dtset%usedmft
 call init_sc_dmft(dtset%dmftbandi,dtset%dmftbandf,dtset%mband,dtset%nband,dtset%nkpt,dtset%nspden, &
& dtset%nspinor,dtset%nsppol,occ,dtset%usedmft,paw_dmft,dtset%usedmft)
!write(std_out,*) "paw_dmft%use_dmft",paw_dmft%use_dmft

!PAW: 1- Initialize values for several arrays unchanged during iterations
!2- Initialize data for LDA+U
!3- Eventually open temporary storage file
 if(psps%usepaw==1) then
!  1-
   if (psp_gencond==1.or.&
&   paw_gencond(1)/=dtset%pawlcutd .or.paw_gencond(2)/=dtset%pawlmix  .or.&
&   paw_gencond(3)/=dtset%pawnphi  .or.paw_gencond(4)/=dtset%pawntheta.or.&
&   paw_gencond(5)/=dtset%pawspnorb.or.paw_gencond(6)/=dtset%pawxcdev) then
     call timab(553,1,tsec)
     diecut_eff=abs(dtset%diecut)*dtset%dilatmx**2
     call pawinit(diecut_eff,psps%indlmn,dtset%pawlcutd,dtset%pawlmix,psps%lmnmax,psps%mpsang,&
&     dtset%pawnphi,dtset%nsym,dtset%pawntheta,psps%ntypat,&
&     pawang,pawrad,dtset%pawspnorb,pawtab,dtset%pawxcdev)
     paw_gencond(1)=dtset%pawlcutd ; paw_gencond(2)=dtset%pawlmix
     paw_gencond(3)=dtset%pawnphi  ; paw_gencond(4)=dtset%pawntheta
     paw_gencond(5)=dtset%pawspnorb; paw_gencond(6)=dtset%pawxcdev
     call timab(553,2,tsec)
   else
     if (pawtab(1)%has_kij  ==1) pawtab(1:psps%ntypat)%has_kij  =2
     if (pawtab(1)%has_nabla==1) pawtab(1:psps%ntypat)%has_nabla=2
   end if
   if (psp_gencond==1.or.nsym_old/=dtset%nsym) then
     call setsymrhoij(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,&
&     rprimd,symrec,pawang%zarot)
     nsym_old=dtset%nsym
   end if
   psps%n1xccc=maxval(pawtab(1:psps%ntypat)%usetcore)
!  2-Initialize and compute data for LDA+U, EXX, or LDA+DMFT
   pawtab(:)%usepawu=0
   pawtab(:)%useexexch=0
   pawtab(:)%exchmix=zero
   if(paw_dmft%use_dmft==1) call print_sc_dmft(paw_dmft,dtset%pawprtvol)
   if (dtset%usepawu>0.or.dtset%useexexch>0.or.paw_dmft%use_dmft>0) then
     call pawpuxinit(dtset%dmatpuopt,dtset%exchmix,dtset%jpawu,dtset%lexexch,dtset%lpawu,&
&     psps%indlmn,psps%lmnmax,psps%ntypat,pawang,dtset%pawprtvol,pawrad,pawtab,dtset%upawu,&
&     dtset%usedmft,dtset%useexexch,dtset%usepawu)
   end if
!  3-Eventually open temporary storage file
   if(dtset%mkmem==0) then
     open(dtfil%unpaw,file=dtfil%fnametmp_paw,form='unformatted',status='unknown')
     rewind(unit=dtfil%unpaw)
   end if
 end if

!DEBUG
!write(std_out,*)'gstate before call of initberry'
!END DEBUG

!###########################################################
!### 11. Initialize (eventually) electron-positron data

!Initialize (eventually) electron-positron data
 nullify (electronpositron)
 if (dtset%positron/=0) then
   call init_electronpositron(dtfil%ireadwf,dtset,electronpositron,mpi_enreg,nfftf,pawrhoij,pawtab)
 end if

!Electric and magnetic field: initialization
 dtefield%has_qijb = 0
 dtefield%has_expibi = 0
 dtefield%has_expibr = 0
 dtefield%has_twdij0 = 0
 dtefield%has_tweijkl = 0
 dtefield%has_Lij = 0
 dtefield%has_Lijr3 = 0
 dtefield%has_rij = 0
 dtefield%usecprj = 0
 call nullify_efield(dtefield)

 if ((dtset%berryopt < 0).or.(dtset%berryopt == 4) &
& .or. (abs(dtset%berryopt)==5)) then
   nullify(pwind,pwnsfac)
   call initberry(dtefield,dtset,gmet,gprimd,kg,&
&   dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,&
&   dtset%natom,dtset%nkpt,npwarr,dtset%nsppol,&
&   dtset%nsym,dtset%ntypat,occ,pawang,pawrad,pawtab,&
&   psps,pwind,pwind_alloc,pwnsfac,rprimd,symrec,&
&   dtset%typat,psps%usepaw,xred)
 else
   pwind_alloc = 1
   ABI_ALLOCATE(pwind,(pwind_alloc,2,3))
   ABI_ALLOCATE(pwnsfac,(2,pwind_alloc))
   pwind(:,:,:)=zero
   pwnsfac(:,:)=zero
 end if

!###########################################################
!### 12. Operations dependent of iscf value

!Get starting charge density : rhor as well as rhog
!Also initialize the kinetic energy density
 if (scf_history%history_size>=0) then
   rhor => scf_history%rhor_last
   taur => scf_history%taur_last
 else
   ABI_ALLOCATE(rhor,(nfftf,dtset%nspden))
   ABI_ALLOCATE(taur,(nfftf,dtset%nspden*dtset%usekden))
 end if
 ABI_ALLOCATE(rhog,(2,nfftf))
 ABI_ALLOCATE(taug,(2,nfftf*dtset%usekden))

 if (has_to_init) then
   if (dtset%iscf>0) then
     if(dtfil%ireadden/=0.and.dtset%positron<=0)then

       rdwr=1;rdwrpaw=psps%usepaw;if(dtfil%ireadwf/=0) rdwrpaw=0
       call ioarr(accessfil,rhor,dtset,results_gs%etotal,fformr,dtfil%fildensin,hdr,&
&       mpi_enreg, nfftf,pawrhoij,rdwr,rdwrpaw,wvl%descr)
       if(dtfil%ireadkden/=0 .and. dtset%usekden==1 )then
         call ioarr(accessfil,taur,dtset,results_gs%etotal,fformr,dtfil%filkdensin,hdr,&
&         mpi_enreg, nfftf,pawrhoij,rdwr,rdwrpaw,wvl%descr)
       end if
       if (rdwrpaw/=0) then
         call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
&         residm,rprimd,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)
       end if
!      Compute up+down rho(G) by fft
       ABI_ALLOCATE(work,(nfftf))
       work(:)=rhor(:,1)
       call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
       if(dtset%usekden==1)then
         work(:)=taur(:,1)
         call fourdp(1,taug,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
       end if
       ABI_DEALLOCATE(work)

     else if(dtfil%ireadwf/=0)then
       izero=0
!      Obtain the charge density from wfs that were read previously
!      Be careful: in PAW, rho does not include the compensation
!      density (to be added in scfcv.F90) !
       call status(0,dtfil%filstat,iexit,level,'call mkrho    ')
!      tim_mkrho=1 ; mpi_enreg%paralbd=0
       tim_mkrho=1
       if (psps%usepaw==1) then
         ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
         ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
!        write(std_out,*) "mkrhogstate"
         call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&         mpi_enreg,npwarr,occ,paw_dmft,phnons,rhowfg,rhowfr,rprimd,tim_mkrho,ucvol,&
&         dtfil%unkg,wffnow,wvl%wfs,wvl%descr)
         call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,rhog,rhowfr,rhor)
         ABI_DEALLOCATE(rhowfg)
         ABI_DEALLOCATE(rhowfr)
       else
         call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&         mpi_enreg,npwarr,occ,paw_dmft,phnons,rhog,rhor,rprimd,tim_mkrho,ucvol,&
&         dtfil%unkg,wffnow,wvl%wfs,wvl%descr)
         if(dtset%usekden==1)then
           call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&           mpi_enreg,npwarr,occ,paw_dmft,phnons,taug,taur,rprimd,tim_mkrho,ucvol,&
&           dtfil%unkg,wffnow,wvl%wfs,wvl%descr,option=1)
         end if

       end if

     else if(dtfil%ireadwf==0.and.dtset%positron/=1)then

!      Crude, but realistic initialisation of the density
!      There is not point to compute it from random wavefunctions
!      except with wavelets.
       call status(0,dtfil%filstat,iexit,level,'call initro   ')
       if (dtset%usewvl == 0) then
         call initro(atindx,dtset%densty,gmet,gsqcut_eff,psps%usepaw,&
&         mgfftf,mpi_enreg,psps%mqgrid_vl,dtset%natom,nattyp,nfftf,&
&         ngfftf,dtset%nspden,psps%ntypat,dtset%paral_kgb,pawtab,ph1df,&
&         psps%qgrid_vl,rhog,rhor,dtset%spinat,ucvol,psps%usepaw,&
&         dtset%ziontypat,dtset%znucl)
!        Update initialized density taking into account jellium slab
         if(dtset%jellslab/=0) then
           option=2
           ABI_ALLOCATE(work,(nfftf))
           call jellium(gmet,gsqcut_eff,mpi_enreg,nfftf,ngfftf,dtset%nspden,&
&           option,dtset%paral_kgb,dtset%slabwsrad,rhog,rhor,rprimd,work,dtset%slabzbeg,dtset%slabzend)
           ABI_DEALLOCATE(work)
         end if ! of usejell
!        Kinetic energy density initialized to zero (used only in metaGGAs ... )
         if(dtset%usekden==1)then
           taur=zero ; taug=zero
         end if
       else if (dtset%usewvl/=0) then
         call wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wvl%wfs, wvl%descr)
       end if

     end if

   else if ((dtset%iscf==-1.or.dtset%iscf==-2.or.dtset%iscf==-3).and.dtset%positron<=0) then

     call status(0,dtfil%filstat,iexit,level,'call ioarr    ')
!    Read rho(r) from a disk file
     rdwr=1;rdwrpaw=psps%usepaw
!    Note : results_gs%etotal is read here,
!    and might serve in the tddft routine, but it is contrary to the
!    intended use of results_gs ...
!    Warning : should check the use of results_gs%e_fermie
!    Warning : should check the use of results_gs%residm
!    One might make them separate variables.

     call ioarr(accessfil,rhor,dtset, results_gs%etotal,fformr,dtfil%fildensin,hdr,&
&     mpi_enreg,nfftf,pawrhoij,rdwr,rdwrpaw,wvl%descr)
     if(dtfil%ireadkden/=0 .and. dtset%usekden==1)then
       call ioarr(accessfil,taur,dtset, results_gs%etotal,fformr,dtfil%filkdensin,hdr,&
&       mpi_enreg,nfftf,pawrhoij,rdwr,rdwrpaw,wvl%descr)
     end if

!    Compute up+down rho(G) by fft
     call status(0,dtfil%filstat,iexit,level,'call fourdp   ')
     ABI_ALLOCATE(work,(nfftf))
     work(:)=rhor(:,1)
     call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
     if(dtset%usekden==1)then
       work(:)=taur(:,1)
       call fourdp(1,taug,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
     end if
     ABI_DEALLOCATE(work)

   else ! Disallowed value for iscf
     write(message,'(a,i12,a)')'  iscf has the disallowed value=',dtset%iscf,'.'
     MSG_BUG(message)
   end if
 end if ! has_to_init

!Debugging : print the different parts of rhor
!MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(dtset%prtvol==-level)then
   write(message,'(a)') '   ir     rhor(ir)     '
   call wrtout(std_out,message,'COLL')
   do ir=1,nfftf
     if(ir<=11 .or. mod(ir,301)==0 )then
       write(message,'(i5,a,es13.6)')ir,' ',rhor(ir,1)
       call wrtout(std_out,message,'COLL')
       if(dtset%nsppol==2)then
         write(message,'(a,es13.6)')'      ',rhor(ir,2)
         call wrtout(std_out,message,'COLL')
       end if
     end if
   end do
 end if

!###########################################################
!### 13. If needed, initialize SCF history variables

!If needed, initialize atomic density in SCF history
 if (scf_history%history_size>0.and.has_to_init) then
!  If rhor is an atomic density, just store it in history
   if (.not.read_wf_or_den) then
     scf_history%atmrho_last(:)=rhor(:,1)
   else
!    If rhor is not an atomic density, has to compute rho_at(r)
     ABI_ALLOCATE(rhowfg,(2,nfftf))
     ABI_ALLOCATE(rhowfr,(nfftf,1))
     ABI_ALLOCATE(spinat_dum,(3,dtset%natom))
     spinat_dum=zero
     call initro(atindx,dtset%densty,gmet,gsqcut_eff,psps%usepaw,mgfftf,mpi_enreg,&
&     psps%mqgrid_vl,dtset%natom,nattyp,nfftf,ngfftf,1,psps%ntypat,dtset%paral_kgb,pawtab,&
&     ph1df,psps%qgrid_vl,rhowfg,rhowfr,spinat_dum,ucvol,&
&     psps%usepaw,dtset%ziontypat,dtset%znucl)
     scf_history%atmrho_last(:)=rhowfr(:,1)
     ABI_DEALLOCATE(rhowfg)
     ABI_DEALLOCATE(rhowfr)
     ABI_DEALLOCATE(spinat_dum)
   end if
 end if

 if ((.not.read_wf_or_den).or.(scf_history%history_size>0.and.has_to_init))  then
   ABI_DEALLOCATE(ph1df)
 end if

 fatvshift=one
 rprimd_orig(:,:)=rprimd

 call status(0,dtfil%filstat,iexit,level,'end gstate(1) ')
#if defined HAVE_DFT_BIGDFT
 if (dtset%usewvl == 1) call wvl_timing(mpi_enreg%me,'INIT','PR')
#endif

 if(dtset%prtvol==-level)then
   write(message,'(a1,a,a1,a,i1,a)') ch10,&
&   ' gstate : before scfcv, move or brdmin ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Check whether exiting was required by the user.
!If found then do not start minimization steps
!At this first call to chkexi, initialize cpus, if it
!is non-zero (which would mean that no action has to be taken)
!Should do this in driver ...
 cpus=dtset%cpus
 if(abs(cpus)>1.0d-5)cpus=cpus+cpui
 openexit=1 ; if(dtset%chkexit==0) openexit=0
 call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
!If immediate exit, and wavefunctions were not read, must zero eigenvalues
 if (iexit/=0) then
   eigen(:)=zero
 end if

 call timab(34,2,tsec)

 if (iexit==0) then

!  ###########################################################
!  ### 14. Move atoms and acell acording to ionmov value


!  Eventually symmetrize atomic coordinates over space group elements:
!  call symzat(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,xred)

!  If atoms are not being moved and U should not be determined,
!  use scfcv directly; else
!  call move, pawuj_drive or brdmin which in turn calls scfcv.

   call timab(35,3,tsec)

   iapp=0

   call scfcv_init(ab_scfcv_in,ab_scfcv_inout,atindx,atindx1,cg,cpus,dtefield,dtfil,&
&   dtpawuj,dtset,ecore,eigen,hdr,iapp,&
&   indsym,initialized,irrzon,kg,mcg,mpi_enreg,nattyp,ndtpawuj,&
&   nfftf,npwarr,occ,pawang,pawfgr,pawrad,&
&   pawrhoij,pawtab,phnons,psps,pwind,pwind_alloc,pwnsfac,&
&   rec_set,resid,results_gs,scf_history,&
&   fatvshift,symrec,taug,taur,wvl,ylm,ylmgr)

   call scfcv_init2(scfcv_args,&
&   ab_scfcv_in,ab_scfcv_inout,&
&   dtset,paw_dmft,wffnew,wffnow)

   write(message,'(a,80a)')ch10,('=',mu=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   if (dtset%ionmov==0) then

!    Should merge this call with the call for dtset%ionmov==4 and 5
     iapp=0
!    mpi_enreg%paralbd=0

     call status(0,dtfil%filstat,iexit,level,'call scfcv    ')

     if (dtset%macro_uj==0) then

       call scfcv_new2(scfcv_args,&
&       electronpositron,&
&       rhog,rhor,&
&       rprimd,&
&       xred,xred_old)

!      call scfcv_new(ab_scfcv_in,ab_scfcv_inout,&
!      &       dtset,electronpositron,&
!      &       paw_dmft,resid,&
!      &       rhog,rhor,&
!      &       rprimd,&
!      &       wffnew,wffnow,&
!      &       xred,xred_old)

     else
!      Conduct determination of U

       call pawuj_drive(ab_scfcv_in,ab_scfcv_inout,&
&       dtset,electronpositron,paw_dmft,&
&       rhog,rhor,&
&       rprimd,&
&       wffnew,wffnow,&
&       xred,xred_old)

     end if

!    ======================================BEGIN===
!    New structure for geometry optimization
!    ==============================================
   else if (dtset%ionmov>50.or.dtset%ionmov<=21) then

     call status(0,dtfil%filstat,iexit,level,'call mover    ')

     call mover(scfcv_args,&
&     ab_xfh,acell,amass,dtfil,&
&     electronpositron,&
&     rhog,rhor,&
&     rprimd,&
&     vel,&
&     xred,xred_old)

!    Compute rprim from rprimd and acell
     do kk=1,3
       do jj=1,3
         rprim(jj,kk)=rprimd(jj,kk)/acell(kk)
       end do
     end do

!    ==============================================
!    New structure for geometry optimization
!    ========================================END===

   else if (dtset%ionmov == 30) then

     call scphon(amass,ab_scfcv_in,ab_scfcv_inout,&
&     dtset,electronpositron,&
&     paw_dmft,&
&     rhog,rhor,&
&     rprimd,&
&     wffnew,wffnow,&
&     xred,xred_old)

   else ! Not an allowed option
     write(message, '(a,i12,a,a)' )&
&     ' Disallowed value for ionmov=',dtset%ionmov,ch10,&
&     ' Allowed values are:1,2,3,4,5,6,7,8,9,10,11,12,13,14,20,21 and 30\n'
     MSG_BUG(message)
   end if

   call timab(35,2,tsec)

!  ###########################################################
!  ### 15. Final operations and output for gstate

!  End of the check of hasty exit
 end if

 call timab(36,3,tsec)

 write(message, '(80a,a,a,a,a)' ) ('=',mu=1,80),ch10,ch10,&
& ' ----iterations are completed or convergence reached----',ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!Mark this GS computation as done
 initialized=1
#if defined HAVE_DFT_BIGDFT
 if (dtset%usewvl == 1) call wvl_timing(mpi_enreg%me,'WFN_OPT','PR')
#endif

!Close the unneeded temporary data files, if any.
!Other files are closed in clnup1.
 if (dtset%mkmem==0) then
   close (unit=dtfil%unkg,status='delete')
   if (psps%useylm==1) close (unit=dtfil%unylm,status='delete')
   if (psps%usepaw==1) close (unit=dtfil%unpaw,status='delete')
   call WffDelete(wffnew,ierr)
 end if

!Will be put here later.
!! ! WVL - maybe compute the tail corrections to energy
!! if (dtset%tl_radius > real(0, dp)) then
!!    ! Store xcart for each atom
!!    allocate(xcart(3, dtset%natom))
!!    call xredxcart(dtset%natom, 1, rprimd, xcart, xred)
!!    ! Use the tails to improve energy precision.
!!    call wvl_tail_corrections(dtset, results_gs%energies, results_gs%etotal, &
!!         & mpi_enreg, occ, psps, vtrial, wvl, xcart)
!!    deallocate(xcart)
!! end if

!Update the header, before using it
!WARNING : There is a problem now (time of writing, 6.0.4, but was in ABINITv5 and had ever been there) to update
!the header with change of rprim. Might be due to the planewave basis set definition.
!Put the original rprimd .
 call hdr_update(bantot,results_gs%etotal,results_gs%energies%e_fermie,hdr,dtset%natom,&
& results_gs%residm,rprimd_orig,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)

!if(dtset%tfkinfunc/=2) then
 call status(0,dtfil%filstat,iexit,level,'call outwf    ')

 if(dtset%nqpt==0)filnam=dtfil%fnameabo_wfk
 if(dtset%nqpt==1)filnam=dtfil%fnameabo_wfq
 call outwf(cg,dtset,eigen,filnam,hdr,kg,dtset%kptns,&
& dtset%mband,mcg,dtset%mkmem,mpi_enreg,dtset%mpw,ab_xfh%mxfh,dtset%natom,&
& dtset%nband,dtset%nkpt,npwarr,dtset%nsppol,dtset%nstep,&
& ab_xfh%nxfh,occ,resid,response,dtfil%unwff2,wffnow,wvl%wfs,wvl%descr,ab_xfh%xfhist)
!end if

 if(dtset%prtwf==2)then
   call outqmc(cg,dtset,eigen,gprimd,hdr,kg,mcg,mpi_enreg,npwarr,occ,psps,results_gs)
 end if

 call status(0,dtfil%filstat,iexit,level,'call clnup1   ')

 call clnup1(acell,dtset,eigen,results_gs%energies%e_fermie,&
& dtfil%fnameabo_dos,dtfil%fnameabo_eig,results_gs%fred,&
& mpi_enreg,nfftf,ngfftf,occ,dtset%optforces,&
& resid,rhor,rprimd,results_gs%vxcavg,xred)

 if ( (dtset%iscf>0 .or. dtset%iscf==-3) .and. dtset%prtstm==0) then
   call status(0,dtfil%filstat,iexit,level,'call prtene   ')
   call prtene(dtset,results_gs%energies,ab_out,psps%usepaw)
 end if

 ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
 doccde=zero

 call bstruct_init(bantot,tmp_Bst,dtset%nelect,doccde,eigen,hdr%istwfk,hdr%kptns,hdr%nband,&
& hdr%nkpt,hdr%npwarr,hdr%nsppol,hdr%nspinor,hdr%tphysel,hdr%tsmear,hdr%occopt,hdr%occ,hdr%wtk)

 tmp_Bst%fermie = results_gs%energies%e_fermie

 ABI_DEALLOCATE(doccde)

 call ReportGap(tmp_BSt,header="Gap info",unit=std_out,mode_paral="COLL")

 call bstruct_clean(tmp_BSt)

!Open the formatted derivative database file, and write the
!preliminary information
!In the // case, only one processor writes the energy and
!the gradients to the DDB

 if ((me==0).and.(dtset%nimage==1).and.((dtset%iscf > 0).or.&
& (dtset%berryopt == -1).or.(dtset%berryopt) == -3)) then

   call status(0,dtfil%filstat,iexit,level,'call ioddb8_ou')
   vrsddb=100401
   dscrpt=' Note : temporary (transfer) database '
   ddbnm=trim(dtfil%filnam_ds(4))//'_DDB'
!  tolwfr must be initialized here, but it is a dummy value
   tolwfr=1.0_dp
   call ioddb8_out (dscrpt,ddbnm,dtset%natom,dtset%mband,&
&   dtset%nkpt,dtset%nsym,psps%ntypat,dtfil%unddb,vrsddb,&
&   acell,dtset%amu,dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&   dtset%intxc,dtset%iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&   dtset%natom,dtset%nband,ngfft,dtset%nkpt,dtset%nspden,dtset%nspinor,&
&   dtset%nsppol,dtset%nsym,psps%ntypat,occ,dtset%occopt,dtset%pawecutdg,&
&   rprim,dtset%sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&   dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&   dtset%typat,dtset%usepaw,dtset%wtk,xred,psps%ziontypat,dtset%znucl)

   if (dtset%iscf > 0) then
     nblok = 2          ! 1st blok = energy, 2nd blok = gradients
   else
     nblok = 1
   end if
   fullinit = 0 ; choice=2
   call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&   psps%lmnmax,nblok,psps%ntypat,dtfil%unddb,pawtab,&
&   psps%pspso,psps%usepaw,psps%useylm,vrsddb)

   mpert = dtset%natom + 6   ; msize = 3*mpert

   ABI_ALLOCATE(ddb_blk,)
   call nullify_ddb_blk(ddb_blk)
!  create a ddb_blk structure with just one blok
   call create_ddb_blk(msize, 1, ddb_blk)

   ddb_blk%flg = 0
   ddb_blk%qpt = zero
   ddb_blk%nrm = one
   ddb_blk%val = zero

!  Write total energy to the DDB
   if (dtset%iscf > 0) then
     ddb_blk%typ(1) = 0
     ddb_blk%val(1,1,1) = results_gs%etotal
     ddb_blk%flg(1,1) = 1
     call write_blok8(ddb_blk,1,choice,dtset%mband,&
&     mpert,msize,dtset%nkpt,dtfil%unddb)
   end if

!  Write gradients to the DDB
   ddb_blk%typ = 4
   ddb_blk%flg = 0
   ddb_blk%val = zero
   indx = 0
   if (dtset%iscf > 0) then
     do iatom = 1, dtset%natom
       do idir = 1, 3
         indx = indx + 1
         ddb_blk%flg(indx,1) = 1
         ddb_blk%val(1,indx,1) = results_gs%fred(idir,iatom)
       end do
     end do
   end if

   indx = 3*dtset%natom + 3
   if ((abs(dtset%berryopt) == 1).or.(abs(dtset%berryopt) == 3)) then
     do idir = 1, 3
       indx = indx + 1
       if (dtset%rfdir(idir) == 1) then
         ddb_blk%flg(indx,1) = 1
         ddb_blk%val(1,indx,1) = results_gs%pel(idir)
       end if
     end do
   end if

   indx = 3*dtset%natom + 6
   if (dtset%iscf > 0) then
     ddb_blk%flg(indx+1:indx+6,1) = 1
     ddb_blk%val(1,indx+1:indx+6,1) = results_gs%strten(1:6)
   end if

   call write_blok8(ddb_blk,1,choice,dtset%mband,&
&   mpert,msize,dtset%nkpt,dtfil%unddb)

   call destroy_ddb_blk(ddb_blk)
   ABI_DEALLOCATE(ddb_blk)

!  Close DDB
   close(dtfil%unddb)

 end if

 if (dtset%nstep>0 .and. dtset%prtstm==0 .and. dtset%positron/=1) then
   call status(0,dtfil%filstat,iexit,level,'call clnup2   ')
   call clnup2(psps%n1xccc,results_gs%fred,results_gs%gresid,&
&   results_gs%grewtn,results_gs%grxc,dtset%iscf,dtset%natom,&
&   dtset%optforces,dtset%optstress,dtset%prtvol,start,&
&   results_gs%strten,results_gs%synlgr,xred)
 end if

 call status(0,dtfil%filstat,iexit,level,'deallocate    ')

!Deallocate arrays
 ABI_DEALLOCATE(amass)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(atindx1)
 ABI_DEALLOCATE(cg)
 ABI_DEALLOCATE(eigen)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(irrzon)
 ABI_DEALLOCATE(npwarr)
 ABI_DEALLOCATE(nattyp)
 ABI_DEALLOCATE(phnons)
 ABI_DEALLOCATE(resid)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(start)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(taug)
 ABI_DEALLOCATE(ab_xfh%xfhist)
 ABI_DEALLOCATE(pawfgr%fintocoa)
 ABI_DEALLOCATE(pawfgr%coatofin)
 if(dtset%tfkinfunc/=2)  then
   ABI_DEALLOCATE(ylm)
   ABI_DEALLOCATE(ylmgr)
 end if
 if (scf_history%history_size<0) then
   if (psps%usepaw==1) then
     call rhoij_free(pawrhoij)
   end if
   ABI_DEALLOCATE(rhor)
   ABI_DEALLOCATE(taur)
   ABI_DEALLOCATE(pawrhoij)
   ABI_DEALLOCATE(xred_old)
 else
   nullify(rhor,taur,pawrhoij,xred_old)
 end if

!RShaltaf: Changed to include SBC
 if (dtset%icoulomb > 0) then
!  Ask to deallocate the kernel part of Poisson's solver
!  Arguments are dummy ones since iaction == 0.
   call psolver_kernel(dtset, 0, kernel_dummy, mpi_enreg, rprimd, wvl%descr)
 end if

!PAW+U
!TODO: this should be replaced by a call to the destroyer function in 66_paw/m_paw_toolbox.F90
 if (dtset%usepawu>0.or.dtset%useexexch>0) then
   do itypat=1,psps%ntypat
     if((dtset%lpawu(itypat)/=-1).or.(dtset%lexexch(itypat)/=-1)) then
       ABI_DEALLOCATE(pawtab(itypat)%lnproju)
       ABI_DEALLOCATE(pawtab(itypat)%phiphjint)
       ABI_DEALLOCATE(pawtab(itypat)%ph0phiint)
       ABI_DEALLOCATE(pawtab(itypat)%zioneff)
     end if
   end do
 end if

!PAW+DMFT
!write(std_out,*) "before destroy_dmft", paw_dmft%use_dmft
 call destroy_sc_dmft(paw_dmft)

!Destroy electronpositron datastructure
 if (dtset%positron/=0) then
   call destroy_electronpositron(electronpositron)
 end if

!Deallocating the basis set.
 if (dtset%usewvl == 1) then
   call wvl_projectors_free(wvl%projectors)
   call wvl_wfs_free(wvl%wfs)
   call wvl_descr_free(wvl%descr)
 else
   if(dtset%tfkinfunc /=2 )then
     ABI_DEALLOCATE(kg)
   end if
 end if

 if ((dtset%berryopt < 0).or.(dtset%berryopt == 4)) then
   ABI_DEALLOCATE(pwind)
   ABI_DEALLOCATE(pwnsfac)
!  deallocate(dtefield%ikpt_dk)
!  deallocate(dtefield%idxkstr)
!  deallocate(dtefield%sflag)
!  deallocate(dtefield%cgindex)
!  deallocate(dtefield%kgindex)
!  deallocate(dtefield%fkgindex)
!  deallocate(dtefield%fkptns)
!  deallocate(dtefield%indkk_f2ibz)
!  deallocate(dtefield%i2fbz)
!  deallocate(dtefield%coord_str)
!  deallocate(dtefield%str_neigh)
!  deallocate(dtefield%strg_neigh)
!  deallocate(dtefield%lmn_size)
   if (mpi_enreg%paral_compil_kpt == 1) then
     ABI_DEALLOCATE(mpi_enreg%kptdstrb)
     if (dtset%berryopt == 4) then
       ABI_DEALLOCATE(mpi_enreg%kptdstrbi)
!      deallocate(dtefield%cgqindex)
!      deallocate(dtefield%nneigh)
       ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
     end if
   end if
   if (associated(mpi_enreg%kpt_loc2ibz_sp))  then
     ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
   end if
   if (associated(mpi_enreg%kpt_loc2fbz_sp)) then
     ABI_DEALLOCATE(mpi_enreg%kpt_loc2fbz_sp)
   end if
   if (associated(mpi_enreg%fmkmem)) then
     ABI_DEALLOCATE(mpi_enreg%fmkmem)
   end if
   if (associated(mpi_enreg%mkmem)) then
     ABI_DEALLOCATE(mpi_enreg%mkmem)
   end if
 else
   ABI_DEALLOCATE(pwind)
   ABI_DEALLOCATE(pwnsfac)
 end if


!if (dtset%berryopt == 4) deallocate(dtefield%smat)
 call destroy_efield(dtefield)

!deallocate Recursion
 if (dtset%userec == 1)  call CleanRec(rec_set)
!Clean the header
 call hdr_clean(hdr)

 if (me == 0 .and. dtset%prtxml == 1) then
!  The dataset given in argument has been treated, then
!  we output its variables.
!  call outvarsXML()
!  gstate() will handle a dataset, so we output the dataSet markup.
   write(ab_xml_out, "(A)") '  </dataSet>'
 end if

!Clean the MPI informations
 if (dtset%usewvl == 0) then
!  Plane-wave case
   call clnmpi_bandfft(mpi_enreg)
   call clnmpi_atom(mpi_enreg)
   call clnmpi_fft(mpi_enreg)
   call clnmpi_gs(mpi_enreg)
 else
!  Wavelet case
   if (associated(mpi_enreg%nscatterarr))  then
     ABI_DEALLOCATE(mpi_enreg%nscatterarr)
   end if
   if (associated(mpi_enreg%ngatherarr ))  then
     ABI_DEALLOCATE(mpi_enreg%ngatherarr)
   end if
 end if

!Eventually clean cuda runtime
#if defined HAVE_GPU_CUDA
 if (dtset%use_gpu_cuda==1) then
   call dealloc_hamilt_gpu(2,dtset%use_gpu_cuda)
 end if
#endif

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(36,2,tsec)
 call timab(32,2,tsec)

 DBG_EXIT("COLL")

end subroutine gstate
!!***
