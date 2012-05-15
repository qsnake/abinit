!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_positron
!! NAME
!! setup_positron
!!
!! FUNCTION
!! Do various initializations for the positron lifetime calculation
!!
!! NOTE
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (GJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ecore=core psp energy (part of total energy) (hartree)
!!  etotal=current value of total energy
!!  forces_needed=if >0 forces are needed
!!  fred(3,natom)=forces in reduced coordinates
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  gmet(3,3)=reciprocal space metric
!!  grewtn(3,natom)=d(Ewald)/d(xred) (hartree)
!!  gsqcut=cutoff value on G**2 for sphere inside fft box
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  ifirst_gs= 0 if we are in a single ground-state calculation
!!     or in the first ground-state calculation of a structural minimization/dynamics
!!  indsym(4,nsym,natom)=index showing transformation of atom labels
!!                       under symmetry operations (computed in symatm)
!!  istep=index of the number of steps in the routine scfcv
!!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
!!  kg(3,mpw*mkmem)=reduced (integer) coordinates of G vecs in basis sphere
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!!  maxfor=maximum absolute value of fcart (forces)
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  n3xccc=dimension of the xccc3d array (0 or nfftf).
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  nhat(nfftf,nspden*usepaw)= -PAW only- compensation density
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  nvresid(nfftf,nspden)=array for the residual of the density/potential
!!  optres=0 if the potential residual has to be used for forces corrections
!!        =1 if the density residual has to be used for forces corrections
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr(natom*usepaw) <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pel(3)=reduced coordinates of the electronic polarization (a. u.)
!!  ph1d(2,3*(2*mgfftf+1)*natom)=one-dimensional structure factor information (fine FFT grid)
!!  ph1dc(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases (coarse FFT grid)
!!  pion(3)=reduced coordinates of the ionic polarization (a. u.)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  stress_needed=if >0 stresses are needed
!!  strsxc(6)=xc correction to stress
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  ucvol=unit cell volume in bohr**3.
!!  vhartr(nfftf)=array for holding Hartree potential
!!  vpsp(nfftf)=array for holding local psp
!!  vxc(nfftf,nspden)=exchange-correlation potential (hartree) in real space
!!  wffnow=unit number for current wf disk file
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  energies <type(energies_type)>=all part of total energy.
!!  cg(2,mcg)=wavefunctions
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  occ(mband*nkpt*nsppol)=occupation number for each band at each k point
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  rhog(2,nfft)=Fourier transform of total electron/positron density
!!  rhor(nfft,nspden)=total electron/positron density (el/bohr**3)
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      cprj_alloc,cprj_copy,cprj_free,energies_copy,energies_init,forstr
!!      fourdp,getcut,hartre,initrhoij,initro,ioarr,leave_new,pawmknhat
!!      rhoij_alloc,rhoij_copy,rhoij_free,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setup_positron(atindx,atindx1,cg,dtefield,dtfil,dtset,ecore,eigen,etotal,electronpositron,&
&          energies,forces_needed,fred,gmet,gprimd,grewtn,gsqcut,hdr,ifirst_gs,indsym,istep,istep_mix,kg,&
&          kxc,maxfor,mcg,mgfft,mpi_enreg,n3xccc,nattyp,nfft,ngfft,nhat,nkxc,npwarr,nvresid,occ,optres,&
&          paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,pel,ph1d,ph1dc,pion,psps,rhog,rhor,&
&          rprimd,stress_needed,strsxc,symrec,ucvol,usexcnhat,vhartr,vpsp,vxc,wffnow,&
&          xccc3d,xred,ylm,ylmgr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_efield
 use m_energies
 use m_wffile
 use m_electronpositron
 use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_positron'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_44_abitypes_defs
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_62_iowfdenpot
 use interfaces_66_paw
 use interfaces_67_common, except_this_one => setup_positron
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: forces_needed,ifirst_gs,istep,mcg,mgfft,n3xccc,nfft,nkxc,optres,stress_needed,usexcnhat
 integer,intent(inout) :: istep_mix
 real(dp),intent(in) :: ecore,etotal,gsqcut,maxfor,ucvol
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(efield_type),intent(in) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(hdr_type),intent(inout) :: hdr
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type), intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),nattyp(dtset%natom),ngfft(18)
 integer,intent(in) :: npwarr(dtset%nkpt),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),grewtn(3,dtset%natom),kxc(nfft,nkxc),pel(3)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom),ph1dc(2,(3*(2*dtset%mgfft+1)*dtset%natom)*dtset%usepaw)
 real(dp),intent(in) :: pion(3),rprimd(3,3),strsxc(6),vhartr(nfft),vpsp(nfft),vxc(nfft,dtset%nspden)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: cg(2,mcg)
 real(dp),intent(inout) :: nhat(nfft,dtset%nspden*dtset%usepaw)
 real(dp),intent(inout) :: nvresid(nfft,dtset%nspden)
 real(dp),intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol),fred(3,dtset%natom)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: rhog(2,nfft),rhor(nfft,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc),xred(3,dtset%natom)
 type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*dtset%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*dtset%usepaw)
 type(pawtab_type),intent(in)  :: pawtab(dtset%ntypat*dtset%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: accessfil,fformr=52,history_level,iatom,iband,icalctype,icalctype0,icg,ifft,ikpt
 integer :: iocc,ireadwf,ispden,isppol,n3xccc0,nocc,optfor,optstr,rdwr,rdwrpaw
 logical,parameter :: always_restart=.false.  ! Set to true to restart by a pure electronic step at each new atomic structure
 logical :: need_scocc,new_calctype
 real(dp) :: boxcut_dum,diffor_dum,ecut_eff,eigtmp,etotal_read,gsqcut_eff,maxfor_dum
 real(dp) :: maxocc,nelect,occlast,occtmp,rhotmp
 character(len=69) :: TypeCalcStrg
 character(len=500) :: message
 type(energies_type) :: energies_tmp
 type(wvl_internal_type) :: wvl
!arrays
 integer,allocatable :: nlmn(:)
 real(dp) :: cgtmp(2)
 real(dp),parameter :: qphon(3)=(/zero,zero,zero/)
 real(dp),allocatable :: favg_dum(:),fcart_dum(:,:),forold_dum(:,:),fred_tmp(:,:)
 real(dp),allocatable :: gresid_dum(:,:),grhf_dum(:,:),grxc_dum(:,:)
 real(dp),allocatable :: nhatgr(:,:,:),rhog_ep(:,:),scocc(:),str_tmp(:),synlgr_dum(:,:)
 type(cprj_type),allocatable :: cprj(:,:)  ! This is temporary waiting for cprj from SCF cycle...
 type(cprj_type),allocatable :: cprj_tmp(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij_tmp(:)

! *************************************************************************

!Compatibility tests
 if (dtset%positron==0) then
   write(message, '(4a)' ) ch10,&
&   ' setup_positron :  BUG -',ch10,&
&   '  Not valid for dtset%positron=0 !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if (istep>1.and.nfft/=electronpositron%nfft) then
   write(message, '(4a)' ) ch10,&
&   ' setup_positron :  BUG -',ch10,&
&   '  Invalid value for nfft !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if (dtset%positron==1) then
   do isppol=1,dtset%nsppol
     do ikpt=1,dtset%nkpt
       if (dtset%nband(ikpt+dtset%nkpt*(isppol-1))/=dtset%nband(1)) then
         write(message, '(4a)' ) ch10,&
&         ' setup_positron :  ERROR -',ch10,&
&         '  dtset%positron needs nband to be the same at each k-point !'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end do
 end if
 if (dtset%positron<=-10.and.dtset%mkmem==0) then
   write(message, '(4a)' ) ch10,&
&   ' setup_positron :  ERROR -',ch10,&
&   '  mkmem=0 not compatible with dtset%positron<=-10 !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!-----------------------------------------------------------------------
!Compute new value for calctype (type of electron-positron calculation)
!-----------------------------------------------------------------------
 icalctype0=electronpositron%calctype

 new_calctype=.false.
 if (dtset%positron==1.or.dtset%positron==2) then
   if (istep==1) new_calctype=.true.
 else if (dtset%positron<0) then
   if (ifirst_gs/=0.and.istep==1.and.(.not.always_restart)) new_calctype=.true.
   if (electronpositron%scf_converged) new_calctype=.true.
 end if

!Comment:
!history_level=-1:  not used
!history_level= 0:  rhor from scratch, rhor_ep from scratch or read
!history_level= 1:  rhor in memory, rhor_ep from scratch or read
!history_level= 2:  rhor_ep <-rhor, rhor from scratch
!history_level= 3:  rhor_ep <-> rhor
!history_level= 4:  rhor in memory, rhor_ep in memory
 history_level=-1
 if (dtset%positron==1.or.dtset%positron==2) then
   if (ifirst_gs==0.and.istep==1) history_level=0
   if (ifirst_gs/=0.and.istep==1) history_level=4
 else if (dtset%positron<0) then
   if (.not.electronpositron%scf_converged) then
     if (ifirst_gs/=0.and.istep==1.and.(.not.always_restart)) history_level=4
   else if (electronpositron%scf_converged) then
     if (icalctype0==0) history_level=2
     if (icalctype0> 0) history_level=3
   end if
 end if

 electronpositron%calctype=icalctype0
 if (dtset%positron==1.or.dtset%positron==2) then
   electronpositron%calctype=dtset%positron
 else if (dtset%positron<0) then
   if (electronpositron%scf_converged) then
     if (icalctype0==0) electronpositron%calctype=1
     if (icalctype0>0 ) electronpositron%calctype=3-electronpositron%calctype
   else if (ifirst_gs/=0.and.istep==1) then
     if (always_restart) then
       electronpositron%calctype=0
     else
       electronpositron%calctype=2
       if (electronpositron%particle==EP_POSITRON) electronpositron%calctype=1
     end if
   end if
 end if

 electronpositron%scf_converged=.false.
 if (new_calctype) electronpositron%istep=electronpositron%istep+1
 if (istep==1) electronpositron%istep=1
 ireadwf=dtfil%ireadwf;if (electronpositron%istep>1) ireadwf=0

!============================================
!The following lines occur only when the type
!of electron-positron calculation changes
!============================================
 if (new_calctype) then

!  Reset some indexes
   if (electronpositron%calctype==0) then
     electronpositron%particle=EP_NOTHING
   else if (electronpositron%calctype==1) then
     electronpositron%particle=EP_ELECTRON
   else if (electronpositron%calctype==2) then
     electronpositron%particle=EP_POSITRON
   end if
   electronpositron%has_pos_ham=mod(electronpositron%calctype,2)
   electronpositron%istep_scf=1
   istep_mix=1

!  -----------------------------------------------------------------------------------------
!  Update forces and stresses
!  If electronpositron%calctype==1: fred_ep/stress_ep are the electronic fred/stress
!  If electronpositron%calctype==2: fred_ep/stress_ep are the positronic fred/stress
!  -----------------------------------------------------------------------------------------
   if (history_level==2.or.history_level==3) then
     optstr=0;optfor=0
     if (associated(electronpositron%stress_ep)) optstr=stress_needed
     if (associated(electronpositron%fred_ep).and.forces_needed==2) optfor=1
     if (optfor>0.or.optstr>0) then
       ABI_ALLOCATE(favg_dum,(3))
       ABI_ALLOCATE(fcart_dum,(3,dtset%natom))
       ABI_ALLOCATE(forold_dum,(3,dtset%natom))
       ABI_ALLOCATE(gresid_dum,(3,dtset%natom))
       ABI_ALLOCATE(grhf_dum,(3,dtset%natom))
       ABI_ALLOCATE(grxc_dum,(3,dtset%natom))
       ABI_ALLOCATE(synlgr_dum,(3,dtset%natom))
       ABI_ALLOCATE(fred_tmp,(3,dtset%natom))
       ABI_ALLOCATE(str_tmp,(6))
       forold_dum=zero;n3xccc0=n3xccc
       icalctype=electronpositron%calctype;electronpositron%calctype=-icalctype0
       if (electronpositron%calctype==0) electronpositron%calctype=-100
       if (electronpositron%calctype==-1) n3xccc0=0  ! Note: if calctype=-1, previous calculation was positron
       call forstr(atindx1,cg,diffor_dum,dtefield,dtset,eigen,electronpositron,energies,favg_dum,fcart_dum,&
&       forold_dum,fred_tmp,gresid_dum,grewtn,grhf_dum,grxc_dum,gsqcut,indsym,&
&       kg,kxc,maxfor_dum,mcg,mgfft,mpi_enreg,n3xccc0,nattyp,nfft,ngfft,nhat,nkxc,npwarr,&
&       dtset%ntypat,nvresid,occ,optfor,optres,paw_ij,pawang,pawfgr,&
&       pawfgrtab,pawrhoij,pawtab,pel,ph1dc,ph1d,pion,psps,rhog,rhor,rprimd,optstr,&
&       strsxc,str_tmp,symrec,synlgr_dum,ucvol,dtfil%unkg,dtfil%unylm,&
&       usexcnhat,vhartr,vpsp,vxc,wffnow,wvl,xccc3d,xred,ylm,ylmgr)
       electronpositron%calctype=icalctype
       if (optfor>0) electronpositron%fred_ep(:,:)=fred_tmp(:,:)
       if (optstr>0) electronpositron%stress_ep(:)=str_tmp(:)
       ABI_DEALLOCATE(favg_dum)
       ABI_DEALLOCATE(fcart_dum)
       ABI_DEALLOCATE(forold_dum)
       ABI_DEALLOCATE(gresid_dum)
       ABI_DEALLOCATE(grhf_dum)
       ABI_DEALLOCATE(grxc_dum)
       ABI_DEALLOCATE(synlgr_dum)
       ABI_DEALLOCATE(fred_tmp)
       ABI_DEALLOCATE(str_tmp)
     end if
     if (optfor==0.and.forces_needed>0.and.associated(electronpositron%fred_ep)) then
       electronpositron%fred_ep(:,:)=fred(:,:)-electronpositron%fred_ep(:,:)
     end if
   end if

!  ----------------------------------------------------------------------------------------------------
!  Initialize/Update densities
!  If electronpositron%calctype==1: rhor is the positronic density, rhor_ep is the electronic density
!  If electronpositron%calctype==2: rhor is the electronic density, rhor_ep is the positronic density
!  ---------------------------------------------------------------------------------------------------
   ABI_ALLOCATE(rhog_ep,(2,nfft))

!  ===== PREVIOUS DENSITY RHOR_EP:
   if (history_level==0.or.history_level==1) then
!    ----- Read from disk
     if (dtset%positron>0) then
       rdwr=1;rdwrpaw=dtset%usepaw;accessfil=dtset%accesswff
       if (dtset%accesswff==IO_MODE_MPI) accessfil=4
       if (dtset%accesswff==IO_MODE_NETCDF) accessfil=1
       call ioarr(accessfil,electronpositron%rhor_ep,dtset,etotal_read,fformr,dtfil%fildensin,&
&       hdr,mpi_enreg,nfft,electronpositron%pawrhoij_ep,rdwr,rdwrpaw,wvl)
       call fourdp(1,rhog_ep,electronpositron%rhor_ep,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
       if (dtset%usepaw==1.and.associated(electronpositron%nhat_ep)) then
         call pawmknhat(occtmp,1,0,0,0,0,gprimd,mpi_enreg,dtset%natom,dtset%natom,nfft,ngfft,0,&
&         dtset%nspden,dtset%ntypat,dtset%paral_kgb,pawang,pawfgrtab,nhatgr,electronpositron%nhat_ep,&
&         electronpositron%pawrhoij_ep,electronpositron%pawrhoij_ep,pawtab,qphon,rprimd,ucvol,xred)
       end if
     end if
!    ----- Electronic from scratch
     if (dtset%positron<0.and.electronpositron%calctype==1) then
       ecut_eff=dtset%pawecutdg*(dtset%dilatmx)**2
       call getcut(boxcut_dum,ecut_eff,gmet,gsqcut_eff,dtset%iboxcut,std_out,qphon,ngfft)
       call initro(atindx,dtset%densty,gmet,gsqcut_eff,dtset%usepaw,mgfft,mpi_enreg,&
&       psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,dtset%nspden,dtset%ntypat,&
&       dtset%paral_kgb,pawtab,ph1d,psps%qgrid_vl,rhog_ep,electronpositron%rhor_ep,&
&       dtset%spinat,ucvol,dtset%usepaw,dtset%ziontypat,dtset%znucl)
       if (dtset%usepaw==1) then
         ABI_ALLOCATE(pawrhoij_tmp,(dtset%natom))
         call initrhoij(electronpositron%pawrhoij_ep(1)%cplex,psps%indlmn,dtset%lexexch,psps%lmnmax,&
&         dtset%lpawu,mpi_enreg,dtset%natom,mpi_enreg%natom,electronpositron%pawrhoij_ep(1)%nspden,&
&         electronpositron%pawrhoij_ep(1)%nspinor,dtset%nsppol,&
&         dtset%ntypat,pawrhoij_tmp,pawtab,dtset%spinat,dtset%typat,&
&         ngrhoij=electronpositron%pawrhoij_ep(1)%ngrhoij,&
&         nlmnmix=electronpositron%pawrhoij_ep(1)%lmnmix_sz,&
&         use_rhoij_=electronpositron%pawrhoij_ep(1)%use_rhoij_,&
&         use_rhoijres=electronpositron%pawrhoij_ep(1)%use_rhoijres)
         if (electronpositron%pawrhoij_ep(1)%lmnmix_sz>0) then
           do iatom=1,dtset%natom
             pawrhoij_tmp(iatom)%kpawmix(:)=electronpositron%pawrhoij_ep(iatom)%kpawmix(:)
           end do
         end if
         call rhoij_copy(pawrhoij_tmp,electronpositron%pawrhoij_ep)
         call rhoij_free(pawrhoij_tmp)
         ABI_DEALLOCATE(pawrhoij_tmp)
         if (associated(electronpositron%nhat_ep)) then
           call pawmknhat(occtmp,1,0,0,0,0,gprimd,mpi_enreg,dtset%natom,dtset%natom,nfft,ngfft,0,&
&           dtset%nspden,dtset%ntypat,dtset%paral_kgb,pawang,pawfgrtab,nhatgr,electronpositron%nhat_ep,&
&           electronpositron%pawrhoij_ep,electronpositron%pawrhoij_ep,pawtab,qphon,rprimd,ucvol,xred)
         end if
       end if
     end if
!    ----- Positronic from scratch
     if (dtset%positron<0.and.electronpositron%calctype==2) then
       electronpositron%rhor_ep(:,1)=one/ucvol
       if (dtset%nspden>=2) electronpositron%rhor_ep(:,2)=half/ucvol
       if (dtset%nspden==4) electronpositron%rhor_ep(:,3:4)=zero
       rhog_ep=zero;rhog_ep(1,1)=one/ucvol
       if (dtset%usepaw==1) then
         do iatom=1,dtset%natom
           electronpositron%pawrhoij_ep(iatom)%rhoijp(:,:)=zero
           electronpositron%pawrhoij_ep(iatom)%nrhoijsel=0
         end do
         if (associated(electronpositron%nhat_ep)) electronpositron%nhat_ep=zero
       end if
     end if
   end if
!  ----- Deduced from rhor in memory
   if (history_level==2) then
     electronpositron%rhor_ep(:,:)=rhor(:,:)
     rhog_ep(:,:)=rhog(:,:)
     if (dtset%usepaw==1) then
       call rhoij_copy(pawrhoij,electronpositron%pawrhoij_ep)
       if (associated(electronpositron%nhat_ep)) electronpositron%nhat_ep(:,:)=nhat(:,:)
     end if
   end if

!  ===== CURRENT DENSITY RHOR:
   if (history_level==0.or.history_level==2) then
     if (ireadwf==0) then
!      ----- Positronic from scratch
       if (electronpositron%calctype==1) then
         rhor(:,1)=one/ucvol
         if (dtset%nspden>=2) rhor(:,2)=half/ucvol
         if (dtset%nspden==4) rhor(:,3:4)=zero
         rhog=zero;rhog(1,1)=one/ucvol
         if (dtset%usepaw==1) then
           do iatom=1,dtset%natom
             pawrhoij(iatom)%rhoijp(:,:)=zero
             pawrhoij(iatom)%nrhoijsel=0
           end do
           nhat(:,:)=zero
         end if
       end if
!      ----- Electronic from scratch
       if (electronpositron%calctype==2) then
         ecut_eff=dtset%pawecutdg*(dtset%dilatmx)**2
         call getcut(boxcut_dum,ecut_eff,gmet,gsqcut_eff,dtset%iboxcut,std_out,qphon,ngfft)
         call initro(atindx,dtset%densty,gmet,gsqcut_eff,dtset%usepaw,mgfft,mpi_enreg,&
&         psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,dtset%nspden,dtset%ntypat,&
&         dtset%paral_kgb,pawtab,ph1d,psps%qgrid_vl,rhog,rhor,dtset%spinat,ucvol,&
&         dtset%usepaw,dtset%ziontypat,dtset%znucl)
         if (dtset%usepaw==1) then
           ABI_ALLOCATE(pawrhoij_tmp,(dtset%natom))
           call initrhoij(pawrhoij(1)%cplex,psps%indlmn,dtset%lexexch,psps%lmnmax,dtset%lpawu,mpi_enreg,dtset%natom,&
&           mpi_enreg%natom,pawrhoij(1)%nspden,pawrhoij(1)%nspinor,dtset%nsppol,&
&           dtset%ntypat,pawrhoij_tmp,pawtab,dtset%spinat,&
&           dtset%typat,ngrhoij=pawrhoij(1)%ngrhoij,nlmnmix=pawrhoij(1)%lmnmix_sz,&
&           use_rhoij_=pawrhoij(1)%use_rhoij_,use_rhoijres=pawrhoij(1)%use_rhoijres)
           do iatom=1,dtset%natom
             pawrhoij_tmp(iatom)%kpawmix(:)=pawrhoij(iatom)%kpawmix(:)
           end do
           call rhoij_copy(pawrhoij_tmp,pawrhoij)
           call rhoij_free(pawrhoij_tmp)
           ABI_DEALLOCATE(pawrhoij_tmp)
           call pawmknhat(occtmp,1,0,0,0,0,gprimd,mpi_enreg,dtset%natom,dtset%natom,nfft,ngfft,0,&
&           dtset%nspden,dtset%ntypat,dtset%paral_kgb,pawang,pawfgrtab,nhatgr,nhat,&
&           pawrhoij,pawrhoij,pawtab,qphon,rprimd,ucvol,xred)
         end if
       end if
     end if
   end if

!  ===== EXCHANGE POSITRONIC AND ELECTRONIC DENSITY (CURRENT AND PREVIOUS)
   if (history_level==3) then
     do ispden=1,dtset%nspden
       do ifft=1,nfft
         rhotmp=rhor(ifft,ispden)
         rhor(ifft,ispden)=electronpositron%rhor_ep(ifft,ispden)
         electronpositron%rhor_ep(ifft,ispden)=rhotmp
       end do
     end do
     rhog_ep(:,:)=rhog
     call fourdp(1,rhog,rhor,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
!    If PAW, exchange "positronic" and "electronic" rhoij
     if (dtset%usepaw==1) then
       ABI_ALLOCATE(pawrhoij_tmp,(dtset%natom))
       call rhoij_alloc(pawrhoij(1)%cplex,hdr%lmn_size(1:dtset%ntypat),pawrhoij(1)%nspden,&
&       pawrhoij(1)%nspinor,pawrhoij(1)%nsppol,pawrhoij_tmp,dtset%typat,mpi_enreg=mpi_enreg,&
&       ngrhoij=pawrhoij(1)%ngrhoij,nlmnmix=pawrhoij(1)%lmnmix_sz,&
&       use_rhoij_=pawrhoij(1)%use_rhoij_,use_rhoijres=pawrhoij(1)%use_rhoijres)
       call rhoij_copy(pawrhoij,pawrhoij_tmp)
       call rhoij_copy(electronpositron%pawrhoij_ep,pawrhoij)
       call rhoij_copy(pawrhoij_tmp,electronpositron%pawrhoij_ep)
       call rhoij_free(pawrhoij_tmp)
       ABI_DEALLOCATE(pawrhoij_tmp)
       if (associated(electronpositron%nhat_ep)) then
         do ispden=1,dtset%nspden
           do ifft=1,nfft
             rhotmp=nhat(ifft,ispden)
             nhat(ifft,ispden)=electronpositron%nhat_ep(ifft,ispden)
             electronpositron%nhat_ep(ifft,ispden)=rhotmp
           end do
         end do
       end if
     end if
   end if

!  ===== COMPUTE HARTREE POTENTIAL ASSOCIATED TO RHOR_EP
   if (history_level==4) then
     call fourdp(1,rhog_ep,electronpositron%rhor_ep,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   end if
   if (history_level/=-1) then
     call hartre(1,gmet,gsqcut,dtset%usepaw,mpi_enreg,nfft,ngfft,dtset%paral_kgb,qphon,rhog_ep,&
&     electronpositron%vha_ep)
     electronpositron%vha_ep=-electronpositron%vha_ep
   else
     electronpositron%vha_ep=zero
   end if
   ABI_DEALLOCATE(rhog_ep)

!  ----------------------------------------------------------------------
!  Initialize/Update energies
!  ----------------------------------------------------------------------
   electronpositron%etotal_prev=etotal
   electronpositron%maxfor_prev=maxfor

!  Inits/exchange news energies
!  Retrieve energy of non-evolving particle(s)
   if (history_level== 0) then
     call energies_init(energies)
     call energies_init(electronpositron%energies_ep)
     if (dtset%positron>0) energies%e0_electronpositron=etotal_read
     if (dtset%positron<0) energies%e0_electronpositron=zero
   else if (history_level== 1) then
     call energies_init(electronpositron%energies_ep)
     if (dtset%positron>0) energies%e0_electronpositron=etotal_read
   else if (history_level== 2) then
     call energies_copy(energies,electronpositron%energies_ep)
     call energies_init(energies)
     energies%e0_electronpositron=electronpositron%e0
   else if (history_level== 3) then
     call energies_copy(electronpositron%energies_ep,energies_tmp)
     call energies_copy(energies,electronpositron%energies_ep)
     call energies_copy(energies_tmp,energies)
     energies%e0_electronpositron=electronpositron%e0
!    else if (history_level== 4) then
   end if

!  Adjust core psps energy
   if (electronpositron%calctype/=1) energies%e_corepsp=ecore/ucvol

!  -----------------------------------------------------------------------------------------
!  Update wavefunctions
!  If electronpositron%calctype==1: cg are the positronic WFs, cg_ep are the electronic WFs
!  If electronpositron%calctype==2: cg are the electronic WFs, cg_ep are the positronic WFs
!  -----------------------------------------------------------------------------------------
   if (electronpositron%dimcg>0.or.electronpositron%dimcprj>0) then

     if (history_level==0.or.history_level==1) then
       electronpositron%cg_ep=zero
     end if

     if (history_level==2) then
       if (electronpositron%dimcg>0) then
         do icg=1,electronpositron%dimcg
           electronpositron%cg_ep(1:2,icg)=cg(1:2,icg)
         end do
       end if
       if (dtset%usepaw==1.and.electronpositron%dimcprj>0) then
         call cprj_copy(cprj,electronpositron%cprj_ep)
       end if
     end if

     if (history_level==3) then
       if (electronpositron%dimcg>0) then
         do icg=1,electronpositron%dimcg
           cgtmp(1:2)=electronpositron%cg_ep(1:2,icg)
           electronpositron%cg_ep(1:2,icg)=cg(1:2,icg)
           cg(1:2,icg)=cgtmp(1:2)
         end do
       end if
       if (dtset%usepaw==1.and.electronpositron%dimcprj>0) then
         ABI_ALLOCATE(nlmn,(dtset%natom))
         ABI_ALLOCATE(cprj_tmp,(dtset%natom,electronpositron%dimcprj))
         do iatom=1,dtset%natom
           nlmn(iatom)=cprj(iatom,1)%nlmn
         end do
         call cprj_alloc(cprj_tmp,cprj(1,1)%ncpgr,nlmn)
         ABI_DEALLOCATE(nlmn)
         call cprj_copy(electronpositron%cprj_ep,cprj_tmp)
         call cprj_copy(cprj,electronpositron%cprj_ep)
         call cprj_copy(cprj_tmp,cprj)
         call cprj_free(cprj_tmp)
         ABI_DEALLOCATE(cprj_tmp)
       end if
     end if

   end if ! dimcg>0 or dimcprj>0

!  -----------------------------------------------------------------------------------------------------------
!  Initialize/Update occupations
!  If electronpositron%calctype==1: occ are the positronic occupations, occ_ep are the electronic occupations
!  If electronpositron%calctype==2: occ are the electronic occupations, occ_ep are the positronic occupations
!  -----------------------------------------------------------------------------------------------------------
!  When needed, precompute electronic occupations with semiconductor occupancies
   need_scocc=.false.
   if (electronpositron%dimocc>0.and.electronpositron%calctype==1.and. &
&   (history_level==0.or.history_level==1)) need_scocc=.true.
   if (electronpositron%calctype==2.and.ireadwf==0.and. &
&   (history_level==0.or.history_level==2.or. &
&   (history_level==3.and.electronpositron%dimocc==0))) need_scocc=.true.
   if (need_scocc) then
     nelect=-dtset%charge
     do iatom=1,dtset%natom
       nelect=nelect+dtset%ziontypat(dtset%typat(iatom))
     end do
     maxocc=two/real(dtset%nsppol*dtset%nspinor,dp)
     nocc=(nelect-tol8)/maxocc + 1
     nocc=min(nocc,dtset%nband(1)*dtset%nsppol)
     occlast=nelect-maxocc*(nocc-1)
     ABI_ALLOCATE(scocc,(dtset%nband(1)*dtset%nsppol))
     scocc=zero
     if (1<nocc)  scocc(1:nocc-1)=maxocc
     if (1<=nocc) scocc(nocc)=occlast
   end if

!  ===== PREVIOUS OCCUPATIONS OCC_EP:
   if (electronpositron%dimocc>0) then
     if (history_level==0.or.history_level==1) then
!      ----- Electronic from scratch
       if (electronpositron%calctype==1) then
!        Initialize electronic occupations with semiconductor occupancies
         do ikpt=1,dtset%nkpt
           do iband=1,dtset%nband(1)
             do isppol=1,dtset%nsppol
               electronpositron%occ_ep(iband+dtset%nband(1)*(ikpt-1+dtset%nkpt*(isppol-1)))=&
&               scocc(isppol+dtset%nsppol*(iband-1))
             end do
           end do
         end do
       end if
!      ----- Positronic from scratch
       if (electronpositron%calctype==1) then
!        Initialize positronic occupations with only one positron (or less)
         electronpositron%occ_ep(:)=zero
         isppol=1;iocc=1
         do ikpt=1,dtset%nkpt
           electronpositron%occ_ep(iocc)=electronpositron%posocc
           iocc=iocc+dtset%nband(ikpt+dtset%nkpt*(isppol-1))
         end do
       end if
     end if
!    ----- Deduced from occ in memory
     if (history_level==2) then
       electronpositron%occ_ep(:)=occ(:)
     end if
   end if ! dimocc>0

!  ===== CURRENT OCCUPATIONS OCC:
   if (history_level==0.or.history_level==2.or.(history_level==3.and.electronpositron%dimocc==0)) then
     if (ireadwf==0) then
!      ----- Positronic from scratch
       if (electronpositron%calctype==1) then
!        Initialize positronic occupations with only one positron (or less)
         occ(:)=zero
         isppol=1;iocc=1
         do ikpt=1,dtset%nkpt
           occ(iocc)=electronpositron%posocc
           iocc=iocc+dtset%nband(ikpt+dtset%nkpt*(isppol-1))
         end do
       end if
!      ----- Electronic from scratch
       if (electronpositron%calctype==2) then
!        Initialize electronic occupations with semiconductor occupancies
         do ikpt=1,dtset%nkpt
           do iband=1,dtset%nband(1)
             do isppol=1,dtset%nsppol
               occ(iband+dtset%nband(1)*(ikpt-1+dtset%nkpt*(isppol-1)))=&
&               scocc(isppol+dtset%nsppol*(iband-1))
             end do
           end do
         end do
       end if
     end if
   end if

!  ===== EXCHANGE POSITRONIC AND ELECTRONIC OCCUPATIONS (CURRENT AND PREVIOUS)
   if (history_level==3.and.electronpositron%dimocc>0) then
     do iocc=1,electronpositron%dimocc
       occtmp=occ(iocc)
       occ(iocc)=electronpositron%occ_ep(iocc)
       electronpositron%occ_ep(iocc)=occtmp
     end do
   end if

   if (need_scocc)  then
     ABI_DEALLOCATE(scocc)
   end if

!  -----------------------------------------------------------------------------------------------------------
!  Initialize/Update eigen energies
!  If electronpositron%calctype==1: eigen are the positronic eigen E, eigen_ep are the electronic eigen E
!  If electronpositron%calctype==2: eigen are the electronic eigen E, eigen_ep are the positronic eigen E
!  -----------------------------------------------------------------------------------------------------------

!  ===== PREVIOUS EIGEN ENERGIES EIGEN_EP:
   if (electronpositron%dimeigen>0) then
     if (history_level==0.or.history_level==1) then
!      ----- Electronic or positronic from scratch
       electronpositron%eigen_ep(:)=zero
     end if
!    ----- Deduced from eigen in memory
     if (history_level==2) then
       electronpositron%eigen_ep(:)=eigen(:)
     end if
   end if ! dimeigen>0

!  ===== CURRENT EIGEN ENERGIES EIGEN:
   if (history_level==0.or.history_level==2.or.(history_level==3.and.electronpositron%dimeigen==0)) then
     if (ireadwf==0) then
!      ----- Electronic or positronic from scratch
       eigen(:)=zero
     end if
   end if

!  ===== EXCHANGE POSITRONIC AND ELECTRONIC EIGEN ENERGIES (CURRENT AND PREVIOUS)
   if (history_level==3.and.electronpositron%dimeigen>0) then
     do iocc=1,electronpositron%dimeigen
       eigtmp=eigen(iocc)
       eigen(iocc)=electronpositron%eigen_ep(iocc)
       electronpositron%eigen_ep(iocc)=eigtmp
     end do
   end if

!  =============================================
 end if  ! the type of e-p calculation changes
!=============================================

!------------------------------------------------------------------
!Write messages
!------------------------------------------------------------------
 if (istep_mix==1.and.dtset%positron/=0) then
!  Log message
   if (electronpositron%calctype==0) then
     write(message, '(4a)' ) ch10,&
&     ' setup_positron :  COMMENT -',ch10,&
&     '  Were are now performing an electronic ground-state calculation...'
   else if (electronpositron%calctype==1) then
     write(message, '(4a)' ) ch10,&
&     ' setup_positron :  COMMENT -',ch10,&
&     '  Were are now performing a positronic ground-state calculation...'
   else if (electronpositron%calctype==2) then
     write(message, '(4a)' ) ch10,&
&     ' setup_positron :  COMMENT -',ch10,&
&     '  Were are now performing an electronic ground-state calculation in presence of a positron...'
   end if
   call wrtout(std_out,message,'COLL')
!  Output message
   if (dtset%positron<0) then
     if (electronpositron%calctype==0) then
       TypeCalcStrg='ELECTRONIC GROUND-STATE CALCULATION'
     else if (electronpositron%calctype==1) then
       TypeCalcStrg='POSITRONIC GROUND-STATE CALCULATION IN PRESENCE OF ELECTRONS AND IONS'
     else if (electronpositron%calctype==2) then
       TypeCalcStrg='ELECTRONIC GROUND-STATE CALCULATION IN PRESENCE OF A POSITRON'
     end if
     if (istep>1) then
       write(message,'(2a,i3,2a)') ch10,'TC-DFT STEP ',electronpositron%istep,' - ',trim(TypeCalcStrg)
     else
       write(message,'(a,i3,2a)') 'TC-DFT STEP ',electronpositron%istep,' - ',trim(TypeCalcStrg)
     end if
     call wrtout(ab_out,message,'COLL')
   end if
 end if

end subroutine setup_positron
!!***
