!{\src2tex{textfont=tt}}
!!****f* ABINIT/newvtr
!! NAME
!! newvtr
!!
!! FUNCTION
!! Compute new trial potential by mixing new and old values.
!! Call prcref to compute preconditioned residual potential and forces,
!! Then, call one of the self-consistency drivers,
!! then update vtrial.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  cg(2,mcg)=updated wavefunctions; if mkmem>=nkpt, these are kept in a disk file.
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  dielinv(2,npwdiel,nspden,npwdiel,nspden)=
!!                              inverse of the dielectric matrix in rec. space
!!  dielstrt=number of the step at which the dielectric preconditioning begins.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | fixmom=input variable that governs fixed moment calculation
!!   | intxc=control xc quadrature
!!   | iprcch= governs the preconditioning of the atomic charges
!!   | iprcel= governs the preconditioning of the potential residual
!!   | iprcfc=governs the preconditioning of the forces
!!   | iscf=( <= 0 =>non-SCF), >0 => SCF)
!!   |  iscf =1 => determination of the largest eigenvalue of the SCF cycle
!!   |  iscf =2 => SCF cycle, simple mixing
!!   |  iscf =3 => SCF cycle, Anderson mixing
!!   |  iscf =4 => SCF cycle, Anderson mixing (order 2)
!!   |  iscf =5 => SCF cycle, CG based on the minimization of the energy
!!   |  iscf =7 => SCF cycle, Pulay mixing
!!   | isecur=level of security of the computation
!!   | ixc=exchange-correlation choice parameter.
!!   | mffmem=governs the number of FFT arrays which are fit in core memory
!!   |          it is either 1, in which case the array f_fftgr is used,
!!   |          or 0, in which case the array f_fftgr_disk is used
!!   | natom=number of atoms
!!   | nspden=number of spin-density components
!!   | occopt=option for occupancies
!!   | paral_kgb=option for (kpt,g vectors,bands) parallelism
!!   | pawoptmix= - PAW only - 1 if the computed residuals include the PAW (rhoij) part
!!   | prtvol=control print volume and debugging
!!   | typat(natom)=integer type for each atom in cell
!!  etotal=the total energy obtained from the input vtrial
!!  fcart(3,natom)=cartesian forces (hartree/bohr)
!!  ffttomix(nfft*(1-nfftmix/nfft))=Index of the points of the FFT (fine) grid on the grid used for mixing (coarse)
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  initialized= if 0, the initialization of the gstate run is not yet finished
!!  ispmix=1 if mixing is done in real space, 2 if mixing is done in reciprocal space
!!  istep= number of the step in the SCF cycle
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only for electronic!
!     dielectric matrix
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mixtofft(nfftmix*(1-nfftmix/nfft))=Index of the points of the FFT grid used for mixing (coarse) on the FFT (fine) grid
!!  moved_atm_inside= if 1, then the preconditioned forces
!!    as well as the preconditioned potential residual must be computed;
!!    otherwise, compute only the preconditioned potential residual.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftmix=dimension of FFT grid used to mix the densities (used in PAW only)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngfftmix(18)=contain all needed information about 3D FFT, for the grid corresponding to nfftmix
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  npwdiel=number of planewaves for dielectric matrix
!!  nstep=number of steps expected in iterations.
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                         Use here rhoij residuals (and gradients)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vhartr(nfft)=array for holding Hartree potential
!!  vnew_mean(nspden)=constrained mean value of the future trial potential (might be
!!    spin-polarized
!!  vpsp(nfft)=array for holding local psp
!!  vresid(nfft,nspden)=array for the residual of the potential
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  dbl_nnsclo=1 if nnsclo has to be doubled to secure the convergence.
!!
!! SIDE EFFECTS
!!  dtn_pc(3,natom)=preconditioned change of atomic position,
!!                                          in reduced coordinates
!!  vtrial(nfft,nspden)= at input, it is the "in" trial potential that gave vresid=(v_out-v_in)
!!       at output, it is an updated "mixed" trial potential
!!  ===== if iprcch==3 .and. moved_atm_inside==1 =====
!!    ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases
!!  ==== if usepaw==1
!!    pawrhoij(natom)%nrhoijsel,rhoijselect,rhoijp= several arrays
!!                containing new values of rhoij (augmentation occupancies)
!!
!! WARNINGS
!! depending on the value of iprcch and moved_atm_inside,
!! the xc potential or the Hxc potential may have been subtracted from vtrial !
!!
!! NOTES
!!  In case of PAW calculations:
!!    Computations are done either on the fine FFT grid or the coarse grid (depending on dtset%pawmixdg)
!!    All variables (nfft,ngfft,mgfft) refer to the fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  ! Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!!  Subtility in PAW and non-collinear magnetism:
!!    Potentials are stored in (up-up,dn-dn,Re[up-dn],Im[up-dn]) format
!!    On-site occupancies (rhoij) are stored in (n,mx,my,mz)
!!    This is compatible provided that the mixing factors for n and m are identical
!!    and that the residual is not a combination of V_res and rhoij_res (pawoptmix=0).
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      ab6_mixing_copy_current_step,ab6_mixing_eval,ab6_mixing_eval_allocate
!!      ab6_mixing_eval_deallocate,ab6_mixing_use_moving_atoms,fourdp,lavnl
!!      leave_new,mean_fftr,metric,prcref_pma,prctfvw1,prctfvw2,timab,wrtout
!!      xcomm_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine newvtr(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,&
     &  dtn_pc,dtset,efermi,etotal,fcart,ffttomix,&
     &  gmet,grhf,gsqcut,&
     &  initialized,ispmix,&
     &  istep,&
     &  kg_diel,kxc,mgfft,mix,mixtofft,&
     &  moved_atm_inside,mpi_enreg,nattyp,nfft,nfftmix,&
     &  nhat,nhatgr,nhatgrdim,&
     &  ngfft,ngfftmix,nkxc,npawmix,npwdiel,&
     &  nstep,ntypat,n1xccc,optres,optxc,&
     &  pawrhoij,&
     &  ph1d,&
     &  psps,rhor,rprimd,susmat,usepaw,&
     &  vhartr,vnew_mean,vpsp,vresid,&
     &  vtrial,vxc,xred,&
     &  atindx1,cg,deltae,&
     &  dtfil,eeig,eigen,ek,enl,kg,&
     &  mcg,nfftf,&
     &  ngfftf,npwarr,n3xccc,occ,optene,&
     &  pawfgr,pawtab,&
     &  resid,rhog,&
     &  usexcnhat,&
     &  wffnow,wvl,&
     &  ylm,xccc3d )

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_wffile
 use m_ab6_mixing
 use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'newvtr'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_68_rsprc, except_this_one => newvtr
!End of the abilint section

 implicit none

!Arguments-------------------------------
  ! WARNING
  ! BEWARE THERE IS TWO DIFFERENT SIZE DECLARED FOR ARRAY NHAT IN RHOTOV AND RHOHXC
  ! THIS MIGHT RESULT IN A BUG
!scalars
 integer,intent(in) :: dielstrt,initialized,ispmix,istep,mcg,mgfft
 integer,intent(in) :: moved_atm_inside,n1xccc,n3xccc,nfft
 integer,intent(in) :: nfftf,nfftmix,nhatgrdim,nkxc,npawmix,npwdiel,nstep
 integer,intent(in) :: ntypat,optene,optres,optxc,usepaw,usexcnhat
 integer,intent(inout) :: dbl_nnsclo
 real(dp),intent(in) :: deltae,efermi,etotal,gsqcut
 real(dp),intent(out) :: eeig,ek,enl
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(ab6_mixing_object),intent(inout) :: mix
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),kg_diel(3,npwdiel)
 integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft)),nattyp(ntypat)
 integer,intent(in) :: ngfft(18),ngfftf(18),ngfftmix(18),npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: dielar(7),eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: fcart(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(in) :: nhat(nfftf,dtset%nspden*usepaw)
 real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(in) :: vhartr(nfft),vnew_mean(dtset%nspden)
 real(dp),intent(in) :: vxc(nfft,dtset%nspden)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(inout), target :: dtn_pc(3,dtset%natom)
 real(dp),intent(inout) :: gmet(3,3)
 real(dp),intent(inout) :: kxc(nfft,nkxc),ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(inout) :: rhog(2,nfftf),vpsp(nfft)
 real(dp),intent(inout), target :: rhor(nfft,dtset%nspden)
 real(dp),intent(inout) :: vresid(nfft,dtset%nspden),vtrial(nfft,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc)
 real(dp),intent(inout), target :: xred(3,dtset%natom)
 real(dp),intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: compute_lavnlr,cplex,dplex,i_vresid1,i_vrespc1
! integer :: i1,i2,i3,ifft2,ifft3,ifft4,ifft5,ii1,ii2,ii3,ii4,ii5
 integer :: iatom,ifft,indx,irhoij,mpi_comm,old_paral_level
 integer :: ispden,jfft,jrhoij,klmn,kmix,n1,n2,n3,nfftot,nselect
 integer :: use_lavnlr,errid,tim_fourdp
 logical :: mpi_summarize, reset
 real(dp) :: dielng,diemix,fact,ucvol,vme,vme_inold,vme_um
! real(dp) :: sig,x,xm,y,ym,z,zm
 character(len=500) :: message
! character(len=fnlen) :: filapp
!arrays
 real(dp),parameter :: identity(4)=(/1.0_dp,1.0_dp,0.0_dp,0.0_dp/)
 real(dp) :: gprimd(3,3),rmet(3,3),tsec(2),vmean(dtset%nspden)
 real(dp) :: vmean_inold(dtset%nspden),vmean_um(dtset%nspden)
!DEBUG
! real(dp),allocatable :: buffer1(:,:,:),buffer2(:,:,:),dvstar(:,:),g2cart(:)
! real(dp),allocatable :: irdiemac(:),irdiemacf1(:),irdiemacf2(:),irdiemacg(:,:)
! real(dp),allocatable :: ldvstar(:,:),lvres(:,:),rdiemac(:),vres(:,:)
! real(dp),allocatable :: rdiemacf1(:),rdiemacf2(:),rdiemacg(:,:)
!ENDDEBUG
 real(dp),allocatable :: lavnlr(:,:),rhoijrespc(:)
 real(dp),allocatable :: rhoijtmp(:,:),vin_old(:,:),vout_unmixed(:,:)
 real(dp),allocatable :: vresid0(:,:),vrespc(:,:),vreswk(:,:),vtrialg(:,:,:)
 real(dp),pointer :: vtrial0(:,:),vpaw(:)

! *************************************************************************

!DEBUG
 write(std_out,*)' newvtr : enter '
 write(std_out,*)' newvtr : ispmix,nfft,nfftmix=',ispmix,nfft,nfftmix
!ENDDEBUG


 call timab(93,1,tsec)
 call timab(901,1,tsec)
 tim_fourdp=8

!Compatibility tests
 if(nfftmix>nfft) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' newvtr : BUG -',ch10,&
&   '  nfftmix>nfft not allowed !'
   call wrtout(std_out,message,'COLL')
   call leave_new('PERS')
 end if
 if(ispmix/=2.and.nfftmix/=nfft) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' newvtr : BUG -',ch10,&
&   '  nfftmix/=nfft allowed only when ispmix=2 !'
   call wrtout(std_out,message,'COLL')
   call leave_new('PERS')
 end if
 if(dtset%usewvl == 1) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' newvtr : BUG -',ch10,&
&   '  dtset%usewvl == 1 not allowed (use wvl_newtr() instead)!'
   call wrtout(std_out,message,'COLL')
   call leave_new('PERS')
 end if
 if(usepaw==1.and.dtset%nspden==4.and.dtset%pawoptmix==1) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' newvtr : ERROR -',ch10,&
&   '  pawoptmix=1 is not compatible with nspden=4 !'
   call wrtout(std_out,message,'COLL')
   call leave_new('PERS')
 end if

 dielng=dielar(2)
 diemix=dielar(4)
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 if (usepaw==1) then
   cplex=pawrhoij(1)%cplex;dplex=cplex-1
 else
   cplex=0;dplex=0
 end if

!DEBUG
!write(std_out,*)' newvtr : enter '
!write(std_out,*)' newvtr : vnew_mean(:)=',vnew_mean(:)
!write(std_out,*)' newvtr : vtrial(1,:)=',vtrial(1,:)
!stop
!ENDDEBUG

!Get size of FFT grid
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!------Treat the mean of potentiel residual

!Special care must be taken with components of the
!potential that are associated with NO density change.
!In general, only the global mean of the potential has
!such an anomalous feature. However, in the spin
!polarized cas with fixed occupancies, also the
!mean of each spin-potential (independently of the other)
!has such a behaviour. The trick is to remove these
!variables before going in the predictive routines,
!then to put them back

!Compute the mean of the old vtrial
 call mean_fftr(vtrial,vmean,mpi_enreg,nfft,nfftot,dtset%nspden)

!When (collinear) spin-polarized and fixed occupation numbers,
!treat separately spin up and spin down.
!Otherwise, use only global mean
 do ispden=1,dtset%nspden
   if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&   abs(dtset%fixmom+99.99_dp)<1.0d-10)then
     vme=(vmean(1)+vmean(2))*half
   else
     vme=vmean(ispden)
   end if
   vtrial(:,ispden)=vtrial(:,ispden)-vme
 end do

 call timab(901,2,tsec)

!Select components of potential to be mixed
 ABI_ALLOCATE(vtrial0,(ispmix*nfftmix,dtset%nspden))
 ABI_ALLOCATE(vresid0,(ispmix*nfftmix,dtset%nspden))
 if (ispmix==1.and.nfft==nfftmix) then
   vtrial0=vtrial;vresid0=vresid
 else if (nfft==nfftmix) then
   do ispden=1,dtset%nspden
     call fourdp(1,vtrial0(:,ispden),vtrial(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,tim_fourdp)
     call fourdp(1,vresid0(:,ispden),vresid(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,tim_fourdp)
   end do
 else
   ABI_ALLOCATE(vtrialg,(2,nfft,dtset%nspden))
   ABI_ALLOCATE(vreswk,(2,nfft))
   do ispden=1,dtset%nspden
     fact=dielar(4);if (ispden>1) fact=dielar(7)
     call fourdp(1,vtrialg(:,:,ispden),vtrial(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,tim_fourdp)
     call fourdp(1,vreswk,vresid(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,tim_fourdp)
     do ifft=1,nfft
       if (ffttomix(ifft)>0) then
         jfft=2*ffttomix(ifft)
         vtrial0(jfft-1,ispden)=vtrialg(1,ifft,ispden)
         vtrial0(jfft  ,ispden)=vtrialg(2,ifft,ispden)
         vresid0(jfft-1,ispden)=vreswk(1,ifft)
         vresid0(jfft  ,ispden)=vreswk(2,ifft)
       else
         vtrialg(:,ifft,ispden)=vtrialg(:,ifft,ispden)+fact*vreswk(:,ifft)
       end if
     end do
   end do
   ABI_DEALLOCATE(vreswk)
 end if

  call timab(902,1,tsec)

!Choice of preconditioner governed by iprcel, iprcch and iprcfc
 ABI_ALLOCATE(vrespc,(ispmix*nfftmix,dtset%nspden))
 ABI_ALLOCATE(vpaw,(npawmix*usepaw))
 if (usepaw==1)  then
   ABI_ALLOCATE(rhoijrespc,(npawmix))
 end if
!in case we are using a tfvw based preconditioner compute the local average of the non local potential
 if(dtset%userid==999) then
   compute_lavnlr = 1
   use_lavnlr=1
 else
   compute_lavnlr = 0
   use_lavnlr=0
 end if
 ABI_ALLOCATE(lavnlr,(dtset%nfft,dtset%nspden*use_lavnlr))
 if(compute_lavnlr == 1) then
   call lavnl(atindx,atindx1,cg,dtfil,dtset,eigen,&
&   kg,lavnlr,mcg,mpi_enreg,&
&   nattyp,&
&   npwarr,dtset%nspinor,&
&   occ,&
&   ph1d,psps,rhor,rprimd,&
&   wffnow,xred,ylm)
 end if

  call timab(902,2,tsec)
  call timab(903,1,tsec)

 call prcref_PMA(atindx,dielar,dielinv,dielstrt,dtn_pc,dtset,fcart,ffttomix,gmet,gsqcut,&
& istep,kg_diel,kxc,lavnlr,mgfft,moved_atm_inside,mpi_enreg,&
& nattyp,nfft,nfftmix,ngfft,ngfftmix,nkxc,npawmix,npwdiel,ntypat,n1xccc,&
& ispmix,0,pawrhoij,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,susmat,&
& vhartr,vpsp,vresid0,vrespc,vxc,xred,&
& deltae,efermi,etotal,nfftf,nhat,nhatgr,nhatgrdim,optene,optxc,&
& pawtab,usexcnhat,use_lavnlr,vtrial,wvl)

  call timab(903,2,tsec)
  call timab(904,1,tsec)

!In case of Thomas-Fermi-von Weizsaecker charge mixing save the old trial potential for
!later use

 if(dtset%iprctfvw /= 0) then !for tfw mixing
   ABI_ALLOCATE(vin_old,(nfft,dtset%nspden))
   ABI_ALLOCATE(vout_unmixed,(nfft,dtset%nspden))
   vin_old(:,:)=vtrial(:,:)
!  save the output potential before mixing for TFW charge mixing correction
!  the utput potential is the sum of the input potential and the preconditionned residual.
   vout_unmixed(:,:)=vtrial(:,:)+vrespc(:,:)
 end if

!------Compute new vtrial and eventual new atomic positions

 i_vresid1=mix%i_vresid(1)
 i_vrespc1=mix%i_vrespc(1)

!Initialise working arrays for the mixing object.
 if (moved_atm_inside == 1) then
   call ab6_mixing_use_moving_atoms(mix, dtset%natom, xred, dtn_pc)
 end if
 call ab6_mixing_eval_allocate(mix, istep)
!Copy current step arrays.
 if (moved_atm_inside == 1) then
   call ab6_mixing_copy_current_step(mix, vresid0, errid, message, &
&   arr_respc = vrespc, arr_atm = grhf)
 else
   call ab6_mixing_copy_current_step(mix, vresid0, errid, message, &
&   arr_respc = vrespc)
 end if
 if (errid /= AB6_NO_ERROR) then
   call wrtout(std_out, message, 'COLL')
   call leave_new('COLL')
 end if
 ABI_DEALLOCATE(vresid0)

!PAW: either use the array f_paw or the array f_paw_disk
 if (usepaw==1) then
   indx=-dplex
   do iatom=1,dtset%natom
     do ispden=1,pawrhoij(iatom)%nspden
       ABI_ALLOCATE(rhoijtmp,(cplex*pawrhoij(iatom)%lmn2_size,1))
       rhoijtmp=zero
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=cplex*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
         rhoijtmp(klmn:klmn+dplex,1)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
         jrhoij=jrhoij+cplex
       end do
       do kmix=1,pawrhoij(iatom)%lmnmix_sz
         indx=indx+cplex;klmn=cplex*pawrhoij(iatom)%kpawmix(kmix)-dplex
         vpaw(indx:indx+dplex)=rhoijtmp(klmn:klmn+dplex,1)-pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,ispden)
         mix%f_paw(indx:indx+dplex,i_vresid1)=pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,ispden)
         mix%f_paw(indx:indx+dplex,i_vrespc1)=rhoijrespc(indx:indx+dplex)
       end do
       ABI_DEALLOCATE(rhoijtmp)
     end do
   end do
 end if

!------Treat the mean of potentiel residual

!Special care must be taken with components of the
!potential that are associated with NO density change.
!In general, only the global mean of the potential has
!such an anomalous feature. However, in the spin
!polarized cas with fixed occupancies, also the
!mean of each spin-potential (independently of the other)
!has such a behaviour. The trick is to remove these
!variables before going in the predictive routines,
!then to put them back

 if(dtset%iprctfvw/=0) then
!  Compute the mean of the old vtrial
   call mean_fftr(vtrial,vmean,mpi_enreg,nfft,nfftot,dtset%nspden)
   call mean_fftr(vin_old,vmean_inold,mpi_enreg,nfft,nfftot,dtset%nspden)
   call mean_fftr(vout_unmixed,vmean_um,mpi_enreg,nfft,nfftot,dtset%nspden)
!  When spin-polarized and fixed occupation numbers,
!  treat separately spin up and spin down.
!  Otherwise, use only global mean
   do ispden=1,dtset%nspden
     if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&     abs(dtset%fixmom+99.99_dp)<1.0d-10)then
       vme=(vmean(1)+vmean(2))*half
       vme_inold=(vmean_inold(1)+vmean_inold(2))*half
       vme_um=(vmean_um(1)+vmean_um(2))*half
     else
       vme=vmean(ispden)
       vme_inold=vmean_inold(ispden)
       vme_um=vmean_um(ispden)
     end if
     vtrial(:,ispden)=vtrial(:,ispden)-vme
     vin_old(:,ispden)=vin_old(:,ispden)-vme_inold
     vout_unmixed(:,ispden)=vout_unmixed(:,ispden)-vme_um
   end do
 else
   call mean_fftr(vtrial,vmean,mpi_enreg,nfft,nfftot,dtset%nspden)
   do ispden=1,dtset%nspden
     if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&     abs(dtset%fixmom+99.99_dp)<1.0d-10)then
       vme=(vmean(1)+vmean(2))*half
     else
       vme=vmean(ispden)
     end if
     vtrial(:,ispden)=vtrial(:,ispden)-vme
   end do
 end if

!------Prediction of the components of the potential associated
!with a density change

!Init mpi_comm
 if(mpi_enreg%paral_compil_fft==1)then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,mpi_comm,spaceComm_bandfft=mpi_enreg%comm_fft)
   mpi_enreg%paral_level=old_paral_level
   mpi_summarize=.true.
 else
   mpi_comm=0
   mpi_summarize=.false.
 end if

 reset = .false.
 if (initialized == 0) reset = .true.
 call ab6_mixing_eval(mix, vtrial0, istep, nfftot, ucvol, &
& mpi_comm, mpi_summarize, errid, message, &
& reset = reset, isecur = dtset%isecur, &
& pawopt = dtset%pawoptmix, pawarr = vpaw, etotal = etotal, potden = rhor)
 if (errid == AB6_ERROR_MIXING_INC_NNSLOOP) then
   dbl_nnsclo = 1
 else if (errid /= AB6_NO_ERROR) then
   call wrtout(std_out, message, 'COLL')
   call leave_new('COLL')
 end if

!PAW: apply a simple mixing to rhoij (this is temporary)
 if(dtset%iscf==5 .or. dtset%iscf==6)then
   if (usepaw==1) then
     indx=-dplex
     do iatom=1,dtset%natom
       ABI_ALLOCATE(rhoijtmp,(cplex*pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
       rhoijtmp=zero
       if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
         do ispden=1,pawrhoij(iatom)%nspden
           do kmix=1,pawrhoij(iatom)%lmnmix_sz
             indx=indx+cplex;klmn=cplex*pawrhoij(iatom)%kpawmix(kmix)-dplex
             rhoijtmp(klmn:klmn+dplex,ispden)=rhoijrespc(indx:indx+dplex) &
&             -pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,ispden)
           end do
         end do
       end if
       do ispden=1,pawrhoij(iatom)%nspden
         jrhoij=1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=cplex*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
           rhoijtmp(klmn:klmn+dplex,ispden)=rhoijtmp(klmn:klmn+dplex,ispden) &
&           +pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
           jrhoij=jrhoij+cplex
         end do
       end do
       nselect=0
       do klmn=1,pawrhoij(iatom)%lmn2_size
         if (any(abs(rhoijtmp(cplex*klmn-dplex:cplex*klmn,:))>tol10)) then
           nselect=nselect+1
           pawrhoij(iatom)%rhoijselect(nselect)=klmn
           do ispden=1,pawrhoij(iatom)%nspden
             pawrhoij(iatom)%rhoijp(cplex*nselect-dplex:cplex*nselect,ispden)=&
&             rhoijtmp(cplex*klmn-dplex:cplex*klmn,ispden)
           end do
         end if
       end do
       pawrhoij(iatom)%nrhoijsel=nselect
       ABI_DEALLOCATE(rhoijtmp)
     end do
   end if
 end if

 if (usepaw==1)  then
   ABI_DEALLOCATE(rhoijrespc)
 end if

!PAW: restore rhoij from compact storage
 if (usepaw==1.and.dtset%iscf/=5.and.dtset%iscf/=6) then
   indx=-dplex
   do iatom=1,dtset%natom
     ABI_ALLOCATE(rhoijtmp,(cplex*pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
     rhoijtmp=zero
     if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
       do ispden=1,pawrhoij(iatom)%nspden
         jrhoij=1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=cplex*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
           rhoijtmp(klmn:klmn+dplex,ispden)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
           jrhoij=jrhoij+cplex
         end do
       end do
     end if
     do ispden=1,pawrhoij(iatom)%nspden
       do kmix=1,pawrhoij(iatom)%lmnmix_sz
         indx=indx+cplex;klmn=cplex*pawrhoij(iatom)%kpawmix(kmix)-dplex
         rhoijtmp(klmn:klmn+dplex,ispden)=vpaw(indx:indx+dplex)
       end do
     end do
     nselect=0
     if (cplex==1) then
       do klmn=1,pawrhoij(iatom)%lmn2_size
         if (any(abs(rhoijtmp(klmn,:))>tol10)) then
           nselect=nselect+1
           pawrhoij(iatom)%rhoijselect(nselect)=klmn
           do ispden=1,pawrhoij(iatom)%nspden
             pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoijtmp(klmn,ispden)
           end do
         end if
       end do
     else
       do klmn=1,pawrhoij(iatom)%lmn2_size
         if (any(abs(rhoijtmp(2*klmn-1:2*klmn,:))>tol10)) then
           nselect=nselect+1
           pawrhoij(iatom)%rhoijselect(nselect)=klmn
           do ispden=1,pawrhoij(iatom)%nspden
             pawrhoij(iatom)%rhoijp(2*nselect-1:2*nselect,ispden)=rhoijtmp(2*klmn-1:2*klmn,ispden)
           end do
         end if
       end do
     end if
     pawrhoij(iatom)%nrhoijsel=nselect
     ABI_DEALLOCATE(rhoijtmp)
   end do
 end if
 ABI_DEALLOCATE(vpaw)

!apply the Thomas--Fermi--von Weizsaecker charge mixing
!to avoid charge sloshing in large system

 if(dtset%iprctfvw /= 0.and.(istep.gt.1).and.(mod(istep,1)==0))  then
   if(dtset%iprctfvw==1) then
     call prctfvw1(atindx,atindx1,cg,deltae,dtfil,dtset,eeig,&
&     efermi,eigen,ek,enl,etotal,dtset%fixmom,gsqcut,&
&     kg,dtset%mband,mcg,mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,nattyp,nfft,nfftf,ngfftf,&
&     nhat,nhatgr,nhatgrdim,&
&     dtset%nkpt,nkxc,npwarr,dtset%nspden,dtset%nspinor,dtset%nsppol,psps%ntypat,n3xccc,occ,dtset%occopt,&
&     optene,optxc,pawfgr,&
&     ph1d,psps,resid,rhog,rhor,rprimd,&
&     usexcnhat,&
&     vin_old,vout_unmixed,vpsp,vtrial,&
&     wffnow,wvl,xccc3d,xred,ylm)
   else if(dtset%iprctfvw==2.and. istep.gt.0) then
     call prctfvw2(atindx,atindx1,cg,dtfil,dtset,eeig,&
&     efermi,eigen,ek,enl,dtset%fixmom,gsqcut,&
&     kg,dtset%mband,mcg,mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,nattyp,nfft,nfftf,ngfftf,&
&     nhat,nhatgr,nhatgrdim,&
&     dtset%nkpt,nkxc,npwarr,dtset%nspden,dtset%nspinor,dtset%nsppol,psps%ntypat,n3xccc,occ,dtset%occopt,&
&     optene,optres,optxc,pawfgr,&
&     ph1d,psps,resid,rhog,rhor,rprimd,&
&     usexcnhat,&
&     vin_old,vpsp,vtrial,&
&     wffnow,wvl,xccc3d,xred,ylm)
   end if
 end if
 ABI_DEALLOCATE(lavnlr)
 ABI_DEALLOCATE(vrespc)

!Eventually write the data on disk and deallocate f_fftgr_disk
 call ab6_mixing_eval_deallocate(mix)

  call timab(904,2,tsec)

!Restore potential
 if (ispmix==1.and.nfft==nfftmix) then
   vtrial=vtrial0
 else if (nfft==nfftmix) then
   do ispden=1,dtset%nspden
     call fourdp(1,vtrial0(:,ispden),vtrial(:,ispden),+1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,tim_fourdp)
   end do
 else
   do ispden=1,dtset%nspden
     do ifft=1,nfftmix
       jfft=mixtofft(ifft)
       vtrialg(1,jfft,ispden)=vtrial0(2*ifft-1,ispden)
       vtrialg(2,jfft,ispden)=vtrial0(2*ifft  ,ispden)
     end do
     call fourdp(1,vtrialg(:,:,ispden),vtrial(:,ispden),+1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,tim_fourdp)
   end do
   ABI_DEALLOCATE(vtrialg)
 end if
 ABI_DEALLOCATE(vtrial0)

  call timab(905,1,tsec)

!------Treat the mean of the potential

!Compute the mean of the new vtrial
 call mean_fftr(vtrial,vmean,mpi_enreg,nfft,nfftot,dtset%nspden)

!Reset the mean of the new vtrial, to the value vnew_mean
!When spin-polarized and fixed occupation numbers,
!treat separately spin up and spin down.
!Otherwise, use only global mean
 do ispden=1,dtset%nspden
   if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&   abs(dtset%fixmom+99.99_dp)<1.0d-10)then
     vme=(vnew_mean(1)+vnew_mean(2)-vmean(1)-vmean(2))*half
   else
     vme=vnew_mean(ispden)-vmean(ispden)
   end if
   vtrial(:,ispden)=vtrial(:,ispden)+vme
 end do

 if(moved_atm_inside==1 .and. istep/=nstep )then
   if(abs(dtset%iprcch)==1.or.abs(dtset%iprcch)==4)then
!    Subtract current local psp, but also vxc (for core charges)
     do ispden=1,dtset%nspden
       vtrial(:,ispden)=vtrial(:,ispden)-vpsp(:)*identity(ispden)-vxc(:,ispden)
     end do
   else if(abs(dtset%iprcch)==2.or.abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6)then
!    Subtract current vpsp+Hxc from vtrial. This should be rationalized later
     do ispden=1,dtset%nspden
       vtrial(:,ispden)=vtrial(:,ispden)-(vpsp(:)+vhartr(:))*identity(ispden)-vxc(:,ispden)
     end do
   end if
 end if

 call timab(905,2,tsec)
 call timab(93,2,tsec)

!DEBUG
!write(std_out,*)' newvtr : exit '
!write(std_out,*)' newvtr : vtrial(1,:)=',vtrial(1,:)
!stop
!ENDDEBUG
!DEBUG
!deallocate(dvstar,ldvstar)
!deallocate(vres,lvres,rdiemac,rdiemacg)
!deallocate(rdiemacf1,rdiemacf2,buffer1,buffer2)
!deallocate(g2cart,irdiemac,irdiemacg,irdiemacf1)
!deallocate(irdiemacf2)
!ENDDEBUG

end subroutine newvtr
!!***
