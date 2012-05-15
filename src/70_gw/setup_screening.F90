!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_screening
!! NAME
!! setup_screening
!!
!! FUNCTION
!!  Initialize the Ep% data type containing the parameters used during the screening calculation.
!!  as well as basic objects describing the BZ sampling .... TODO list to be completed
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! ikss_fname=Name of the input KSS file.
!! acell(3)=length scales of primitive translations (Bohr).
!! rprim(3,3)=dimensionless real space primitive translations.
!! ngfftf(18)=Contain all needed information about the 3D FFT for densities and potentials.
!!
!! OUTPUT
!! ngfft_gw(18)=Contain all needed information about the 3D FFT for the oscillator strengths.
!!  See ~abinit/doc/input_variables/vargs.htm#ngfft
!! Ltg_q(:)<Little_group>,=
!! Ep<Epsilonm1_parameters>=Parameters for the screening calculation. 
!!  Most part of it is Initialized and checked.  
!! Hdr_kss type(Hdr_type)=Header of the KSS file.
!! Cryst<Crystal_structure>=Definition of the unit cell and its symmetries.
!! Kmesh<BZ_mesh_type>=Structure defining the k-point sampling (wavefunctions).
!! Qmesh<BZ_mesh_type>=Structure defining the q-point sampling (screening)
!! Gsph_wfn<Gvectors_type>=Structure defining the G-sphere for the wavefunctions (not k-dependent).
!! Gsph_epsG0<Gvectors_type>=The G-sphere for the screening, enlarged to take into account for umklapps.
!! Psps <Pseudopotential_type)>=Info on pseudopotential, only for consistency check of the KSS file 
!! Vcp <type vcoul_t> datatype gathering information on the coulombian cutoff technique
!! comm=MPI communicator.
!!
!! SIDE EFFECTS
!! Dtset<Dataset_type>=All input variables for this dataset.
!!  %ecutwfn, %npwwfn, %nshwfn,
!!  %ecuteps, %npweps, %nsheps
!!   might be redefinend in setshells in order to close the shell.
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      bstruct_init,destroy_gsphere,find_qmesh,get_bz_item,get_ng0sh,hdr_check
!!      hdr_init,hdr_update,hdr_vs_dtset,init_crystal_from_hdr,init_gsphere
!!      init_kmesh,initmpi_seq,merge_and_sort_kg,metric,mkrdim
!!      nullify_epsilonm1_parameters,nullify_little_group,print_bz_mesh
!!      print_crystal,rhoij_alloc,rhoij_copy,rhoij_free,setmesh,setshells
!!      setup_little_group,testkss,vcoul_init,wfk_read_ene,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setup_screening(codvsn,acell,rprim,ngfftf,ikss_fname,Dtset,Psps,Pawtab,&
& ngfft_gw,Hdr_kss,Hdr_out,Cryst,Kmesh,Qmesh,KS_BSt,Ltg_q,Gsph_epsG0,Gsph_wfn,Vcp,Ep,comm)

 use m_profiling
    
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_errors

 use m_gwdefs,        only : GW_TOLQ0, GW_TOLQ, GW_Q0_DEFAULT, czero_gw, epsilonm1_parameters, nullify_epsilonm1_parameters,&
&                            gw_uses_wfk_file      
 use m_geometry,      only : normv
 use m_header,        only : hdr_init
 use m_crystal,       only : print_crystal, crystal_structure
 use m_crystal_io,    only : init_crystal_from_hdr
 use m_bz_mesh,       only : bz_mesh_type, init_kmesh, get_ng0sh, print_bz_mesh, find_qmesh, get_BZ_item,&
&                            little_group, setup_little_group, nullify_little_group, make_mesh, destroy_BZ_mesh_type
 use m_ebands,        only : bstruct_init
 use m_vcoul,         only : vcoul_t, vcoul_init
 use m_fft_mesh,      only : setmesh
 use m_gsphere,       only : gvectors_type, init_gsphere, merge_and_sort_kg, destroy_gsphere
 use m_io_kss,        only : testkss

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_screening'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_70_gw, except_this_one => setup_screening
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=6),intent(in) :: codvsn
 character(len=fnlen),intent(in) :: ikss_fname
 type(Dataset_type),intent(inout) :: Dtset !INOUT is due to setshells
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
 type(Epsilonm1_parameters),intent(out) :: Ep 
 type(Hdr_type),intent(out) :: Hdr_kss,Hdr_out
 type(Bandstructure_type),intent(out) :: KS_BSt
 type(BZ_mesh_type),intent(out) :: Kmesh,Qmesh
 type(Crystal_structure),intent(out) :: Cryst
 type(Gvectors_type),intent(out) :: Gsph_epsG0,Gsph_wfn
 type(vcoul_t),intent(out) :: Vcp
!arrays
 integer,intent(in) :: ngfftf(18)
 integer,intent(out) :: ngfft_gw(18)
 real(dp),intent(in) :: acell(3),rprim(3,3)
 type(Little_group),pointer :: Ltg_q(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: NOMEGAGAUSS=30,NOMEGAREAL=201
 integer :: bantot,ib,ibtot,ikibz,iq,iqp,isppol,itypat,pinv
 integer :: jj,mod10,nsym_kss,nbnds_kss,ng_kss,iqbz,isym,iq_ibz,itim
 integer :: mpsang_kss,timrev,mg0sh,use_umklp,ierr
 integer :: npwepG0,nshepspG0,method,enforce_sym,nfftgw_tot !,spin,band,ik_ibz,test_ng
 integer :: pertcase_,restart,restartpaw,ii,istart,iend
 real(dp),parameter :: OMEGAERMAX=100.0/Ha_eV
 real(dp) :: ecutepspG0,ucvol,domegareal
 logical :: remove_inv,ltest,only_one_kpt,found,is_static,has_q0
 character(len=500) :: msg      
 character(len=fnlen) :: wfk_fname
 type(MPI_type) :: MPI_enreg_seq
 type(wvl_internal_type) :: wvl
!arrays
 integer :: ng0sh_opt(3)
 integer,allocatable :: npwarr(:),nlmn(:)
 integer,pointer :: gvec_p(:,:),shlim(:) !,test_gvec_p(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),qtmp(3),sq(3),qbz(3)
 real(dp),pointer :: energies_p(:,:,:)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:)
 !type(Gvectors_type) :: Gsph_chi0
 
! *************************************************************************

 DBG_ENTER('COLL')

 ! === Check for calculations that are not implemented ===
 ABI_CHECK(Dtset%mkmem/=0,'mkmem=0 not yet implemented.')
 ltest=ALL(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol)==Dtset%nband(1))
 ABI_CHECK(ltest,'Dtset%nband must be constant')

 !* Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)

 call mkrdim(acell,rprim,rprimd)  
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol) 

 ! === Set up basic parameters of the calculation ===
 call nullify_epsilonm1_parameters(Ep) 

 Ep%gwcalctyp            =Dtset%gwcalctyp
 Ep%plasmon_pole_model   =.TRUE.  
 Ep%analytic_continuation=.FALSE. 
 Ep%contour_deformation  =.FALSE.

 mod10=MOD(Ep%gwcalctyp,10)
 if (mod10/=0.and.mod10/=8) Ep%plasmon_pole_model   =.FALSE.
 if (mod10==1)              Ep%analytic_continuation=.TRUE.
 if (mod10==2.or.mod10==9)  Ep%contour_deformation  =.TRUE.
 is_static=(mod10==5.or.mod10==6.or.mod10==7)

 Ep%nbnds  =Dtset%nband(1) 
 Ep%symchi =Dtset%symchi
 Ep%inclvkb=Dtset%inclvkb; if (Dtset%usepaw/=0) Ep%inclvkb=0
 Ep%zcut   =Dtset%zcut   

 write(msg,'(2a,i4,2a,f10.6,a)')ch10,&
&  ' GW calculation type              = ',Ep%gwcalctyp,ch10,&
&  ' zcut to avoid poles in chi0 [eV] = ',Ep%zcut*Ha_eV,ch10
 call wrtout(std_out,msg,'COLL')
 !
 Ep%awtr  =Dtset%awtr 
 Ep%npwe  =Dtset%npweps 
 Ep%npwwfn=Dtset%npwwfn 
 Ep%npwvec=MAX(Ep%npwe,Ep%npwwfn)

 timrev = 2 ! This information is not reported in the header
            ! 1 --> do not use time-reversal symmetry 
            ! 2 --> take advantage of time-reversal symmetry
 if (timrev==1.and.Dtset%awtr/=0) then
   msg = "awtr/=0 cannot be used when time-reversal symmetry doesn't hold"
   MSG_ERROR(msg)
 end if

#if 0
!BEGIN GWWFK
 ! this piece of code will replace the setup of the G-sphere for Espilon
 pinv=1 !; if (timrev==2) pinv=-1 ! if time-reversal can be used, 1 otherwise
 nullify(test_gvec_p,shlim)
 call merge_and_sort_kg(Dtset%nkpt,Dtset%kptns,Dtset%ecuteps,Dtset%nsym,&
&    pinv,Dtset%symrel,gprimd,test_gvec_p,Dtset%prtvol,shlim_p=shlim)
 test_ng=SIZE(test_gvec_p,DIM=2)
 ABI_DEALLOCATE(shlim)
!END GWWFK
#endif
 !
 ! === Read parameters from (WFK|KSS) and verifify them === 
 if (gw_uses_wfk_file) then
 !if (.FALSE. .and. gw_uses_wfk_file) then ! Temporary hack so that G sphere is always read from the KSS file.
   wfk_fname = ikss_fname; ii=LEN_TRIM(wfk_fname)
   wfk_fname(ii-2:ii) = "WFK"

   MSG_WARNING("Using WFK file"//TRIM(wfk_fname))

   nbnds_kss=-1
   call wfk_read_ene(wfk_fname,Dtset%accesswff,nbnds_kss,energies_p,Hdr_kss,Dtset%prtvol,comm) 

   nsym_kss  = Hdr_kss%nsym
   mpsang_kss=Psps%mpsang

   pinv=1 !; if (timrev==2) pinv=-1 ! if time-reversal can be used, 1 otherwise
   nullify(gvec_p,shlim)
   call merge_and_sort_kg(Hdr_kss%nkpt,Hdr_kss%kptns,Hdr_kss%ecut_eff,Hdr_kss%nsym,&
&    pinv,Hdr_kss%symrel,gprimd,gvec_p,Dtset%prtvol,shlim_p=shlim)

   ng_kss=SIZE(gvec_p,DIM=2)
   ABI_DEALLOCATE(shlim)

 else
   call testkss(ikss_fname,Dtset%accesswff,nsym_kss,&
&    nbnds_kss,ng_kss,mpsang_kss,gvec_p,energies_p,Hdr_kss,comm)
 end if


#if 0
!BEGIN GWWFK
 jj = ng_kss
 if (ng_kss/=test_ng) then
   jj = MIN(ng_kss,test_ng)
!   MSG_WARNING("ng_kss/=test_ng")
!   write(std_out,*)ng_kss,test_ng
 end if

 if (ANY(test_gvec_p(:,1:jj)/=gvec_p(:,1:jj))) then
   MSG_WARNING("ANY(test_gvec/=gvec_p)")
   do ii=1,jj
!     if (ANY(gvec_p(:,ii)-test_gvec_p(:,ii)/=0)) &
!&      write(std_out,*)"gvec_p, test_gvec",ii,"/",jj,gvec_p(:,ii),test_gvec_p(:,ii)
   end do
   !MSG_ERROR("Cannot continue")
 end if

 !do spin=1,Hdr_kss%nsppol
 !  do ik_ibz=1,Hdr_kss%nkpt
 !     write(std_out,*)energies_p(1:nbnds_kss,ik_ibz,spin)
 !  end do
 !end do
!END GWWFK

 ABI_DEALLOCATE(test_gvec_p)
#endif

 ltest=(Psps%mpsang==mpsang_kss)
 ABI_CHECK(ltest,'Psps%mpsang does not agree with mpsang read from KSS file')

 ! === Get important dimension from Hdr_kss ===
 ! * Check also the consistency btw Hdr_kss and Dtset.
 Ep%nsppol=Hdr_kss%nsppol
 Ep%nkibz =Hdr_kss%nkpt

 call hdr_vs_dtset(Hdr_kss,Dtset) 
 remove_inv=(nsym_kss/=Hdr_kss%nsym) 
 
 call init_crystal_from_hdr(Cryst,Hdr_kss,timrev,remove_inv)
 call print_crystal(Cryst,mode_paral='COLL')

 if (Ep%npwvec>ng_kss) then
   Ep%npwvec=ng_kss
   if (Ep%npwwfn> ng_kss) Ep%npwwfn=ng_kss
   if (Ep%npwe  > ng_kss) Ep%npwe  =ng_kss
   write(msg,'(3a,3(a,i6,a))')ch10,&
&    ' Number of G-vectors found less then required. Calculation will proceed with ',ch10,&
&    '  npwvec = ',Ep%npwvec,ch10,&
&    '  npweps = ',Ep%npwe  ,ch10,&
&    '  npwwfn = ',Ep%npwwfn,ch10
   MSG_WARNING(msg)
 end if

 if (Ep%nbnds>nbnds_kss) then
   Ep%nbnds=nbnds_kss
   Dtset%nband(:)=nbnds_kss
   write(msg,'(4a,i4,a)')ch10,&
&    ' Number of bands found less then required. ',ch10,&
&    ' Calculation will proceed with nbnds = ',nbnds_kss,ch10
   MSG_WARNING(msg)
 end if

 ! === Create basic data types for the calculation ===
 ! * Kmesh defines the k-point sampling for the wavefunctions.
 ! * Qmesh defines the q-point sampling for chi0, all possible differences k1-k2 reduced to the IBZ. 
 ! TODO Kmesh%bz should be [-half,half[ but this modification will be painful!

 call init_kmesh(Kmesh,Cryst,Ep%nkibz,Hdr_kss%kptns,Dtset%kptopt,wrap_1zone=.FALSE.)
 !call init_kmesh(Kmesh,Cryst,Ep%nkibz,Hdr_kss%kptns,Dtset%kptopt,wrap_1zone=.TRUE.)

 call print_BZ_mesh(Kmesh,"K-mesh for the wavefunctions",std_out,Dtset%prtvol,"COLL")
 call print_BZ_mesh(Kmesh,"K-mesh for the wavefunctions",ab_out, 0,           "COLL")
 !
 ! === Find Q-mesh, and do setup for long wavelength limit ===
 ! * Stop if a nonzero umklapp is needed to reconstruct the BZ. In this case, indeed, 
 !   epsilon^-1(Sq) should be symmetrized in csigme using a different expression (G-G_o is needed)
 !
 call find_qmesh(Qmesh,Cryst,Kmesh)

 call print_BZ_mesh(Qmesh,"Q-mesh for the screening function",std_out,Dtset%prtvol,"COLL")
 call print_BZ_mesh(Qmesh,"Q-mesh for the screening function",ab_out ,0           ,"COLL")

 do iqbz=1,Qmesh%nbz
   call get_BZ_item(Qmesh,iqbz,qbz,iq_ibz,isym,itim)
   sq = (3-2*itim)*MATMUL(Cryst%symrec(:,:,isym),Qmesh%ibz(:,iq_ibz))
   if (ANY(ABS(qbz-sq )>1.0d-4)) then
     !write(std_out,*) sq,Qmesh%bz(:,iqbz) 
     write(msg,'(a,3f6.3,a,3f6.3,2a,9i3,a,i2,2a)')&
&      ' qpoint ',qbz,' is the symmetric of ',Qmesh%ibz(:,iq_ibz),ch10,&
&      ' through operation ',Cryst%symrec(:,:,isym),' and itim ',itim,ch10,&
&      ' however a non zero umklapp G_o vector is required and this is not yet allowed'
     MSG_ERROR(msg)
   end if
 end do 

 if (Dtset%gw_nqlwl==0) then
   Ep%nqlwl=1 
   ABI_ALLOCATE(Ep%qlwl,(3,Ep%nqlwl))
   Ep%qlwl(:,1)=GW_Q0_DEFAULT ! Use default shift 0.000010, 0.000020, 0.000030
 else 
   Ep%nqlwl=Dtset%gw_nqlwl
   ABI_ALLOCATE(Ep%qlwl,(3,Ep%nqlwl))
   Ep%qlwl(:,:)=Dtset%gw_qlwl(:,1:Ep%nqlwl)
   ABI_CHECK(Ep%nqlwl==1,"nqlwl/=1 not coded")
 end if
 !write(std_out,*)" Using qlwl = ",Ep%qlwl

 ! === Find optimal value for G-sphere enlargment due to oscillator matrix elements ===
 mg0sh=50
 call get_ng0sh(Kmesh%nbz,Kmesh%bz,Qmesh%nibz,Qmesh%ibz,Kmesh%nbz,Kmesh%bz,gmet,GW_TOLQ0,mg0sh,ng0sh_opt)
 Ep%mG0(:)=ng0sh_opt(:)
 !Ep%mG0(:)=(/3,3,3/) 

 ! === In case of symmetrization, find the little group of the q"s ===
 ! * For the long-wavelength limit we consider a small but finite q. However the oscillators are 
 !  evaluated setting q==0. Thus it is possible to take advantage of symmetries also when q --> 0.
 ! * Here we calculate the enlargement of the G-sphere, npwepG0, needed to account for umklapps.
 ! TODO Switch on use_umklp, write all this stuff to ab_out

 Ep%npwepG0=Ep%npwe
 ABI_ALLOCATE(Ltg_q,(Qmesh%nibz))
 call nullify_little_group(Ltg_q) 

 do iq=1,Qmesh%nibz 
   qtmp(:)=Qmesh%ibz(:,iq); if (normv(qtmp,gmet,'G')<GW_TOLQ0) qtmp(:)=zero; use_umklp=0
   call setup_little_group(qtmp,Kmesh,Cryst,use_umklp,Ltg_q(iq),Ep%npwe,gvec=gvec_p)
 end do

 if (Ep%symchi/=0) then
   ecutepspG0=MAXVAL(Ltg_q(:)%max_kin_gmG0)+tol6; npwepG0=0; nshepspG0=0
   write(std_out,*)" Due to umklapp processes : ecutepspg0= ",ecutepspG0
   call setshells(ecutepspG0,npwepG0,nshepspG0,Cryst%nsym,gmet,gprimd,Cryst%symrel,'eps_pG0',Cryst%ucvol)
   Ep%npwepG0=npwepG0
 end if

 if (Ep%npwepG0>Ep%npwvec) then
   write(msg,'(3a,i5,a,i5)')&
&    ' npwepG0 > npwvec, decrease npweps or increase npwwfn. ',ch10,&
&    ' npwepG0 = ',Ep%npwepG0,' npwvec = ',Ep%npwvec
   MSG_ERROR(msg)
 end if
 !
 ! =======================================================================
 ! ==== Setup of the FFT mesh used for the oscillator matrix elements ====
 ! =======================================================================
 ! * ngfft_gw(7:18) is the same as Dtset%ngfft(7:18), initialized before entering setup_screening. 
 !   Here we just redefine ngfft_gw(1:6) according to the following options:
 !
 !    method==0 ==> FFT grid read from __fft.in__ (only for debugging purpose)
 !    method==1 ==> normal FFT grid 
 !    method==2 ==> slightly augmented FFT grid to calculate exactly rho_tw_g (see setmesh.F90)
 !    method==3 ==> doubled FFT grid, to treat exactly the convolution defining the density,
 !      Useful in sigma if ppmodel=[2,3,4] since rho(G-Gp) or to calculate matrix elements of v_Hxc.
 !
 !    enforce_sym==1 ==> enforce a direct space FFT mesh compatible with all symmetries operation
 !    enforce_sym==0 ==> Find the smallest FFT grid compatibile with the library, do not care about symmetries
 !
 ngfft_gw(1:18)=Dtset%ngfft(1:18); method=2
 if (Dtset%fftgw==00 .or. Dtset%fftgw==01) method=0
 if (Dtset%fftgw==10 .or. Dtset%fftgw==11) method=1
 if (Dtset%fftgw==20 .or. Dtset%fftgw==21) method=2
 if (Dtset%fftgw==30 .or. Dtset%fftgw==31) method=3
 enforce_sym=MOD(Dtset%fftgw,10) 

 ! npwepG0 to account for umklapps.
 call setmesh(gmet,gvec_p,ngfft_gw,Ep%npwvec,Ep%npwepG0,Ep%npwwfn,nfftgw_tot,method,Ep%mG0,Cryst,enforce_sym)

 ! === Create structure describing the G-sphere used for chi0/espilon and Wfns ===
 ! * The cutoff is >= ecuteps to allow for umklapp
 ! MG, Fabien, I modified setup_G_rotation to speed up the loops, 
 ! sincerely I dont like only_one_kpt because some entries in Gsphere are not correctly filled!
 only_one_kpt=(Kmesh%nbz==1) 
 call init_gsphere(Gsph_epsG0,only_one_kpt,Cryst,Ep%npwepG0,gvec=gvec_p)
 call init_gsphere(Gsph_wfn  ,only_one_kpt,Cryst,Ep%npwvec ,gvec=gvec_p)

 ABI_DEALLOCATE(gvec_p)

 ! TODO use this when WFK files are used.
 ierr=0
#if 0
 call init_gsphere(Gsph_chi0,only_one_kpt,Cryst,0,ecut=ecutepspG0)


 ii = MIN(Gsph_chi0%ng,Gsph_epsG0%ng) 
 do jj=1,ii
   if (ANY(Gsph_epsG0%gvec(:,jj)/=Gsph_chi0%gvec(:,jj))) then
     ierr=ierr+1
     write(std_out,*)" Gsph_epsG0, Gsph_chi0 ",jj,"/",ii,Gsph_epsG0%gvec(:,jj),Gsph_chi0%gvec(:,jj)
   end if
 end do

 ABI_CHECK(ierr==0,"Mismatch between Gsph_chi0 and Gsph_epsG0")
 call destroy_gsphere(Gsph_chi0)
#endif

 ! FIXME this wont work if nqptdm/=0
 call vcoul_init(Vcp,Gsph_epsG0,Qmesh,Kmesh,Dtset%rcut,Dtset%icutcoul,Dtset%vcutgeo,Ep%npwe,Ep%nqlwl,&
&  Ep%qlwl,Cryst%rprimd,ngfftf,comm)

#if 0
 ! Using the random q for the optical limit is one of the reasons
 ! why sigma breaks the initial energy degeneracies.
 Vcp%i_sz=zero
 Vcp%vc_sqrt(1,:)=czero
 Vcp%vcqlwl_sqrt(1,:)=czero
#endif

 ! Value of scissor energy
 Ep%soenergy=Dtset%soenergy

 ! === Define the frequency mesh for epsilon according to the method used ===
 Ep%nomegaei=1 
 Ep%nomegaer=1; if (is_static) Ep%nomegaer=0
 Ep%nomegaec=0
 Ep%omegaermax=zero

 ! === For ppmodels 2,3,4, only omega=0 is needed ===
 if (Ep%plasmon_pole_model.and.Dtset%nfreqre==1.and.Dtset%nfreqim==0) then
   Ep%nomegaer=1; Ep%nomegaei=0
   write(msg,'(7a)')ch10,&
&    ' The inverse dielectric matrix will be calculated on zero frequency only',ch10,&
&    ' please note that the calculated epsilon^-1 cannot be used ',ch10,&
&    ' to calculate QP corrections using plasmonpole model 1',ch10
   call wrtout(std_out,msg,'COLL') 
   call wrtout(ab_out,msg,'COLL')
 end if 
 !
 ! === Max number of omega along the imaginary axis ===
 if (Ep%analytic_continuation.or.Ep%contour_deformation) then
   Ep%nomegaei=Dtset%nfreqim
   if (Dtset%cd_custom_imfrqs/=0) then
     write(msg,'(a)')' Custom imaginary grid specified. Assuming experienced user.'
     MSG_WARNING(msg)
     Ep%nomegaei=Dtset%cd_custom_imfrqs
   end if
   if (Ep%nomegaei==-1) then 
     Ep%nomegaei=NOMEGAGAUSS
     write(msg,'(a,i5)')' Number of imaginary frequencies set to default= ',NOMEGAGAUSS
     MSG_WARNING(msg)
   end if
   if (Ep%nomegaei==0) then 
     write(msg,'(a)')' nfreqim = 0 ! Assuming experienced user merging several frequency calculations.'
     MSG_WARNING(msg)
   end if
 end if

 ! === Range and total number of real frequencies ===
 Ep%omegaermin = zero
 if (Ep%contour_deformation) then
   Ep%nomegaer=Dtset%nfreqre; Ep%omegaermin=Dtset%freqremin; Ep%omegaermax=Dtset%freqremax
   if (Dtset%cd_use_tangrid==1) then
     Ep%omegaermax=Dtset%cd_max_freq
     write(msg,'(a)')' Tangent transfom grid will be used '
     MSG_WARNING(msg)
   end if
   if (Ep%nomegaer==-1) then 
     Ep%nomegaer=NOMEGAREAL
     write(msg,'(a,i5)')' Number of real frequencies set to default= ',NOMEGAREAL
     MSG_WARNING(msg)
   end if
   if (Ep%nomegaer==0) then 
     write(msg,'(a)')'  nfreqre = 0 ! Assuming experienced user merging several frequency calculations.'
     MSG_WARNING(msg)
   end if
   if (ABS(Ep%omegaermin)<TOL16) then 
     Ep%omegaermin=zero
     write(msg,'(a,f8.4)')' Min real frequency set to default [Ha] = ',Ep%omegaermin
     MSG_WARNING(msg)
   end if
   if (Ep%omegaermin>Ep%omegaermax) then
     MSG_ERROR(' freqremin > freqremax !')
   end if
   if (Ep%omegaermax<TOL16) then 
     Ep%omegaermax=OMEGAERMAX
     write(msg,'(a,f8.4)')' Max real frequency set to default [Ha] = ',OMEGAERMAX
     MSG_WARNING(msg)
   end if 
   ! Check if a subset of the frequencies is to be used
   if (Dtset%cd_subset_freq(1)/=0) then
     istart = Dtset%cd_subset_freq(1)
     iend   = Dtset%cd_subset_freq(2)
     if (istart>iend.or.istart<0.or.iend<0) then
       MSG_ERROR(' check indices of cd_subset_freq!')
     end if
     write(msg,'(2(a,i0))')' Using cd_subset_freq to only do freq. from ',istart,' to ',iend
     MSG_WARNING(msg)
     ! Reset the numbers
     if (Dtset%cd_use_tangrid/=1) then ! Normal equidistant grid 
       Ep%nomegaer = iend-istart+1
       domegareal=(Ep%omegaermax-Ep%omegaermin)/(Ep%nomegaer-1)
       Ep%omegaermin = Ep%omegaermin+(istart-1)*domegareal
       Ep%omegaermax = Ep%omegaermin+(iend-1)*domegareal
     else
       Ep%nomegaer = iend-istart+1
     end if
   end if
 end if

 ! Check full grid calculations
 if (Dtset%cd_full_grid/=0) then
   MSG_WARNING(' FULL GRID IN COMPLEX PLANE CALCULATED.')
   MSG_WARNING(' YOU MIGHT NOT BE ABLE TO USE SCREENING FILES!')
   if (Dtset%cd_subset_freq(1)/=0) then
     MSG_ERROR(' cd_subset_freq cannot be used with cd_full_grid!')
   end if
   Ep%nomegaec = Ep%nomegaei*(Ep%nomegaer-1)
 end if

 Ep%nomega=Ep%nomegaer+Ep%nomegaei+Ep%nomegaec ! Total number of frequencies.

 ! ==== Setup of the spectral method ====
 Ep%spmeth  =Dtset%spmeth 
 Ep%nomegasf=Dtset%nomegasf 
 Ep%spsmear =Dtset%spbroad

 if (Ep%spmeth/=0) then 
   write(msg,'(2a,i3,2a,i8)')ch10,&
&    ' setup_screening : using spectral method = ',Ep%spmeth,ch10,&
&    '  Number of frequencies for imaginary part= ',Ep%nomegasf
   call wrtout(std_out,msg,'COLL')
   if (Ep%spmeth==2) then 
     write(msg,'(a,f8.5,a)')' gaussian broadening = ',Ep%spsmear*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
   end if  
 end if

 Ep%nI=1; Ep%nJ=1 
 if (Dtset%nspinor==2) then 
   if (Dtset%usepaw==1.and.Dtset%pawspnorb>0) then
     Ep%nI=1; Ep%nJ=4 
   end if
   ! For spin-spin interaction
   ! Ep%nI=4; Ep%nJ=4 
   if (Ep%npwepG0/=Ep%npwe) STOP "If spinor npwepG0 must be == npwe"
   if (Ep%symchi/=0) STOP "symchi/0 and spinor not available"
 end if
 !
 ! === Enable the calculations of chi0 on user-specified q-points ===
 Ep%nqibz=Qmesh%nibz
 ABI_ALLOCATE(Ep%qibz,(3,Ep%nqibz))
 Ep%qibz(:,:)=Qmesh%ibz(:,:)

 Ep%nqcalc=Ep%nqibz
 if (Dtset%nqptdm>0) Ep%nqcalc=Dtset%nqptdm

 ABI_ALLOCATE(Ep%qcalc,(3,Ep%nqcalc))
 if (Ep%nqcalc/=Ep%nqibz) then
   write(msg,'(6a)')ch10,&
&    ' Dielectric matrix will be calculated only for some ',ch10,&
&    ' selected q points provided by the user through the input variables ',ch10,&
&    ' nqptdm and qptdm'
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   ltest=(Ep%nqcalc<=Qmesh%nibz)
   msg = 'nqptdm should not exceed the number of q points in the IBZ'
   ABI_CHECK(ltest,msg)
   Ep%qcalc(:,:)=Dtset%qptdm(:,1:Ep%nqcalc)
   do iq=1,Ep%nqcalc ! * Check whether the q-points provided are correct.
     found=.FALSE.
     do iqp=1,Qmesh%nibz
       qtmp(:)=Ep%qcalc(:,iq)-Qmesh%ibz(:,iqp) 
       found=(normv(qtmp,gmet,'G')<GW_TOLQ) 
       if (found) EXIT
     end do
     msg = 'One or more points specified by Dtset%qptdm do not satisfy q=k1-k2'
     ABI_CHECK(found,msg)
   end do 
 else 
   Ep%qcalc(:,:)=Ep%qibz(:,:)
 end if 

 ! To write the SCR header correctly, with heads and wings, we have 
 ! to make sure that q==0, if present, is the first q-point in the list.
 !$has_q0=(ANY(normv(Ep%qcalc(:,:),gmet,'G')<GW_TOLQ0)) !commented to avoid problems with sunstudio12
 has_q0=.FALSE.
 do iq=1,Ep%nqcalc
   if (normv(Ep%qcalc(:,iq),gmet,'G')<GW_TOLQ0) then
     has_q0=.TRUE.; EXIT
   end if
 end do

 if (has_q0.and.normv(Ep%qcalc(:,1),gmet,'G')>=GW_TOLQ0) then
   write(msg,'(5a)')&
&    ' The list of q-points to be calculated contains the Gamma point, ',ch10,&
&    ' however Gamma is not the first point in the list. ' ,ch10,&
&    ' Please, change your input file accordingly. '
   MSG_ERROR(msg)
 end if

 ! === Initialize the band structure datatype ===
 ! * Copy KSS energies and occupations up to Ep%nbnds==Dtset%nband(:)
 ! TODO Recheck symmorphy and inversion

 bantot=SUM(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol))

 ABI_ALLOCATE(doccde,(bantot))
 ABI_ALLOCATE(eigen,(bantot))
 ABI_ALLOCATE(occfact,(bantot))
 doccde(:)=zero; eigen(:)=zero; occfact(:)=zero 

 jj=0; ibtot=0
 do isppol=1,Dtset%nsppol
   do ikibz=1,Dtset%nkpt
     do ib=1,Hdr_kss%nband(ikibz+Dtset%nkpt*(isppol-1))
       ibtot=ibtot+1
       if (ib<=Ep%nbnds) then 
         jj=jj+1
         occfact(jj)=Hdr_kss%occ(ibtot)
         eigen  (jj)=energies_p(ib,ikibz,isppol)
       end if
     end do
   end do
 end do
 ABI_DEALLOCATE(energies_p)

 ! Make sure that Dtset%wtk==Kmesh%wt due to the dirty treatment of 
 ! the symmetry operations in the old GW code (symmorphy and inversion) 
 ltest=(ALL(ABS(Dtset%wtk(1:Kmesh%nibz)-Kmesh%wt(1:Kmesh%nibz))<tol6))
 ABI_CHECK(ltest,'Mismatch between Dtset%wtk and Kmesh%wt')

 ABI_ALLOCATE(npwarr,(Hdr_kss%nkpt))
 npwarr(:)=Ep%npwwfn 

 call bstruct_init(bantot,KS_BSt,Dtset%nelect,doccde,eigen,Dtset%istwfk,Kmesh%ibz,Dtset%nband,&
& Kmesh%nibz,npwarr,Dtset%nsppol,Dtset%nspinor,Dtset%tphysel,Dtset%tsmear,Dtset%occopt,occfact,Kmesh%wt) 

 ! TODO modify outkss in order to calculate the eigenvalues also if NSCF calculation.
 ! this fails simply because in case of NSCF occ  are zero
 !ltest=(ALL(ABS(occfact-KS_BSt%occ)<1.d-2)) 
 !call assert(ltest,'difference in occfact')
 !write(std_out,*)MAXVAL(ABS(occfact(:)-KS_BSt%occ(:))) 

 !TODO call update_occ here
 !$call update_occ(KS_BSt,fixmom,stmbias,Dtset%prtvol)

 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(eigen)
 ABI_DEALLOCATE(npwarr)

 ! === Initialize abinit header for the screening part ===
 pertcase_=0
 call hdr_init(KS_BSt,codvsn,Dtset,Hdr_out,Pawtab,pertcase_,Psps,wvl)

 ! === Get Pawrhoij from the header of the KSS file ===
 ABI_ALLOCATE(Pawrhoij,(Cryst%natom*Dtset%usepaw))
 if (Dtset%usepaw==1) then
   ABI_ALLOCATE(nlmn,(Cryst%ntypat))
   do itypat=1,Cryst%ntypat
     nlmn(itypat)=Pawtab(itypat)%lmn_size
   end do
   call rhoij_alloc(1,nlmn,Dtset%nspden,Dtset%nspinor,Dtset%nsppol,Pawrhoij,Cryst%typat,MPI_enreg=MPI_enreg_seq)
   ABI_DEALLOCATE(nlmn)
   call rhoij_copy(Hdr_kss%Pawrhoij,Pawrhoij,MPI_enreg=MPI_enreg_seq)
 end if

 call hdr_update(bantot,1.0d20,1.0d20,Hdr_out,Cryst%natom,1.0d20,&
&  Cryst%rprimd,occfact,MPI_enreg_seq,Pawrhoij,Dtset%usepaw,Cryst%xred)

! call hdr_update(bantot,1.0d20,1.0d20,Hdr_out,Cryst%natom,1.0d20,&
!& Cryst%rprimd,KS_BSt%occ,MPI_enreg_seq,Pawrhoij,Dtset%usepaw,Cryst%xred)

 ! This is just to do a check, the file format is wrong!
 call hdr_check(1002,1002,Hdr_out,Hdr_kss,'COLL',restart,restartpaw)

 ABI_DEALLOCATE(occfact)
 call rhoij_free(Pawrhoij)
 ABI_DEALLOCATE(Pawrhoij)
 !
 ! ==== Setup of extrapolar technique ====
 Ep%gwcomp   = Dtset%gwcomp
 Ep%gwencomp = Dtset%gwencomp

 if (Ep%gwcomp==1) then
  write(msg,'(a,f8.2,a)')' Using the completeness correction with gwencomp ',Ep%gwencomp*Ha_eV,' [eV] '
  call wrtout(std_out,msg,'COLL')
 end if

 ! === Final compatibility tests ===
 if (ANY(KS_BSt%istwfk/=1)) then
   MSG_WARNING('istwfk/=1 is still under development')
 end if

 ltest=(KS_BSt%mband==Ep%nbnds.and.ALL(KS_BSt%nband==Ep%nbnds))
 ABI_CHECK(ltest,'BUG in definition of KS_BSt%nband')

 if (Ep%gwcomp==1.and.Ep%spmeth>0) then
   msg = "Hilbert transform and extrapolar method are not compatible"
   MSG_ERROR(msg)
 end if

 DBG_EXIT('COLL')

end subroutine setup_screening
!!***
