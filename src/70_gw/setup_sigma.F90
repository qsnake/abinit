!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_sigma
!! NAME
!! setup_sigma
!!
!! FUNCTION
!!  Initialize the data type containing parameters for a sigma calculation.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! Dtfil<type(datafiles_type)>=variables related to files
!! rprim(3,3)=dimensionless real space primitive translations
!! ngfft(18)=information on the (fine) FFT grid used for the density.
!! Psps <Pseudopotential_type)>=Info on pseudopotential, only for consistency check of the KSS file 
!!
!! OUTPUT
!! Sigp<Sigma_parameters>=Parameters governing the self-energy calculation.
!! Kmesh <BZ_mesh_type>=Structure describing the k-point sampling.
!! Qmesh <BZ_mesh_type>=Structure describing the q-point sampling.
!! Cryst<Crystal_structure>=Info on unit cell and symmetries.
!! Gsph_Max<Gvectors_type>=Info on the G-sphere
!! Gsph_c<Gvectors_type>=Info on the G-sphere for W and Sigma_c
!! Hdr_kss<hdr_type>=The header of the KSS file
!! Hdr_out<hdr_type>=The header to be used for the results of sigma calculations.
!! Vcp<vcoul_t>= Datatype gathering information on the coulombian interaction and the cutoff technique.
!! Er<Epsilonm1_results>=Datatype storing data used to construct the screening (partially Initialized in OUTPUT)
!! KS_BSt<Bandstructure_type>=The KS energies and occupation factors.
!! gwc_ngfft(18), gwx_ngfft(18)= FFT meshes for the oscillator strengths used for the correlated and the
!!   exchange part of the self-energy, respectively.
!! comm=MPI communicator.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      destroy_sigijtab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setup_sigma(codvsn,acell,rprim,ngfftf,Dtset,Dtfil,Psps,Pawtab,&
& gwx_ngfft,gwc_ngfft,Hdr_kss,Hdr_out,Cryst,Kmesh,Qmesh,KS_BSt,Gsph_Max,Gsph_c,Vcp,Er,Sigp,comm)

 use m_profiling
    
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_errors

 use m_gwdefs,        only : GW_Q0_DEFAULT, SIG_GW_AC, nullify_sigma_parameters, sigma_parameters, sigma_is_herm, sigma_needs_w,&
&                            gw_uses_wfk_file
 use m_fstrings,      only : basename
 use m_header,        only : hdr_init
 use m_crystal,       only : print_crystal, idx_spatial_inversion, crystal_structure
 use m_crystal_io,    only : init_crystal_from_hdr
 use m_bz_mesh,       only : bz_mesh_type, init_kmesh, has_BZ_item, isamek, get_ng0sh, print_BZ_mesh,&
&                            get_bz_item, has_IBZ_item, find_qmesh
 use m_ebands,        only : bstruct_init, enclose_degbands
 use m_vcoul,         only : vcoul_t, vcoul_init
 use m_fft_mesh,      only : setmesh
 use m_gsphere,       only : gvectors_type, init_gsphere, merge_and_sort_kg
 use m_screening,     only : nullify_epsilonm1_results, init_er_from_file, epsilonm1_results
 use m_io_kss,        only : testkss
 use m_io_screening,  only : print_scrhdr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_sigma'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_70_gw, except_this_one => setup_sigma
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=6),intent(in) :: codvsn
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(inout) :: Dtset
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
 type(Sigma_parameters),intent(out) :: Sigp
 type(Epsilonm1_results),intent(out) :: Er
 type(Bandstructure_type),intent(out) :: KS_BSt
 type(BZ_mesh_type),intent(out) :: Kmesh,Qmesh
 type(Crystal_structure),intent(out) :: Cryst
 type(Gvectors_type),intent(out) :: Gsph_Max,Gsph_c
 type(Hdr_type),intent(out) :: Hdr_kss,Hdr_out
 type(vcoul_t),intent(out) :: Vcp
!arrays
 integer,intent(in) :: ngfftf(18)
 integer,intent(out) :: gwc_ngfft(18),gwx_ngfft(18)
 real(dp),intent(in) :: acell(3),rprim(3,3)

!Local variables-------------------------------
!scalars
 integer :: bantot,enforce_sym,ib,ibtot,ii,ikcalc,ikibz,io,isppol,itypat,jj,method
 integer :: mg0sh,mod10,mqmem,mpsang_kss,nbnds_kss,ikcalc2bz,ierr,ng
 integer :: gwc_nfftot,gwx_nfftot,ng_kss,nsym_kss,nqlwl !,test_ng
 integer :: pertcase_,timrev !restart,restartpaw,
 integer :: iqbz,isym,iq_ibz,itim,ic,pinv,ig1,ig2,ng_sigx,spin !,ik_ibz,band
 real(dp),parameter :: OMEGASIMIN=0.01d0
 real(dp) :: domegas,domegasi,ucvol,tol_enedif
 logical,parameter :: linear_imag_mesh=.TRUE.
 logical :: ltest,only_one_kpt,remove_inv,found,ltmp,changed,silent
 character(len=500) :: msg                   
 character(len=fnlen) :: fname,fcore,string
 character(len=fnlen) :: wfk_fname
 type(MPI_type) :: MPI_enreg_seq
 type(wvl_internal_type) :: wvl
!arrays
 integer :: ng0sh_opt(3),G0(3),q_umklp(3),gswp(3)
 integer,allocatable :: npwarr(:),nlmn(:)
 integer,pointer :: gvec_p(:,:),gsphere_sigx_p(:,:),shlim(:)
 !integer,pointer :: test_gvec_p(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),sq(3),q_bz(3)
 real(dp) :: gamma_point(3,1)
 real(dp),pointer :: energies_p(:,:,:)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:),qlwl(:,:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:)
 
! *************************************************************************
 
 DBG_ENTER('COLL')

 ! === Check for calculations that are not implemented ===
 ltest=ALL(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol)==Dtset%nband(1))
 ABI_CHECK(ltest,'Dtset%nband must be constant')

 !* Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)

 ! === Nullify pointers in the output structures ===
 call nullify_sigma_parameters(Sigp)
 !
 ! === Basic parameters ===
 Sigp%ppmodel    = Dtset%ppmodel
 Sigp%gwcalctyp  = Dtset%gwcalctyp
 Sigp%nbnds      = Dtset%nband(1) 
 Sigp%symsigma   = Dtset%symsigma
 Sigp%zcut       = Dtset%zcut
 Sigp%soenergy   = Dtset%soenergy
 timrev=  2 ! This information is not reported in the header
            ! 1 => do not use time-reversal symmetry 
            ! 2 => take advantage of time-reversal symmetry
 !
 ! === For HF, SEX or COHSEX use Hybertsen-Louie PPM (only $\omega=0$) ===
 ! * Use fake screening for HF.
 ! FIXME Why, we should not redefine Sigp%ppmodel
 mod10=MOD(Sigp%gwcalctyp,10)
 if (mod10==5.or.mod10==6.or.mod10==7) Sigp%ppmodel=2
 if (mod10<5.and.MOD(Sigp%gwcalctyp,1)/=1) then ! * One shot GW (PPM or contour deformation).
   if (Dtset%nomegasrd==1) then ! avoid division by zero!
     Sigp%nomegasrd  =1 
     Sigp%maxomega4sd=zero 
     Sigp%deltae     =zero
   else
     Sigp%nomegasrd   = Dtset%nomegasrd 
     Sigp%maxomega4sd = Dtset%omegasrdmax
     Sigp%deltae     = (2*Sigp%maxomega4sd)/(Sigp%nomegasrd-1)
   endif
 else ! * For AC no need to evaluate derivative by finite differences.
   Sigp%nomegasrd  =1 
   Sigp%maxomega4sd=zero 
   Sigp%deltae     =zero
 end if
 !
 !=== For analytic continuation define the number of imaginary frequencies for Sigma ===
 ! * Tests show than more than 12 freqs in the Pade approximant worsen the results!
 Sigp%nomegasi=0

 if (mod10==1) then 
   Sigp%nomegasi  =Dtset%nomegasi
   Sigp%omegasimax=Dtset%omegasimax 
   Sigp%omegasimin=OMEGASIMIN 
   write(msg,'(4a,i3,2(2a,f8.3),a)')ch10,&
&    ' Parameters for analytic continuation : ',ch10,&
&    '  number of imaginary frequencies for sigma =  ',Sigp%nomegasi,ch10,&
&    '  min frequency for sigma on imag axis [eV] =  ',Sigp%omegasimin*Ha_eV,ch10,&
&    '  max frequency for sigma on imag axis [eV] =  ',Sigp%omegasimax*Ha_eV,ch10
   call wrtout(std_out,msg,'COLL')

   !TODO this should not be done here but in init_sigma_results
   ABI_ALLOCATE(Sigp%omegasi,(Sigp%nomegasi))

   if (linear_imag_mesh) then  ! * Linear mesh along the imaginary axis.
     domegasi=Sigp%omegasimax/(Sigp%nomegasi-1)
     do io=1,Sigp%nomegasi
       Sigp%omegasi(io)=CMPLX(zero,(io-1)*domegasi)
     end do
   else ! * Logarithmic mesh along the imaginary axis.
     MSG_ERROR("AC + log mesh not implemented")
     !domegasi=(Sigp%omegasimax/Sigp%omegasimin)**(one/(Sigp%nomegasi-1))
     !Sigp%omegasi(1)=czero; ldi=domegasi
     !do io=2,Sigp%nomegasi
     ! omega(io)=CMPLX(zero,ldi*Sigp%omegasimin)
     ! Sigp%omegasi(io)=ldi*domegasi
     !end do
   end if
   
   write(msg,'(4a)')ch10,&
&    ' setup_sigma : calculating Sigma(iw)',&
&    ' at imaginary frequencies [eV] (Fermi Level set to 0) ',ch10
   call wrtout(std_out,msg,'COLL') 
   call wrtout(ab_out,msg,'COLL')
   do io=1,Sigp%nomegasi
     write(msg,'(2(f10.3,2x))')Sigp%omegasi(io)*Ha_eV
     call wrtout(std_out,msg,'COLL') 
     call wrtout(ab_out,msg,'COLL')
   end do

   ltest=(Sigp%omegasimax>0.1d-4.and.Sigp%nomegasi>0)
   ABI_CHECK(ltest,'Wrong value of omegasimax or nomegasi')

   if (Sigp%gwcalctyp/=1) then ! only one shot GW is allowed for AC.
     MSG_ERROR("SC-GW with analytic continuation is not coded")
   end if 
 end if 

 if (Sigp%symsigma/=0.and.Sigp%gwcalctyp>=20) then
   msg = "SC-GW with symmetries is still under development. Use at your own risk!"
   MSG_WARNING(msg)
 end if
 !
 ! === Setup parameters for Spectral function ===

 if (Dtset%gw_custom_freqsp/=0) then
   Sigp%nomegasr = Dtset%gw_custom_freqsp
   write(msg,'(a)')' Custom grid for spectral function specified. Assuming experienced user.'
   MSG_WARNING(msg)
   if (Dtset%gw_custom_freqsp/=0) then
     Dtset%nfreqsp = Dtset%gw_custom_freqsp
     write(msg,'(a)')' nfreqsp has been set to the same number as gw_custom_freqsp'
     MSG_WARNING(msg)
   end if
 else
   Sigp%nomegasr = Dtset%nfreqsp
   Sigp%minomega_r=Dtset%freqspmin
   Sigp%maxomega_r=Dtset%freqspmax
 end if

 if (Sigp%nomegasr>0) then
   if (Dtset%gw_custom_freqsp==0) then
     ! Check
     if (Sigp%minomega_r.GE.Sigp%maxomega_r) then
       MSG_ERROR(' freqspmin must be smaller than freqspmax!')
     end if
     domegas=(Sigp%maxomega_r-Sigp%minomega_r)/(Sigp%nomegasr-1)
     !TODO this should be moved to Sr% and done in init_sigma_results
     ABI_ALLOCATE(Sigp%omega_r,(Sigp%nomegasr))
     do io=1,Sigp%nomegasr
       Sigp%omega_r(io) = CMPLX(Sigp%minomega_r + domegas*(io-1),zero)
     end do
     write(msg,'(4a,i8,3(2a,f8.3),a)')ch10,&
&      ' Parameters for the calculation of the spectral function : ',ch10,&
&      '  Number of points    = ',Sigp%nomegasr,ch10,&
&      '  Min frequency  [eV] = ',Sigp%minomega_r*Ha_eV,ch10,&
&      '  Max frequency  [eV] = ',Sigp%maxomega_r*Ha_eV,ch10,&
&      '  Frequency step [eV] = ',domegas*Ha_eV,ch10
     call wrtout(std_out,msg,'COLL')
   else
     Sigp%minomega_r = MINVAL(Dtset%gw_freqsp(:))
     Sigp%maxomega_r = MAXVAL(Dtset%gw_freqsp(:))
     !TODO this should be moved to Sr% and done in init_sigma_results
     ABI_ALLOCATE(Sigp%omega_r,(Sigp%nomegasr))
     do io=1,Sigp%nomegasr
       Sigp%omega_r(io) = CMPLX(Dtset%gw_freqsp(io),zero)
     end do
     write(msg,'(4a,i8,2(2a,f8.3),3a)')ch10,&
&      ' Parameters for the calculation of the spectral function : ',ch10,&
&      '  Number of points    = ',Sigp%nomegasr,ch10,&
&      '  Min frequency  [eV] = ',Sigp%minomega_r*Ha_eV,ch10,&
&      '  Max frequency  [eV] = ',Sigp%maxomega_r*Ha_eV,ch10,&
&      '  A custom set of frequencies is used! See the input file for values.',ch10
     call wrtout(std_out,msg,'COLL')
   end if
 else
   !In indefo all these quantities are set to zero
   !Sigp%nomegasr=1 
   !allocate(Sigp%omega_r(Sigp%nomegasr))
   !Sigp%omega_r(1)=0
 end if
 !
 ! === Dimensional primitive translations rprimd (from input), gprimd, metrics and unit cell volume ===
 call mkrdim(acell,rprim,rprimd)  
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol) 
 !
 Sigp%npwwfn=Dtset%npwwfn 
 Sigp%npwx  =Dtset%npwsigx 
 !
 ! === Read parameters of the KSS, verifify them and retrieve all G-vectors ===


#if 0
!BEGIN GWWFK
 ! this piece of code will replace the setup of the G-sphere for Espilon
 pinv=1 !; if (timrev==2) pinv=-1 ! if time-reversal can be used, 1 otherwise
 nullify(test_gvec_p,shlim)
 call merge_and_sort_kg(Dtset%nkpt,Dtset%kptns,Dtset%ecut,Dtset%nsym,&
&    pinv,Dtset%symrel,gprimd,test_gvec_p,Dtset%prtvol,shlim_p=shlim)
 test_ng=SIZE(test_gvec_p,DIM=2)
 ABI_DEALLOCATE(shlim)
!END GWWFK
#endif

 if (gw_uses_wfk_file) then
 !if (.FALSE. .and. gw_uses_wfk_file) then  !Temporary hack so that the G sphere is always read from the KSS file.
   wfk_fname = Dtfil%fnameabi_kss; ii=LEN_TRIM(wfk_fname)
   wfk_fname(ii-2:ii) = "WFK"
   MSG_WARNING("Using WFK file"//TRIM(wfk_fname))

   nbnds_kss=-1
   call wfk_read_ene(wfk_fname,Dtset%accesswff,nbnds_kss,energies_p,Hdr_kss,Dtset%prtvol,comm) 

   nsym_kss = Hdr_kss%nsym
   !nbnds_kss=MINVAL(Hdr_kss%nband)
   mpsang_kss=Psps%mpsang

   pinv=1 !; if (timrev==2) pinv=-1 ! if time-reversal can be used, 1 otherwise
   nullify(gvec_p,shlim)
   call merge_and_sort_kg(Hdr_kss%nkpt,Hdr_kss%kptns,Hdr_kss%ecut_eff,Hdr_kss%nsym,&
&    pinv,Hdr_kss%symrel,gprimd,gvec_p,Dtset%prtvol,shlim_p=shlim)
   ng_kss=SIZE(gvec_p,DIM=2)
   
   ABI_DEALLOCATE(shlim)
 else
   call testkss(Dtfil%fnameabi_kss,Dtset%accesswff,nsym_kss,nbnds_kss,&
&    ng_kss,mpsang_kss,gvec_p,energies_p,Hdr_kss,comm) 
 end if

!BEGIN GWWFK
#if 0
 jj = ng_kss
 if (ng_kss/=test_ng) then
   jj = MIN(ng_kss,test_ng)
   !MSG_WARNING("ng_kss/=test_ng")
   !write(std_out,*)ng_kss,test_ng
 end if

 if (ANY(test_gvec_p(:,1:jj)/=gvec_p(:,1:jj))) then
   !MSG_WARNING("ANY(test_gvec/=gvec_p)")
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
 !ABI_CHECK(ltest,'Psps%mpsang/=mpsang_kss') Removed because it triggers a bug on fock_xlf

 ! === Get important dimensions from the KSS header ===
 ! * Check also the consistency btw Hdr_kss and Dtset.
 Sigp%nsppol =Hdr_kss%nsppol
 Sigp%nspinor=Hdr_kss%nspinor 
 Sigp%nsig_ab=Hdr_kss%nspinor**2  !FIXME Is it useful calculating only diagonal terms?
 call hdr_vs_dtset(Hdr_kss,Dtset) 
                                                                                     
 ! === Check input ===
 if (Sigp%ppmodel==3.or.Sigp%ppmodel==4) then 
   if (Sigp%gwcalctyp>=10) then 
     write(msg,'(a,i3,a)')' The ppmodel chosen and gwcalctyp ',Dtset%gwcalctyp,' are not compatible. '
     MSG_ERROR(msg)
   end if
   if (Sigp%nspinor==2) then 
     write(msg,'(a,i3,a)')' The ppmodel chosen and nspinor ',Sigp%nspinor,' are not compatible. '
     MSG_ERROR(msg)
   end if
 end if 

 ! === Create crystal_structure data type ===
 remove_inv=(nsym_kss/=Hdr_kss%nsym) 

 call init_crystal_from_hdr(Cryst,Hdr_kss,timrev,remove_inv)
 call print_crystal(Cryst)

 !Sigp%npwvec=MAX(Sigp%npwwfn,Sigp%npwx)

 if (Sigp%npwwfn>ng_kss) then ! cannot use more G"s for the wfs than those stored on file
   Sigp%npwwfn    =ng_kss 
   Dtset%npwwfn =ng_kss
   write(msg,'(2a,(a,i8,a))')&
&    ' Number of G-vectors for WFS found in the KSS file is less than required',ch10,&
&    '  calculation will proceed with npwwfn  = ',Sigp%npwwfn,ch10
   MSG_WARNING(msg)
 end if

 if (Sigp%npwx>ng_kss) then ! Have to recalcuate the (large) sphere for Sigma_x.

  pinv=1; if (remove_inv.and.Cryst%timrev==2) pinv=-1
  gamma_point(:,1) = (/zero,zero,zero/); nullify(gsphere_sigx_p)

  call merge_and_sort_kg(1,gamma_point,Dtset%ecutsigx,Cryst%nsym,pinv,Cryst%symrel,&
&  Cryst%gprimd,gsphere_sigx_p,Dtset%prtvol)

!  call merge_and_sort_kg(Hdr_kss%nkpt,Hdr_kss%kptns,Dtset%ecutsigx,Cryst%nsym,pinv,Cryst%symrel,&
!&  Cryst%gprimd,gsphere_sigx_p,Dtset%prtvol) !this call overestimates npwsigx due the the k-list

  ng_sigx=SIZE(gsphere_sigx_p,DIM=2)
  Sigp%npwx     = ng_sigx
  Dtset%npwsigx = ng_sigx

  write(msg,'(2a,(a,i8,a))')&
&   ' Number of G-vectors for Sigma_x found in the KSS file is less than required',ch10,&
&   '  calculation will proceed with npwsigx = ',Sigp%npwx,ch10
  MSG_WARNING(msg)

  ltest = (Sigp%npwx >= ng_kss)
  ABI_CHECK(ltest,"Sigp%npwx<ng_kss!")
  !
  ! gvec_p is contained in gsphere_sigx_p but the ordering might differ as we have used 
  ! a different cutoff value in merge_and_sort_kg. Thus swap the elements in gsphere_sigx_p 
  ! in order to have gvec_p(:,1:ng_kss) == gsphere_sigx_p(:,1:ng_kss)

  do ig1=1,ng_kss

   if (ANY(gvec_p(:,ig1)/=gsphere_sigx_p(:,ig1))) then
    gswp = gsphere_sigx_p(:,ig1)

    ig2=ig1; found=.FALSE.
    do while (ig2<ng_sigx.and..not.found)
     ig2=ig2+1
     found = ALL(gvec_p(:,ig1)==gsphere_sigx_p(:,ig2))
    end do

    if (.not.found) MSG_ERROR("gvec_p(:,ig1) not in gsphere_sigx_p")
    gsphere_sigx_p(:,ig1) = gvec_p(:,ig1)
    gsphere_sigx_p(:,ig2) = gswp
   end if

  end do

  !ltest = ALL(gvec_p(:,1:ng_kss)==gsphere_sigx_p(:,1:ng_kss))
  !ABI_CHECK(ltest,"")

  ! * Fill gvec_p with larger sphere.
  ABI_DEALLOCATE(gvec_p)
  ABI_ALLOCATE(gvec_p,(3,Sigp%npwx))
  gvec_p = gsphere_sigx_p
  ABI_DEALLOCATE(gsphere_sigx_p)

 end if

 Sigp%npwvec=MAX(Sigp%npwwfn,Sigp%npwx)

 if (Sigp%nbnds>nbnds_kss) then
   Sigp%nbnds    =nbnds_kss 
   Dtset%nband(:)=nbnds_kss
   Dtset%mband   =MAXVAL(Dtset%nband)
   write(msg,'(3a,i4,a)')&
&    ' Number of bands found less then required',ch10,&
&    ' calculation will proceed with nbnds = ',nbnds_kss,ch10
   MSG_WARNING(msg)
 end if

 !==== Set up of the k-points and tables in the whole BZ ===
 call init_kmesh(Kmesh,Cryst,Hdr_kss%nkpt,Hdr_kss%kptns,Dtset%kptopt,wrap_1zone=.FALSE.)
 !call init_kmesh(Kmesh,Cryst,Hdr_kss%nkpt,Hdr_kss%kptns,Dtset%kptopt,wrap_1zone=.TRUE.)

 call print_BZ_mesh(Kmesh,"K-mesh for the wavefunctions",std_out,Dtset%prtvol,"COLL")
 call print_BZ_mesh(Kmesh,"K-mesh for the wavefunctions",ab_out, 0,           "COLL")

 ! === Initialize the band structure datatype ===
 ! * Copy KSS energies and occupations up to Sigp%nbnds==Dtset%nband(:)
 ! TODO Recheck symmorphy and inversion

 bantot=SUM(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol))
 ABI_ALLOCATE(doccde,(bantot))
 ABI_ALLOCATE(eigen,(bantot))
 ABI_ALLOCATE(occfact,(bantot))
 doccde(:)=zero; eigen(:)=zero; occfact(:)=zero 

 jj=0; ibtot=0
 do isppol=1,Dtset%nsppol
   do ikibz=1,Dtset%nkpt
     do ib=1,Hdr_kss%nband(ikibz+(isppol-1)*Dtset%nkpt)
       ibtot=ibtot+1
       if (ib<=Sigp%nbnds) then 
         jj=jj+1
         occfact(jj)=Hdr_kss%occ(ibtot)
         eigen  (jj)=energies_p(ib,ikibz,isppol)
       end if
     end do
   end do
 end do
 ABI_DEALLOCATE(energies_p)
 !
 ! * Make sure that Dtset%wtk==Kmesh%wt due to the dirty treatment of 
 !   symmetry operations in the old GW code (symmorphy and inversion) 
 ltest=(ALL(ABS(Dtset%wtk(1:Kmesh%nibz)-Kmesh%wt(1:Kmesh%nibz))<tol6))
 ABI_CHECK(ltest,'Mismatch between Dtset%wtk and Kmesh%wt')

 ABI_ALLOCATE(npwarr,(Dtset%nkpt))
 npwarr(:)=Sigp%npwwfn

 call bstruct_init(bantot,KS_BSt,Dtset%nelect,doccde,eigen,Dtset%istwfk,Kmesh%ibz,Dtset%nband,&
&  Kmesh%nibz,npwarr,Dtset%nsppol,Dtset%nspinor,Dtset%tphysel,Dtset%tsmear,Dtset%occopt,occfact,Kmesh%wt) 

 !TODO call update_occ here
 !$call update_occ(KS_BSt,fixmom,stmbias,Dtset%prtvol)

 ! this fails simply because in case of NSCF occ are zero
 ! TODO outkss should calculate occ factors in case of NSCF run. 
 !ltest=(ALL(ABS(occfact-KS_BSt%occ)<tol6)) 
 !call assert(ltest,'difference in occfact')
 !write(std_out,*)MAXVAL(ABS(occfact(:)-KS_BSt%occ(:))) 

 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(eigen)
 ABI_DEALLOCATE(npwarr)

 ! === Create Sigma header === 
 ! TODO Fix problems with symmorphy and k-points
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
 !call hdr_check(1002,1002,Hdr_out,Hdr_kss,'COLL',restart,restartpaw)

 ABI_DEALLOCATE(occfact)
 call rhoij_free(Pawrhoij)
 ABI_DEALLOCATE(Pawrhoij)
 !
 ! ===========================================================
 ! ==== Setup of k-points and bands for the GW corrections ====
 ! ===========================================================
 ! * maxbdgw and minbdgw are the Max and min band index for GW corrections over k-points. 
 !   They are used to dimension wfr_gw and calculate the matrix elements.
 !
 if (Dtset%nkptgw==0) then
   ! * If not precised, calculate all k-points in the IBZ and all bands.
   !   This convention is particularly useful for self-consistent HF.
   Dtset%nkptgw=Kmesh%nibz
   Sigp%nkptgw =Dtset%nkptgw
   ABI_ALLOCATE(Sigp%kptgw,(3,Sigp%nkptgw))
   ABI_ALLOCATE(Sigp%kptgw2bz,(Sigp%nkptgw))
   ABI_ALLOCATE(Sigp%minbnd,(Sigp%nkptgw,Sigp%nsppol))
   ABI_ALLOCATE(Sigp%maxbnd,(Sigp%nkptgw,Sigp%nsppol))
   Sigp%kptgw(:,:)=Kmesh%ibz(:,:) 
   Sigp%minbnd=1;          Sigp%minbdgw=MINVAL(Sigp%minbnd) 
   Sigp%maxbnd=Sigp%nbnds; Sigp%maxbdgw=MAXVAL(Sigp%maxbnd)

 else ! * Treat only the k-points and bands specified in the input file.
   Sigp%nkptgw=Dtset%nkptgw
   ABI_ALLOCATE(Sigp%kptgw,(3,Sigp%nkptgw))
   ABI_ALLOCATE(Sigp%kptgw2bz,(Sigp%nkptgw))
   ABI_ALLOCATE(Sigp%minbnd,(Sigp%nkptgw,Sigp%nsppol))
   ABI_ALLOCATE(Sigp%maxbnd,(Sigp%nkptgw,Sigp%nsppol))
   do ii=1,3
     do ikcalc=1,Sigp%nkptgw
       Sigp%kptgw(ii,ikcalc)=Dtset%kptgw(ii,ikcalc)
     end do
   end do
   
   do spin=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw
       if (Dtset%bdgw(2,ikcalc,spin)>Sigp%nbnds) then
         write(std_out,'(a,2i0,2(a,i0),2a,i0)')&
&          " For (k,s) ",ikcalc,spin," bdgw= ",Dtset%bdgw(2,ikcalc,spin), " > nbnds=",Sigp%nbnds,ch10,&
&          " Calculation will continue with bdgw =",Sigp%nbnds
         Dtset%bdgw(2,ikcalc,spin)=Sigp%nbnds
         MSG_COMMENT(msg)
       end if
     end do
   end do
   !
   ! Make sure that all the degenerate states are included.
   ! * We will have to average the GW corrections over degenerate states if symsigma=1 is used.
   ! * KS states belonging to the same irreducible representation should be included in the basis set used for SCGW.
   if (Sigp%symsigma/=0 .or. Sigp%gwcalctyp>=10) then
     tol_enedif = 0.001/Ha_eV
     do isppol=1,Sigp%nsppol
       do ikcalc=1,Sigp%nkptgw

         if (has_IBZ_item(Kmesh,Sigp%kptgw(:,ikcalc),ikibz,G0)) then
         !if (has_BZ_item(Kmesh,Sigp%kptgw(:,ikcalc),ikcalc2bz,G0)
           call enclose_degbands(KS_BSt,ikibz,isppol,Dtset%bdgw(1,ikcalc,isppol),Dtset%bdgw(2,ikcalc,isppol),changed,tol_enedif)
           if (changed) then 
             write(msg,'(2(a,i0),2a,2(1x,i0))')&
&              " Not all the degenerate states at ikcalc= ",ikcalc,", spin= ",isppol,ch10,&
&              " were included in the bdgw set. bdgw has been changed to: ",Dtset%bdgw(:,ikcalc,isppol) 
             MSG_COMMENT(msg)
           end if
         else
           write(msg,'(a,3(f6.3,1x),a)')' k-point ',Sigp%kptgw(:,ikcalc),' not in the IBZ'
           MSG_ERROR(msg)
         end if

       end do
     end do
   end if

   do spin=1,Sigp%nsppol 
     Sigp%minbnd(:,spin)=Dtset%bdgw(1,:,spin) 
     Sigp%maxbnd(:,spin)=Dtset%bdgw(2,:,spin) 
   end do

   Sigp%minbdgw=MINVAL(Sigp%minbnd) 
   Sigp%maxbdgw=MAXVAL(Sigp%maxbnd) 
 end if
 !
 !=== Check if the k-points are in the BZ ===
 !FB TODO Honestly the code is not able to treat k-points, which are not in the IBZ.
 !This extension should require to change the code in different places.
 !Therefore, one should by now prevent the user from calculating sigma for a k-point not in the IBZ.
 !
 do ikcalc=1,Sigp%nkptgw
   if (has_BZ_item(Kmesh,Sigp%kptgw(:,ikcalc),ikcalc2bz,G0)) then
     !found = has_IBZ_item(Kmesh,Sigp%kptgw(:,ikcalc),ikcalc2bz,G0)
     Sigp%kptgw2bz(ikcalc) = ikcalc2bz
   else 
     write(msg,'(a,3(f6.3,1x),a)')' k-point ',Sigp%kptgw(:,ikcalc),' not in the set of kbz'
     MSG_ERROR(msg)
   end if
 end do
 !
 ! Check if there are duplicated k-point in Sigp%
 do ii=1,Sigp%nkptgw
   do jj=ii+1,Sigp%nkptgw
     if (isamek(Sigp%kptgw(:,ii),Sigp%kptgw(:,jj),G0)) then
       write(msg,'(5a)')&
&        ' kptgw contains duplicated k-points. This is not allowed since ',ch10,&
&        ' the QP corrections for this k-point will be calculated more than once. ',ch10,& 
&        ' Check your input file. '
       MSG_ERROR(msg)
     end if
   end do
 end do
 !
 ! Warn the user if SCGW run and not all the k-points are included.
 if (Sigp%gwcalctyp>=10 .and. Sigp%nkptgw/=Hdr_kss%nkpt) then
   write(msg,'(3a,2(a,i0),2a)')ch10,&
&    " COMMENT: In a self-consistent GW run, the QP corrections should be calculated for all the k-points of the KSS file ",ch10,&
&    " but nkptgw= ",Sigp%nkptgw," and KSS nkpt= ",Hdr_kss%nkpt,ch10,&
&    " Assuming expert user. Execution will continue. "
   MSG_COMMENT(msg)
   call wrtout(ab_out,msg,"COLL")
 end if
 !
 ! Setup of the table used in the case of SCGW on wavefunctions to reduce the number 
 ! of elements <i,kgw,s|\Sigma|j,kgw,s> that have to be calculated. No use of symmetries, except for Hermiticity.
 !
 call sigma_tables(Sigp,Kmesh)

 ! === Read external file and initialize basic dimension of Er% ===
 ! TODO use mqmem as input variable instead of gwmem

 ! === If required, use a matrix for $\Sigma_c$ which is smaller than that stored on file ===
 ! * By default the entire matrix is read and used,
 ! * Define consistently npweps, nsheps, and ecuteps for \Sigma_c according the input
 if (Dtset%npweps>0.or.Dtset%ecuteps>0.or.Dtset%nsheps>0) then
   if      (Dtset%npweps>0) then; Dtset%ecuteps=zero; Dtset%nsheps=0;
   else if (Dtset%nsheps>0) then; Dtset%ecuteps=zero; Dtset%npweps=0;
   else                         ; Dtset%npweps=0;     Dtset%nsheps=0; end if
   call setshells(Dtset%ecuteps,Dtset%npweps,Dtset%nsheps,Dtset%nsym,gmet,gprimd,Dtset%symrel,'eps',ucvol)
 end if

 mqmem=0; if (Dtset%gwmem/10==1) mqmem=1

 if (Dtset%getscr/=0.or.Dtset%irdscr/=0) then
   fname=Dtfil%fnameabi_scr
 else if (Dtset%getsuscep/=0.or.Dtset%irdsuscep/=0) then
   fname=Dtfil%fnameabi_sus
 else 
   fname=Dtfil%fnameabi_scr
   !FIXME this has to be cleaned, in tgw2_3 Dtset%get* and Dtset%ird* are  not defined
   !ABI_DIE("getsuscep or irdsuscep are not defined")
 end if
 !
 ! === Setup of q-mesh in the whole BZ ===
 ! * Stop if a nonzero umklapp is needed to reconstruct the BZ. In this case, indeed, 
 !   epsilon^-1(Sq) should be symmetrized in csigme using a different expression (G-G_o is needed)
 !
 call nullify_epsilonm1_results(Er)

 if (sigma_needs_w(Sigp)) then
   call init_Er_from_file(Er,fname,mqmem,Dtset%npweps,Dtset%accesswff,comm)

   Sigp%npwc=Er%npwe          
   if (Sigp%npwc>Sigp%npwx) then 
     Sigp%npwc=Sigp%npwx
     MSG_COMMENT("Found npw_correlation > npw_exchange, Imposing npwc=npwx")
     ! There is a good reason for doing so, see csigme.F90 and the size of the arrays 
     ! rhotwgp and rhotwgp: we need to define a max size and we opt for Sigp%npwx.
   end if 
   Er%npwe=Sigp%npwc
   Dtset%npweps=Er%npwe

   call init_kmesh(Qmesh,Cryst,Er%nqibz,Er%qibz,Dtset%kptopt)
 else
   Er%npwe     =1 
   Sigp%npwc   =1
   Dtset%npweps=1
   call find_qmesh(Qmesh,Cryst,Kmesh)
   ABI_ALLOCATE(Er%gvec,(3,1))
   Er%gvec(:,1) = (/0,0,0/)
 end if

 call print_BZ_mesh(Qmesh,"Q-mesh for screening function",std_out,Dtset%prtvol,"COLL")
 call print_BZ_mesh(Qmesh,"Q-mesh for screening function",ab_out ,0           ,"COLL")

 do iqbz=1,Qmesh%nbz
   call get_BZ_item(Qmesh,iqbz,q_bz,iq_ibz,isym,itim,umklp=q_umklp)

   if (ANY(q_umklp/=0)) then
     sq = (3-2*itim)*MATMUL(Cryst%symrec(:,:,isym),Qmesh%ibz(:,iq_ibz))
     write(std_out,*) sq,Qmesh%bz(:,iqbz) 
     write(msg,'(a,3f6.3,a,3f6.3,2a,9i3,a,i2,2a)')&
&      ' qpoint ',Qmesh%bz(:,iqbz),' is the symmetric of ',Qmesh%ibz(:,iq_ibz),ch10,&
&      ' through operation ',Cryst%symrec(:,:,isym),' and itim ',itim,ch10,&
&      ' however a non zero umklapp G_o vector is required and this is not yet allowed'
     MSG_ERROR(msg)
   end if
 end do 
 !
 ! === Find optimal value for G-sphere enlargment due to oscillator matrix elements ===
 ! * Here I have to be sure that Qmesh%bz is always inside the BZ, not always true size bz is buggy
 ! * -one is used because we loop over all the possibile differences, unlike screening
 mg0sh=5
 call get_ng0sh(Sigp%nkptgw,Sigp%kptgw,Kmesh%nbz,Kmesh%bz,Qmesh%nbz,Qmesh%bz,gmet,-one,mg0sh,ng0sh_opt)
 Sigp%mG0(:)=ng0sh_opt(:)

 ! === Make biggest G-sphere of Sigp%npwvec vectors ===
 only_one_kpt=(Kmesh%nbz==1)
 call init_gsphere(Gsph_Max,only_one_kpt,Cryst,Sigp%npwvec,gvec=gvec_p)

 ! G-sphere for W and Sigma_c is initialized from the SCR file.
 call init_gsphere(Gsph_c,only_one_kpt,Cryst,Sigp%npwc,gvec=Er%gvec)
 ABI_DEALLOCATE(gvec_p)

!BEGINDEBUG
 ! Make sure that the two G-spheres are equivalent.
 ierr=0
#if 1
 if (sigma_needs_w(Sigp)) then
   ng = MIN(Gsph_c%ng,Er%npwe)
   do ig1=1,ng
     if (ANY(Gsph_c%gvec(:,ig1)/=Er%gvec(:,ig1))) then
       ierr=ierr+1
       write(std_out,*)" Gsph_c, Er ",ig1,"/",ng,Gsph_c%gvec(:,ig1),Er%gvec(:,ig1)
     end if
   end do
   ABI_CHECK(ierr==0,"Mismatch between Gsph_chi0 and Gsph_epsG0")
 end if
#endif
!ENDDEBUG
 !
 ! === Get Fourier components of the Coulombian for all q-points in the IBZ ===
 ! * If required, use a cutoff in the interaction 
 ! * Pcv%vc_sqrt contains Vc^{-1/2}
 ! * Setup also the analytical calculation of the q->0 component
 ! FIXME recheck ngfftf since I got different charge outside the cutoff region

 if (Dtset%gw_nqlwl==0) then
   nqlwl=1 
   ABI_ALLOCATE(qlwl,(3,nqlwl))
   qlwl(:,1)= GW_Q0_DEFAULT
 else 
   nqlwl=Dtset%gw_nqlwl
   ABI_ALLOCATE(qlwl,(3,nqlwl))
   qlwl(:,:)=Dtset%gw_qlwl(:,1:nqlwl)
 end if

 call vcoul_init(Vcp,Gsph_Max,Qmesh,Kmesh,Dtset%rcut,Dtset%icutcoul,Dtset%vcutgeo,Sigp%npwx,nqlwl,qlwl,Cryst%rprimd,ngfftf,comm)

#if 0
 ! Using the random q for the optical limit is one of the reasons
 ! why sigma breaks the initial energy degeneracies.
 Vcp%i_sz=zero
 Vcp%vc_sqrt(1,:)=czero
 Vcp%vcqlwl_sqrt(1,:)=czero
#endif

 ABI_DEALLOCATE(qlwl)

 ! === Setup of the FFT mesh for the oscilator strengths === 
 ! * gwc_ngfft(7:18)==Dtset%ngfft(7:18) which is initialized before entering screening.
 ! * Here we redefine gwc_ngfft(1:6) according to the following options :
 !
 ! method==0 --> FFT grid read from fft.in (debugging purpose)
 ! method==1 --> Normal FFT mesh
 ! method==2 --> Slightly augmented FFT grid to calculate exactly rho_tw_g (see setmesh.F90)
 ! method==3 --> Doubled FFT grid, same as the the FFT for the density,
 !
 ! enforce_sym==1 ==> Enforce a FFT mesh compatible with all the symmetry operation and FFT library
 ! enforce_sym==0 ==> Find the smallest FFT grid compatbile with the library, do not care about symmetries
 !
 gwc_ngfft(1:18)=Dtset%ngfft(1:18)
 gwx_ngfft(1:18)=Dtset%ngfft(1:18)
 method=2
 if (Dtset%fftgw==00 .or. Dtset%fftgw==01) method=0
 if (Dtset%fftgw==10 .or. Dtset%fftgw==11) method=1
 if (Dtset%fftgw==20 .or. Dtset%fftgw==21) method=2
 if (Dtset%fftgw==30 .or. Dtset%fftgw==31) method=3
 enforce_sym=MOD(Dtset%fftgw,10) 

 ! FFT mesh for sigma_x.
 call setmesh(gmet,Gsph_Max%gvec,gwx_ngfft,Sigp%npwvec,Sigp%npwx,Sigp%npwwfn,&
&  gwx_nfftot,method,Sigp%mG0,Cryst,enforce_sym)

 ! FFT mesh for sigma_c.
 call setmesh(gmet,Gsph_Max%gvec,gwc_ngfft,Sigp%npwvec,Sigp%npwc,Sigp%npwwfn,&
&  gwc_nfftot,method,Sigp%mG0,Cryst,enforce_sym,silent)

 ! ======================================================================
 ! ==== Check for presence of files with core orbitals, for PAW only ====
 ! ======================================================================
 Sigp%use_sigxcore=0
 if (Dtset%usepaw==1.and.Dtset%gw_sigxcore==1) then
   ii = 0
   do itypat=1,Cryst%ntypat
     string = Psps%filpsp(itypat)
     fcore = "CORE_"//TRIM(basename(string))
     ic = INDEX (TRIM(string), "/" , back=.TRUE.)
     if (ic>0 .and. ic<LEN_TRIM(string)) then ! string defines a path, prepend path to fcore
       fcore = Psps%filpsp(itypat)(1:ic)//TRIM(fcore)
     end if
     inquire(file=fcore,exist=ltmp)
     if (ltmp) then 
       ii=ii+1
     else 
       msg=" HF decoupling is required but could not find file: "//TRIM(fcore)
       MSG_WARNING(msg)
     end if
   end do

   Sigp%use_sigxcore=1
   if (ii/=Cryst%ntypat) then 
     MSG_ERROR("Files with core orbitals not found")
   end if
 end if ! PAW+HF decoupling
 !
 ! ==============================
 ! ==== Extrapolar technique ====
 ! ==============================
 Sigp%gwcomp   = Dtset%gwcomp
 Sigp%gwencomp = Dtset%gwencomp

 if (Sigp%gwcomp==1) then
   write(msg,'(6a,e11.4,a)')ch10,&
&    ' Using the extrapolar approximation to accelerate convergence',ch10,&
&    ' with respect to the number of bands included',ch10,&
&    ' with gwencomp: ',Sigp%gwencomp*Ha_eV,' [eV]'
   call wrtout(std_out,msg,'COLL')
 end if
 !
 ! ===================================
 ! ==== Final compatibility tests ====
 ! ===================================
 if (ANY(KS_BSt%istwfk/=1)) then
   MSG_WARNING('istwfk/=1 is still under development')
 end if

 ltest=(KS_BSt%mband==Sigp%nbnds.and.ALL(KS_BSt%nband==Sigp%nbnds))
 ABI_CHECK(ltest,'BUG in definition of KS_BSt%nband')

 ! FIXME
 if (Dtset%symsigma/=0 .and. Sigp%nomegasr/=0) then
   if (idx_spatial_inversion(Cryst) == 0) then 
     write(msg,'(5a)')' setup_sigma : BUG :',ch10,&
&      ' It is not possible to use symsigma/=0 to calculate the spectral function ',ch10,&
&      ' when the system does not have the spatial inversion. Please use symsigma=0 '
     MSG_WARNING(msg)
   end if
 end if

 if (mod10==SIG_GW_AC) then
   if (Sigp%gwcalctyp/=1) &
      MSG_ERROR("Self-consistency with AC not implemented")
   if (Sigp%gwcomp==1) &
      MSG_ERROR("AC with extrapolar technique not implemented")
 end if

 if (Sigp%nspinor==2) then
   ABI_CHECK(Sigp%symsigma==0,'symsigma=1 and nspinor=2 not implemented')
 end if

 DBG_EXIT('COLL')

end subroutine setup_sigma
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/sigma_tables
!! NAME
!! sigma_tables
!!
!! FUNCTION
!!
!! INPUTS
!! Kmesh <BZ_mesh_type>=Structure describing the k-point sampling.
!! [Bnd_sym(Kmesh%nibz,Sigp%nsppol)] <type(Bands_Symmetries)>
!!
!! SiDE EFFECTS
!! Sigp<Sigma_parameters>=This routine initializes the tables:
!!   %Sigcij_tab
!!   %Sigxij_tab
!!  that are used to select the matrix elements of the self-energy that have to be calculated.
!!
!! PARENTS
!!      setup_sigma,sigma
!!
!! CHILDREN
!!      destroy_sigijtab
!!
!! SOURCE


subroutine sigma_tables(Sigp,Kmesh,Bnd_sym)

 use m_profiling
    
 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use m_errors

 use m_bz_mesh,     only : bz_mesh_type
 use m_bands_sym,   only : bands_symmetries, bsym_failed

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigma_tables'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Sigma_parameters),intent(inout) :: Sigp
 type(BZ_mesh_type),intent(in) :: Kmesh
!arrays
 type(Bands_Symmetries),optional,intent(in) :: Bnd_sym(Kmesh%nibz,Sigp%nsppol)

!Local variables-------------------------------
!scalars
 integer :: spin,ikcalc,ik_ibz,bmin,bmax,bcol,brow
 integer :: ii,idx_x,idx_c,irr_idx1,irr_idx2
!arrays
 integer,allocatable :: sigc_bidx(:),sigx_bidx(:)
 logical :: use_sym_at(Kmesh%nibz,Sigp%nsppol)

! *************************************************************************

 ! Recreate the Sig_ij tables taking advantage of the classification of the bands.
 if (associated(Sigp%Sigxij_tab)) call destroy_sigijtab(Sigp%Sigxij_tab)
 if (associated(Sigp%Sigcij_tab)) call destroy_sigijtab(Sigp%Sigcij_tab)

 ABI_ALLOCATE(Sigp%Sigcij_tab,(Sigp%nkptgw,Sigp%nsppol))
 ABI_ALLOCATE(Sigp%Sigxij_tab,(Sigp%nkptgw,Sigp%nsppol))

 use_sym_at=.FALSE.
 if (PRESENT(Bnd_sym)) then
   do spin=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw
      ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
      use_sym_at(ik_ibz,spin) = ( .not.bsym_failed(Bnd_sym(ik_ibz,spin)) )
     end do
   end do
 end if

 do spin=1,Sigp%nsppol
   do ikcalc=1,Sigp%nkptgw
     ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))

     if (use_sym_at(ik_ibz,spin)) then
       if (Sigp%gwcalctyp<20) then  
         MSG_ERROR("You should not be here!")
       end if

       bmin=Sigp%minbnd(ikcalc,spin); bmax=Sigp%maxbnd(ikcalc,spin)
       ABI_ALLOCATE(Sigp%Sigxij_tab(ikcalc,spin)%col,(bmin:bmax))
       ABI_ALLOCATE(Sigp%Sigcij_tab(ikcalc,spin)%col,(bmin:bmax))

       do bcol=bmin,bmax
         ABI_ALLOCATE(sigc_bidx,(bmax-bmin+1))
         ABI_ALLOCATE(sigx_bidx,(bmax-bmin+1))

         if (Bnd_sym(ik_ibz,spin)%err_status/=0) then   ! Band classification failed.
           sigc_bidx = (/(ii,ii=bmin,bmax)/) 
           idx_c = bmax-bmin+1
           sigx_bidx = (/(ii,ii=bmin,bcol)/) ! Hermitian 
           idx_x = bcol-bmin+1
         else
           irr_idx2 = Bnd_sym(ik_ibz,spin)%b2irrep(bcol)
           idx_c = 0
           do brow=bmin,bmax
             irr_idx1 = Bnd_sym(ik_ibz,spin)%b2irrep(brow)
             if (sigma_is_herm(Sigp).and.bcol<brow) CYCLE  ! Only the upper triangle for HF, SEX, or COHSEX.
             if (irr_idx1 == irr_idx2) then ! same character, add this row to the list.
               idx_c = idx_c +1
               sigc_bidx(idx_c) = brow
             end if
           end do
           idx_x = 0
           do brow=bmin,bcol
             irr_idx1 = Bnd_sym(ik_ibz,spin)%b2irrep(brow) 
             if (bcol<brow) CYCLE  ! Sig_x is always Hermitian.
             if (irr_idx1 == irr_idx2) then ! same character, add this row to the list.
               idx_x = idx_x +1
               sigx_bidx(idx_x) = brow
             end if
           end do
         end if
         !
         ! Table for Sigma_x matrix elements taking into account symmetries of the bands.
         ABI_ALLOCATE(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx,(idx_x))

         Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%size1= idx_x
         Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(:) = sigx_bidx(1:idx_x)
         !write(std_out,*)" Sigxij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol
         !write(std_out,*)" size: ",idx_x,(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(ii),ii=1,idx_x)
         !
         ! Table for Sigma_c matrix elements taking into account symmetries of the bands.
         ABI_ALLOCATE(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx,(idx_c))

         Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%size1= idx_c
         Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(:) = sigc_bidx(1:idx_c)
         !write(std_out,*)" Sigcij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol
         !write(std_out,*)" size: ",idx_c,(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(ii), ii=1,idx_c)

         ABI_DEALLOCATE(sigx_bidx)
         ABI_DEALLOCATE(sigc_bidx)
       end do ! bcol

     else  ! Symmetries cannot be used for this (k,s).

       bmin=Sigp%minbnd(ikcalc,spin); bmax=Sigp%maxbnd(ikcalc,spin)
       ABI_ALLOCATE(Sigp%Sigcij_tab(ikcalc,spin)%col,(bmin:bmax))
       ABI_ALLOCATE(Sigp%Sigxij_tab(ikcalc,spin)%col,(bmin:bmax))

       if (Sigp%gwcalctyp<20) then  ! QP wavefunctions == KS, therefore only diagonal elements are calculated.
         do bcol=bmin,bmax
           ABI_ALLOCATE(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx,(1:1))
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%size1= 1
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(1) = bcol
           ABI_ALLOCATE(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx,(1:1))
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%size1= 1
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(1) = bcol
         end do
       else  
         ! Use QP wavefunctions, Sigma_ij matrix is sparse but we have to classify the states in sigma.
         ! The only thing we can do here is filling the entire matrix taking advantage of Hermiticity (if any).
         do bcol=bmin,bmax
           ABI_ALLOCATE(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx,(bcol-bmin+1))
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%size1= bcol-bmin+1
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(:) = (/(ii,ii=bmin,bcol)/) ! Sigma_x is Hermitian.
           !write(std_out,*)"Sigxij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol,Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(:)

           ABI_ALLOCATE(sigc_bidx,(bmax-bmin+1))
           idx_c = 0
           do brow=bmin,bmax
             if (sigma_is_herm(Sigp).and.bcol<brow) CYCLE  ! Only the upper triangle of Sigc_ij is needed (SEX, COHSEX). 
             idx_c = idx_c +1
             sigc_bidx(idx_c) = brow
           end do
           ABI_ALLOCATE(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx,(idx_c))
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%size1= idx_c  
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(:) = sigc_bidx(1:idx_c)
           ABI_DEALLOCATE(sigc_bidx)
           !write(std_out,*)"Sigcij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol,Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(:)
         end do
       end if
     end if

   end do !ikcalc
 end do !spin

end subroutine sigma_tables
!!***
