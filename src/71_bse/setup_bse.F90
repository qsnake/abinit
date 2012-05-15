!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_bse
!! NAME
!!  setup_bse
!!
!! FUNCTION
!!  This routine performs the initialization of basic objects and quantities used for Bethe-Salpeter calculations.
!!  In particular the excparam data type that defines the parameters of the calculation is completely
!!  initialized starting from the content of Dtset and the parameters read from the external KSS and SCR (SUSC) file.
!!
!! COPYRIGHT
!! Copyright (C) 1992-2009 EXC group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida)
!! Copyright (C) 2009-2012 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! codvsn=Code version
!! ngfft_gw(18)=Information about 3D FFT for density and potentials, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! acell(3)=Length scales of primitive translations (bohr)
!! rprim(3,3)=Dimensionless real space primitive translations.
!! Dtset<dataset_type>=All input variables for this dataset.
!!  Some of them might be redefined here TODO
!! Dtfil=filenames and unit numbers used in abinit.
!! Psps <pseudopotential_type>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat*Dtset%usepaw)<pawtab_type>=PAW tabulated starting data
!!
!! OUTPUT
!! Cryst<crystal_structure>=Info on the crystalline Structure. 
!! Kmesh<BZ_mesh_type>=Structure defining the k-sampling for the wavefunctions.
!! Qmesh<BZ_mesh_type>=Structure defining the q-sampling for the symmetrized inverse dielectric matrix.
!! Gsph_Max<Gvectors_type=Data type gathering info on the G-sphere for wave functions and e^{-1}, 
!! KS_BSt<Bandstructure_type>=The KS band structure (energies, occupancies, k-weights...)
!! Vcp<vcoul_t>=Structure gathering information on the Coulomb interaction in reciprocal space,
!!   including a possible cutoff in real space.
!! ngfft_osc(18)=Contain all needed information about the 3D FFT for the oscillator matrix elements.
!!   See ~abinit/doc/input_variables/vargs.htm#ngfft
!! Bsp<excparam>=Basic parameters defining the Bethe-Salpeter run. Completely initialed in output.
!! Hdr_kss<Hdr_type>=The header of the KSS file.
!! Hdr_bse<Hdr_type>=Local header initialized from the parameters used for the Bethe-Salpeter calculation.
!! BS_files<excfiles>=Files used in the calculation.
!! w_file=File name used to construct W. Set to ABI_NOFILE if no external file is used.
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      apply_scissor,bsp_calctype2str,bstruct_init,copy_bandstructure
!!      find_qmesh,free_scrhdr,get_bz_item,get_ng0sh,hdr_check,hdr_init
!!      hdr_update,hdr_vs_dtset,init_crystal_from_hdr,init_gsphere,init_kmesh
!!      init_transitions,initmpi_seq,make_mesh,matrginv,metric,mkrdim
!!      nullify_bs_parameters,print_bandstructure,print_bs_files
!!      print_bs_parameters,print_bz_mesh,print_crystal,print_gsphere
!!      print_ngfft,rdgw,reportgap,rhoij_alloc,rhoij_copy,rhoij_free,scr_hdr_io
!!      scrhdr_comm,setmesh,testkss,update_occ,vcoul_init,wrtout,xbarrier_mpi
!!      xcast_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setup_bse(codvsn,acell,rprim,ngfftf,ngfft_osc,Dtset,Dtfil,BS_files,Psps,Pawtab,BSp,&
& Cryst,Kmesh,Qmesh,KS_BSt,QP_bst,Hdr_kss,Gsph_Max,Gsph_c,Vcp,Hdr_bse,w_fname,comm,Wvl)

 use m_profiling

 use defs_basis
 use m_bs_defs
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_errors
 use m_timer
 use m_xmpi

 use m_gwdefs,        only : GW_Q0_DEFAULT
 use m_fstrings,      only : toupper
 use m_io_tools,      only : file_exist, get_unit
 use m_geometry,      only : normv
 use m_header,        only : hdr_init, hdr_get_nelect_byocc
 use m_crystal,       only : print_crystal, idx_spatial_inversion, crystal_structure
 use m_crystal_io,    only : init_crystal_from_hdr
 use m_bz_mesh,       only : bz_mesh_type, init_kmesh, get_ng0sh, print_BZ_mesh, get_BZ_item, find_qmesh, make_mesh
 use m_ebands,        only : bstruct_init, print_bandstructure, copy_bandstructure, bstruct_clean, &
&                            update_occ, get_valence_idx, apply_scissor, ReportGap
 use m_vcoul,         only : vcoul_t, vcoul_init
 use m_fft_mesh,      only : setmesh, print_ngfft
 use m_gsphere,       only : gvectors_type, init_gsphere, print_gsphere
 use m_io_screening,  only : print_scrhdr, ScrHdr_type, free_scrhdr, scr_hdr_io, scrhdr_comm
 use m_io_kss,        only : testkss
 use m_qparticles,    only : rdgw

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_bse'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=6),intent(in) :: codvsn
 character(len=fnlen),intent(out) :: w_fname
 type(dataset_type),intent(inout) :: Dtset 
 type(datafiles_type),intent(in) :: Dtfil
 type(pseudopotential_type),intent(in) :: Psps
 type(excparam),intent(inout) :: Bsp
 type(hdr_type),intent(out) :: Hdr_kss,Hdr_bse
 type(crystal_structure),intent(out) :: Cryst
 type(BZ_mesh_type),intent(out) :: Kmesh,Qmesh
 type(Gvectors_type),intent(out) :: Gsph_Max,Gsph_c
 type(Bandstructure_type),intent(out) :: KS_BSt,QP_Bst
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
 type(vcoul_t),intent(out) :: Vcp
 type(excfiles),intent(out) :: BS_files
 type(wvl_internal_type), intent(in) :: Wvl
!arrays
 integer,intent(in) :: ngfftf(18)
 integer,intent(out) :: ngfft_osc(18)
 real(dp),intent(in) :: acell(3),rprim(3,3)

!Local variables ------------------------------
!scalars
 integer,parameter :: pertcase0=0
 integer :: bantot,enforce_sym,ib,ibtot,ik_ibz,isppol,itypat,jj,method
 integer :: mpsang_kss,nbnds_kss,ng_kss,nsym_kss,io,istat
 integer :: nfftot_osc,spin,hexc_size,nqlwl,iq,ish,jsh,osl
 integer :: restart,restartpaw,timrev,iq_bz,isym,iq_ibz,itim,ig,mg0sh
 integer :: my_rank,master,w_unt,ios,rdwr,fform,npwe_file,ierr
 integer :: first_dig,second_dig,it
 real(dp) :: ucvol,nelect_hdr,qnorm,mvec1,mvec2
 logical :: ltest,only_one_kpt,remove_inv
 character(len=500) :: msg
 character(len=fnlen) :: gw_fname,test_file
 type(MPI_type) :: MPI_enreg_seq
 type(ScrHdr_type) :: Hscr
!arrays
 integer :: ng0sh_opt(3),val_idx(Dtset%nsppol)
 integer,allocatable :: npwarr(:),nlmn(:)
 integer,allocatable :: shlimtmp(:) 
 integer,allocatable :: val_indeces(:,:)
 integer,pointer :: gvec_p(:,:)
 real(dp) :: gg(3),qpt_bz(3)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),sq(3)
 real(dp) :: qred2cart(3,3),qcart2red(3,3)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:),qlwl(:,:)
 real(dp),allocatable :: igwene(:,:,:)
 real(dp),pointer :: energies_p(:,:,:)
 complex(dpc),allocatable :: gw_energy(:,:,:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:)

!************************************************************************

 DBG_ENTER("COLL")

 my_rank = xcomm_rank(comm)
 master  = 0

 call nullify_bs_parameters(BSp)

 ! === Check for calculations that are not implemented ===
 ltest=ALL(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol)==Dtset%nband(1))
 ABI_CHECK(ltest,'Dtset%nband must be constant')
 ABI_CHECK(Dtset%nspinor==1,"nspinor==2 not coded")

 !* Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)

 ! === Dimensional primitive translations rprimd (from input), gprimd, metrics and unit cell volume ===
 call mkrdim(acell,rprim,rprimd)  
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol) 

 ! === Define consistently npw, nsh, and ecut for wavefunctions and W ===
 !call setshells(Dtset%ecutwfn,Dtset%npwwfn,Dtset%nshwfn,Dtset%nsym,gmet,gprimd,Dtset%symrel,'wfn',ucvol)
 !call setshells(Dtset%ecuteps,Dtset%npweps,Dtset%nsheps,Dtset%nsym,gmet,gprimd,Dtset%symrel,'eps',ucvol)

 call testkss(Dtfil%fnameabi_kss,Dtset%accesswff,nsym_kss,nbnds_kss,&
&  ng_kss,mpsang_kss,gvec_p,energies_p,Hdr_kss,comm) 

 nelect_hdr = hdr_get_nelect_byocc(Hdr_kss)

 ! === Create crystal_structure data type ===
 remove_inv=(nsym_kss/=Hdr_kss%nsym) 
 timrev=  2 ! This information is not reported in the header
            ! 1 => do not use time-reversal symmetry 
            ! 2 => take advantage of time-reversal symmetry
                                                            
 call init_crystal_from_hdr(Cryst,Hdr_kss,timrev,remove_inv)
 call print_crystal(Cryst)

! here we reorder the shell, because most likely are badly written in outkss.
 call wrtout(std_out,' Reordering shells',"COLL")

 ABI_ALLOCATE(shlimtmp,(ng_kss))
 do ig=1,ng_kss
   shlimtmp(ig)=ig
 end do

 jsh   = 1
 mvec1 = zero
 do ish=2,ng_kss
   osl = shlimtmp(ish-1) + 1
   gg(:)  = gvec_p(:, osl)
   mvec2 = normv(gvec_p(:,osl),Cryst%gmet,"G")
!  write(std_out,*) ish, mvec2, gg

   if (abs(mvec2 - mvec1) > tol5) then
     jsh=jsh + 1
     mvec1 = mvec2
   end if
    
   shlimtmp(jsh) = shlimtmp(ish)
 end do
 BSp%nsh = jsh

 ABI_ALLOCATE(BSp%shlim,(BSp%nsh))
 BSp%shlim(:)=shlimtmp(1:BSp%nsh)
 ABI_DEALLOCATE(shlimtmp)
 
 ltest=(Psps%mpsang==mpsang_kss) 
 ABI_CHECK(ltest,'Psps%mpsang/=mpsang_kss')

 call hdr_vs_dtset(Hdr_kss,Dtset) 
 !
 ! Setup of the k-point list and symmetry tables in the  BZ -----------------------------------
 if (Dtset%chksymbreak==0) then
   MSG_WARNING("Calling make_mesh")
   call make_mesh(Kmesh,Cryst,Dtset%kptopt,Dtset%kptrlatt,Dtset%nshiftk,Dtset%shiftk,break_symmetry=.TRUE.)
   ! TODO
   !Check if kibz from KSS file corresponds to the one returned by make_mesh.
 else
   call init_kmesh(Kmesh,Cryst,Hdr_kss%nkpt,Hdr_kss%kptns,Dtset%kptopt)
 end if

 call print_BZ_mesh(Kmesh,"K-mesh for the wavefunctions",std_out,Dtset%prtvol,"COLL")
 call print_BZ_mesh(Kmesh,"K-mesh for the wavefunctions",ab_out, 0,           "COLL")

 nqlwl = 0

 w_fname = ABI_NOFILE
 if (Dtset%getscr/=0.or.Dtset%irdscr/=0) then
   w_fname=Dtfil%fnameabi_scr
 else if (Dtset%getsuscep/=0.or.Dtset%irdsuscep/=0) then
   w_fname=Dtfil%fnameabi_sus
   MSG_ERROR("(get|ird)suscep not implemented")
 end if

 if (w_fname /= ABI_NOFILE) then ! Read dimensions from the external file
   !
   if (my_rank==master) then ! Read npw and nqlwl from file.
     call wrtout(std_out,' Testing file: '//TRIM(w_fname),"COLL")
     w_unt = get_unit()
     open(unit=w_unt,file=w_fname,status='old',form='unformatted',iostat=ios)
     msg = ' Opening file '//TRIM(w_fname)//' as old '
     ABI_CHECK(ios==0,msg)

     rdwr=5 ! Read the header.
     call scr_hdr_io(fform,rdwr,w_unt,xmpi_self,0,Dtset%accesswff,Hscr)
     close(w_unt)
     if (Dtset%prtvol>0) then ! Echo the header.
       rdwr=4
       call scr_hdr_io(fform,rdwr,std_out,xmpi_self,0,Dtset%accesswff,Hscr)
     end if

     npwe_file = Hscr%npwe ! Have to change %npweps if it was larger than dim on disk.
     nqlwl     = Hscr%nqlwl

     if (Dtset%npweps>npwe_file) then
       write(msg,'(2(a,i0),2a,i0)')&
&       " The number of G-vectors stored on file (",npwe_file,") is smaller than Dtset%npweps = ",Dtset%npweps,ch10,&
&       " Calculation will proceed with the maximum available set, npwe_file = ",npwe_file
       MSG_WARNING(msg)
       Dtset%npweps = npwe_file
     else  if (Dtset%npweps<npwe_file) then
       write(msg,'(2(a,i0),2a,i0)')&
&       " The number of G-vectors stored on file (",npwe_file,") is larger than Dtset%npweps = ",Dtset%npweps,ch10,&
&       " Calculation will proceed with Dtset%npweps = ",Dtset%npweps
       MSG_COMMENT(msg)
     end if
   end if
   
   call xbarrier_mpi(comm)
   call scrhdr_comm(Hscr,master,my_rank,comm)
   call xcast_mpi(Dtset%npweps,master,comm,ierr)
   call xcast_mpi(nqlwl,master,comm,ierr)

   if (nqlwl>0) then
     ABI_ALLOCATE(qlwl,(3,nqlwl))
     qlwl = Hscr%qlwl
   end if
   !
   ! Init Qmesh from the SCR file.
   call init_kmesh(Qmesh,Cryst,Hscr%nqibz,Hscr%qibz,Dtset%kptopt)

   call free_scrhdr(Hscr)
 else 
   ! Init Qmesh from the K-mesh reported in the KSS file.
   call find_qmesh(Qmesh,Cryst,Kmesh)
 end if

 if (nqlwl==0) then 
   nqlwl=1
   ABI_ALLOCATE(qlwl,(3,nqlwl))
   qlwl(:,nqlwl)= GW_Q0_DEFAULT
   write(msg,'(3a,i2,a,3f9.6)')&
&    " The Header of the screening file does not contain the list of q-point for the optical limit ",ch10,&
&    " Using nqlwl= ",nqlwl," and qlwl = ",qlwl(:,1)  
   MSG_COMMENT(msg)
 end if
 !write(std_out,*)"nqlwl and qlwl for Coulomb singularity and e^-1",nqlwl,qlwl

 ! === Setup of q-mesh in the whole BZ ===
 ! * Stop if a nonzero umklapp is needed to reconstruct the BZ. In this case, indeed, 
 !   epsilon^-1(Sq) should be symmetrized in csigme using a different expression (G-G_o is needed)
 !
 call print_BZ_mesh(Qmesh,"Q-mesh for the screening function",std_out,Dtset%prtvol,"COLL")
 call print_BZ_mesh(Qmesh,"Q-mesh for the screening function",ab_out ,0           ,"COLL")

 do iq_bz=1,Qmesh%nbz
   call get_BZ_item(Qmesh,iq_bz,qpt_bz,iq_ibz,isym,itim)
   sq = (3-2*itim)*MATMUL(Cryst%symrec(:,:,isym),Qmesh%ibz(:,iq_ibz))
   if (ANY(ABS(Qmesh%bz(:,iq_bz)-sq )>1.0d-4)) then
     write(std_out,*) sq,Qmesh%bz(:,iq_bz) 
     write(msg,'(a,3f6.3,a,3f6.3,2a,9i3,a,i2,2a)')&
&      ' qpoint ',Qmesh%bz(:,iq_bz),' is the symmetric of ',Qmesh%ibz(:,iq_ibz),ch10,&
&      ' through operation ',Cryst%symrec(:,:,isym),' and itim ',itim,ch10,&
&      ' however a non zero umklapp G_o vector is required and this is not yet allowed'
     MSG_ERROR(msg)
   end if
 end do 

 BSp%algorithm = Dtset%bs_algorithm
 BSp%nstates   = Dtset%bs_nstates
 Bsp%nsppol    = Dtset%nsppol
 Bsp%hayd_term = Dtset%bs_hayd_term
 !
 ! Define the algorithm for solving the BSE.
 if (BSp%algorithm == BSE_ALGO_HAYDOCK) then
   BSp%niter       = Dtset%bs_haydock_niter
   BSp%haydock_tol = Dtset%bs_haydock_tol
 else if (BSp%algorithm == BSE_ALGO_CG) then
   ! FIXME For the time being use an hardcoded value.
   ! TODO change name in Dtset%
   BSp%niter       = 100   !Dtset%bs_haydock_niter
   BSp%cg_tolwfr   = tol12 !Dtset%bs_haydock_tol
   BSp%nline       = Dtset%nline
   BSp%nbdbuf      = 10    !Dtset%nbdbuf
   BSp%nstates     = Dtset%bs_nstates
   MSG_WARNING("Check CG setup")
 else
   !BSp%niter       = 0
   !BSp%tol_iter    = HUGE(one)
 end if
 !
 ! Setup of the frequency mesh for the absorption spectrum.
 BSp%broad  = Dtset%zcut

 Bsp%omegai = Dtset%bs_freq_mesh(1)
 Bsp%omegae = Dtset%bs_freq_mesh(2)
 Bsp%domega = Dtset%bs_freq_mesh(3)

 ! The frequency mesh (including the complex imaginary shift)
 BSp%nomega = (BSp%omegae - BSp%omegai)/BSp%domega + 1
 ABI_ALLOCATE(BSp%omega,(BSp%nomega))
 do io=1,BSp%nomega
   BSp%omega(io) = (BSp%omegai + (io-1)*BSp%domega)  + j_dpc*BSp%broad
 end do

 ! Possible cutoff on the transitions.
 BSp%ircut     = Dtset%bs_eh_cutoff(1)
 BSp%uvcut     = Dtset%bs_eh_cutoff(2)
 BSp%stripecut = Dtset%bs_eh_cutoff(3) 
 !
 ! Shall we include Local field effects?
 SELECT CASE (Dtset%bs_exchange_term)
 CASE (0,1)
   BSp%exchange_term = Dtset%bs_exchange_term
 CASE DEFAULT 
   write(msg,'(a,i0)')" Wrong bs_exchange_term: ",Dtset%bs_exchange_term
   MSG_ERROR(msg)
 END SELECT
 !
 !
 ! Treatment of the off-diagonal coupling block.
 SELECT CASE (Dtset%bs_coupling)
 CASE (0)
   BSp%use_coupling = 0
   msg = 'RESONANT ONLY CALCULATION'
 CASE (1)
   BSp%use_coupling = 1
   msg = ' RESONANT+COUPLING CALCULATION '
 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong bs_coupling: ",Dtset%bs_coupling
   MSG_ERROR(msg)
 END SELECT 
 call wrtout(std_out,msg,"COLL")
 !
 ! Treatment of the Coulomb term.
 BSp%use_diagonal_Wgg = .FALSE.
 Bsp%use_coulomb_term = .TRUE.
 BSp%eps_inf=zero; Bsp%mdlf_type=0

 first_dig =MOD(Dtset%bs_coulomb_term,10)
 second_dig=Dtset%bs_coulomb_term/10

 Bsp%wtype = second_dig
 SELECT CASE (second_dig)
 CASE (BSE_WTYPE_NONE)
   call wrtout(std_out,"Coulomb term won't be calculated","COLL")
   Bsp%use_coulomb_term = .FALSE.
 CASE (BSE_WTYPE_FROM_SCR)
   call wrtout(std_out,"W is read from an external SCR file","COLL")
   Bsp%use_coulomb_term = .TRUE.
 CASE (BSE_WTYPE_FROM_MDL)
   call wrtout(std_out,"W is approximated with the model dielectric function","COLL")
   MSG_WARNING(" Model dielectric function is still under testing")
   Bsp%use_coulomb_term = .TRUE.
   BSp%mdlf_type=1
   BSp%eps_inf = 26
 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong second digit in bs_coulomb_term: ",Dtset%bs_coulomb_term
   MSG_ERROR(msg)
 END SELECT
 ! 
 ! Diagonal approximation or full matrix?
 BSp%use_diagonal_Wgg = .TRUE.
 if (Bsp%wtype /= BSE_WTYPE_NONE) then 
   SELECT CASE (first_dig)
   CASE (0)
     call wrtout(std_out,"Using diagonal approximation W_GG","COLL")
     BSp%use_diagonal_Wgg = .TRUE.
   CASE (1)
     call wrtout(std_out,"Using full W_GG' matrix ","COLL")
     BSp%use_diagonal_Wgg = .FALSE.
   CASE DEFAULT
     write(msg,'(a,i0)')" Wrong first digit in bs_coulomb_term: ",Dtset%bs_coulomb_term
     MSG_ERROR(msg)
   END SELECT
 end if
 !
 ! Dimensions and parameters of the calculation.
 ! TODO one should add npwx as well
 BSp%npweps=Dtset%npweps
 BSp%npwwfn=Dtset%npwwfn

 BSp%lomo  = Dtset%bs_loband
 if (Bsp%lomo > Dtset%nelect/2) then ! correct only for unpolarized semiconductors 
   write(msg,'(a,i0,a,f8.3)') " Bsp%lomo = ",Bsp%lomo," cannot be greater than nelec/2 = ",Dtset%nelect/2
   MSG_ERROR(msg)
 end if
 !
 ! ==============================================
 ! ==== Setup of the q for the optical limit ====
 ! ==============================================
 Bsp%inclvkb = Dtset%inclvkb

 qred2cart = two_pi*Cryst%gprimd 
 qcart2red = qred2cart
 call matrginv(qcart2red,3,3)

 if (Dtset%gw_nqlwl==0) then
   BSp%nq = 6
   ABI_ALLOCATE(BSp%q,(3,BSp%nq))
   BSp%q(:,1) = (/one,zero,zero/)  ! (100)
   BSp%q(:,2) = (/zero,one,zero/)  ! (010)
   BSp%q(:,3) = (/zero,zero,one/)  ! (001)
   BSp%q(:,4) = MATMUL(qcart2red,(/one,zero,zero/)) ! (x)
   BSp%q(:,5) = MATMUL(qcart2red,(/zero,one,zero/)) ! (y)
   BSp%q(:,6) = MATMUL(qcart2red,(/zero,zero,one/)) ! (z)
 else
   BSp%nq = Dtset%gw_nqlwl
   ABI_ALLOCATE(BSp%q,(3,BSp%nq))
   BSp%q = Dtset%gw_qlwl
 end if

 do iq=1,BSp%nq ! normalization 
   qnorm = normv(BSp%q(:,iq),Cryst%gmet,"G") 
   BSp%q(:,iq) = BSp%q(:,iq)/qnorm 
 end do
 !
 ! ======================================================
 ! === Define the flags defining the calculation type ===
 ! ======================================================
 Bsp%calc_type = Dtset%bs_calctype
 
 BSp%soenergy = zero ! Shall we use the scissors operator to open the gap?
 if (ABS(Dtset%soenergy)>tol6) BSp%soenergy = Dtset%soenergy

!now test input parameters from input and kss file and assume some defaults
!
! TODO Add the possibility of using a randomly shifted k-mesh with nsym>1.
! so that densities and potentials are correctly symmetrized but 
! the list of the k-point in the IBZ is not expanded.
!MG NOW Have to copy some values to BSp%
 BSp%nkibz  = Kmesh%nibz  !We might allow for a smaller number of points....
 BSp%npwwfn = Dtset%npwwfn
 BSp%npweps = Dtset%npweps 

 !TODO add new dim for exchange part and consider the possibility of 
 !having npwsigx > npwwfn (see setup_sigma).

 if (BSp%npwwfn>ng_kss) then
   BSp%npwwfn=ng_kss 
   write(msg,'(2a,(a,i8,a))')&
&    ' Number of G-vectors in KSS file found less than required',ch10,&
&    '  calculation will proceed with npwwfn  = ',BSp%npwwfn,ch10
   MSG_WARNING(msg)
 end if

!test if shells are closed
 do ish=1,BSp%nsh
   if (BSp%npweps >  BSp%shlim(ish)) cycle
   if (BSp%npweps == BSp%shlim(ish)) then
     EXIT
   else  ! that means npweps is in between two closed shells
     BSp%npweps=BSp%shlim(ish-1)
     write(msg,'(3a,i0)')&
&      ' npweps is not a closed shell',ch10,&
&      ' calculation continues with npweps =', BSp%shlim(ish-1)
     MSG_COMMENT(msg)
     EXIT
   end if
 end do

 do ish=1,BSp%nsh
   if (BSp%npwwfn  > BSp%shlim(ish)) cycle
   if (BSp%npwwfn == BSp%shlim(ish)) then
     EXIT
   else  ! that means npwwfn is in between two closed shells
     BSp%npwwfn=BSp%shlim(ish-1)
     write(msg,'(a,i0)')&
&      ' npwwfn is not a closed shell. Calculation continues with npweps =', BSp%shlim(ish-1)
     MSG_COMMENT(msg)
     EXIT 
   end if
 end do
     
!calculate NPWVEC as the biggest between npweps and npwwfn.
!MG RECHECK this part.
 BSp%npwvec=MAX(BSp%npwwfn,BSp%npweps)

 if (ng_kss < BSp%npwvec) then
   write(msg,'(a,i0,a)')" KSS file contains only ", ng_kss, " planewaves"
   MSG_WARNING(msg)
   BSp%npwvec = ng_kss
   if (BSp%npwwfn > ng_kss) then
     BSp%npwwfn = ng_kss
     write(msg,'(a,i0)') ' Calculation will be done with npwwfn = ',BSp%npwwfn
     call wrtout(std_out,msg,"COLL")
   end if
   if (BSp%npweps > ng_kss) then
     BSp%npweps = ng_kss
     write(msg,'(a,i0)') ' Calculation will be done with npweps = ',BSp%npweps
     call wrtout(std_out,msg,"COLL")
   end if
 end if

 if (ABS(Dtset%nelect-nelect_hdr)>tol6) then
   write(msg,'(2(a,f8.2))')&
&   " File contains ", nelect_hdr," electrons but nelect initialized from input is ",Dtset%nelect
   MSG_ERROR(msg)
 end if 

 if (nbnds_kss < Dtset%nband(1)) then
   write(msg,'(2(a,i0),3a,i0)')&
&    ' KSS file contains only ', nbnds_kss,' levels instead of ',Dtset%nband(1),' required;',ch10,&
&    ' The calculation will be done with nbands= ',nbnds_kss
   MSG_WARNING(msg)
   Dtset%nband(:) = nbnds_kss
 end if

 BSp%nbnds = Dtset%nband(1) ! TODO Note the change in the meaning of input variables

 if (BSp%nbnds<=Dtset%nelect/2) then
   write(msg,'(2a,2(a,i0))')&
&    ' BSp%nbnds cannot be smaller than homo ',ch10,&
&    ' while BSp%nbnds = ',BSp%nbnds,' and Dtset%nelect= ',Dtset%nelect
   MSG_ERROR(msg)
 end if

 ! === Make biggest G-sphere of BSP%npwvec vectors ===
 only_one_kpt=.FALSE. !(Kmesh%nbz==1)
 call init_gsphere(Gsph_Max,only_one_kpt,Cryst,BSp%npwvec,gvec=gvec_p)

 call print_gsphere(Gsph_Max,unit=std_out,prtvol=Dtset%prtvol)

 ! G-sphere for W and Sigma_c is initialized from the SCR file.
 call init_gsphere(Gsph_c,only_one_kpt,Cryst,BSp%npweps,gvec=gvec_p)

 ABI_DEALLOCATE(gvec_p)

 call vcoul_init(Vcp,Gsph_Max,Qmesh,Kmesh,Dtset%rcut,Dtset%icutcoul,Dtset%vcutgeo,BSp%npwvec,nqlwl,qlwl,Cryst%rprimd,ngfftf,comm)
 
 ABI_DEALLOCATE(qlwl)

 bantot=SUM(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol))                                                       
 ABI_ALLOCATE(doccde,(bantot))
 ABI_ALLOCATE(eigen,(bantot))
 ABI_ALLOCATE(occfact,(bantot))
 doccde=zero; eigen=zero; occfact=zero 

 jj=0; ibtot=0
 do isppol=1,Dtset%nsppol
   do ik_ibz=1,Dtset%nkpt
     do ib=1,Hdr_kss%nband(ik_ibz+(isppol-1)*Dtset%nkpt)
       ibtot=ibtot+1
       if (ib<=BSP%nbnds) then  
         jj=jj+1
         occfact(jj)=Hdr_kss%occ(ibtot)
         eigen  (jj)=energies_p(ib,ik_ibz,isppol)
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
 npwarr=BSP%npwwfn

 call bstruct_init(bantot,KS_BSt,Dtset%nelect,doccde,eigen,Dtset%istwfk,Kmesh%ibz,Dtset%nband,&
&  Kmesh%nibz,npwarr,Dtset%nsppol,Dtset%nspinor,Dtset%tphysel,Dtset%tsmear,Dtset%occopt,occfact,Kmesh%wt) 

 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(eigen)
 ABI_DEALLOCATE(npwarr)

 !TODO Occupancies are zero if NSCF. One should calculate the occupancies from the energies when 
 ! the occupation scheme for semiconductors is used.
 call update_occ(KS_BSt,Dtset%fixmom,prtvol=Dtset%prtvol)

 call print_bandstructure(KS_BSt,"Band structure read from the KSS file",unit=std_out,prtvol=Dtset%prtvol)

 call ReportGap(KS_BSt,header=" KS band structure",unit=std_out,mode_paral="COLL")

 ABI_ALLOCATE(val_indeces,(KS_BSt%nkpt,KS_BSt%nsppol))

 val_indeces = get_valence_idx(KS_BSt)

 do spin=1,KS_BSt%nsppol
   val_idx(spin) = val_indeces(1,spin)
   write(std_out,*)" spin : ",spin," val_idx ",val_idx(spin)
   if ( ANY(val_indeces(1,spin) /= val_indeces(:,spin)) ) then
     MSG_ERROR(" BSE code does not support metals") 
   end if
 end do

 ABI_DEALLOCATE(val_indeces)
 !
 ! === Create the BSE header === 
 call hdr_init(KS_BSt,codvsn,Dtset,Hdr_bse,Pawtab,pertcase0,Psps,wvl)

 ! === Get Pawrhoij from the header of the KSS file ===
 ABI_ALLOCATE(Pawrhoij,(Cryst%natom*Dtset%usepaw))
 if (Dtset%usepaw==1) then
   ABI_ALLOCATE(nlmn,(Cryst%ntypat))
   do itypat=1,Cryst%ntypat
     nlmn(itypat)=Pawtab(itypat)%lmn_size
   end do
   call rhoij_alloc(1,nlmn,Dtset%nspden,Dtset%nspinor,Dtset%nsppol,Pawrhoij,Cryst%typat)
   ABI_DEALLOCATE(nlmn)
   call rhoij_copy(Hdr_kss%Pawrhoij,Pawrhoij)
 end if

 call hdr_update(bantot,1.0d20,1.0d20,Hdr_bse,Cryst%natom,1.0d20,&
&  Cryst%rprimd,occfact,MPI_enreg_seq,Pawrhoij,Dtset%usepaw,Cryst%xred)

! call hdr_update(bantot,1.0d20,1.0d20,Hdr_bse,Cryst%natom,1.0d20,&
!& Cryst%rprimd,KS_BSt%occ,Pawrhoij,Dtset%usepaw,Cryst%xred)

 ABI_DEALLOCATE(occfact)

 ! This is just to do a check, the file format is wrong!
 call hdr_check(1002,1002,Hdr_bse,Hdr_kss,'COLL',restart,restartpaw)

 if (Dtset%usepaw==1) call rhoij_free(Pawrhoij) 
 ABI_DEALLOCATE(Pawrhoij)

 ! === Find optimal value for G-sphere enlargment due to oscillator matrix elements ===
 ! * Here I have to be sure that Qmesh%bz is always inside the BZ, not always true size bz is buggy
 ! * -one is used because we loop over all the possibile differences, unlike screening
 mg0sh=5
 call get_ng0sh(Kmesh%nbz,Kmesh%bz,Kmesh%nbz,Kmesh%bz,Qmesh%nbz,Qmesh%bz,Cryst%gmet,-one,mg0sh,ng0sh_opt)
 BSp%mg0(:)=ng0sh_opt(:)

 ! === Setup of the FFT mesh for the oscilator strengths === 
 ! * ngfft_osc(7:18)==Dtset%ngfft(7:18) which is initialized before entering screening.
 ! * Here we redefine ngfft_osc(1:6) according to the following options :
 !
 ! method==0 --> FFT grid read from fft.in (debugging purpose)
 ! method==1 --> Normal FFT mesh
 ! method==2 --> Slightly augmented FFT grid to calculate exactly rho_tw_g (see setmesh.F90)
 ! method==3 --> Doubled FFT grid, same as the the FFT for the density,
 !
 ! enforce_sym==1 ==> Enforce a FFT mesh compatible with all the symmetry operation and FFT library
 ! enforce_sym==0 ==> Find the smallest FFT grid compatbile with the library, do not care about symmetries
 !
 ngfft_osc(1:18)=Dtset%ngfft(1:18); method=2
 if (Dtset%fftgw==00 .or. Dtset%fftgw==01) method=0
 if (Dtset%fftgw==10 .or. Dtset%fftgw==11) method=1
 if (Dtset%fftgw==20 .or. Dtset%fftgw==21) method=2
 if (Dtset%fftgw==30 .or. Dtset%fftgw==31) method=3
 enforce_sym=MOD(Dtset%fftgw,10) 
 
 call setmesh(gmet,Gsph_Max%gvec,ngfft_osc,BSp%npwvec,BSp%npweps,BSp%npwwfn,nfftot_osc,method,BSp%mg0,Cryst,enforce_sym)
 nfftot_osc=PRODUCT(ngfft_osc(1:3))

 call print_ngfft(ngfft_osc,"FFT mesh for oscillator matrix elements",std_out,"COLL",prtvol=Dtset%prtvol)
 !
 ! BSp%homo gives the highest occupied band
 BSp%homo  = val_idx(1)

 ! TODO generalize the code to account for this unlikely case.
 if (Dtset%nsppol==2) then 
   ABI_CHECK(BSp%homo == val_idx(2),"Different valence indeces for spin up and down")
 end if

 BSp%lumo  = BSp%homo + 1
 BSp%humo  = BSp%nbnds

 BSp%nbndv = BSp%homo  - BSp%lomo + 1
 BSp%nbndc = BSp%nbnds - BSp%homo

 BSp%nkbz = Kmesh%nbz

 call copy_bandstructure(KS_BSt,QP_bst)
 ABI_ALLOCATE(igwene,(QP_bst%mband,QP_bst%nkpt,QP_bst%nsppol))
 igwene=zero

 call bsp_calctype2str(Bsp,msg)
 call wrtout(std_out,"Calculation type: "//TRIM(msg))

 SELECT CASE (Bsp%calc_type)
 CASE (BSE_HTYPE_RPA_KS)
   if (ABS(BSp%soenergy)>tol6) then
     write(msg,'(a,f8.2,a)')' Applying a scissors operator energy= ',BSp%soenergy*Ha_eV," [eV] on top of the KS energies."
     call wrtout(std_out,msg,"COLL")
     call apply_scissor(QP_BSt,BSp%soenergy)
   else
     write(msg,'(a,f8.2,a)')' Using KS energies since soenergy= ',BSp%soenergy*Ha_eV," [eV]."
     call wrtout(std_out,msg,"COLL")
   end if
   !
 CASE (BSE_HTYPE_RPA_QPENE)
   ! Read _GW files with the corrections TODO here I should introduce variable getgw
  
   gw_fname=TRIM(Dtfil%filnam_ds(4))//'_GW'
   gw_fname="__in.gw__"
   if (.not.file_exist(gw_fname)) then
     msg = " File "//TRIM(gw_fname)//" not found. Aborting now"
     MSG_ERROR(msg)
   end if
                                                                                                                 
   call rdgw(QP_Bst,gw_fname,igwene,extrapolate=.FALSE.) ! here gwenergy is real
                                                                                                                 
   do isppol=1,Dtset%nsppol
     write(std_out,*) ' k       GW energies [eV]'
     do ik_ibz=1,Kmesh%nibz
       write(std_out,'(i3,7x,10f7.2/50(10x,10f7.2/))')ik_ibz,(QP_bst%eig(ib,ik_ibz,isppol)*Ha_eV,ib=1,BSp%nbnds)
     end do
     write(std_out,*) ' k       Im GW energies [eV]'
     do ik_ibz=1,Kmesh%nibz
       write(std_out,'(i3,7x,10f7.2/50(10x,10f7.2/))')ik_ibz,(igwene(ib,ik_ibz,isppol)*Ha_eV,ib=1,BSp%nbnds)
     end do
   end do
   !
   ! If required apply the scissors operator on top of the QP bands structure (!)
   if (ABS(BSp%soenergy)>tol6) then
     write(msg,'(a,f8.2,a)')' Applying a scissors operator ',BSp%soenergy*Ha_eV," [eV] on top of the QP energies!"
     MSG_COMMENT(msg)
     call apply_scissor(QP_BSt,BSp%soenergy)
   end if
 CASE (BSE_HTYPE_RPA_QP)
   MSG_ERROR("Not implemented error!")
 CASE DEFAULT
   write(msg,'(a,i0)')"Unknown value for Bsp%calc_type= ",Bsp%calc_type
   MSG_ERROR(msg)
 END SELECT

 call ReportGap(QP_BSt,header=" QP band structure",unit=std_out,mode_paral="COLL")

 ! Transitions are ALWAYS ordered in c-v-k mode with k being the slowest index.
 ! FIXME: linewidths not coded.
 ABI_ALLOCATE(gw_energy,(BSp%nbnds,Kmesh%nibz,Dtset%nsppol))
 istat = ABI_ALLOC_STAT
 gw_energy = QP_BSt%eig

 BSp%have_complex_ene = ANY(igwene > tol16)

 ABI_CHECK(BSp%nbndv==BSp%homo-BSp%lomo+1,"Wrong nbndv")
 ABI_CHECK(BSp%nbndc==BSp%humo-BSp%lumo+1,"Wrong nbndc")

 ABI_ALLOCATE(Bsp%nreh,(Dtset%nsppol))
 Bsp%nreh=zero

 call init_transitions(BSp%Trans,BSp%lomo,BSp%humo,BSp%ircut,Bsp%uvcut,BSp%nkbz,Bsp%nbnds,Bsp%nkibz,Dtset%nsppol,Dtset%nspinor,&
&  gw_energy,QP_BSt%occ,Kmesh%tab,Bsp%nreh)

 ABI_DEALLOCATE(gw_energy)
 ABI_DEALLOCATE(igwene)

 do spin=1,Dtset%nsppol
   write(msg,'(a,i2,a,i0)')" For spin: ",spin,' the number of resonant e-h transitions is: ',BSp%nreh(spin)
   call wrtout(std_out,msg,"COLL")
 end do
                                                                                                            
 if (ANY(Bsp%nreh/=Bsp%nreh(1))) then
   write(msg,'(a,(i0))')" BSE code does not support different number of transitions for the two spin channels",Bsp%nreh
   MSG_ERROR(msg)
 end if
 !
 ! Create transition table vcks2t
 ABI_ALLOCATE(Bsp%vcks2t,(BSp%lomo:BSp%homo,BSp%lumo:BSp%humo,BSp%nkbz,Dtset%nsppol))
 Bsp%vcks2t = 0
                                                                                                            
 do spin=1,Dtset%nsppol
   do it=1,BSp%nreh(spin)
     BSp%vcks2t(BSp%Trans(it,spin)%v,BSp%Trans(it,spin)%c,BSp%Trans(it,spin)%k,spin) = it
   end do 
 end do

 hexc_size = SUM(Bsp%nreh); if (Bsp%use_coupling>0) hexc_size=2*hexc_size
 if (Bsp%nstates<=0) then
   Bsp%nstates=hexc_size
 else
   if (Bsp%nstates>hexc_size) then
      Bsp%nstates=hexc_size
      write(msg,'(2(a,i0),2a)')&
&      " Since the total size of excitonic Hamiltonian ",hexc_size," is smaller than Bsp%nstates ",Bsp%nstates,ch10,&
&      " the number of excitonic states nstates has been modified"
     MSG_WARNING(msg)
   end if
 end if

 msg=' Fundamental parameters for the solution of the Bethe-Salpeter equation:'
 call print_bs_parameters(BSp,unit=std_out,header=msg,mode_paral="COLL",prtvol=Dtset%prtvol) 
 call print_bs_parameters(BSp,unit=ab_out, header=msg,mode_paral="COLL") 

 if ( ANY (Cryst%symrec(:,:,1) /= RESHAPE ( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )) .or. &
&     ANY( ABS(Cryst%tnons(:,1)) > tol6) ) then
   write(msg,'(3a,9i2,2a,3f6.3,2a)')&
&    " The first symmetry operation should be the Identity with zero tnons while ",ch10,&
&    " symrec(:,:,1) = ",Cryst%symrec(:,:,1),ch10,&
&    " tnons(:,1)    = ",Cryst%tnons(:,1),ch10,&
&    " This is not allowed, sym_rhotwgq0 should be changed."
   MSG_ERROR(msg)
 end if
 !
 ! Prefix for generic output files.
 BS_files%out_basename = TRIM(Dtfil%filnam_ds(4))
 !
 ! Search for files to restart from. 
 if (Dtset%gethaydock/=0 .or. Dtset%irdhaydock/=0) then
   BS_files%in_haydock_basename = TRIM(Dtfil%fnameabi_haydock)
 end if

 test_file = Dtfil%fnameabi_bsham_reso
 if (file_exist(test_file)) then
   BS_files%in_hreso = test_file
 else 
   BS_files%out_hreso = TRIM(Dtfil%filnam_ds(4))//'_BSR'
 end if

 test_file = Dtfil%fnameabi_bsham_coup
 if ( file_exist(test_file) ) then
   BS_files%in_hcoup = test_file
 else 
   BS_files%out_hcoup = TRIM(Dtfil%filnam_ds(4))//'_BSC'
 end if
 !
 ! in_eig is the name of the input file with eigenvalues and eigenvectors 
 ! constructed from getbseig or irdbseig. out_eig is the name of the output file
 ! produced by this dataset. in_eig_exists checks for the presence of the input file.
 !
 if ( file_exist(Dtfil%fnameabi_bseig) ) then
   BS_files%in_eig = Dtfil%fnameabi_bseig
 else 
   BS_files%out_eig = TRIM(BS_files%out_basename)//"_BSEIG"
 end if

 call print_bs_files(BS_files,unit=std_out)
 !
 ! ==========================================================
 ! ==== Final check on the parameters of the calculation ====
 ! ==========================================================
 if ( Bsp%use_coupling>0 .and. ALL(Bsp%algorithm/=(/BSE_ALGO_DDIAGO, BSE_ALGO_HAYDOCK/)) ) then
   msg = " Resonant+Coupling is only available with the direct diagonalization or the haydock method."
   MSG_ERROR(msg)
 end if

 DBG_EXIT("COLL")

end subroutine setup_bse
!!***
