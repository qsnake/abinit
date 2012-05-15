!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_bs_defs
!! NAME
!!  m_bs_defs
!!
!! FUNCTION
!!  This module defines basic structures used for Bethe-Salpeter calculations.
!!
!! COPYRIGHT
!! Copyright (C) 1992-2009 EXC and ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_bs_defs

 use m_profiling

 use defs_basis 
 use m_errors 

 use m_io_tools,  only : get_unit

 implicit none

 private 
!!***

! Algorithm used to solve the BS problem.
 integer,public,parameter :: BSE_ALGO_DDIAGO  =1  ! Direct diagonalization.
 integer,public,parameter :: BSE_ALGO_HAYDOCK =2  ! Haydock recursion method.
 integer,public,parameter :: BSE_ALGO_CG      =3  ! Iterative diagonalization via CG method.

! Flags defining the content of the files used to restart the Haydock method.
 integer,public,parameter :: BSE_HAYD_IMEPS   =1
 integer,public,parameter :: BSE_HAYD_DOS     =2

! Approximations for the excitonic Hamiltonian.
 integer,public,parameter :: BSE_HTYPE_RPA_KS    =1  ! Use KS results to construct the RPA polarizability.
 integer,public,parameter :: BSE_HTYPE_RPA_QPENE =2  ! Use KS orbitals and QP energies to construct the RPA polarizability.
 integer,public,parameter :: BSE_HTYPE_RPA_QP    =3  ! Use QP orbitals and energies to construct the RPA polarizability.

! Flags for the treatment of W(G1,G2)
 integer,public,parameter :: BSE_WTYPE_NONE       =0 ! Coulomb term not included
 integer,public,parameter :: BSE_WTYPE_FROM_SCR   =1 ! W is read from a SCR file
 integer,public,parameter :: BSE_WTYPE_FROM_MDL   =2 ! W is approximated with a model dielectric function.

!$! Treatment of W(omega)
!! integer,public,parameter :: BSE_WFREQ_STATIC=1
!! integer,public,parameter :: BSE_WFREQ_PPM   =2
!! integer,public,parameter :: BSE_WFREQ_FULL  =3

 character(len=fnlen),public,parameter :: BSE_NOFILE="None"

!----------------------------------------------------------------------

!!****t* m_bs_defs/transition
!! NAME
!!  transition
!!
!! FUNCTION
!!  The transition derived data type is used to store the correspondence 
!!  between the transition index and the set of quantum numbers (ik_bz,v,c) 
!!  The energy of the transition is stored as well.
!!
!! SOURCE                                                                                   

 type,public :: transition
   integer :: k=0      ! Index of the k-point in the BZ
   integer :: v=0      ! Valence band index.
   integer :: c=0      ! Conduction band index.
   complex(dpc) :: en  ! Transition energy 
 end type transition
!!***

 public :: init_transitions     ! Main creation method.
 public :: repr_trans           ! Returns a string representing the transition or a couple of transitions.

 interface repr_trans 
   module procedure repr_1trans
   module procedure repr_2trans
 end interface repr_trans 

!----------------------------------------------------------------------

!!****t* m_bs_defs/excparam
!! NAME
!!  excparam
!!
!! FUNCTION
!!  The excparam derived data type contains the parameters controlling the BS calculation. 
!!
!! SOURCE

type,public :: excparam

!scalars
  integer :: algorithm         ! Algorithm used. 
  integer :: calc_type         ! Calculation type (see Dtset%bs_calc_type).
  integer :: hayd_term         ! Option for the terminator used in the Haydock solver.
  integer :: use_coupling      ! Include off-diagonal block coupling resonant and anti-resonant transitions.
  integer :: exchange_term     ! Include the exchange term in the BS Hamiltonian.
  integer :: inclvkb
  integer :: mdlf_type         ! Model dielectric function type. 
  integer :: nline             ! Number of line minimization used for CG minimization.
  integer :: nbdbuf            ! Number of states in the buffer that will be excluded from the convergence check.
  integer :: nstates           ! Number of states that will be considered in the CG minimization.
  integer :: npweps            ! No. of G in the Screening.
  integer :: npwwfn            ! No. of G for wave functions.
  !$integer :: npwx            ! No. of G for the exchange part.
  integer :: npwvec            ! MAX between npwwfn and npweps
  integer :: nsh               ! Number of shells corresponding to npwvec=MAX(npwwfn,npweps)
  integer :: nbnds             ! Total number of bands considered.
  integer :: nbndv             ! No. of valence states treated (homo-lomo+1)
  integer :: nbndc             ! No. of conduction states (humo-lumo+1)
  integer :: lomo,homo         ! Lowest and highest occupied orbital considered.
  integer :: lumo,humo         ! Lowest and highest unoccupied orbital considered.
  integer :: niter             ! No. of iterations for (Haydock|CG).
  integer :: nkibz, nkbz       ! No. of k-points in the IBZ and BZ (resp.)
  integer :: nomega            ! No. of frequencies for epsilon.
  integer :: nq                ! Number of "small" q for optical limit.
  integer :: nsppol            ! Number of independent spin polarizations.
  integer :: wtype             ! Option used for dealing with W (see BSE_WTYPE_) flags

  real(dp) :: eps_inf          ! Electronic dielectric constant used for the model dielectric function.
  real(dp) :: soenergy         ! Scissors energy (used if it absolute value is > tol6)
  real(dp) :: omegai           ! First omega for epsilon.
  real(dp) :: omegae           ! Last omega for epsilon (defaults to 10eV)
  real(dp) :: domega           ! Step of the frequency mesh.
  real(dp) :: broad            ! Lorentzian Broadening.
  real(dp) :: ircut            ! Infrared cutoff for transitions
  real(dp) :: uvcut            ! Ultraviolet cutoff for transitions.
  real(dp) :: stripecut        ! Cutoff on (t,t") transitions.
  real(dp) :: haydock_tol(2)   ! Tolerance for stopping the Haydock algorithm.
  real(dp) :: cg_tolwfr        ! Tolerance for stopping the CG algorithm 

  logical :: use_diagonal_Wgg  ! Use diagonal approximation for Wgg.
  logical :: use_coulomb_term  ! Include W term in the BS Hamiltonian.
  logical :: have_complex_ene  ! .TRUE. if energies have a non-zero imaginary part.

!arrays
  integer :: mg0(3)            ! For each reduced direction gives the max G0 component
                               ! to account for umklapp processes
  integer,pointer :: nreh(:)
  ! nreh(nsppol)
  ! Number of resonant electron-hole transitions for each spin.

  integer,pointer :: vcks2t(:,:,:,:)   SET2NULL
  ! vcks2t(v,c,ik_bz,spin) gives the transition index associated to (v,c,kbz,spin)

  integer,pointer :: shlim(:)   SET2NULL

  real(dp),pointer :: q(:,:)           ! Q-points for optical limit (reduced coordinates).

  complex(dpc),pointer :: omega(:)  SET2NULL
  ! omega(nomega)
  ! Frequency mesh for epsilon (including the complex imaginary shift)

  type(transition),pointer :: Trans(:,:) SET2NULL

end type excparam
!!***

 public :: nullify_bs_parameters 
 public :: destroy_bs_parameters
 public :: print_bs_parameters
 public :: bsp_calctype2str

!!****t* m_bs_defs/excfiles
!! NAME
!!  excfiles
!!
!! FUNCTION
!!  The excfiles derived data type contains file names and unit numbers used to store 
!!  temporary or final results of the Bethe-Salpeter calculation. 
!!
!! SOURCE

type,public :: excfiles

  character(len=fnlen) :: in_hreso = BSE_NOFILE   
  ! Name of the input file with the resonant part of the Hamiltonian (Hermitian).

  character(len=fnlen) :: out_hreso = BSE_NOFILE
  ! Name of the output file with the resonant part of the Hamiltonian (Hermitian).

  character(len=fnlen) :: in_hcoup = BSE_NOFILE   
  ! Name of the input file with the coupling part of the Hamiltonian (Symmetric).
                                                                                   
  character(len=fnlen) :: out_hcoup = BSE_NOFILE
  ! Name of the output file with the coupling part of the Hamiltonian (Symmetric).

  character(len=fnlen) :: in_eig = BSE_NOFILE
  ! Name of the input file with the eigenvalues and the eigenvectors of the Hamiltonian.

  character(len=fnlen) :: out_eig = BSE_NOFILE
  ! Name of the output file with the eigenvalues and the eigenvectors of the Hamiltonian.

  character(len=fnlen) :: in_haydock_basename = BSE_NOFILE
  ! Name of the input file used to restart Haydock algorithm. 

  character(len=fnlen) :: out_basename = BSE_NOFILE
  ! Prefix to be used for other output files.

end type excfiles
!!***

public :: print_bs_files    ! Printout of the excfiles data type.


CONTAINS  !========================================================================================================
!!***

!!****f* m_bs_defs/nullify_bs_parameters
!! NAME
!!  nullify_bs_parameters
!!
!! FUNCTION
!!  Set all pointers in a structure of type excparam to NULL.
!!
!! SIDE EFFECTS
!!  Bsp<excparam>=All pointers are initialized to NULL.
!!
!! PARENTS
!!      setup_bse
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine nullify_bs_parameters(BSp)
    
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_bs_parameters'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(excparam),intent(inout) :: BSp

!************************************************************************

 !@excparam
 nullify(BSp%q)

 nullify(Bsp%nreh)
 nullify(Bsp%vcks2t)
 nullify(Bsp%shlim)
 nullify(Bsp%omega)

 nullify(Bsp%Trans)

end subroutine nullify_bs_parameters
!!***

!----------------------------------------------------------------------

!!****f* m_bs_defs/destroy_bs_parameters
!! NAME
!!  destroy_bs_parameters
!!
!! FUNCTION
!!  Free all memory allocated in a structure of type excparam
!!
!! SIDE EFFECTS
!!  Bsp<excparam>=All associated pointers are deallocated.
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_bs_parameters(BSp)
    
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_bs_parameters'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(excparam),intent(inout) :: BSp

!************************************************************************
 
 !@excparam
 if (associated(BSp%q))    then
   ABI_DEALLOCATE(BSp%q)
 end if

 if (associated(Bsp%nreh))    then
   ABI_DEALLOCATE(Bsp%nreh)
 end if
 if (associated(Bsp%vcks2t))  then
   ABI_DEALLOCATE(Bsp%vcks2t)
 end if
 if (associated(Bsp%shlim))   then
   ABI_DEALLOCATE(Bsp%shlim)
 end if
 if (associated(Bsp%omega))   then
   ABI_DEALLOCATE(Bsp%omega)
 end if

 if (associated(Bsp%Trans))  then
   ABI_DEALLOCATE(Bsp%Trans)
 end if

end subroutine destroy_bs_parameters
!!***

!----------------------------------------------------------------------

!!****f* m_bs_defs/print_bs_parameters
!! NAME
!!  print_bs_parameters
!!
!! FUNCTION
!!  Printout of the parameters used for the BS calculation.
!!
!! INPUTS
!!  p<excparam>=Datatype storing the parameters of the Bethe-Salpeter calculation.
!!
!! OUTPUT
!!  Only printing.
!!
!! PARENTS
!!      setup_bse
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_bs_parameters(BSp,header,unit,mode_paral,prtvol) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_bs_parameters'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral 
 character(len=*),optional,intent(in) :: header
 type(excparam),intent(inout) :: BSp

!Local variables ------------------------------
!scalars
 integer :: my_unt,my_prtvol,iq,ii
 character(len=4) :: my_mode
 character(len=500) :: msg      

! ********************************************************************* 

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol 
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Parameters of the Bethe-Salpeter  run ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 select case (Bsp%algorithm)
 case (BSE_ALGO_DDIAGO)
   msg  = " Algorithm: Direct diagonalization."
 case (BSE_ALGO_HAYDOCK)
   msg  = " Algorithm: Haydock technique."
 case (BSE_ALGO_CG)
   msg  = " Algorithm: Conjugate gradient."
 case default
   msg  = " Algorithm: Unknown!."
 end select
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(4(a,i5,a),a,i5)')&
&  ' Dimension of the v, W matrices,  npweps  = ',BSp%npweps,ch10,&
&  ' Cutoff for the wavefunctions,    npwwfn  = ',BSp%npwwfn,ch10,&
&  ' Highest empty band included,     nband   = ',BSp%nbnds,ch10, &
&  ' Lowest occupied band included,   loband  = ',BSp%lomo,ch10,  &
&  ' Number of k-points in the IBZ,   nkibz   = ',BSp%nkibz
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(3(a,f6.2,a),a,f6.2)')&
&  ' Minimum frequency [eV]           Emin    = ',BSp%omegai*Ha_eV,ch10,&
&  ' Maximum frequency [eV]           Emax    = ',BSp%omegae*Ha_eV,ch10,&
&  ' Frequency step [eV]              dE      = ',BSp%domega*Ha_eV,ch10,&
&  ' Lorentzian broadening [eV]       eta     = ',BSp%broad*Ha_eV
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(5(a,i5,a),a,i5)')&
&  ' Lowest occupied state                   ',BSp%lomo,ch10,&
&  ' Highest occupied state                  ',BSp%homo,ch10,&
&  ' Lowest unoccupied state                 ',BSp%lumo,ch10,&
&  ' Highest unoccupied state                ',BSp%nbnds,ch10,&
&  ' Number of valence bands                 ',BSp%nbndv,ch10,&
&  ' Number of conduction bands              ',BSp%nbndc
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(a,i0)')&
&  " Number of resonant transitions   nreh(1) =  ",BSp%nreh(1)
 call wrtout(my_unt,msg,my_mode)
 if (BSp%nsppol==2) then
   write(msg,'(a,i0)')&
&    " Number of resonant transitions   nreh(2) =  ",BSp%nreh(2)
   call wrtout(my_unt,msg,my_mode)
 end if
 !
 ! Calculation type
 call bsp_calctype2str(Bsp,msg)
 call wrtout(my_unt,msg,my_mode)

 if (ABS(Bsp%soenergy)>tol6) then
   write(msg,'(a,f5.2)')&
&   " Scissors operator energy [eV] =         ",Bsp%soenergy*Ha_eV
   call wrtout(my_unt,msg,my_mode)
 end if

 msg=' Local fields effects (v term) excluded'
 if (BSp%exchange_term>0) msg=' Local fields effects (v term) included'
 call wrtout(my_unt,msg,my_mode)

 msg=' Excitonic effects (W term) excluded'
 if (BSp%use_coulomb_term) msg=' Excitonic effects (W term) included'
 call wrtout(my_unt,msg,my_mode)

 if (BSp%use_coulomb_term) then 
   msg=" Full W_GG' included"
   if (BSp%use_diagonal_Wgg) msg=' Only diagonal term W_GG included'
   call wrtout(my_unt,msg,my_mode)
   if (BSp%wtype==BSE_WTYPE_FROM_SCR) call wrtout(my_unt," W is read from an external SCR file",my_mode)
   if (BSp%wtype==BSE_WTYPE_FROM_MDL) call wrtout(my_unt," W is approximated with the model dielectric function",my_mode)
 end if

 msg=' Resonant-only calculation (Hermitian case)'
 if (BSp%use_coupling>0) msg=' Resonant + Coupling calculation'
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(a)')ch10
 call wrtout(my_unt,msg,my_mode)

 call wrtout(my_unt,' Calculating epsilon_Macro(q-->0,w), along the following directions:',my_mode)
 do iq=1,BSp%nq
   write(msg,'(a,3f10.6,2a)')' q = (',(BSp%q(ii,iq),ii=1,3),  ') [r.l.u.]'
   call wrtout(my_unt,msg,my_mode)
 end do

 !TODO 
 !Add file sizes and size of the buffer used for the matrix.

end subroutine print_bs_parameters
!!***

!----------------------------------------------------------------------


!!****f* m_bs_defs/bsp_calctype2str
!! NAME
!!  bsp_calctype2str
!!
!! FUNCTION
!!  Returns a string with the calculation type.
!!
!! PARENTS
!!      exc_spectra,m_bs_defs,setup_bse
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine bsp_calctype2str(BSp,str)
    
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bsp_calctype2str'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=500),intent(out) :: str
 type(excparam),intent(in) :: BSp

!************************************************************************

 SELECT CASE (Bsp%calc_type)
 CASE (BSE_HTYPE_RPA_KS)
   str = " RPA L0 with KS energies and KS wavefunctions"
 CASE (BSE_HTYPE_RPA_QPENE)
   str = " RPA L0 with QP energies and KS wavefunctions"
 CASE (BSE_HTYPE_RPA_QP)
   str = " RPA L0 with QP energies and QP wavefunctions"
 CASE DEFAULT
   str = " Unknown"
 END SELECT

end subroutine bsp_calctype2str
!!***

!----------------------------------------------------------------------

!!****f* m_bs_defs/init_transitions
!! NAME
!!  init_transitions
!!
!! FUNCTION
!!  Main creation method for the transition structured datatype.
!!
!! INPUTS
!!  lomo
!!  humo
!!  ir_cut,uv_cut
!!  nkbz
!!  nbnds
!!  nkibz
!!  nsppol
!!  nspinor
!!  gw_energy
!!  occ
!!  ktab
!!
!! OUTPUT
!!  nreh(nsppol)=Number of resonant transitions for each spin.
!!
!! SIDE EFFECTS
!!  Trans(:,:)
!!    input:  pointer to null()
!!    output: Trans(max_nreh,nsppol) stores the correspondence t -> (band,kbz,spin) and the transition energy.
!!
!! PARENTS
!!      m_shexc,setup_bse
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine init_transitions(Trans,lomo,humo,ir_cut,uv_cut,nkbz,nbnds,nkibz,nsppol,nspinor,gw_energy,occ,ktab,nreh) 
                                                                        
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_transitions'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lomo,humo
 integer,intent(in) :: nkbz,nbnds,nkibz,nsppol,nspinor
 real(dp),intent(in) :: ir_cut,uv_cut
 type(transition),pointer :: Trans(:,:)
!arrays
 integer,intent(in) :: ktab(nkbz) 
 integer,intent(out) :: nreh(nsppol)
 real(dp),intent(in) :: occ(nbnds,nkibz,nsppol)
 complex(dpc),intent(in) :: gw_energy(nbnds,nkibz,nsppol)

!Local variables ------------------------------
!scalars
 integer :: spin,it,ik_bz,ik_ibz,iv,ic,istat,max_occ,sweep,max_nreh
 real(dp) :: ene_t, delta_f
 complex(dpc) :: cplx_enet
 logical :: add_transition

!************************************************************************
      
 ! Find transitions                                                  
 max_occ=two/(nsppol*nspinor)
 nreh=zero
 !
 ! sweep=1 calculats the number of resonants transitions taking into 
 !         account a possible energy cutoff.
 ! sweep=2 initializes the tables describing the e-h transition.
 !
 do sweep=1,2
   !
   if (sweep==2) then
     max_nreh = MAXVAL(nreh)
     ABI_ALLOCATE(Trans,(max_nreh,nsppol))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0," out of memory in Trans ")
   end if
   !
   do spin=1,nsppol
     it=0
     do ik_bz=1,nkbz 
       ik_ibz=ktab(ik_bz) 
       !
       do iv=lomo,humo
         do ic=lomo,humo
           delta_f   = ( occ(ic,ik_ibz,spin)-occ(iv,ik_ibz,spin) ) / max_occ
           cplx_enet = gw_energy(ic,ik_ibz,spin)-gw_energy(iv,ik_ibz,spin)
           ene_t = DBLE(cplx_enet)

           add_transition =                       &
&             (ene_t > tol12) .and.               &  ! Resonant transition.
&             ( ABS(delta_f) > tol12) .and.       &  ! c-v transition.
&             (ene_t < uv_cut .and. ene_t > ir_cut)  ! Energy cutoff.

           if (add_transition) it = it + 1 
           if (add_transition.and.sweep==2) then 
             Trans(it,spin)%k  = ik_bz 
             Trans(it,spin)%v  = iv 
             Trans(it,spin)%c  = ic 
             Trans(it,spin)%en = cplx_enet
           end if

         end do 
       end do 
     end do ! ik_bz
     if (sweep==1) nreh(spin) = it 
   end do ! spin
   !
 end do ! sweep

end subroutine init_transitions                                           
!!***

!----------------------------------------------------------------------

!!****f* m_bs_defs/repr_1trans
!! NAME
!!  repr_1trans
!!
!! FUNCTION
!!  Returns a string with info on the (k,v,c) transition.
!!
!! INPUTS
!!  Trans<transition>=structure datatype containing indececes and info on the optical transition.
!!  [prtvol]=Verbosity level. Defaults to 0.
!!
!! OUTPUT
!!  str(len=500)=The string representing the transition.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function repr_1trans(Trans,prtvol) result(str)
 
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'repr_1trans'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol
 character(len=500) :: str
 type(transition),intent(in) :: Trans 

!Local variables ------------------------------
!scalars
 integer :: my_prtvol

!************************************************************************

 my_prtvol=0; if (PRESENT(prtvol)) my_prtvol=prtvol

 if (my_prtvol==0) then
   write(str,'(3(a,i3))')" k= ",Trans%k," v= ",Trans%v," c= ",Trans%c
 else 
   write(str,'(3(a,i3),a,2f6.2)')" k= ",Trans%k," v= ",Trans%v," c= ",Trans%c," ene= ",Trans%en*Ha_eV
 end if

end function repr_1trans
!!***

!----------------------------------------------------------------------

!!****f* m_bs_defs/repr_2trans
!! NAME
!!  repr_2trans
!!
!! FUNCTION
!!  Returns a string with info on two transitions
!!
!! INPUTS
!!  Trans1, Trans2<transition>=structure datatypes containing indececes and info on the optical transitions
!!  [prtvol]=Verbosity level. Defaults to 0.
!!
!! OUTPUT
!!  str(len=500)=The string representing the transition.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function repr_2trans(Trans1,Trans2,prtvol) result(str)
 
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'repr_2trans'
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol
 character(len=500) :: str
 type(transition),intent(in) :: Trans1,Trans2

!Local variables ------------------------------
!scalars
 integer :: my_prtvol

!************************************************************************

 my_prtvol=0; if (PRESENT(prtvol)) my_prtvol=prtvol

 str = repr_1trans(Trans1,my_prtvol)

 str = TRIM(str)//"|"//TRIM(repr_1trans(Trans2,my_prtvol))

end function repr_2trans
!!***

!----------------------------------------------------------------------

!!****f* m_bs_defs/print_bs_files
!! NAME
!!  print_bs_files
!!
!! FUNCTION
!!  Printout of the content of the excfiles structure.
!!
!! INPUTS
!!  BS_files<excfiles>=An object of type excfile storing the filenames used in the Bethe-Salpeter code.
!!  [unit]=Unit number for output
!!  [prtvol]=Verbosity level
!!  [mode_paral]=Either "COLL" or "PERS"
!!  [header]=String to be printed as header for additional info.
!!
!! OUTPUT
!!  Only printing.
!!
!! PARENTS
!!      setup_bse
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_bs_files(BS_files,header,unit,mode_paral,prtvol)                                   
 
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_bs_files'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none 

!Arguments ------------------------------------
!scalars
 type(excfiles),intent(in) :: BS_files
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral 
 character(len=*),optional,intent(in) :: header
!arrays

!Local variables ------------------------------
!scalars
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg      
! ********************************************************************* 

 !@excfiles
 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol 
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral
                                                                    
 msg=' ==== Files used for the Bethe-Salpeter calculation  ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 if (BS_files%in_hreso /= BSE_NOFILE) then
   msg = " Resonant block will be read from: "//TRIM(BS_files%in_hreso)
   call wrtout(my_unt,msg,my_mode)
 end if

 if (BS_files%in_hcoup /= BSE_NOFILE) then
   msg = " Coupling block will be read from: "//TRIM(BS_files%in_hcoup)
   call wrtout(my_unt,msg,my_mode)
 end if

 if (BS_files%in_eig /= BSE_NOFILE) then
   msg = " BS eigenstates will be read from: "//TRIM(BS_files%in_eig)
   call wrtout(my_unt,msg,my_mode)
 end if

 if (BS_files%in_haydock_basename /= BSE_NOFILE) then
   msg = " Haydock restart files have basename : "//TRIM(BS_files%in_haydock_basename)
   call wrtout(my_unt,msg,my_mode)
 end if

end subroutine print_bs_files
!!***

!----------------------------------------------------------------------

END MODULE m_bs_defs
!!***
