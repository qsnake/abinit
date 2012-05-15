#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"  

module m_gamma

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon
 use m_kptrank
 use m_io_gkk
 use m_errors
                                                                       
 !use m_kptrank,  only : kptrank_type, nullify_kptrank, destroy_kptrank
 use m_fstrings,       only : toupper
 use m_numeric_tools,  only : arth
 use m_io_tools,       only : get_unit
 use m_crystal,        only : crystal_structure
 use m_bz_mesh,        only : isamek, make_path
 use m_header,         only : hdr_clean

 implicit none

 private 
!!***

!----------------------------------------------------------------------

!!****t* m_gamma/gamma_t
!! NAME
!! gamma_t
!!
!! FUNCTION
!! Contains the phonon linewidths and provides methods for performing the interpolation in q-space.
!!
!! SOURCE

 type,public :: gamma_t

  integer :: natom
  ! Number of atoms per unit cell.
  
  integer :: nbranch
  ! Number of phonon branches i.e. 3*natom.
  
  integer :: nsppol
  ! Number of independent spin polarizations.
  
  integer :: nqibz
  ! Number of q-points in the IBZ.
  
  integer :: nqbz
  ! Number of q-points in the BZ.
  
  integer :: ep_scalprod
  !
  integer :: nrpt
  ! Number of points in the real space representation of the gamma matrices.

  integer :: qptrlatt(3,3)

  integer,pointer :: qirredtofull(:)     SET2NULL
  ! qirredtofull(nqibz)
  ! Mapping qibz ==> qbz

  integer,pointer :: qpttoqpt(:,:,:)      SET2NULL
  ! qpttoqpt(2,Cryst%nsym,nqbz)

  real(dp) :: gprim(3,3)
  ! Needed for Fourier interpolation.
  ! NOTE: gprim (not gprimd) is used for all FT interpolations,
  ! to be consistent with the dimensions of the rpt, which come from anaddb.

  real(dp),pointer :: qibz(:,:)  SET2NULL
  ! qibz(3,nqibz)
  ! Reduced coordinates of the q-points in the IBZ.

  real(dp),pointer :: qbz(:,:)   SET2NULL
  ! qbz(3,nqbz)
  ! Reduced coordinates of the q-points in the BZ.

  real(dp),pointer :: rpt(:,:)   SET2NULL
  ! rpt(3,nrpt)
  !  Reduced coordinates ***in terms of rprim*** of the lattice points used 
  !  for the Fourier transform of the phonon linewidths.

  real(dp),pointer :: wghatm(:,:,:) SET2NULL
  ! wghatm(natom,natom,nrpt)
  ! Weights used in the FT of the phonon linewidths.

  real(dp),pointer :: matij_qibz(:,:,:,:,:)  SET2NULL
  ! matij_qibz(2,nbranch,nbranch,nqibz,nsppol)) in reduced coordinates for each q-point in the IBZ.
  ! matii_qibz {\tau'\alpha',\tau\alpha} =
  !   <psi_{k+q,ib2} | H(1)_{\tau'\alpha'} | psi_{k,ib1}>*  \cdot
  !   <psi_{k+q,ib2} | H(1)_{\tau \alpha } | psi_{k,ib1}>

  !NOTE: choice to put nsppol before or after nqbz is a bit arbitrary
  !   abinit uses nband,nkpt,nsppol, but here for convenience nkpt_phon,nsppol,nqbz
  !   as interpolation is on qpt
  !
  !MG: I think that nsppol should be the last dimension set call to ftgam.

  real(dp),pointer :: gamma_qpt(:,:,:,:)  SET2NULL
  ! gamma_qpt(2,nbranch**2,nsppol,nqbz)
  ! gamma matrices integrated over kpoints coeff and bands: still depends on qpt and spin.
                                                                                
  real(dp),pointer :: gamma_rpt(:,:,:,:)  SET2NULL
  ! gamma_rpt(2,nbranch**2,nsppol,nrpt)
  ! gamma matrices in real space.

 end type gamma_t

 public :: gamma_nullify        ! Initialized all pointers to null().
 public :: gamma_free           ! Destructor method.
 public :: gamma_init           ! Creation method.
 public :: gamma_interp         ! Interpolates the phonon linewidths.
 public :: gamma_interp_setup   ! Precalculates the internal tables used for the interpolation in q-space.
 public :: gamma_linwid         ! Calculates linewidths along a given q-path.

!!***

!----------------------------------------------------------------------

!!****t* m_gamma/a2f_t
!! NAME
!! a2f_t
!!
!! FUNCTION
!! Object storing the Eliashberg function a2F(w).
!!
!! SOURCE

 type,public :: a2f_t

  integer :: nomega
  ! Number of frequency points in a2f(w).

  integer :: nsppol
  ! Number of independent spin polarizations.

  real(dp) :: omega_min,omega_max
  ! min and Max frequency.

  !?? real(dp) :: omega_step
  ! Step of the linear mesh 

  real(dp) :: a2fsmear
  ! Gaussian broadening used to approximated the Dirac distribution.

  integer :: nqshift
  ! Number of shifts in the q-mesh used for calculating a2f(w).

  integer :: qptrlatt(3,3)
  ! The q-mesh used for calculating a2f(w).

  integer,pointer :: qshift(:,:)  SET2NULL
  ! qshift(3,nqshift)
  ! The shifts used to generate the q-mesh.

  real(dp),pointer :: n0(:)   SET2NULL
  ! n0(nsppol)
  ! Electronic DOS at the Fermi level.

  real(dp),pointer :: omega(:)  SET2NULL
  ! omega(nomega)
  ! Frequency mesh (linear).

  real(dp),pointer :: a2f(:,:)  SET2NULL
  ! a2f(nomega,nsppol)
  ! Eliashberg function.
 end type a2f_t
 !
 public :: a2f_nullify   ! Set all pointers to null().
 public :: a2f_free      ! Free the memory allocated in the structure.
 public :: a2f_init      ! Calculates the FS averaged alpha^2F(w) function.
 !public :: a2f_moment    ! Returns moments of alpha^2F(w).
 !public :: a2f_logmoment ! Returns log moment of alpha^2F(w).
 public :: a2f_dump      ! Dumps alpha^2F(w) on an external file.
!!***

CONTAINS  !=========================================================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_gamma/gamma_nullify
!! NAME
!! gamma_nullify
!!
!! FUNCTION
!!  Set all pointers to null()
!!
!! SIDE EFFECTS
!!  Gam<gamma_t>=All pointers are set to null.
!!
!! PARENTS
!!      elphon,m_gamma
!!
!! CHILDREN
!!
!! SOURCE

subroutine gamma_nullify(Gam)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gamma_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(gamma_t),intent(inout) :: Gam

! *************************************************************************
 
 !@gamma_t
 !
 !integer
 nullify(Gam%qirredtofull)
 nullify(Gam%qpttoqpt)
 !
 !real
 nullify(Gam%qibz)
 nullify(Gam%qbz )

 nullify(Gam%rpt   )
 nullify(Gam%wghatm)

 nullify(Gam%matij_qibz)
 nullify(Gam%gamma_qpt)
 nullify(Gam%gamma_rpt)

end subroutine gamma_nullify

!!***

!----------------------------------------------------------------------

!!****f* m_gamma/gamma_free
!! NAME
!! gamma_free
!!
!! FUNCTION
!!  Free the dynamic memory.
!!
!! SIDE EFFECTS
!!  Gam<gamma_t>=All associated pointers are deallocated.
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine gamma_free(Gam)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gamma_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(gamma_t),intent(inout) :: Gam

! *************************************************************************
 
 !@gamma_t
 !
 !integer 
 if (associated(Gam%qirredtofull))   then
   ABI_DEALLOCATE(Gam%qirredtofull)
 end if
 if (associated(Gam%qpttoqpt    ))   then
   ABI_DEALLOCATE(Gam%qpttoqpt)
 end if
 !
 !real
 if (associated(Gam%qibz))   then
   ABI_DEALLOCATE(Gam%qibz)
 end if
 if (associated(Gam%qbz ))   then
   ABI_DEALLOCATE(Gam%qbz)
 end if

 if (associated(Gam%rpt   ))   then
   ABI_DEALLOCATE(Gam%rpt)
 end if
 if (associated(Gam%wghatm))   then
   ABI_DEALLOCATE(Gam%wghatm)
 end if

 if (associated(Gam%matij_qibz))   then
   ABI_DEALLOCATE(Gam%matij_qibz)
 end if
 if (associated(Gam%gamma_qpt ))   then
   ABI_DEALLOCATE(Gam%gamma_qpt)
 end if
 if (associated(Gam%gamma_rpt ))   then
   ABI_DEALLOCATE(Gam%gamma_rpt)
 end if

end subroutine gamma_free

!!***

!----------------------------------------------------------------------

!!****f* m_gamma/gamma_init
!! NAME
!! gamma_init
!!
!! FUNCTION
!!  Creation method for the gamma_t datatype.
!!
!! INPUTS
!! nrpt
!! gkk_fname
!! Cryst<crystal_structure>
!! Bst<bandstructure_type>
!! elph_ds<elph_type>
!! phon_ds<phon_type>
!! FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
!! qptrlatt(3,3)
!! gprim(3,3)
!! rpt(3,nrpt)
!! wghatm(Cryst%natom,Cryst%natom,nrpt)
!!
!! OUTPUT
!! Gam<gamma_t>
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine gamma_init(Gam,gkk_fname,elph_ds,Cryst,Bst,Phon_ds,FSfullpqtofull,gprim,qptrlatt,nrpt,rpt,wghatm)

 use defs_basis
 use m_header

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gamma_init'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nrpt
 character(len=fnlen),intent(in) :: gkk_fname
 type(crystal_structure),intent(in) :: Cryst
 type(bandstructure_type),intent(in) :: Bst
 type(elph_type),intent(in) :: elph_ds
 type(phon_type),intent(inout) :: phon_ds
 type(gamma_t),intent(out) :: Gam
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 integer,intent(in) :: qptrlatt(3,3)
 real(dp),intent(in) :: gprim(3,3)
 real(dp),intent(in) :: rpt(3,nrpt)
 real(dp),intent(in) :: wghatm(Cryst%natom,Cryst%natom,nrpt)

!Local variables-------------------------------
!scalars
 integer :: nsppol,nbranch,nFSband,minFSband,ikpt_phon,nqibz,nqbz 
 integer :: ib1,ib2,ii,ipert1,iqptfull,spin,qtimrev
 integer :: iq_ibz,iq_bz,nkwant,pert,pertcase,iqpt_fullbz,ikpt_phonq,ibeff,istat
 integer :: isym,symrankkpt_phon
 real(dp) :: sd1,sd2,lambda_tot !e1mef,e2mef,
 character(len=500) :: msg
 logical :: found
 type(hdr_type) :: GS_Hdr
 type(gkkfd_t) :: Fd
 type(kptrank_type) :: Qrank
!arrays
 integer :: symq(4,2,Cryst%nsym),g0(3)
 integer,allocatable :: gkk_flag(:,:,:,:)
 real(dp) :: qirr(3),qfull(3),tmp_kpt(3)
 real(dp) :: lambda(elph_ds%nsppol)
 real(dp) :: displ_cart(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: displ_red(2,3*Cryst%natom,3*Cryst%natom),eigval(3*Cryst%natom)
 real(dp) :: eigvec(2,3*Cryst%natom,3*Cryst%natom),phfrq(3*Cryst%natom)
 real(dp),allocatable :: buf_h1(:,:,:,:)
 real(dp),pointer :: kwanted(:,:),qibz(:,:),qbz(:,:)
 real(dp),allocatable :: h1_mat_el(:,:,:,:,:),h1_mat_el_sq(:,:,:,:,:)
 real(dp),allocatable :: accum_mat(:,:,:,:)
 real(dp),allocatable :: accum_mat2(:,:,:,:)
 real(dp),allocatable :: gkq_sum_bands(:,:,:)
 real(dp),allocatable :: gam_now(:,:,:) 

! *************************************************************************

 MSG_WARNING("Entering NEW_GKK section")

 ABI_UNUSED(Bst%bantot)

 !@gamma_t
 call gamma_nullify(Gam)

 call hdr_nullify(GS_Hdr)

 nsppol    = elph_ds%nsppol
 !nbranch   = elph_ds%nbranch
 nbranch   =3*Cryst%natom
 nFSband   = elph_ds%nFSband
 minFSband = elph_ds%minFSband

 nqibz     = elph_ds%nqptirred
 nqbz      = elph_ds%nqpt_full

 nkwant  =  elph_ds%k_phon%nkpt
 kwanted => elph_ds%k_phon%kpt

 call gkkfd_init(Fd,gkk_fname,GS_Hdr)

 ABI_CHECK(nqibz==Fd%nqibz,"nqibz/=Fd%nqibz")

 do iq_ibz=1,Fd%nqibz
   iqpt_fullbz = elph_ds%qirredtofull(iq_ibz)
   if (.not.isamek(elph_ds%qpt_full(:,iqpt_fullbz),Fd%qibz(:,iq_ibz),g0)) then
     write(std_out,*)elph_ds%qpt_full(:,iqpt_fullbz),Fd%qibz(:,iq_ibz)
     MSG_BUG("qpt_full(:,iqpt_fullbz) /= Fd%qibz(:,iq_ibz)")
   end if
 end do

 qibz => Fd%qibz
 qbz  => elph_ds%qpt_full
 !
 ! Set basic dimensions.
 Gam%natom       = Cryst%natom
 Gam%nbranch     = nbranch
 Gam%nsppol      = nsppol
 Gam%nqibz       = nqibz
 Gam%nqbz        = nqbz
 Gam%ep_scalprod = elph_ds%ep_scalprod
 Gam%nrpt        = nrpt

 Gam%qptrlatt    = qptrlatt
 Gam%gprim       = gprim
 !
 ABI_ALLOCATE(Gam%qibz,(3,nqibz))
 Gam%qibz = qibz
 ABI_ALLOCATE(Gam%qbz ,(3,nqbz))
 Gam%qbz  = qbz

 ABI_ALLOCATE(Gam%rpt,(3,nrpt))
 Gam%rpt = rpt
 ABI_ALLOCATE(Gam%wghatm,(Gam%natom,Gam%natom,nrpt))
 Gam%wghatm = wghatm

 ! TODO
 ABI_ALLOCATE(Gam%qirredtofull,(nqbz))
 Gam%qirredtofull=0

 do iq_ibz=1,nqibz
   qirr = qibz(:,iq_ibz)
   found = .FALSE.
   do iq_bz=1,nqbz
     qfull = qbz(:,iq_bz)
     if ( isamek(qirr,qfull,g0) ) then
       found = .TRUE.
       Gam%qirredtofull(iq_ibz) = iq_bz
       EXIT
     end if
   end do
   ABI_CHECK(found,"qirr not in BZ!")
 end do
 !
 ! TODO See mkqptequiv
 ABI_ALLOCATE(Gam%qpttoqpt,(2,Cryst%nsym,nqbz))
 Gam%qpttoqpt = -1

 call mkkptrank(qbz,nqbz,Qrank)

 do iq_bz=1,nqbz
   do isym=1,Cryst%nsym
     tmp_kpt(:) =  Cryst%symrec(:,1,isym)*qbz(1,iq_bz) &
&                + Cryst%symrec(:,2,isym)*qbz(2,iq_bz) &
&                + Cryst%symrec(:,3,isym)*qbz(3,iq_bz)

     call get_rank_1kpt(tmp_kpt,symrankkpt_phon,Qrank)
     if (Qrank%invrank(symrankkpt_phon) == -1) then
       msg = "looks like no kpoint equiv to q by symmetry without time reversal!!!"
       MSG_ERROR(msg)
     end if
     Gam%qpttoqpt(1,isym,Qrank%invrank(symrankkpt_phon)) = iq_bz

     tmp_kpt = -tmp_kpt
     call get_rank_1kpt (tmp_kpt,symrankkpt_phon,Qrank)
     if (Qrank%invrank(symrankkpt_phon) == -1) then
       msg = ' mkqptequiv : Error : looks like no kpoint equiv to q by symmetry with time reversal!!!'
       MSG_ERROR(msg)
     end if
     Gam%qpttoqpt(2,isym,Qrank%invrank(symrankkpt_phon)) = iq_bz
   end do
 end do

 call destroy_kptrank(Qrank)
 !
 ! Allocate matrices in the IBZ.
 ABI_ALLOCATE(Gam%matij_qibz,(2,nbranch,nbranch,nqibz,nsppol))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in %matij_qibz")
 Gam%matij_qibz = zero

 ABI_ALLOCATE(gkk_flag,(nbranch,nbranch,nkwant,nsppol))

 ABI_ALLOCATE(h1_mat_el   ,(2, nFSband**2, nbranch,   nkwant, nsppol))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory h1_mat_el")

 ABI_ALLOCATE(h1_mat_el_sq,(2, nFSband**2, nbranch**2,nkwant, nsppol))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory h1_mat_el_sq")

 do iq_ibz=1,Fd%nqibz
   iqptfull = elph_ds%qirredtofull(iq_ibz)

   write (msg,'(a,i5,a,3es16.8)')&
&   ' NEW_gkk : full zone qpt number ',iqptfull,' is ',elph_ds%qpt_full(:,iqptfull)
   call wrtout(std_out,msg,'COLL')

   h1_mat_el    = HUGE(zero)
   h1_mat_el_sq = HUGE(zero)
   !h1_mat_el    = zero
   !h1_mat_el_sq = zero
   gkk_flag     = -1
   !
   ! Examine the symmetries of the q wavevector
   ! these will be used to complete the perturbations for other atoms and idir
   call symq3(Cryst%nsym,Fd%qibz(:,iq_ibz),symq,Cryst%symrec,qtimrev,prtvol=0)
   !
   ! Determine dynamical matrix, phonon frequencies and displacement vector for qpoint
   call inpphon(displ_cart,eigval,eigvec,phfrq,phon_ds,Fd%qibz(:,iq_ibz))
   !
   ! Displacement vectors in reduced coordinates
   call phdispl_cart2red(Cryst%natom,Cryst%gprimd,displ_cart,displ_red)
   !
   do pert=1,Fd%npert_q(iq_ibz)
     pertcase = Fd%pertcase_q(pert,iq_ibz)
     ABI_ALLOCATE(buf_h1,(2,nFSband**2,nkwant,nsppol))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,"out of memory buf_h1")      

     call gkkfd_read_h1me(Fd,iq_ibz,pertcase,nFSband,minFSband,nkwant,kwanted,Cryst,buf_h1)

     !do ii=1,nkwant 
     !  do ibb=1,nFSband**2
     !    write(445,*)pertcase,ibb,ii,buf_h1(:,ibb,ii,1)
     !  end do
     !end do

     gkk_flag(pertcase,pertcase,:,:) = 3 ! read from GKK file.

     ! Save values
     do spin=1,nsppol
       h1_mat_el(:,:,pertcase,:,spin) = buf_h1(:,:,:,spin)
     end do

     ABI_DEALLOCATE(buf_h1)
   end do
   !
   ! ========================================================================
   !  Now use more general symops to complete the other equivalent
   !  perturbations: the kpoints are also shuffled by these symops
   !  afterwards h1_mat_el_sq contains gamma_\tau\alpha,\tau'\alpha' in reduced coordinates
   !    
   !  gamma_{\tau'\alpha',\tau\alpha} =
   !    <psi_{k+q,ib2} | H(1)_{\tau'\alpha'} | psi_{k,ib1}>*  \cdot
   !    <psi_{k+q,ib2} | H(1)_{\tau \alpha } | psi_{k,ib1}>
   ! ========================================================================
   !
   !write(555,*)h1_mat_el

   call completeperts(Cryst,nbranch,nFSband,nkwant,nsppol,&
&    gkk_flag,h1_mat_el,h1_mat_el_sq,Fd%qibz(:,iq_ibz),symq,qtimrev)

   !write(777,*)h1_mat_el_sq

   ABI_ALLOCATE(accum_mat ,(2,nbranch,nbranch,nsppol))
   accum_mat =zero
   ABI_ALLOCATE(accum_mat2,(2,nbranch,nbranch,nsppol))
   accum_mat2=zero

if (.TRUE.) then
   call nmsq_pure_gkk_sumFS(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,h1_mat_el_sq,iq_ibz)
else
   iqpt_fullbz = elph_ds%qirredtofull(iq_ibz)
   ABI_ALLOCATE(gkq_sum_bands,(2,nbranch,nbranch))
   !
   !accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
   !accum_mat2 is used to store the phonon-linewidhts before interpolation
   ! FIXME: 
   !  BE careful here, since it wont' work if kwanted /= elph_ds%k_phon%kpt
   ! due to the tables and the way used to access the eigevalues.
   do spin=1,nsppol
     do ikpt_phon=1,nkwant
       !
       ! The index of k+q in the BZ.
       ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)

       !% ikptgs = irredtoGS_phon(elph_ds%k_phon%full2irr(1,ikpt_phon))
       !
       ! gkq_sum_bands = 
       !   \sum_{ib1,ib2} <k+q| H^{(1)}_{q,\tau_i,\alpha_i} |k> \cdot <k| H^{(1)}_{q,\tau_j,\alpha_j}|k+q>
       !
       ! where ibranch = (\tau_i,\alpha_i) and  jbranch = (\tau_j,\alpha_j).
       gkq_sum_bands(:,:,:) = zero

       do ib1=1,nFSband
         !%% e1mef = Bst%eig(ikptgs,minFSband-1+ib1,spin) - Bst%fermie 
         !%% sd1 = gaussian(e2mef,sigma)
         sd1 = elph_ds%k_phon%wtk(ib1,ikpt_phon,spin)      !  weights for distance from the fermi surface

         do ib2=1,nFSband
           !%% e2mef = Bst%eig(ikptgs,minFSband-1+ib2,spin) - Bst%fermie 
           !%% sd2 = gaussian(e2mef,sigma)
           sd2 = elph_ds%k_phon%wtk(ib2,ikpt_phonq,spin)  !  weights for distance from the fermi surface
           ibeff=ib2+(ib1-1)*nFSband

           gkq_sum_bands = gkq_sum_bands + &
&            sd1*sd2*pi * reshape(h1_mat_el_sq(:,ibeff,:,ikpt_phon,spin),(/2,nbranch,nbranch/))
         end do !ib2
       end do !ib1
       !
       ! gamma matrix contribution in reduced coordinates (ie interpolatable form)
       ! The sum over Fermi surface bands is done here, and fed into (ib1,ib2)=(1,1)
       h1_mat_el_sq(:,1,:,ikpt_phon,spin) = reshape(gkq_sum_bands,(/2,nbranch**2/))

       accum_mat(:,:,:,spin) = accum_mat(:,:,:,spin) + gkq_sum_bands(:,:,:)
     end do ! kpt_phon
   end do ! spin

   ABI_DEALLOCATE(gkq_sum_bands)
   !
   !MG20060603
   ! scalar product wit displ_red to calculate the ph lwdth before interpolation (stored in accum_mat2)
   ABI_ALLOCATE(gam_now,(2,nbranch,nbranch))
   !allocate(zgemm_tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch))

   do spin=1,nsppol
     !zgemm_tmp_mat = accum_mat(:,:,:,spin)
     !
     !call gam_mult_displ(nbranch, displ_red, zgemm_tmp_mat, gam_now)
     call gam_mult_displ(nbranch, displ_red, accum_mat(:,:,:,spin), gam_now)

     do ipert1=1,nbranch
       accum_mat2(1,ipert1,ipert1,spin) = accum_mat2(1,ipert1,ipert1,spin) + gam_now(1,ipert1,ipert1)
     end do
     !
   end do

   ABI_DEALLOCATE(gam_now)
   !deallocate(zgemm_tmp_mat)
endif

   !MG: values without the good prefactor
   accum_mat = accum_mat * elph_ds%occ_factor/nkwant 
   !
   ! Save results in Gam%
   do spin=1,nsppol
     Gam%matij_qibz(:,:,:,iq_ibz,spin) = accum_mat(:,:,:,spin)
   end do

   !MG: accum_mat2 contains the line-widhts before the Fourier interpolation
   accum_mat2 = accum_mat2 * elph_ds%occ_factor/nkwant 

   !MG20060531i
   !write e-ph quantities before Fourier interpolation
   !save e-ph values in the temporary array qdata that will be copied into elph_ds%qgrid_data

    write (msg,'(4a,3es16.6,63a)')ch10,                    &
&    ' NEW_GKK : Phonon linewidths before interpolation ',ch10,&
&    ' Q point = ',Fd%qibz(:,iq_ibz),ch10,('=',ii=1,60),ch10,  &
&    ' Mode          Frequency (Ha)  Linewidth (Ha)  Lambda '
    call wrtout(std_out,msg,'COLL')

    lambda_tot = zero
    do spin=1,nsppol
      do ii=1,nbranch
        lambda(spin)=zero
        if (abs(phfrq(ii)) > tol10) then ! The tolerance factor is somehow arbitrary
          lambda(spin)=accum_mat2(1,ii,ii,spin)/ (pi*elph_ds%n0(spin)*phfrq(ii)**2)
        end if
        lambda_tot=lambda_tot+lambda(spin)
        write(msg,'(i8,es20.6,2es16.6)' )ii,phfrq(ii),accum_mat2(1,ii,ii,spin),lambda(spin)
        call wrtout(std_out,msg,'COLL')
        !save values
        !qdata(ii,spin,1)=phfrq(ii)
        !qdata(ii,spin,2)=accum_mat2(1,ii,ii,spin)
        !qdata(ii,spin,3)=lambda(spin)
      end do !loop over branch
    end do !loop over sppol

    ! normalize for number of spins
    lambda_tot = lambda_tot / nsppol

    write(msg,'(61a,44x,es16.6,62a)' )('=',ii=1,60),ch10,lambda_tot,ch10,('=',ii=1,60),ch10
    call wrtout(std_out,msg,'COLL')
!   ENDMG20060531

!   immediately calculate linewidths:
    write(std_out,*) 'summed my_accum_mat = '
    write(std_out,'(3(2E18.6,1x))') accum_mat(:,:,:,1)
    write(std_out,*) 'summed my_accum_mat2 = '
    write(std_out,'(3(2E18.6,1x))')  (accum_mat2(:,ii,ii,1),ii=1,nbranch)
    write(std_out,*) 'displ_red  = '
    write(std_out,'(3(2E18.6,1x))') displ_red

   ABI_DEALLOCATE(accum_mat)
   ABI_DEALLOCATE(accum_mat2)
 end do ! iq_ibz

 ABI_DEALLOCATE(gkk_flag)
 ABI_DEALLOCATE(h1_mat_el)
 ABI_DEALLOCATE(h1_mat_el_sq)

 call gkkfd_free(Fd)
 call hdr_clean(GS_Hdr)

end subroutine gamma_init
!!***

!----------------------------------------------------------------------

!!****f* m_gamma/gamma_interp_setup
!! NAME
!! gamma_interp_setup
!!
!! FUNCTION
!!  This routines performs the (allocation|deallocation) of the internal tables
!!  used to interpolate the linewidths in q-space
!!
!! INPUTS
!!  mode = 
!!    "INIT" to allocate and compute the internal tables.
!!    "FREE" to deallocate the internal tables.
!!
!! SIDE EFFECTS 
!!  Gam<gamma_t>= Gam%gamma_qpt and Gam%gamma_rpt, depending on mode. 
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine gamma_interp_setup(Gam,Cryst,mode)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gamma_interp_setup'
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: mode
 type(gamma_t),intent(inout) :: Gam
 type(crystal_structure),intent(in) :: Cryst
!arrays

!Local variables-------------------------------
!scalars
 integer :: natom,nbranch,nsppol,ep_scalprod,iq_bz,iq_ibz
 integer :: nqibz,nqbz,istat,nrpt,spin,qtor 
 !character(len=500) :: msg
 !arrays
 real(dp),pointer :: qbz(:,:)

! *************************************************************************
 
 !@gamma_t
 natom       = Cryst%natom 
 nbranch     = Gam%nbranch
 nsppol      = Gam%nsppol
 nqibz       = Gam%nqibz
 nqbz        = Gam%nqbz
 nrpt        = Gam%nrpt
 ep_scalprod = Gam%ep_scalprod

 qbz         => Gam%qbz !elph_ds%qpt_full)

 !?? call integrate_gamma(elph_ds,FSfullpqtofull,nrpt)

 select case (toupper(mode))

 case ("INIT")
   if (.not.associated(Gam%gamma_qpt)) then
     ABI_ALLOCATE(Gam%gamma_qpt,(2,nbranch**2,nsppol,nqbz))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,'out of memory in Gam%gamma_qpt')
     Gam%gamma_qpt = zero

     do iq_ibz=1,nqibz
       iq_bz = Gam%qirredtofull(iq_ibz)
       do spin=1,nsppol
         Gam%gamma_qpt(:,:,spin,iq_bz) = RESHAPE(Gam%matij_qibz(:,:,:,iq_ibz,spin), (/2,nbranch**2/))
       end do
     end do
     !
     ! Complete the gamma_qpt in the full BZ.
     !?? if (elph_ds%symgkq ==1) then
     call complete_gamma(Cryst,nbranch,nsppol,nqibz,nqbz,ep_scalprod,Gam%qirredtofull,Gam%qpttoqpt,Gam%gamma_qpt)
     !end if
   end if
   !
   ! Now FT to real space too
   ! NOTE: gprim (not gprimd) is used for all FT interpolations,
   ! to be consistent with the dimensions of the rpt, which come from anaddb.
   ! TODO: this is needed only if FT is used, no when the linear interpolation is employed.
   if (.not.associated(Gam%gamma_rpt)) then
     ABI_ALLOCATE(Gam%gamma_rpt,(2,nbranch**2,nsppol,nrpt))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,'out of memory in Gam%gamma_rpt')
     Gam%gamma_rpt = zero

     qtor = 1 ! q --> r
     do spin=1,nsppol
       call ftgam(Gam%wghatm,Gam%gamma_qpt(:,:,spin,:),Gam%gamma_rpt(:,:,spin,:),Gam%gprim,natom,nqbz,Gam%nrpt,qtor,Gam%rpt,qbz)
     end do
   end if
   !
 CASE ("FREE")
   if (associated(Gam%gamma_qpt))   then
     ABI_DEALLOCATE(Gam%gamma_qpt)
   end if
   if (associated(Gam%gamma_rpt))   then
     ABI_DEALLOCATE(Gam%gamma_rpt)
   end if
   !
 CASE DEFAULT
   MSG_BUG("Wrong mode "//TRIM(mode))
 END SELECT

end subroutine gamma_interp_setup

!!***

!----------------------------------------------------------------------

!!****f* m_gamma/gamma_interp
!! NAME
!! gamma_interp
!!
!! FUNCTION
!!  Interpolates the linewidths at a given q-point.
!!
!! INPUTS
!!  alter_int_gam,spin
!!  Gam<gamma_t>
!!  Cryst<crystal_structure>
!!  qpt(3)
!!  displ_cart(2,3*Cryst%natom,3*Cryst%natom)
!!
!! OUTPUT
!!  eigval(Gam%nbranch)
!!
!! NOTES
!!  This routine assumes that the internal tables gamma_qpt and gamma_rpt have been already computed. 
!!  by calling gamma_interp_setup.
!!
!! PARENTS
!!      m_gamma
!!
!! CHILDREN
!!
!! SOURCE

subroutine gamma_interp(Gam,Cryst,alter_int_gam,spin,qpt,displ_cart,eigval)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gamma_interp'
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: alter_int_gam,spin
 type(gamma_t),intent(inout) :: Gam
 type(crystal_structure),intent(in) :: Cryst
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(in) :: displ_cart(2,3*Cryst%natom,3*Cryst%natom)
 real(dp),intent(out) :: eigval(Gam%nbranch)

!Local variables-------------------------------
!scalars
 integer :: natom,nbranch,nsppol,ep_scalprod
 integer :: nqibz,nqbz,nrpt,ibranch,jbranch,qtor
 real(dp) :: diagerr
 character(len=500) :: msg
 !arrays
 real(dp),pointer :: qbz(:,:)
 real(dp) :: displ_red(2,Gam%nbranch,Gam%nbranch)
 real(dp) :: gam_now(2,Gam%nbranch**2)
 real(dp) :: imeigval(Gam%nbranch)
 real(dp) :: pheigvec(2*Gam%nbranch**2) 
 real(dp) :: tmp_gam1(2,Gam%nbranch,Gam%nbranch)
 real(dp) :: tmp_gam2(2,Gam%nbranch,Gam%nbranch)

! *************************************************************************
 
 !@gamma_t
 natom       = Cryst%natom 
 nbranch     = Gam%nbranch
 nsppol      = Gam%nsppol
 nqibz       = Gam%nqibz
 nqbz        = Gam%nqbz
 nrpt        = Gam%nrpt
 ep_scalprod = Gam%ep_scalprod
 qbz         => Gam%qbz !elph_ds%qpt_full)

 ! Taken from mkph_linwid

 ! This reduced version of ftgkk supposes the kpoints have been integrated
 ! in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.
 if (alter_int_gam == 0) then
   qtor = 0
   call ftgam(Gam%wghatm,gam_now,Gam%gamma_rpt(:,:,spin,:),Gam%gprim,natom,1,Gam%nrpt,qtor,Gam%rpt,qpt)
 else if (alter_int_gam == 1) then
   call lin_interpq_gam(Gam%gamma_qpt,nbranch,nqbz,nsppol,gam_now,spin,Gam%qptrlatt,qpt)
 end if
 !
 ! If the matrices do not contain the scalar product with the displ_cart vectors yet do it now.
 eigval = zero
 SELECT CASE (ep_scalprod)
 CASE (0)
   call phdispl_cart2red(natom,Cryst%gprimd,displ_cart,displ_red)

   tmp_gam2 = reshape (gam_now, (/2,nbranch,nbranch/))
   call gam_mult_displ(nbranch, displ_red, tmp_gam2, tmp_gam1)

   do jbranch=1,nbranch
     eigval(jbranch)   = tmp_gam1(1, jbranch, jbranch)
     imeigval(jbranch) = tmp_gam1(2, jbranch, jbranch)
     !
     if (abs(imeigval(jbranch)) > tol8) then
       write (msg,'(a,i0,a,es16.8)')' non-zero imaginary part for branch= ',jbranch,', imeigval= ',imeigval(jbranch)
       MSG_WARNING(msg)
     end if
     !
   end do

 CASE (1)
   !
   ! Diagonalize gamma matrix at qpoint (complex matrix).
   ! MJV NOTE: gam_now is recast implicitly here to matrix 
   call ZGEMM('N','N',nbranch,nbranch,nbranch,cone,gam_now, nbranch,pheigvec,nbranch,czero,tmp_gam1,nbranch)

   call ZGEMM('C','N',nbranch,nbranch,nbranch,cone,pheigvec,nbranch,tmp_gam1 ,nbranch,czero,tmp_gam2,nbranch)

   diagerr = zero
   do ibranch=1,nbranch

     eigval(ibranch) = tmp_gam2(1,ibranch,ibranch)

     do jbranch=1,ibranch-1
       diagerr = diagerr + abs(tmp_gam2(1,jbranch,ibranch))+abs(tmp_gam2(2,jbranch,ibranch))
     end do
     do jbranch=ibranch+1,nbranch
       diagerr = diagerr + abs(tmp_gam2(1,jbranch,ibranch))+abs(tmp_gam2(2,jbranch,ibranch))
     end do
     diagerr = diagerr + abs(tmp_gam2(2,ibranch,ibranch))
   end do

   if (diagerr > tol12) then
     write (msg,'(a,es14.6)')' Numerical error in diagonalization of gamma with phon eigenvectors: ',diagerr
     MSG_WARNING(msg)
   end if

 CASE DEFAULT
   write (msg,'(a,i0)')' Wrong value for ep_scalprod= ',ep_scalprod
   MSG_BUG(msg)
 END SELECT 
 ! end taken from mkph_linwid

end subroutine gamma_interp

!!***

!----------------------------------------------------------------------

!!****f* m_gamma/gamma_linwid
!! NAME
!! gamma_linwid
!!
!! FUNCTION
!!  Calculates the phonon linewidths on a trajectory in q space
!!
!! INPUTS
!!  Cryst<crystal_structure>=Info on the unit cell and symmetries.
!!  alter_int_gam = (1 or 0) use of alternative method to integrate over FS
!!  elph_ds = datastructure with phonon matrix elements
!!  nqpath = dimension of qpath_vertices
!!  phon_ds = datastructure with interatomic force constants
!!  qpath_vertices = vertices of reciprocal space trajectory
!!
!! OUTPUT
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine gamma_linwid(Gam,Cryst,alter_int_gam,elph_ds,nqpath,phon_ds,qpath_vertices)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gamma_linwid'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: alter_int_gam,nqpath
 type(crystal_structure),intent(in) :: Cryst
 type(gamma_t),intent(inout) :: Gam
 type(elph_type),intent(inout) :: elph_ds
 type(phon_type),intent(inout) :: phon_ds
!arrays
 real(dp),intent(in) :: qpath_vertices(3,nqpath)

!Local variables-------------------------------
!scalars
 integer :: ibranch,natom,ii,indx,ios,ipoint,nbranch,nqbz,nsppol
 integer :: spin,jbranch,qtor,unit_bs,unit_lambda,unit_lwd,npt_tot,nrpt
 real(dp) :: diagerr,res
 character(len=500) :: msg
 character(len=fnlen) :: fname,base_name
!arrays
 integer :: ndiv(nqpath-1)
 integer, allocatable :: indxprtqpt(:)
 real(dp) :: gprim(3,3)
 real(dp) :: displ_cart(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: displ_red(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: eigval(3*Cryst%natom)
 real(dp) :: gam_now(2,(3*Cryst%natom)**2)
 real(dp) :: imeigval(3*Cryst%natom)
 real(dp) :: lambda(3*Cryst%natom),pheigval(3*Cryst%natom)
 real(dp) :: pheigvec(2*3*Cryst%natom*3*Cryst%natom),phfrq_tmp(3*Cryst%natom)
 real(dp) :: qpt(3),redkpt(3)
 real(dp) :: tmpgam1(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: tmpgam2(2,3*Cryst%natom,3*Cryst%natom)
 real(dp),pointer :: qpath(:,:),rpt(:,:),wghatm(:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 natom     = Cryst%natom
 nbranch   = elph_ds%nbranch 
 nsppol    = elph_ds%nsppol
 base_name = "NEW" !elph_ds%elph_base_name
 gprim     = Gam%gprim
 nrpt      = Gam%nrpt

 rpt       => Gam%rpt
 wghatm    => Gam%wghatm

 if (.not.associated(Gam%gamma_qpt) .or. .not.associated(Gam%gamma_rpt)) then
   MSG_ERROR("Gam%gamma_qpt or Gam%gamma_rpt are not associated")
 end if
 !
 ! Define the q-path along which ph linwid will be interpolated.
 nullify(qpath)
 call make_path(nqpath,qpath_vertices,Cryst%gmet,'G',20,ndiv,npt_tot,qpath)
 ABI_ALLOCATE(indxprtqpt,(npt_tot))
 indxprtqpt = 0
 !
 ! =======================================
 ! === Open _LWD file and write header ===
 ! =======================================
 unit_lwd=get_unit()
 fname=trim(base_name) // '_LWD'
 open(unit=unit_lwd,file=fname,status='unknown',iostat=ios)
 ABI_CHECK(ios==0,"Opening "//trim(fname))

 write(unit_lwd,'(a)')       '#'
 write(unit_lwd,'(a)')       '# ABINIT package : Phonon linewidth file'
 write(unit_lwd,'(a)')       '#'
 write(unit_lwd,'(a,i10,a)') '#  Phonon linewidths calculated on ',npt_tot,' points along the qpath'
 write(unit_lwd,'(a)')       '#  Description of the Q-path :'
 write(unit_lwd, '(a,i10)')  '#  Number of line segments = ',nqpath-1
 write(unit_lwd,'(a)')       '#  Vertices of the Q-path and corresponding index = '
 indx=1
 indxprtqpt(1) = 1
 indxprtqpt(npt_tot) = 1
 do ii=1,nqpath
   write (unit_lwd,'(a,3(e16.6,1x),i8)')'#  ',qpath_vertices(:,ii),indx
   if (ii<nqpath) then
     indx=indx+ndiv(ii)
     indxprtqpt(indx) = 1
   end if
 end do
 write (unit_lwd,'(a)')'#'
 !
 ! =======================================
 ! === Open _BST file and write header ===
 ! =======================================
 unit_bs=get_unit()
 fname=trim(base_name) // '_BST'
 open(unit=unit_bs,file=fname,status='unknown',iostat=ios)
 ABI_CHECK(ios==0,"opening "//trim(fname))

 write(unit_bs, '(a)')      '#'
 write(unit_bs, '(a)')      '# ABINIT package : Phonon band structure file'
 write(unit_bs, '(a)')      '#'
 write(unit_bs, '(a,i10,a)')'# Phonon BS calculated on ', npt_tot,' points along the qpath'
 write(unit_bs, '(a,i10)')  '# Number of line segments = ', nqpath-1
 indx=1
 do ii=1,nqpath
   write (unit_bs,'(a,3(E16.6,1x),i8)')'#  ',qpath_vertices(:,ii),indx
   if (ii<nqpath) indx=indx+ndiv(ii)
 end do
 write (unit_bs,'(a)')'#'

!MG20060606
!==========================================================
!open _LAMBDA file and write header
!contains \omega(q,n) and \lambda(q,n) and can be plotted using xmgrace
!==========================================================
 unit_lambda=get_unit()
 fname=trim(base_name) // '_LAMBDA'
 open(unit=unit_lambda,file=fname,status='unknown',iostat=ios)
 ABI_CHECK(ios==0,"opening "//trim(fname))

 write(unit_lambda,'(a)')      '#'
 write(unit_lambda,'(a)')      '# ABINIT package : Lambda file'
 write(unit_lambda,'(a)')      '#'
 write(unit_lambda,'(a,i10,a)')'#  Lambda(q,nu) calculated on ',npt_tot,' Q-points'
 write(unit_lambda,'(a)')      '# Description of the Q-path :'
 write(unit_lambda,'(a,i10)')  '# Number of line segments = ',nqpath-1
 write(unit_lambda,'(a)')      '# Vertices of the Q-path and corresponding index = '

 indx=1
 do ii=1,nqpath
   write (unit_lambda,'(a,3(E16.6,1x),i8)')'#  ',qpath_vertices(:,ii),indx
   if (ii<nqpath) indx=indx+ndiv(ii)
 end do
 write(unit_lambda,'(a)')'#'
 write(unit_lambda,'(a)')'# index frequency lambda(q,n) frequency lambda(q,n) .... lambda_tot'
 write(unit_lambda,'(a)')'#'

!initialize the maximum phonon frequency
 elph_ds%omega_min = zero
 elph_ds%omega_max = zero

 write(std_out,*) ' gamma_linwid : shape(elph_ds%gamma_qpt) = ',shape(elph_ds%gamma_qpt)
 nqbz =  SIZE(elph_ds%gamma_qpt,DIM=4)
 write(std_out,*) " nqbz =  SIZE(elph_ds%gamma_qpt,DIM=4) = ",nqbz
 !
 ! Big do loop over spin polarizations
 ! could put in locally, so phonon stuff is not done twice...
 !
 do spin=1,nsppol
   indx=1
   !
   ! Output to the main output file
   write(msg,'(a,a)')ch10,&
&    ' Output of the linewidths for the first point of each segment. Linewidths are given in Hartree.'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write(std_out,*) ' gamma_linwid : elph_ds%ep_scalprod = ', elph_ds%ep_scalprod
   !
   ! Interpolation along specified path in q space
   do ipoint=1,npt_tot
     !
     ! Get qpoint along the path from qpath_vertices
     qpt(:) = qpath(:,ipoint)

     call wrap2_pmhalf(qpt(1),redkpt(1),res)
     call wrap2_pmhalf(qpt(2),redkpt(2),res)
     call wrap2_pmhalf(qpt(3),redkpt(3),res)
     qpt(:) = redkpt(:)
     !
     ! Get phonon frequencies and eigenvectors.
     call inpphon(displ_cart,pheigval,pheigvec,phfrq_tmp,phon_ds,qpt)

     !if (ipoint==1) then 
     !  write(std_out,*)"interp for qpt ",qpt
     !  write(std_out,*)"pheigvec :",pheigvec
     !  write(std_out,*)"phfrq_tmp: ",phfrq_tmp
     !  write(std_out,*)"dipl_cart: ",displ_cart
     !end if

if (.TRUE.) then
     !if (ipoint==1) qpt(:) = tol6 !zero
     call gamma_interp(Gam,Cryst,alter_int_gam,spin,qpt,displ_cart,eigval)

!BEGIN DEBUG
     if ( ANY( ABS( Gam%gamma_qpt - elph_ds%gamma_qpt) > tol16 )) then
       MSG_ERROR("Gam%gamma_qpt /= elph_ds%gamma_qpt")
     end if

     if (ANY( ABS( Gam%gamma_rpt - elph_ds%gamma_rpt) > tol16 )) then
       MSG_ERROR("Gam%gamma_rpt /= elph_ds%gamma_rpt")
     end if
!END DEBUG

else
     !
     ! This reduced version of ftgkk supposes the kpoints have been integrated
     ! in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.
     if (alter_int_gam == 0) then
       qtor = 0 !real space to q space
       call ftgam(wghatm,gam_now,elph_ds%gamma_rpt(:,:,spin,:),gprim,natom,1,nrpt,qtor,rpt,qpt)
     else if (alter_int_gam == 1) then
       call lin_interpq_gam(elph_ds%gamma_qpt,nbranch,nqbz,nsppol,gam_now,spin,Gam%qptrlatt,qpt)
     end if
     !    
     ! If the matrices do not contain the scalar product with the displ_cart vectors yet do it now
     if (elph_ds%ep_scalprod == 0) then

       call phdispl_cart2red(natom,Cryst%gprimd,displ_cart,displ_red)

       tmpgam2 = reshape (gam_now, (/2,nbranch,nbranch/))
       call gam_mult_displ(nbranch, displ_red, tmpgam2, tmpgam1)

       do jbranch=1,nbranch
         eigval(jbranch)   = tmpgam1(1, jbranch, jbranch)
         imeigval(jbranch) = tmpgam1(2, jbranch, jbranch)

         if (abs(imeigval(jbranch)) > tol8) then
           write (msg,'(a,i0,a,es16.8)')' imaginary values for branch = ',jbranch,' imeigval = ',imeigval(jbranch)
           MSG_WARNING(msg)
         end if
       end do

     else if (elph_ds%ep_scalprod == 1) then
       !
       ! Diagonalize gamma matrix at qpoint (complex matrix).
       ! MJV NOTE: gam_now is recast implicitly here to matrix 
       call ZGEMM( 'N', 'N', 3*natom, 3*natom, 3*natom, cone, gam_now, 3*natom, pheigvec, 3*natom, czero, tmpgam1, 3*natom)

       call ZGEMM( 'C', 'N', 3*natom, 3*natom, 3*natom, cone, pheigvec, 3*natom, tmpgam1, 3*natom, czero, tmpgam2, 3*natom)

       diagerr = zero
       do ibranch=1,nbranch

         eigval(ibranch) = tmpgam2(1,ibranch,ibranch)

         do jbranch=1,ibranch-1
           diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))+abs(tmpgam2(2,jbranch,ibranch))
         end do
         do jbranch=ibranch+1,nbranch
           diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))+abs(tmpgam2(2,jbranch,ibranch))
         end do
         diagerr = diagerr + abs(tmpgam2(2,ibranch,ibranch))
       end do

       if (diagerr > tol12) then
         write (msg,'(a,es14.6)')' Numerical error in diagonalization of gamma with phon eigenvectors: ', diagerr
         MSG_WARNING(msg)
       end if

     else
       write (msg,'(a,i0)')' Wrong value for elph_ds%ep_scalprod = ',elph_ds%ep_scalprod
       MSG_BUG(msg)
     end if ! end elph_ds%ep_scalprod if
end if
     !
     ! ==========================================================
     ! write data to files for each q point
     ! ==========================================================
     write(unit_lwd,'(i5)', advance='no') indx
     write(unit_lwd,'(18E16.5)',advance='no') (eigval(ii),ii=1,nbranch)
     write(unit_lwd,*)

     ! only print phonon BS for spin 1: independent of electron spins
     if (spin==1) then
       write(unit_bs,'(i5)', advance='no') indx
       write(unit_bs,'(18E16.5)',advance='no') phfrq_tmp
       write(unit_bs,*)
     end if

     write(unit_lambda,'(i5)', advance='no') indx
     do ii=1,nbranch
       lambda(ii)=zero
       if (abs(phfrq_tmp(ii)) > tol10) lambda(ii)=eigval(ii)/(pi*elph_ds%n0(spin)*phfrq_tmp(ii)**2)
       write(unit_lambda,'(18es16.8)',advance='no')phfrq_tmp(ii),lambda(ii)
     end do
     write(unit_lambda,'(es16.8)',advance='no') sum(lambda)
     write(unit_lambda,*)
     !
     ! MG NOTE: I wrote a piece of code to output all these quantities using units
     ! chosen by the user, maybe in version 5.2?
     ! In this version the output of lambda(q,\nu) has been added
     !
     ! Output to the main output file, for first point in segment
     if (indxprtqpt(ipoint)==1) then
       write(msg,'(a,a,3es16.6,a,i4,a,a)')ch10,&
&       ' Q point =',qpt(:),'   spin = ',spin,ch10,&
&       ' Mode number    Frequency (Ha)  Linewidth (Ha)  Lambda(q,n)'
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       do ii=1,nbranch
         write(msg,'(i8,es20.6,2es16.6)' )ii,phfrq_tmp(ii),eigval(ii),lambda(ii)
         call wrtout(std_out,msg,'COLL')
         call wrtout(ab_out,msg,'COLL')
       end do
     end if
     !
     ! Find max/min phonon frequency along path chosen
     ! presumed to be representative of full BZ to within 10 percent
     elph_ds%omega_min = min(elph_ds%omega_min,1.1_dp*phfrq_tmp(1))
     elph_ds%omega_max = max(elph_ds%omega_max,1.1_dp*phfrq_tmp(nbranch))

     indx = indx+1
   end do ! end ipoint do
   !
   ! Add blank lines to output files between sppol
   write(msg,'(a)' ) ''
   call wrtout(unit_lwd,msg,'COLL')
   call wrtout(unit_lambda,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end do ! spin

 close(unit_lwd)
 close(unit_bs)
 close(unit_lambda)

 ABI_DEALLOCATE(qpath)
 ABI_DEALLOCATE(indxprtqpt)

 DBG_EXIT("COLL")

end subroutine gamma_linwid
!!***

!----------------------------------------------------------------------

!!****f* m_gamma/a2f_nullify
!! NAME
!! a2f_nullify
!!
!! FUNCTION
!!  Initialize the pointers to null()
!!
!! SIDE EFFECTS
!!  A2f<a2f_t>=Structure storing the Eliashberg function a2F.      
!!
!! PARENTS
!!      m_gamma
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2f_nullify(A2f)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'a2f_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(a2f_t),intent(inout) :: A2f

! *********************************************************************

 ! @a2f_t
 ! integer
 nullify(A2f%qshift)
 !
 ! real
 nullify(A2f%n0)
 nullify(A2f%omega)
 nullify(A2f%a2f)

end subroutine a2f_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_gamma/a2f_free
!! NAME
!! a2f_free
!!
!! FUNCTION
!!  Free the memory allocated in A2f
!!
!! SIDE EFFECTS
!!  A2f<a2f_t>=Structure storing the Eliashberg function a2F.      
!!
!! OUTPUT
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2f_free(A2f)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'a2f_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(a2f_t),intent(inout) :: A2f

! *********************************************************************

 ! @a2f_t
 ! integer
 if (associated(A2f%qshift))  then
   ABI_DEALLOCATE(A2f%qshift)
 end if
 !
 ! real
 if (associated(A2f%n0)   )   then
   ABI_DEALLOCATE(A2f%n0)
 end if
 if (associated(A2f%omega))   then
   ABI_DEALLOCATE(A2f%omega)
 end if
 if (associated(A2f%a2f)  )   then
   ABI_DEALLOCATE(A2f%a2f)
 end if

end subroutine a2f_free
!!***

!!***

!----------------------------------------------------------------------

!!****f* m_gamma/a2f_init
!! NAME
!! a2f_init
!!
!! FUNCTION
!!  Calculates the FS averaged alpha^2F function
!!
!! INPUTS
!!  Cryst<crystal_structure>=Info on the unit cell.
!!  Gam<gamma_t>=Structure storing the phonon linewidths.     
!!  Phon_ds<phon_type>= datastructure with interatomic force constants to interpolate phonons.
!!  nomega=Number of points in the linear mesh used for a2F(w)
!!  a2fsmear=Gaussian broadening used to approximate the Dirac delta.
!!  n0(Gam%nsppol)=Density of states at the Fermi level.
!!  alter_int_gam=Option for the interpolation of the phonon linewidths.
!!  qptrlatt(3,3)=Defines the Q-mesh used for interpolating the phonon linewidths (see also nqshift and qshift).
!!  nqshift=Number of shifts used to generated the Q-mesh.
!!  qshift(3,nqshift)=The shifts.
!!  [qptopt]=Controls the generation of the q-points. If not specified, the routine takes fully into account 
!!    the symmetries of the system to generate the q points in the Irreducible Brillouin Zone i.e. qptopt=1
!!    Other values of qptopt can be used for debugging purpose.
!!
!! OUTPUT
!!  A2f<a2f_t>=Structure storing the Eliashberg function a2F(w).
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2f_init(A2f,Cryst,Gam,Phon_ds,nomega,a2fsmear,n0,alter_int_gam,qptrlatt,nqshift,qshift,&
&  qptopt)  ! optional

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'a2f_init'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
 use interfaces_56_recipspace
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: alter_int_gam,nomega,nqshift
 integer,intent(in),optional :: qptopt
 real(dp),intent(in) :: a2fsmear
 type(gamma_t),intent(inout) :: Gam
 type(crystal_structure),intent(in) :: Cryst
 type(phon_type),intent(inout) :: Phon_ds
 type(a2f_t),intent(out) :: A2f
!arrays
 integer, intent(in) :: qptrlatt(3,3) 
 real(dp),intent(in) :: qshift(3,nqshift)
 real(dp),intent(in) :: n0(Gam%nsppol)

!Local variables -------------------------
!scalars
 integer,parameter :: iout0=0,chksymbreak0=0,iscf2=2,qptopt1=1
 integer :: my_qptopt,iq_ibz,nqibz,nqpt_computed
 integer :: my_nqshift,ibranch,iw,nbranch,nsppol,spin
 real(dp) :: a2fprefactor,gaussfactor,gaussprefactor,gaussval,lambda_2,lambda_3,lambda_4,lambda_5
 real(dp) :: lambda_iso,omega,omega_log,xx,domega,qptrlen,omega_min,omega_max
 character(len=500) :: msg
!arrays
 integer,parameter :: vacuum0(3)=(/0,0,0/)
 integer :: my_qptrlatt(3,3) 
 real(dp),allocatable :: my_qshift(:,:)
 real(dp) :: displ_cart(2,Gam%nbranch,Gam%nbranch)
 real(dp) :: eigval(Gam%nbranch),pheigval(Gam%nbranch)
 real(dp) :: pheigvec(2*Gam%nbranch**2),phfrq(Gam%nbranch)
 real(dp) :: tmp_a2f(nomega),a2f_1d(nomega)
 real(dp),allocatable :: qibz(:,:),wtq(:)
 real(dp),allocatable :: a2f1mom(:),a2f2mom(:),a2f3mom(:),a2f4mom(:)
 real(dp),allocatable :: a2f_1mom(:),a2f_1mom_int(:),a2flogmom(:),a2flogmom_int(:)
 real(dp),allocatable :: freq(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 nbranch   =  Gam%nbranch
 nsppol    =  Gam%nsppol
 !
 ! Generate the q-mesh finding the IBZ and the corresponding weights.
 my_qptopt   = qptopt1; if (PRESENT(qptopt)) my_qptopt = qptopt
 my_qptrlatt = qptrlatt
 my_nqshift  = nqshift ! Be careful as getkgrid expects shiftk(3,8).
 ABI_CHECK(my_nqshift>0.and.my_nqshift<=8,"nqshift must be in [1,8]") 

 ABI_ALLOCATE(my_qshift,(3,8))
 my_qshift=zero; my_qshift(:,1:nqshift) = qshift(:,:)
 !
 ! First call to getkgrid to obtain nqibz.
 call getkgrid(chksymbreak0,iout0,iscf2,qibz,my_qptopt,my_qptrlatt,qptrlen,&
&  Cryst%nsym,0,nqibz,my_nqshift,Cryst%nsym,Cryst%rprimd,my_qshift,Cryst%symafm,Cryst%symrel,vacuum0,wtq)
 !
 ! Recall getkgrid to get qibz and wtq.
 ABI_ALLOCATE(qibz,(3,nqibz))
 ABI_ALLOCATE(wtq,(nqibz))

 call getkgrid(chksymbreak0,iout0,iscf2,qibz,my_qptopt,my_qptrlatt,qptrlen,&
&  Cryst%nsym,nqibz,nqpt_computed,my_nqshift,Cryst%nsym,Cryst%rprimd,my_qshift,Cryst%symafm,Cryst%symrel,vacuum0,wtq)

 !do iq_ibz=1,nqibz
 !  write(std_out,*)"wtq",wtq(iq_ibz)
 !end do
 !write(std_out,*)"SUM(wtq)",SUM(wtq)

 ! Store quantities that cannot be easily (and safely) calculated if we only know the IBZ.
 A2f%nqshift  = my_nqshift
 A2f%qptrlatt = my_qptrlatt

 ABI_ALLOCATE(A2f%qshift,(3,my_nqshift))
 A2f%qshift=my_qshift(:,1:my_nqshift)
 ABI_DEALLOCATE(my_qshift)
 !
 ! Calculate and store the phonon frequencies in the IBZ so that we have the min and Max frequency for the mesh.
 ABI_ALLOCATE(freq,(nbranch,nqibz))
 do iq_ibz=1,nqibz
   call inpphon(displ_cart,pheigval,pheigvec,freq(:,iq_ibz),Phon_ds,qibz(:,iq_ibz))
 end do 
 !
 ! Min and max frequency for the mesh.
 omega_min = MINVAL(freq(1,:))
 omega_max = MAXVAL(freq(nbranch,:))

 omega_min = omega_min - 0.1*ABS(omega_min)
 omega_max = omega_max + 0.1*ABS(omega_max)

 !@a2f_t
 ! Initialization of the a2f_t structure.

 domega = (omega_max-omega_min)/(nomega-one)
 !elph_ds%domega  = domega  ! MG Why do we need to store domega in elph_ds?
 !omega_min       = elph_ds%omega_min
 !omega_max       = elph_ds%omega_max

 call a2f_nullify(A2f)
 A2f%nomega    = nomega
 A2f%nsppol    = nsppol
 A2f%a2fsmear  = a2fsmear
 A2f%omega_min = omega_min
 A2f%omega_max = omega_max

 ABI_ALLOCATE(A2f%n0,(nsppol))
 A2f%n0=n0
 ABI_ALLOCATE(A2f%a2f,(nomega,nsppol))
 A2f%a2f=zero 
 ABI_ALLOCATE(A2f%omega,(nomega))
 !
 ! Build linear mesh.
 A2f%omega = arth(omega_min,domega,nomega)
 !
 ! TODO: it seems there's a factor 2 missing here!
 gaussprefactor = sqrt(piinv) / a2fsmear
 gaussfactor = one / a2fsmear
 !
 do spin=1,nsppol
   a2f_1d = zero
   !
   ! Loop over qpoints in the IBZ 
   do iq_ibz=1,nqibz

     call inpphon(displ_cart,pheigval,pheigvec,phfrq,Phon_ds,qibz(:,iq_ibz))
     !phfrq = freq(:,iq_ibz)

     call gamma_interp(Gam,Cryst,alter_int_gam,spin,qibz(:,iq_ibz),displ_cart,eigval)
     !eigval = one
     !
     ! Add all contributions from the phonon modes at this qpoint to a2f (note that unstable modes are included).
     do ibranch=1,nbranch
       if (ABS(phfrq(ibranch)) < tol10) then
         a2fprefactor= zero
       else
         a2fprefactor = eigval(ibranch)/(two_pi*ABS(phfrq(ibranch))*n0(spin))
       end if

       tmp_a2f(:) = zero
       do iw=1,nomega
         xx = (A2f%omega(iw)-phfrq(ibranch))*gaussfactor
         gaussval = gaussprefactor*EXP(-xx*xx)

         tmp_a2f(iw) = tmp_a2f(iw) + gaussval*a2fprefactor
         !% tmp_a2f(iw) = tmp_a2f(iw) + gaussian(xx,a2fsmear)
       end do

       a2f_1d(:) = a2f_1d(:) + tmp_a2f(:) * wtq(iq_ibz)
     end do ! ibranch
   end do ! iq_ibz
   !
   ! 1/nkpt factor for the integration weights.
   !a2f_1d     = a2f_1d/nkpt

   A2f%a2f(:,spin) = a2f_1d

   write(std_out,*) ' from a2f_init: for spin ', spin
   write(std_out,'(a,i2,a,e16.6)') '# The DOS at Fermi level for spin ',spin,' is ',n0(spin)
   !
   ! Do isotropic calculation of lambda and output lambda, Tc(MacMillan)
   !
   ABI_ALLOCATE(a2f_1mom,(nomega))
   ABI_ALLOCATE(a2f_1mom_int,(nomega))
   ABI_ALLOCATE(a2f1mom,(nomega))
   ABI_ALLOCATE(a2f2mom,(nomega))
   ABI_ALLOCATE(a2f3mom,(nomega))
   ABI_ALLOCATE(a2f4mom,(nomega))

   a2f_1mom=zero; a2f_1mom_int=zero
   a2f1mom =zero; a2f2mom     =zero
   a2f3mom =zero; a2f4mom     =zero
   
   omega = omega_min
   do iw=1,nomega
     if (ABS(omega) > tol10) then
       a2f_1mom(iw) = two*a2f_1d(iw)/ABS(omega)   ! first inverse moment of alpha2F
       a2f1mom(iw)  = two*a2f_1d(iw)*ABS(omega)   ! first positive moment of alpha2F
       a2f2mom(iw)  =     a2f1mom(iw)*ABS(omega)  ! second positive moment of alpha2F. Factor of 2 is included in a2f1mom recursively
       a2f3mom(iw)  =     a2f2mom(iw)*ABS(omega)  ! third positive moment of alpha2F
       a2f4mom(iw)  =     a2f3mom(iw)*ABS(omega)  ! fourth positive moment of alpha2F
     end if
     omega=omega + domega
   end do
   !  
   ! From Allen PRL 59 1460
   !  \lambda <\omega^n> = 2 \int_0^{\infty} d\omega [\alpha^2F / \omega] \omega^n
   !  
   call simpson_int(nomega,domega,a2f_1mom,a2f_1mom_int)
   lambda_iso = a2f_1mom_int(nomega)

   call simpson_int(nomega,domega,a2f1mom,a2f_1mom_int)
   lambda_2 = a2f_1mom_int(nomega)

   call simpson_int(nomega,domega,a2f2mom,a2f_1mom_int)
   lambda_3 = a2f_1mom_int(nomega)

   call simpson_int(nomega,domega,a2f3mom,a2f_1mom_int)
   lambda_4 = a2f_1mom_int(nomega)

   call simpson_int(nomega,domega,a2f4mom,a2f_1mom_int)
   lambda_5 = a2f_1mom_int(nomega)

   ABI_DEALLOCATE(a2f_1mom)
   ABI_DEALLOCATE(a2f_1mom_int)
   ABI_DEALLOCATE(a2f1mom)
   ABI_DEALLOCATE(a2f2mom)
   ABI_DEALLOCATE(a2f3mom)
   ABI_DEALLOCATE(a2f4mom)

   write(std_out,*)' a2f_init: elphon coupling lambdas for spin = ',spin
   write(std_out,*)' a2f_init: isotropic lambda', lambda_iso
   write(std_out,*)' a2f_init: positive moments of alpha2F:'
   write(std_out,*)' lambda <omega^2> = ',lambda_2
   write(std_out,*)' lambda <omega^3> = ',lambda_3
   write(std_out,*)' lambda <omega^4> = ',lambda_4
   write(std_out,*)' lambda <omega^5> = ',lambda_5
   !
   ! Get log moment of alpha^2F.
   ABI_ALLOCATE(a2flogmom,(nomega))
   ABI_ALLOCATE(a2flogmom_int,(nomega))
   a2flogmom = zero

   omega = omega_min
   do iw=1,nomega
     if (ABS(omega) > tol10) then
       a2flogmom(iw) = (two/lambda_iso)*a2f_1d(iw)*LOG(ABS(omega))/ABS(omega)
     end if
     omega=omega + domega
   end do
   call simpson_int(nomega,domega,a2flogmom,a2flogmom_int)
   omega_log = EXP(a2flogmom_int(nomega))

   ABI_DEALLOCATE(a2flogmom)
   ABI_DEALLOCATE(a2flogmom_int)
   
   !tc_macmill = omega_log/1.2_dp * EXP((-1.04_dp*(one+lambda_iso)) / (lambda_iso-mustar*(one+0.62_dp*lambda_iso)))

   if (nsppol > 1) then
     write(msg,'(3a)') ch10,&
&     ' Warning : some of the following quantities should be integrated over spin', ch10
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if

   write(msg,'(3a)') ch10,&
&   ' Superconductivity : isotropic evaluation of parameters from electron-phonon coupling.',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write(msg,'(a,2es16.6)' )' a2f_init: isotropic lambda = ',lambda_iso,a2f_moment(A2f,spin,0)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write(msg,'(a,es16.6)' )' a2f_init: lambda <omega^2> = ',lambda_2
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write(msg,'(a,es16.6)' )' a2f_init: lambda <omega^3> = ',lambda_3
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write(msg,'(a,es16.6)' )' a2f_init: lambda <omega^4> = ',lambda_4
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write(msg,'(a,es16.6)' )' a2f_init: lambda <omega^5> = ',lambda_5
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write(msg,'(a,es16.6,a,es16.6,a)' )' a2f_init: omegalog  = ',omega_log,' (Ha) ', omega_log/kb_HaK, ' (Kelvin) '
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   !write(msg,'(a,es16.6,a,es16.6,a)')' a2f_init: MacMillan Tc = ',tc_macmill,' (Ha) ', tc_macmill/kb_HaK, ' (Kelvin) '
 end do ! spin

 ABI_DEALLOCATE(freq)
 ABI_DEALLOCATE(qibz)
 ABI_DEALLOCATE(wtq)

 DBG_ENTER("COLL")

end subroutine a2f_init
!!***

!----------------------------------------------------------------------

!!****f* m_gamma/a2f_moment
!! NAME
!! a2f_moment
!!
!! FUNCTION
!!  This routine calculates 
!!     2 \int_0^{\infty} d\omega [\alpha^2F / \omega] \omega^n
!!
!! INPUTS
!!  A2f<a2f_t>=Structure storing the Eliashberg function.
!!  spin=The spin components 
!!  nn
!!
!! OUTPUT
!!  a2f_moment = 2 \int_0^{\infty} d\omega [\alpha^2F / \omega] \omega^n
!!  [out_int(w)] = 2 \int_0^w d\omega [\alpha^2F / \omega] \omega^n
!!
!! NOTES
!!  From Allen PRL 59 1460 (See also Grimvall, Eq 6.72 page 175)
!!  \lambda <\omega^n> = 2 \int_0^{\infty} d\omega [\alpha^2F / \omega] \omega^n
!!  [out_int(A2f%nomega)]
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function a2f_moment(A2f,spin,nn,out_int)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'a2f_moment'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin,nn
 real(dp) :: a2f_moment
 type(a2f_t),intent(in) :: A2f
!arrays
 real(dp),intent(out),optional :: out_int(A2f%nomega)

!Local variables -------------------------
!scalars
 integer :: iw,nomega
 real(dp) :: domega,omg,omg_nm1
!arrays
 real(dp),allocatable :: ff(:),int_ff(:)

! *********************************************************************

 nomega = A2f%nomega
 domega = (A2f%omega_max - A2f%omega_min)/(nomega-one)
 !
 ! Construct the integrand function.
 ! From Allen PRL 59 1460
 !  \lambda <\omega^n> = 2 \int_0^{\infty} d\omega [\alpha^2F / \omega] \omega^n
 ABI_ALLOCATE(ff,(nomega))
 ABI_ALLOCATE(int_ff,(nomega))
 ff=zero; int_ff=zero

 if (nn-1>=0) then
   do iw=1,nomega
     !omg = A2f%omega(iw)
     omg = ABS(A2f%omega(iw)) ! FIXME: this trick is needed to reproduce the automatic tests!
     omg_nm1 = omg**(nn-1)
     ff(iw) = two* A2f%a2f(iw,spin) * omg_nm1
   end do
 else 
   do iw=1,nomega
     !omg = A2f%omega(iw)
     omg = ABS(A2f%omega(iw)) ! FIXME: this trick is needed to reproduce the automatic tests!
     omg_nm1 = zero 
     if (ABS(omg) > tol10) omg_nm1 = omg**(nn-1)
     ff(iw) = two * A2f%a2f(iw,spin) * omg_nm1
   end do
 end if
 !  
 ! Integration with Simpson rule on a linear mesh.
 call simpson_int(nomega,domega,ff,int_ff)

 a2f_moment = int_ff(nomega)
 if (PRESENT(out_int)) out_int = int_ff

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(int_ff)

end function a2f_moment
!!***

!----------------------------------------------------------------------

!!****f* m_gamma/a2f_logmoment
!! NAME
!! a2f_logmoment
!!
!! FUNCTION
!!
!! INPUTS
!!  A2f<a2f_t>=Structure storing the Eliashberg function.
!!  spin=The spin components 
!!  nn
!!
!! OUTPUT
!!  a2f_logmoment = 
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function a2f_logmoment(A2f,spin,nn)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'a2f_logmoment'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin,nn
 real(dp) :: a2f_logmoment
 type(a2f_t),intent(in) :: A2f
!arrays

!Local variables -------------------------
!scalars
 integer :: iw,nomega
 real(dp) :: domega,omg,lambda_iso
!arrays
 real(dp),allocatable :: a2flogmom(:),a2flogmom_int(:)

! *********************************************************************

 ABI_UNUSED(nn)

 nomega = A2f%nomega
 domega = (A2f%omega_max - A2f%omega_min)/(nomega-one)
 !
 ! Get log moment of alpha^2F.
 ABI_ALLOCATE(a2flogmom,(nomega))
 ABI_ALLOCATE(a2flogmom_int,(nomega))
 a2flogmom = zero

 lambda_iso = a2f_moment(A2f,spin,0)
                                                                            
 do iw=1,nomega
   omg = A2f%omega(iw)
   if (ABS(omg) > tol10) then
     a2flogmom(iw) = (two/lambda_iso) * A2f%a2f(iw,spin) * LOG(ABS(omg))/ABS(omg)
   end if
 end do

 call simpson_int(nomega,domega,a2flogmom,a2flogmom_int)
 a2f_logmoment = EXP(a2flogmom_int(nomega))

 ABI_DEALLOCATE(a2flogmom)
 ABI_DEALLOCATE(a2flogmom_int)

end function a2f_logmoment
!!***

!----------------------------------------------------------------------

!!****f* m_gamma/a2f_dump
!! NAME
!! a2f_dump
!!
!! FUNCTION
!!
!! INPUTS
!!  A2f<a2f_t>=Container storing the Eliashberg functions.
!!  fname=Filename for output.
!!
!! OUTPUT
!!  Output is written on file.
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine a2f_dump(A2f,fname)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'a2f_dump'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: fname
 type(a2f_t),intent(in) :: A2f

!Local variables -------------------------
!scalars
 integer :: iw,ios,spin,unt
 !real(dp) :: omega_log,lambda_2,lambda_3,lambda_4,lambda_5
 !character(len=500) :: msg

! *********************************************************************
 !
 ! Open the file.
 unt = get_unit()
 open(unit=unt,file=fname,status='unknown',iostat=ios)
 ABI_CHECK(ios==0,"Opening "//TRIM(fname))
 !
 ! Output the a2f header.
 write(unt,'(a)')                '#'
 write(unt,'(a)')                '# ABINIT package : a2f file'
 write(unt,'(a)')                '#'
 write(unt,'(a)')                '# a2f function integrated over the FS. omega in a.u.'
 !write(unt,'(a,i10)')            '#  number of q-points integrated over : ',nkpt
 write(unt,'(a,i0)')             '#  number of energy points : ',A2f%nomega
 write(unt,'(a,e16.6,a,e16.6,a)')'#  between omega_min = ',A2f%omega_min,' Ha and omega_max = ',A2f%omega_max,' Ha'
 write(unt,'(a,e16.6)')          '#  and the smearing width for gaussians is ',A2f%a2fsmear
 !write(unt,'(a)')                '#'
 !write(unt,'(a)')                '#'
 ! TODO Do isotropic calculation of lambda and output lambda, Tc(MacMillan)
 !
 do spin=1,A2f%nsppol
   !
   write(unt,'(a,i2,a,e16.6)')'# The DOS at Fermi level for spin ',spin,' is ',A2f%n0(spin)
   write(unt,'(a)')           '#'
   !
   do iw=1,A2f%nomega
     write(unt,*) A2f%omega(iw), A2f%a2f(iw,spin)
   end do
   write(unt,*)
   !
 end do

 close(unt)

end subroutine a2f_dump

end module m_gamma
!!***
