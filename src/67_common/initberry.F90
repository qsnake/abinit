!{\src2tex{textfont=tt}}
!!****f* ABINIT/initberry
!! NAME
!! initberry
!!
!! FUNCTION
!! Initialization of Berryphase calculation of the polarization, the
!! ddk and the response of an insulator to a homogenous electric field.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVeithen).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables in this dataset
!!  gmet(3,3) = reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3) = primitive translations in recip space
!!  kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!!  mband = maximum number of bands
!!  mkmem = maximum number of k-points in core memory
!!  mpw = maximum number of plane waves
!!  natom = number of atoms in unit cell
!!  nkpt = number of k points
!!  npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!!  nsppol = 1 for unpolarized, 2 for spin-polarized
!!  nsym = number of symmetry operations
!!  ntypat = number of types of atoms in unit cell
!!  occ(mband*nkpt*nsppol) = occup number for each band at each k point
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3) = dimensional primitive vectors
!!  symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!    reciprocal space primitive translations
!!  typat = typat(natom) list of atom types
!!  usepaw = flag for PAW (1 PAW, 0 NCPP)
!!  xred(3,natom) = location of atoms in reduced units
!!
!! OUTPUT
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations
!!  pwind(pwind_alloc,2,3) = array used to compute the overlap matrix smat
!!                         between k-points k and k +- dk where dk is
!!                         parallel to the direction idir
!!    jpw = pwind(ipw,ifor,idir)
!!      * ipw = index of plane wave vector G for a given k-point k
!!      * ifor = 1: k + dk
!!               2: k - dk
!!      * idir = direction of the polarization/ddk calculation [dk(idir)
!!               is the only non-zero element of dk(:)]
!!      * jpw = index of plane wave vector G (+dG) at k +- dk
!!              where dG is a shift of one reciprocal lattice vector
!!              (required to close the strings of k-points using the
!!               periodic gauge condition)
!!    In case a G-vector of the basis sphere of plane waves at k
!!    does not belong to the basis sphere of plane waves at k+dk, jpw = 0.
!!   pwind_alloc = first dimension of pwind and pwnsfac
!!   pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!
!! SIDE EFFECTS
!!  mpi_enreg = informations about MPI parallelization
!!    kptdstrb(nproc,nneighbour,fmkmem_max*nsppol) : Array required
!!      by berryphase_new.f for MPI // over k-points. Defined
!!      for k-points in the fBZ
!!    kptdstrbi(nproc,nneighbour,mkmem_max*nsppol) : Same as kptdstrb
!!      but for k-points in the iBZ. Used by vtorho.f
!!           nproc = number of cpus
!!           nneighbour = number of neighbours for each k-point (= 6)
!!
!! TO DO
!!
!! NOTES
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      cprj_alloc,expibi,kpgsph,leave_new,lij,listkk,pawtwdij_1,pawtwdij_2a
!!      pawtwdij_2b,pawtwdij_2c,pawtwdij_2d,pawtwdij_2e,pawtwdij_2f,qijb_bk
!!      qijb_kk,set_twind,setsymrhoij,smpbz,symatm,wrtout,xcomm_world,xmax_mpi
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initberry(dtefield,dtset,gmet,gprimd,kg,mband,&
&              mkmem,mpi_enreg,mpw,natom,nkpt,npwarr,nsppol,&
&              nsym,ntypat,occ,pawang,pawrad,pawtab,psps,&
&              pwind,pwind_alloc,pwnsfac,&
&              rprimd,symrec,typat,usepaw,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initberry'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
 use interfaces_66_paw
 use interfaces_66_wfs
 use interfaces_67_common, except_this_one => initberry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,natom,nkpt,nsppol,nsym,ntypat,usepaw
 integer,intent(out) :: pwind_alloc
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset
 type(efield_type),intent(out) :: dtefield
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
 integer,intent(in) :: symrec(3,3,nsym),typat(natom)
 integer,pointer :: pwind(:,:,:)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 real(dp),pointer :: pwnsfac(:,:)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: count,exchn2n3d,flag,flag_kpt,fnkpt_computed,iatom,iband,icg,icprj
 integer :: idir,idum,idum1,ierr,ifor,ihcg,ikg,ikg1,ikpt,ikpt1,ikpt1f
 integer :: ikpt1i,ikpt2,ikpt_loc,ikptf,ikpti,ikstr,index,ineigh,ipw,ipwnsfac
 integer :: isppol,istr,istwf_k,isym,isym1,itrs,itypat,iunmark,jpw,klmn,lmax,lmn2_size_max
 integer :: me,me_g0,mkmem_,my_nspinor,nband_k,nband_occ_k,ncpgr,nkstr,nproc,npw_k,npw_k1,spaceComm
 integer :: option, brav, mkpt, nkptlatt
 integer :: jstr,ii,jj
 integer :: dk_flag, coord1, coord2
 integer :: mult
!integer :: imult

 real(dp) :: c1,ecut_eff,eg,eg_ev,rdum
 real(dp) :: dist_, max_dist, last_dist, dist
 character(len=500) :: message
 logical :: fieldflag,use_symrec
!arrays
 integer :: dg(3),iadum(3),iadum1(3),neigh(6)
 integer,allocatable :: buffer(:),dimlmn(:),kg1_k(:,:),kpt_mark(:)
 real(dp) :: diffk(3),dk(3),dum33(3,3),eg_dir(3)
 real(dp) :: kpt1(3)
 real(dp) :: delta_str3(2), dstr(2),dk_str(2,2,3)
 real(dp),allocatable :: spkpt(:,:)
! real(dp),allocatable :: dist_str2(:,:)

! *************************************************************************

!DEBUG
 write(std_out,*)'initberry: enter'
!entry_count=entry_count+1
!write(std_out,*)' initberry count=',entry_count
!call leave_new('COLL')
!ENDDEBUG

!save the current value of berryopt
 dtefield%berryopt = dtset%berryopt

!----------------------------------------------------------------------------
!-------------------- Obtain k-point grid in the full BZ --------------------
!----------------------------------------------------------------------------

 if(dtset%kptopt==1 .or. dtset%kptopt==2 .or. dtset%kptopt==4)then
!  Compute the number of k points in the G-space unit cell
   nkptlatt=dtset%kptrlatt(1,1)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,3) &
&   +dtset%kptrlatt(1,2)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,1) &
&   +dtset%kptrlatt(1,3)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,2) &
&   -dtset%kptrlatt(1,2)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,3) &
&   -dtset%kptrlatt(1,3)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,1) &
&   -dtset%kptrlatt(1,1)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,2)

!  Call smpbz to obtain the list of k-point in the full BZ - without symmetry reduction
   option = 0
   brav = 1
   mkpt=nkptlatt*dtset%nshiftk
   ABI_ALLOCATE(spkpt,(3,mkpt))
   call smpbz(1,ab_out,dtset%kptrlatt,mkpt,fnkpt_computed,dtset%nshiftk,option,dtset%shiftk,spkpt)
   dtefield%fnkpt = fnkpt_computed
   ABI_ALLOCATE(dtefield%fkptns,(3,dtefield%fnkpt))
   dtefield%fkptns(:,:)=spkpt(:,1:dtefield%fnkpt)
   ABI_DEALLOCATE(spkpt)
 else if(dtset%kptopt==3.or.dtset%kptopt==0)then
   dtefield%fnkpt=nkpt
   ABI_ALLOCATE(dtefield%fkptns,(3,dtefield%fnkpt))
   dtefield%fkptns(1:3,1:dtefield%fnkpt)=dtset%kpt(1:3,1:dtefield%fnkpt)
   if(dtset%kptopt==0)then
     write(message,'(10a)') ch10,&
&     ' initberry : WARNING -',ch10,&
&     '  you have defined manually the k-point grid with kptopt = 0',ch10,&
&     '  the berry phase calculation works only with a regular k-points grid,',ch10,&
&     '  abinit doesn''t check if your grid is regular...'
     call wrtout(std_out,message,'PERS')
   end if
 end if

!call listkk to get mapping from FBZ to IBZ
 rdum=1.0d-5  ! cutoff distance to decide when two k points match
 ABI_ALLOCATE(dtefield%indkk_f2ibz,(dtefield%fnkpt,6))

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)

!ji: The following may need modification in the future
!**** no spin-polarization doubling ; allow use of time reversal symmetry ****

!Here is original call
!
!call listkk(rdum,gmet,dtefield%indkk_f2ibz,dtset%kptns,dtefield%fkptns,nkpt,&
!& dtefield%fnkpt,dtset%nsym,1,dtset%symafm,dtset%symrel,1)

 use_symrec = .TRUE.
 call listkk(rdum,gmet,dtefield%indkk_f2ibz,dtset%kptns,dtefield%fkptns,nkpt,&
& dtefield%fnkpt,dtset%nsym,1,dtset%symafm,symrec,1,use_symrec)


!Construct i2fbz and f2ibz
 ABI_ALLOCATE(dtefield%i2fbz,(nkpt))
 idum=0
 do ikpt=1,dtefield%fnkpt
   if (dtefield%indkk_f2ibz(ikpt,2)==1 .and. &
&   dtefield%indkk_f2ibz(ikpt,6) == 0 .and. &
&   maxval(abs(dtefield%indkk_f2ibz(ikpt,3:5))) == 0 ) then
     dtefield%i2fbz(dtefield%indkk_f2ibz(ikpt,1))=ikpt
     idum=idum+1
   end if
 end do
 if (idum/=nkpt)then
   write(message,'(a,a,a,a)')ch10,&
&   ' initberry: ERROR - ',ch10,&
&   '   Found wrong number of k-points in IBZ'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!set a flag to indicate finite electric fields or magnetization calcs
 fieldflag = .FALSE.
 if ( (dtset%berryopt==4) .or. (abs(dtset%berryopt)==5) ) fieldflag = .TRUE.

!----------------------------------------------------------------------------
!------------- Allocate PAW space if necessary ------------------------------
!----------------------------------------------------------------------------

 if (usepaw == 1) then

   dtefield%natom = natom
   dtefield%usepaw = usepaw

   ABI_ALLOCATE(dtefield%lmn_size,(ntypat))
   ABI_ALLOCATE(dtefield%lmn2_size,(ntypat))
   do itypat = 1, ntypat
     dtefield%lmn_size(itypat) = pawtab(itypat)%lmn_size
     dtefield%lmn2_size(itypat) = pawtab(itypat)%lmn2_size
   end do

   lmn2_size_max = psps%lmnmax*(psps%lmnmax+1)/2

   ABI_ALLOCATE(dtefield%qijb_kk,(2,lmn2_size_max,natom,3))
   ABI_ALLOCATE(dtefield%expibi,(2,natom,9))
   dtefield%has_expibi = 1
   if (abs(dtset%berryopt) == 5) then
     ABI_ALLOCATE(dtefield%qijb_bk,(2,lmn2_size_max,natom,6))
     ABI_ALLOCATE(dtefield%twdij0,(2,psps%lmnmax,psps%lmnmax,natom,24))
     ABI_ALLOCATE(dtefield%twdij,(2,psps%lmnmax,psps%lmnmax,natom,24))
     dtefield%has_twdij0 = 1
     ABI_ALLOCATE(dtefield%tweijkl,(2,lmn2_size_max,lmn2_size_max,natom,6))
     dtefield%has_tweijkl = 1
   end if
   dtefield%has_qijb = 1

   if (dtset%berryopt==4 .and. dtefield%has_rij==0) then
     lmn2_size_max = psps%lmnmax*(psps%lmnmax+1)/2
     ABI_ALLOCATE(dtefield%rij,(lmn2_size_max,ntypat,3))
     dtefield%has_rij = 1
   end if

   if (abs(dtset%berryopt)==5 .and. &
   (dtefield%has_Lij==0 .or. dtefield%has_Lijr3==0)) then
     lmn2_size_max = psps%lmnmax*(psps%lmnmax+1)/2
     ABI_ALLOCATE(dtefield%Lij,(2,lmn2_size_max,ntypat,3))
     ABI_ALLOCATE(dtefield%Lijr3,(2,lmn2_size_max,ntypat,3))
     dtefield%has_Lij = 1
     dtefield%has_Lijr3 = 1
   end if

   if ( fieldflag .and. dtefield%usecprj == 0) then
     ABI_ALLOCATE(dimlmn,(natom))
     do iatom=1,natom
       dimlmn(iatom)=pawtab(typat(iatom))%lmn_size
     end do
!    allocate space for cprj at kpts in BZ (IBZ or FBZ)
     ABI_ALLOCATE(dtefield%cprj,(natom, mband*dtset%nkpt*nsppol))
!    write(std_out,*) "initberry alloc of cprj ", shape(dtefield%cprj)
     ncpgr = 3 ! allows for gradients wrt atom positions, may have to add strains later
     call cprj_alloc(dtefield%cprj,ncpgr,dimlmn)
     dtefield%usecprj = 1
     ABI_DEALLOCATE(dimlmn)
   end if

   ABI_ALLOCATE(dtefield%cprjindex,(nkpt,nsppol))
   dtefield%cprjindex(:,:) = 0

   if (dtset%kptopt /= 3) then
     ABI_ALLOCATE(dtefield%atom_indsym,(4,nsym,natom))
     call symatm(dtefield%atom_indsym,natom,nsym,symrec,dtset%tnons,tol8,typat,xred)
     lmax = psps%mpsang - 1
     ABI_ALLOCATE(dtefield%zarot,(2*lmax+1,2*lmax+1,lmax+1,nsym))
     call setsymrhoij(gprimd,lmax,nsym,1,rprimd,symrec,dtefield%zarot)
     dtefield%nsym = nsym
     dtefield%lmax = lmax
     dtefield%lmnmax = psps%lmnmax
   end if

   if (abs(dtset%berryopt) == 5) then
     ABI_ALLOCATE(dtefield%mag_local_k,(3,nkpt*nsppol))
     ABI_ALLOCATE(dtefield%mag_k,(2,nkpt*nsppol,3))
   end if

 end if

!------------------------------------------------------------------------------
!------------------- Compute variables related to MPI // ----------------------
!------------------------------------------------------------------------------

 call xcomm_world(mpi_enreg,spaceComm,myrank=me,mysize=nproc)

 if (mpi_enreg%paral_compil_kpt == 0) then  ! no MPI //

   dtefield%fmkmem = dtefield%fnkpt
   dtefield%fmkmem_max = dtefield%fnkpt
   dtefield%mkmem_max = nkpt

!  we allocate also those
   ABI_ALLOCATE(mpi_enreg%kpt_loc2fbz_sp,(0:nproc-1,1:dtefield%fmkmem_max*nsppol, 1:2))
   ABI_ALLOCATE(mpi_enreg%mkmem,(0:nproc-1))
   ABI_ALLOCATE(mpi_enreg%fmkmem,(0:nproc-1))
   mpi_enreg%kpt_loc2fbz_sp(:,:,:) = 0
   mpi_enreg%fmkmem(:) = 0
   mpi_enreg%mkmem(:) = 0

   ABI_ALLOCATE(mpi_enreg%kpt_loc2ibz_sp,(0:nproc-1,1:dtefield%mkmem_max*nsppol, 1:2))
   mpi_enreg%kpt_loc2ibz_sp(:,:,:) = 0

 else    ! MPI //

!  Number of k-points in the FBZ for each cpu

   dtefield%fmkmem = 0
   do ikpt = 1, dtefield%fnkpt
     ikpti = dtefield%indkk_f2ibz(ikpt,1)
     nband_k = dtset%nband(ikpti)
     if (minval(abs(mpi_enreg%proc_distrb(ikpti,1:nband_k,1:nsppol)-me))==0) &
&     dtefield%fmkmem = dtefield%fmkmem + 1
   end do

!  Maximum value of mkmem and fmkmem
   call xmax_mpi(dtefield%fmkmem,dtefield%fmkmem_max,spaceComm,ierr)

   mkmem_ = mkmem   ! I have to use the dummy variable mkmem_ because
!  mkmem is declared as intent(in) while the first
!  argument of xmax_mpi must be intent(inout)

   call xmax_mpi(mkmem_,dtefield%mkmem_max,spaceComm,ierr)

   ABI_ALLOCATE(mpi_enreg%kptdstrb,(nproc,6,dtefield%fmkmem_max*nsppol*2))
   mpi_enreg%kptdstrb(:,:,:) = 0

   ABI_ALLOCATE(mpi_enreg%kpt_loc2fbz_sp,(0:nproc-1,1:dtefield%fmkmem_max*nsppol, 1:2))
   ABI_ALLOCATE(mpi_enreg%mkmem,(0:nproc-1))
   ABI_ALLOCATE(mpi_enreg%fmkmem,(0:nproc-1))
   mpi_enreg%kpt_loc2fbz_sp(:,:,:) = 0
   mpi_enreg%fmkmem(:) = 0
   mpi_enreg%mkmem(:) = 0


   if (fieldflag) then
     ABI_ALLOCATE(mpi_enreg%kptdstrbi,(nproc,6,dtefield%mkmem_max*nsppol*2))
     ABI_ALLOCATE(dtefield%cgqindex,(3,6,nkpt*nsppol))
     ABI_ALLOCATE(dtefield%nneigh,(nkpt))
     mpi_enreg%kptdstrbi(:,:,:) = 0
     dtefield%cgqindex(:,:,:) = 0 ; dtefield%nneigh(:) = 0

   end if

   ABI_ALLOCATE(mpi_enreg%kpt_loc2ibz_sp,(0:nproc-1,1:dtefield%mkmem_max*nsppol, 1:2))
   mpi_enreg%kpt_loc2ibz_sp(:,:,:) = 0

 end if

 pwind_alloc = mpw*dtefield%fmkmem_max
 ABI_ALLOCATE(pwind,(pwind_alloc,2,3))
 ABI_ALLOCATE(pwnsfac,(2,pwind_alloc))


!------------------------------------------------------------------------------
!---------------------- Compute efield_type variables -------------------------
!------------------------------------------------------------------------------

!Initialization of efield_type variables
 mult=dtset%useria+1
 dtefield%efield_dot(:) = zero
 dtefield%dkvecs(:,:) = zero
 dtefield%maxnstr = 0    ; dtefield%maxnkstr  = 0
 dtefield%nstr(:) = 0    ; dtefield%nkstr(:) = 0
 ABI_ALLOCATE(dtefield%ikpt_dk,(dtefield%fnkpt,2,3))
 ABI_ALLOCATE(dtefield%cgindex,(nkpt,nsppol))
 ABI_ALLOCATE(dtefield%kgindex,(nkpt))
 ABI_ALLOCATE(dtefield%fkgindex,(dtefield%fnkpt))
 dtefield%ikpt_dk(:,:,:) = 0
 dtefield%cgindex(:,:) = 0
 dtefield%nband_occ = 0
 dtefield%kgindex(:) = 0
 dtefield%fkgindex(:) = 0

 if (fieldflag) then
   dtset%rfdir(1:3) = 1
 end if


!Compute spin degeneracy
 if (nsppol == 1 .and. dtset%nspinor == 1) then
   dtefield%sdeg = two
 else if (nsppol == 2 .or. my_nspinor == 2) then
   dtefield%sdeg = one
 end if

!Compute the number of occupied bands and check that
!it is the same for each k-point

 index = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt

     nband_occ_k = 0
     nband_k = dtset%nband(ikpt + (isppol - 1)*nkpt)

     do iband = 1, nband_k
       index = index + 1
       if (abs(occ(index) - dtefield%sdeg) < tol8) nband_occ_k = nband_occ_k + 1
     end do

     if (fieldflag) then
       if (nband_k /= nband_occ_k) then
         write(message,'(a,a,a,a,a,a)')ch10,&
&         ' initberry: ERROR - ',ch10,&
&         '  In a finite electric field or magnetization calculation, nband must be equal ',&
&         ch10,&
&         '  to the number of valence bands.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end if

     if ((ikpt > 1).or.(isppol > 1)) then
       if (dtefield%nband_occ /= nband_occ_k) then
         write(message,'(a,a,a,a)')ch10,&
&         ' initberry: ERROR - ',ch10,&
&         '   The number of valence bands is not the same for every k-point'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     else
       dtefield%nband_occ = nband_occ_k
     end if

   end do                ! close loop over ikpt
 end do                ! close loop over isppol

 if (fieldflag) then
   ABI_ALLOCATE(dtefield%smat,(2,dtefield%nband_occ,dtefield%nband_occ,nkpt*nsppol,2,3))

   dtefield%smat(:,:,:,:,:,:) = zero
 end if

 if (abs(dtset%berryopt)==5) then
   ABI_ALLOCATE(dtefield%twh,(2,dtefield%nband_occ,dtefield%nband_occ,nkpt*nsppol,24))
   dtefield%twh(:,:,:,:,:) = zero
   ABI_ALLOCATE(dtefield%chern_k,(2,nkpt*nsppol,3))
   dtefield%chern_k(:,:,:) = zero

   ABI_ALLOCATE(dtefield%emat,(2,dtefield%nband_occ,nkpt*nsppol))

 end if 

 if (dtset%berryopt == 5) then
   ihcg = 0
   do ikpt = 1, nkpt
     ihcg = ihcg + npwarr(ikpt)*dtefield%nband_occ
   end do
   dtefield%mhcg = ihcg
   ABI_ALLOCATE(dtefield%hcg,(2,dtefield%mhcg,24))
   dtefield%hcg(:,:,:) = zero

!  I think that dtefield%cgindex (defined below) will also give location
!  in hcg

 end if

 ABI_ALLOCATE(dtefield%sflag,(dtefield%nband_occ,nkpt*nsppol,2,3))
 dtefield%sflag(:,:,:,:) = 0

!Compute the location of each wavefunction

 icg = 0
 icprj = 0
!ikg = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt

     nband_k = dtset%nband(ikpt + (isppol-1)*nkpt)

     if (mpi_enreg%paral_compil_kpt == 1) then
       if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) then
         cycle
       end if
     end if

     dtefield%cgindex(ikpt,isppol) = icg
     npw_k = npwarr(ikpt)
     icg = icg + npw_k*nband_k

     if (usepaw == 1) then
       dtefield%cprjindex(ikpt,isppol) = icprj
       icprj = icprj + nband_k
     end if

   end do
 end do

 ikg = 0
 do ikpt = 1, nkpt
   if (mpi_enreg%paral_compil_kpt == 1) then
     if ((minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,1)-me))/=0).and.&
&     (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,nsppol)-me))/=0)) then
       cycle
     end if
   end if
   npw_k = npwarr(ikpt)
   dtefield%kgindex(ikpt) = ikg
   ikg = ikg + npw_k
 end do


!Compute the reciprocal lattice coordinates of the electric field
 if (dtset%berryopt == 4) then

   dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
   dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
   dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

   write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&   ' initberry: Reciprocal lattice coordinates of the electric field',ch10,&
&   '  efield_dot(1:3) = ',dtefield%efield_dot(1:3),ch10
   call wrtout(std_out,message,'COLL')

 end if

!Store magnetic field
 if (dtset%berryopt == 5) dtefield%bfield(:)=dtset%bfield(:)

!------------------------------------------------------------------------------
!---------------------- Build the strings of k-points -------------------------
!------------------------------------------------------------------------------

 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

!    Compute dk(:), the vector between a k-point and its nearest
!    neighbour along the direction idir

     dk(:) = zero
     dk(idir) = 1._dp   ! 1 mean there is no other k-point un the direction idir
     do ikpt = 2, dtefield%fnkpt
       diffk(:) = abs(dtefield%fkptns(:,ikpt) - dtefield%fkptns(:,1))
       if ((diffk(1) < dk(1)+tol8).and.(diffk(2) < dk(2)+tol8).and.&
&       (diffk(3) < dk(3)+tol8)) dk(:) = diffk(:)
     end do
     dtefield%dkvecs(:,idir) = dk(:)
!    DEBUG
!    write(std_out,*)' initberry : idir, dk', idir, dk
!    ENDDEBUG

!    For each k point, find k_prim such that k_prim= k + dk mod(G)
!    where G is a vector of the reciprocal lattice

     do ikpt = 1, dtefield%fnkpt

!      First: k + dk
       do ikpt1 = 1, dtefield%fnkpt
         diffk(:) = abs(dtefield%fkptns(:,ikpt1) - &
&         dtefield%fkptns(:,ikpt) - dk(:))
         if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           dtefield%ikpt_dk(ikpt,1,idir) = ikpt1
           exit
         end if
       end do

!      Second: k - dk
       do ikpt1 = 1, dtefield%fnkpt
         diffk(:) = abs(dtefield%fkptns(:,ikpt1) - &
&         dtefield%fkptns(:,ikpt) + dk(:))
         if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           dtefield%ikpt_dk(ikpt,2,idir) = ikpt1
           exit
         end if
       end do

     end do     ! ikpt

!    Find the string length, starting from k point 1
!    (all strings must have the same number of points)

     nkstr = 1
     ikpt1 = 1
     do ikpt = 1, dtefield%fnkpt
       ikpt1 = dtefield%ikpt_dk(ikpt1,1,idir)
       if (ikpt1 == 1) exit
       nkstr = nkstr + 1
     end do

!    Check that the string length is a divisor of nkpt
     if(mod(dtefield%fnkpt,nkstr) /= 0) then
       write(message,'(a,a,a,a,i5,a,i7)')ch10,&
&       ' berryphase: BUG -',ch10,&
&       '  The string length = ',nkstr,&
&       ', is not a divisor of fnkpt =',dtefield%fnkpt
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     dtefield%nkstr(idir) = nkstr
     dtefield%nstr(idir)  = dtefield%fnkpt/nkstr

   end if      ! dtset%rfdir(idir) == 1

   write(message,'(a,i1,a,i3,a,i3)')&
&   '  initberry: for direction ',idir,', nkstr = ',dtefield%nkstr(idir),&
&   ', nstr = ',dtefield%nstr(idir)
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

 end do     ! close loop over idir

 dtefield%maxnstr  = maxval(dtefield%nstr(:))
 dtefield%maxnkstr = maxval(dtefield%nkstr(:))
 ABI_ALLOCATE(dtefield%idxkstr,(dtefield%maxnkstr,dtefield%maxnstr,3))
 dtefield%idxkstr(:,:,:) = 0

!for the geometry of the string space :
 ABI_ALLOCATE(dtefield%coord_str,(2,dtefield%maxnstr,3))
 ABI_ALLOCATE(dtefield%str_neigh,(-2:2,dtefield%maxnstr,3))
 ABI_ALLOCATE(dtefield%strg_neigh,(-2:2,dtefield%maxnstr,2,3))
 dtefield%coord_str(:,:,:) = 0.d0
 dtefield%str_neigh(:,:,:)=0
 dtefield%strg_neigh(:,:,:,:)=0
 dtefield%gmet_str(:,:,:)=0.d0


!Build the different strings

 ABI_ALLOCATE(kpt_mark,(dtefield%fnkpt))
 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

     iunmark = 1
     kpt_mark(:) = 0
     do istr = 1, dtefield%nstr(idir)

       do while(kpt_mark(iunmark) /= 0)
         iunmark = iunmark + 1
       end do
       dtefield%idxkstr(1,istr,idir) = iunmark
       kpt_mark(iunmark) = 1
       do ikstr = 2, dtefield%nkstr(idir)
         ikpt1 = dtefield%idxkstr(ikstr-1,istr,idir)
         ikpt2 = dtefield%ikpt_dk(ikpt1,1,idir)
         dtefield%idxkstr(ikstr,istr,idir) = ikpt2
         kpt_mark(ikpt2) = 1
       end do

     end do    ! istr

!    compute distance between strings
!    compute the metric matrix of the strings space in the direction idir
     do ii = 1,3
       do jj = 1,3
         if (ii<idir.and.jj<idir) dtefield%gmet_str(ii  ,jj  ,idir) = &
&         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
         if (ii<idir.and.jj>idir) dtefield%gmet_str(ii  ,jj-1,idir) = &
&         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
         if (ii>idir.and.jj<idir) dtefield%gmet_str(ii-1,jj  ,idir) = &
&         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
         if (ii>idir.and.jj>idir) dtefield%gmet_str(ii-1,jj-1,idir) = &
&         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
       end do
     end do
!    DEBUG
!    write(std_out,*)'gmet'
!    do ii=1,3
!    write(std_out,*)gmet(ii,:)
!    end do
!    write(std_out,*)'gmet_str'
!    do ii=1,2
!    write(std_out,*)dtefield%gmet_str(ii,:,idir)
!    end do
!    ENDDEBUG
     do istr = 1, dtefield%nstr(idir)
       do ii = 1,3
         if (ii<idir) dtefield%coord_str(ii,istr,idir)=dtefield%fkptns(ii,dtefield%idxkstr(1,istr,idir))
         if (ii>idir) dtefield%coord_str(ii-1,istr,idir)=dtefield%fkptns(ii,dtefield%idxkstr(1,istr,idir))
       end do
     end do

!    the following is very similar to getshell
     dist_ = 0._dp
     do ii = 1,2
       dist_ = dist_ + dtefield%gmet_str(ii,ii,idir)
     end do
     max_dist = 2._dp * dist_ * 2._dp

     dk_str(:,:,idir) = 0._dp
     last_dist = 0._dp
!    ishell = 0
!    dtefield%str_neigh(:,:,:) = 0
     dk_flag = 0
     do while (dk_flag /= 2)
!      Advance shell counter
!      ishell = ishell + 1

!      Search the smallest distance between two strings
       dist = max_dist
       do istr = 1,dtefield%nstr(idir)
         delta_str3(:) = dtefield%coord_str(:,1,idir) - dtefield%coord_str(:,istr,idir)
         do coord1 = -1,1  !two loop to search also on the border of the BZ
           do coord2 = -1,1
             dist_ = 0._dp
             dstr(:) = delta_str3(:) - nint(delta_str3(:))
             dstr(1) = dstr(1) + real(coord1,dp)
             dstr(2) = dstr(2) + real(coord2,dp)
             do ii = 1,2
               do jj = 1,2
                 dist_ = dist_ + dstr(ii)*dtefield%gmet_str(ii,jj,idir)*dstr(jj)
               end do
             end do
             if ((dist_ < dist).and.(dist_ - last_dist > tol8)) then
               dist = dist_
             end if
           end do
         end do
       end do

       last_dist = dist

!      search the connecting vectors for that distance
       do istr = 1,dtefield%nstr(idir)
         delta_str3(:) = dtefield%coord_str(:,istr,idir) - dtefield%coord_str(:,1,idir)
         do coord1 = -1,1
           do coord2 = -1,1
             dist_ = 0._dp
             dstr(:) = delta_str3(:) - nint(delta_str3(:))
             dstr(1) = dstr(1) + real(coord1,dp)
             dstr(2) = dstr(2) + real(coord2,dp)
             do ii = 1,2
               do jj = 1,2
                 dist_ = dist_ + dstr(ii)*dtefield%gmet_str(ii,jj,idir)*dstr(jj)
               end do
             end do
             if (abs(dist_ - dist) < tol8) then
               if (dk_flag == 0) then
                 dk_str(:,1,idir) = dstr(:)
                 dk_flag = 1
!                DEBUG
!                write(std_out,'(a,i4,2e15.4)')'1st connect', istr, dstr
!                ENDDEBUG
               elseif (dk_str(1,1,idir)*dstr(2)-dk_str(2,1,idir)*dstr(1) > tol8) then
                 dk_str(:,2,idir) = dstr(:)
                 dk_flag = 2
!                DEBUG
!                write(std_out,'(a,i4,2e15.4)')'2nd connect', istr, dstr
!                ENDDEBUG
                 exit
               end if
             end if
           end do
           if (dk_flag == 2) exit
         end do
         if (dk_flag == 2) exit
       end do

     end do ! do while

!    search the two neighbours for each string
     do istr = 1,dtefield%nstr(idir)
       dtefield%str_neigh(0,istr,idir) = istr
       dtefield%strg_neigh(0,istr,:,idir) = 0
       do jstr = 1,dtefield%nstr(idir)
         delta_str3(:) = dtefield%coord_str(:,jstr,idir) - dtefield%coord_str(:,istr,idir)
         do coord1 = -1,1
           do coord2 = -1,1
             dist_ = 0._dp
             dstr(:) = delta_str3(:) - nint(delta_str3(:))
             dstr(1) = dstr(1) + real(coord1,dp)
             dstr(2) = dstr(2) + real(coord2,dp)
             do ii = 1,2
               if (sum(abs(dstr(:)-dk_str(:,ii,idir)))<tol8) then
                 dtefield%str_neigh(ii,istr,idir) = jstr
                 dtefield%strg_neigh(ii,istr,1,idir) = coord1
                 dtefield%strg_neigh(ii,istr,2,idir) = coord2
               elseif (sum(abs(dstr(:)+dk_str(:,ii,idir)))<tol8) then
                 dtefield%str_neigh(-ii,istr,idir) = jstr
                 dtefield%strg_neigh(-ii,istr,1,idir) = coord1
                 dtefield%strg_neigh(-ii,istr,2,idir) = coord2
               end if
             end do
           end do
         end do
       end do
     end do

!    DEBUG
!    write(std_out,'(a,e15.4,e15.4,e15.4,e15.4)')'dk_str',dk_str(1,1,idir),dk_str(2,1,idir),dk_str(1,2,idir),dk_str(2,2,idir)
!    write(std_out,*)'istr, neigh1, strg(1,:), neigh2, strg(2,:),neigh-1, strg(-1,:), neigh-2, strg(-2,:)'
!    do istr=1,dtefield%nstr(idir)
!    write(std_out,'(13i4)')istr, &
!    &       dtefield%str_neigh(1,istr,idir), dtefield%strg_neigh(1,istr,:,idir),&
!    &       dtefield%str_neigh(2,istr,idir), dtefield%strg_neigh(2,istr,:,idir),&
!    &       dtefield%str_neigh(-1,istr,idir), dtefield%strg_neigh(-1,istr,:,idir),&
!    &       dtefield%str_neigh(-2,istr,idir), dtefield%strg_neigh(-2,istr,:,idir)
!    end do
!    ENDDEBUG


   end if         ! rfdir(idir) == 1

 end do           ! close loop over idir

 ABI_DEALLOCATE(kpt_mark)

!------------------------------------------------------------------------------
!------------ Compute PAW on-site terms if necessary --------------------------
!------------------------------------------------------------------------------

 if (usepaw == 1 .and. dtefield%has_expibi == 1) then
   call expibi(dtefield,gprimd,natom,rprimd,xred)
 end if

 if (usepaw == 1 .and. dtefield%has_qijb == 1) then
   call qijb_kk(dtefield,gprimd,natom,ntypat,pawang,pawrad,pawtab,typat)
   if (abs(dtset%berryopt)==5) then
     call set_twind(dtefield)
     call qijb_bk(dtefield,gprimd,natom,ntypat,pawang,pawrad,pawtab,typat)
   end if
 end if

 if (usepaw == 1 .and. dtefield%has_Lij == 1 .and. dtefield%has_Lijr3 == 1 ) then
   call Lij(dtefield,ntypat,pawrad,pawtab,psps)
 end if

 if (usepaw == 1 .and. dtefield%has_rij == 1) then
   c1=sqrt(four_pi/three)
   do itypat = 1, ntypat
     do klmn = 1, pawtab(itypat)%lmn2_size
       dtefield%rij(klmn,itypat,1) = c1*pawtab(itypat)%qijl(4,klmn) ! S_{1,1} ~ x
       dtefield%rij(klmn,itypat,2) = c1*pawtab(itypat)%qijl(2,klmn) ! S_{1,-1} ~ y
       dtefield%rij(klmn,itypat,3) = c1*pawtab(itypat)%qijl(3,klmn) ! S_{1,0} ~ z
     end do ! end loop over klmn
   end do ! end loop over itypat
   dtefield%has_rij = 2
 end if !

 if (usepaw == 1 .and. dtefield%has_twdij0==1 ) then
!  all three terms checked ok at dk = 0 
   write(std_out,'(a)')' making pawtwdij0 terms for magnetic field ... '

   dtefield%twdij0(:,:,:,:,:) = zero

   call pawtwdij_1(dtefield,gprimd,natom,ntypat,pawrad,pawtab,psps,typat)
   call pawtwdij_2b(dtefield,gprimd,natom,ntypat,pawrad,pawtab,psps,typat)
   call pawtwdij_2e(dtefield,gprimd,natom,ntypat,pawrad,pawtab,psps,typat)
   dtefield%has_twdij0 = 2

   write(std_out,'(a)')' done making pawtwdij0 terms for magnetic field ... '

 end if

 if (usepaw == 1 .and. dtefield%has_tweijkl == 1) then
!  all four terms checked perfect at dk = 0 
   write(std_out,'(a)')' making pawtweijkl terms for magnetic field ... '

   dtefield%tweijkl(:,:,:,:,:) = zero
   
   call pawtwdij_2a(dtefield,gprimd,natom,ntypat,pawrad,pawtab,typat)
   call pawtwdij_2c(dtefield,gprimd,natom,ntypat,pawrad,pawtab,typat)
   call pawtwdij_2d(dtefield,gprimd,natom,ntypat,pawrad,pawtab,typat)
   call pawtwdij_2f(dtefield,gprimd,natom,ntypat,pawrad,pawtab,typat)

   dtefield%has_tweijkl = 2

   write(std_out,'(a)')' done making pawtweijkl terms for magnetic field ... '

 end if


!------------------------------------------------------------------------------
!------------ Build the array pwind that is needed to compute the -------------
!------------ overlap matrices at k +- dk                         -------------
!------------------------------------------------------------------------------

 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 exchn2n3d = 0 ; istwf_k = 1 ; ikg1 = 0
 pwind(:,:,:) = 0
 pwnsfac(1,:) = 1.0_dp
 pwnsfac(2,:) = 0.0_dp
 ABI_ALLOCATE(kg1_k,(3,mpw))

 ipwnsfac = 0

 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

     dk(:) = dtefield%dkvecs(:,idir)

     do ifor = 1, 2

       if (ifor == 2) dk(:) = -1._dp*dk(:)

!      Build pwind and kgindex
!      NOTE: The array kgindex is important for parallel execution.
!      In case nsppol = 2, it may happent that a particular processor
!      treats k-points at different spin polarizations.
!      In this case, it is not possible to address the elements of
!      pwind correctly without making use of the kgindex array.

       ikg = 0 ; ikpt_loc = 0 ; isppol = 1
       do ikpt = 1, dtefield%fnkpt

         ikpti = dtefield%indkk_f2ibz(ikpt,1)
         nband_k = dtset%nband(ikpti)
         ikpt1f = dtefield%ikpt_dk(ikpt,ifor,idir)
         ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

         if (mpi_enreg%paral_compil_kpt == 1) then

           if ((minval(abs(mpi_enreg%proc_distrb(ikpti,1:nband_k,1)-me))/=0).and.&
&           (minval(abs(mpi_enreg%proc_distrb(ikpti,1:nband_k,nsppol)-me))/=0)) then
!            if (minval(abs(mpi_enreg%proc_distrb(ikpti,1:nband_k,1:nsppol)-me))/=0) then
             cycle
           end if

           ikpt_loc = ikpt_loc + 1

         end if

!        Build basis sphere of plane waves for the nearest neighbour of
!        the k-point (important for MPI //)

         kg1_k(:,:) = 0
         kpt1(:) = dtset%kptns(:,ikpt1i)
         call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg1_k,kpt1,&
&         1,mpi_enreg,mpw,npw_k1)
         me_g0=mpi_enreg%me_g0


!        ji: fkgindex is defined here !
         dtefield%fkgindex(ikpt) = ikg

!        
!        Deal with symmetry transformations
!        

!        bra k-point k(b) and IBZ k-point kIBZ(b) related by
!        k(b) = alpha(b) S(b)^t kIBZ(b) + G(b)
!        where alpha(b), S(b) and G(b) are given by indkk_f2ibz
!        
!        For the ket k-point:
!        k(k) = alpha(k) S(k)^t kIBZ(k) + G(k) - GBZ(k)
!        where GBZ(k) takes k(k) to the BZ
!        

         isym  = dtefield%indkk_f2ibz(ikpt,2)
         isym1 = dtefield%indkk_f2ibz(ikpt1f,2)

!        Construct transformed G vector that enters the matching condition:
!        alpha(k) S(k)^{t,-1} ( -G(b) - GBZ(k) + G(k) )

         dg(:) = -dtefield%indkk_f2ibz(ikpt,3:5) &
&         -nint(-dtefield%fkptns(:,ikpt) - dk(:) - tol10 + &
&         dtefield%fkptns(:,ikpt1f)) &
&         +dtefield%indkk_f2ibz(ikpt1f,3:5)

!        old code
!        iadum(:)=0
!        do idum=1,3
!        iadum(:)=iadum(:)+ symrec(:,idum,isym1)*dg(idum)
!        end do

!        new code
         iadum(:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),dg(:))

         dg(:) = iadum(:)

         if ( dtefield%indkk_f2ibz(ikpt1f,6) == 1 ) dg(:) = -dg(:)

!        Construct S(k)^{t,-1} S(b)^{t}

         dum33(:,:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),symrec(:,:,isym))

!        Construct alpha(k) alpha(b)

         if (dtefield%indkk_f2ibz(ikpt,6) == dtefield%indkk_f2ibz(ikpt1f,6)) then
           itrs=0
         else
           itrs=1
         end if


         npw_k  = npwarr(ikpti)
!        npw_k1 = npwarr(ikpt1i)

!        loop over bra G vectors
         do ipw = 1, npw_k

!          NOTE: the bra G vector is taken for the sym-related IBZ k point,
!          not for the FBZ k point
           iadum(:) = kg(:,dtefield%kgindex(ikpti) + ipw)

!          Store non-symmorphic operation phase factor exp[i2\pi \alpha G \cdot t]

           if ( ipwnsfac == 0 ) then
!            old code
             rdum=0.0_dp
             do idum=1,3
               rdum=rdum+dble(iadum(idum))*dtset%tnons(idum,isym)
             end do
             rdum=two_pi*rdum
             if ( dtefield%indkk_f2ibz(ikpt,6) == 1 ) rdum=-rdum
             pwnsfac(1,ikg+ipw) = cos(rdum)
             pwnsfac(2,ikg+ipw) = sin(rdum)
!            
!            new code
!            rdum = DOT_PRODUCT(dble(iadum(:)),dtset%tnons(:,isym))
!            rdum= two_pi*rdum
!            if ( dtefield%indkk_f2ibz(ikpt,6) == 1 ) rdum=-rdum
!            pwnsfac(1,ikg+ipw) = cos(rdum)
!            pwnsfac(2,ikg+ipw) = sin(rdum)

           end if

!          to determine r.l.v. matchings, we transformed the bra vector
!          Rotation
           iadum1(:)=0
           do idum1=1,3
             iadum1(:)=iadum1(:)+dum33(:,idum1)*iadum(idum1)
           end do
           iadum(:)=iadum1(:)
!          Time reversal
           if (itrs==1) iadum(:)=-iadum(:)
!          Translation
           iadum(:) = iadum(:) + dg(:)

           do jpw = 1, npw_k1
             iadum1(1:3) = kg1_k(1:3,jpw)
             if ( (iadum(1) == iadum1(1)).and. &
&             (iadum(2) == iadum1(2)).and. &
&             (iadum(3) == iadum1(3)) ) then
               pwind(ikg + ipw,ifor,idir) = jpw
!              write(std_out,'(a,2x,3i4,2x,i4)') 'Found !:',iadum1(:),jpw
               exit
             end if
           end do
         end do

         ikg  = ikg + npw_k

       end do    ! close loop over ikpt

       ipwnsfac = 1

     end do    ! close loop over ifor

   end if      ! rfdir(idir) == 1

 end do        ! close loop over idir

!Build mpi_enreg%kptdstrb
!array required to communicat the WFs between cpus in berryphase_new.f
!(MPI // over k-points)

 if (mpi_enreg%paral_compil_kpt == 1) then
   do idir = 1, 3
     if (dtset%rfdir(idir) == 1) then
       do ifor = 1, 2

         ikpt_loc = 0
         do isppol = 1, nsppol

           do ikpt = 1, dtefield%fnkpt

             ikpti = dtefield%indkk_f2ibz(ikpt,1)
             nband_k = dtset%nband(ikpti)
             ikpt1f = dtefield%ikpt_dk(ikpt,ifor,idir)
             ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

             if (minval(abs(mpi_enreg%proc_distrb(ikpti,1:nband_k,isppol)-me))/=0) then
               cycle
             end if

             ikpt_loc = ikpt_loc + 1
             mpi_enreg%kptdstrb(me + 1,ifor+2*(idir-1),ikpt_loc) = &
&             ikpt1i + (isppol - 1)*nkpt

             mpi_enreg%kptdstrb(me+1,ifor+2*(idir-1),&
&             ikpt_loc+dtefield%fmkmem_max*nsppol) = &
&             ikpt1f + (isppol - 1)*dtefield%fnkpt

           end do   ! ikpt
         end do     ! isppol
       end do       ! ifor
     end if         ! dtset%rfdir(idir) == 1
   end do           ! idir
 end if             ! mpi_enreg%paral_compil_kpt == 1

!build mpi_enreg%kpt_loc2fbz_sp and mpi_enreg%fmkmem
 ikpt_loc = 0
 do isppol = 1, nsppol
   do ikpt = 1, dtefield%fnkpt

     ikpti = dtefield%indkk_f2ibz(ikpt,1)
     nband_k = dtset%nband(ikpti)

     if(mpi_enreg%paral_compil_kpt == 1)then
       if (minval(abs(mpi_enreg%proc_distrb(ikpti,1:nband_k,isppol)-me))/=0) cycle
     end if

     ikpt_loc = ikpt_loc + 1

     mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 1) = ikpt
     mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 2) = isppol

   end do
 end do


!parallel case only :
!build mpi_enreg%kpt_loc2ibz_sp, dtefield%cgqindex and dtefield%nneigh
 if ((fieldflag).and.(mpi_enreg%paral_compil_kpt == 1)) then
   ikpt_loc = 0
   do isppol = 1, nsppol
     do ikpt = 1, nkpt

       ikptf = dtefield%i2fbz(ikpt)
       nband_k = dtset%nband(ikpti)

       neigh(:) = 0 ; icg = 0 ; ikg = 0 ; flag_kpt = 0; icprj = 0
       do idir=1, 3

!        skip idir values for which efield_dot(idir) = 0
         if (abs(dtefield%efield_dot(idir)) < tol12 .and. dtset%berryopt == 4) cycle

         do ifor = 1, 2

           flag = 0

           ikpt1f = dtefield%ikpt_dk(ikptf,ifor,idir)
           ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

           dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ikg
           ikg = ikg + npwarr(ikpt1i)

!          check if this neighbour is also a previous neighbour
           do ineigh = 1, (ifor+2*(idir-1))
             if (neigh(ineigh) == ikpt1i) then
               flag = 1
               dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ineigh
               dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
&               dtefield%cgqindex(2,ineigh,ikpt+(isppol-1)*nkpt)
               exit
             end if
           end do
!          create the cgqindex of the neighbour if necessary
           if (flag == 0) then
             neigh(ifor+2*(idir-1)) = ikpt1i
             dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
&             ifor+2*(idir-1)
             dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = icg
             if (isppol == 1) dtefield%nneigh(ikpt) = dtefield%nneigh(ikpt) + 1
             icg = icg + npwarr(ikpt1i)*nband_k
           end if

!          if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) then
!          !ikpt is not one of my kpt_loc
!          cycle
!          end if

         end do !ifor
       end do !idir

       if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))==0) then
!        ikpt is one of my kpt_loc
         ikpt_loc = ikpt_loc + 1

         mpi_enreg%kpt_loc2ibz_sp(me, ikpt_loc, 1) = ikpt
         mpi_enreg%kpt_loc2ibz_sp(me, ikpt_loc, 2) = isppol
       end if

     end do !ikpt

   end do !isppol

 end if !mpi_enreg%paral_compil_kpt == 1

!should be temporary
!unassigned mpi_enreg%kpt_loc2fbz_sp are empty ; inform other cpu (there are better ways...)
 mpi_enreg%fmkmem(me) = dtefield%fmkmem
 mpi_enreg%mkmem(me) = mkmem
!do ii=ikpt_loc+1,dtefield%fmkmem_max
!mpi_enreg%kpt_loc2fbz_sp(me, ii, 1) = -1
!end do


!Build mpi_enreg%kptdstrbi
!(same as mpi_enreg%kptdstrb but for k-points in the iBZ),
!dtefield%cgqindex and dtefield%nneigh

 if ((fieldflag).and.(mpi_enreg%paral_compil_kpt == 1)) then

   ikpt_loc = 1
   do isppol = 1, nsppol
     do ikpt = 1, nkpt

       nband_k = dtset%nband(ikpt)
       ikptf = dtefield%i2fbz(ikpt)

       neigh(:) = 0 ; icg = 0 ; ikg = 0 ; flag_kpt = 0; icprj = 0
       do idir = 1, 3

!        skip idir values for which efield_dot(idir) = 0
         if (abs(dtefield%efield_dot(idir)) < tol12 .and. dtset%berryopt == 4) cycle

         do ifor = 1, 2

           ikpt1f = dtefield%ikpt_dk(ikptf,ifor,idir)
           ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

!          dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ikg
           ikg = ikg + npwarr(ikpt1i)

           flag = 0
           do ineigh = 1, (ifor+2*(idir-1))
             if (neigh(ineigh) == ikpt1i) then
               flag = 1
!              dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ineigh
!              dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
!              &               dtefield%cgqindex(2,ineigh,ikpt+(isppol-1)*nkpt)
               exit
             end if
           end do
           if (flag == 0) then
!            neigh(ifor+2*(idir-1)) = ikpt1i
!            dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
!            &             ifor+2*(idir-1)
!            dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = icg
!            if (isppol == 1) dtefield%nneigh(ikpt) = dtefield%nneigh(ikpt) + 1
!            icg = icg + npwarr(ikpt1i)*dtset%nspinor*nband_k
           end if

           if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) then
             cycle
           end if

           mpi_enreg%kptdstrbi(me + 1,ifor+2*(idir-1),&
&           ikpt_loc + dtefield%mkmem_max*nsppol) = &
&           ikpt1f + (isppol - 1)*dtefield%fnkpt

           flag_kpt = 1
           if (flag == 0) then
             mpi_enreg%kptdstrbi(me + 1,ifor+2*(idir-1),ikpt_loc) = &
&             ikpt1i + (isppol - 1)*nkpt
           end if

!          MVeithen: the if condition allows to avoid that the same wavefunction
!          is send several times to a particular cpu

         end do    ! ifor
       end do    ! idir

       if (flag_kpt == 1) ikpt_loc = ikpt_loc + 1

     end do    ! ikpt
   end do    ! isppol

 end if   ! fieldflag and paral_compil_kpt

 if (mpi_enreg%paral_compil_kpt == 1) then

   count = nproc*6*dtefield%fmkmem_max*nsppol*2
   ABI_ALLOCATE(buffer,(count))
   buffer(:) = reshape(mpi_enreg%kptdstrb(:,:,:),(/count/))
   call xsum_mpi(buffer,spaceComm,ierr)
   mpi_enreg%kptdstrb(:,:,:) = reshape(buffer(:),&
&   (/nproc,6,dtefield%fmkmem_max*nsppol*2/))
   ABI_DEALLOCATE(buffer)

   count = nproc*dtefield%fmkmem_max*nsppol*2
   ABI_ALLOCATE(buffer,(count))
   buffer(:) = reshape(mpi_enreg%kpt_loc2fbz_sp(:,:,:),(/count/))
   call xsum_mpi(buffer,spaceComm,ierr)
   mpi_enreg%kpt_loc2fbz_sp(0:nproc-1,:,:) = reshape(buffer(:),&
&   (/nproc,dtefield%fmkmem_max*nsppol,2/))
   ABI_DEALLOCATE(buffer)

   call xsum_mpi(mpi_enreg%fmkmem,spaceComm,ierr)

   if (fieldflag) then
     count = nproc*6*dtefield%mkmem_max*nsppol*2
     ABI_ALLOCATE(buffer,(count))
     buffer(:) = reshape(mpi_enreg%kptdstrbi(:,:,:),(/count/))
     call xsum_mpi(buffer,spaceComm,ierr)
     mpi_enreg%kptdstrbi(:,:,:) = reshape(buffer(:),&
&     (/nproc,6,dtefield%mkmem_max*nsppol*2/))
     ABI_DEALLOCATE(buffer)


     count = nproc*nsppol*dtefield%mkmem_max*2
     ABI_ALLOCATE(buffer,(count))
     buffer(:) = reshape(mpi_enreg%kpt_loc2ibz_sp(:,:,:),(/count/))
     call xsum_mpi(buffer,spaceComm,ierr)
     mpi_enreg%kpt_loc2ibz_sp(0:nproc-1,:,:) = reshape(buffer(:),&
&     (/nproc,nsppol*dtefield%mkmem_max,2/))
     ABI_DEALLOCATE(buffer)

     call xsum_mpi(mpi_enreg%mkmem,spaceComm,ierr)
   end if

 end if



!------------------------------------------------------------------------------
!------------------------ Estimate critical field -----------------------------
!------------------------------------------------------------------------------

!Compute the minimal value of the bandgap required to be below
!the critical field as defined by the relation
!| E_i*a_i | < E_g/n_i

 if (dtset%berryopt == 4) then

   do idir = 1, 3
     eg_dir(idir) = abs(dtefield%efield_dot(idir))*dtefield%nkstr(idir)
   end do
   eg = maxval(eg_dir)
   eg_ev = eg*Ha_eV

   write(message,'(a,a,a,a,a,a,a,a,f7.2,a,a)')ch10,&
&   ' initberry: COMMENT - ',ch10,&
&   '  As a rough estimate,',ch10,&
&   '  to be below the critical field, the bandgap of your system',ch10,&
&   '  should be larger than ',eg_ev,' eV.',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

 end if

 ABI_DEALLOCATE(kg1_k)
!!! deallocate(dum_fwtk) !! by MM

!DEBUG
!write(std_out,*)'initberry: exit'
!ENDDEBUG

end subroutine initberry
!!***
