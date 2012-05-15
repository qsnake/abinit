!{\src2tex{textfont=tt}}
!!****f* ABINIT/elphon
!!
!! NAME
!! elphon
!!
!! FUNCTION
!! This routine extracts the electron phonon coupling matrix
!! elements and calculates related properties - Tc, phonon linewidths...
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   anaddb_dtset=dataset with input variables
!!     anaddb_dtset%a2fsmear = smearing for alpha2F function
!!     anaddb_dtset%brav = type of Bravais lattice
!!     anaddb_dtset%dipdip  =dipole dipole interaction flag
!!     anaddb_dtset%elphsmear = smearing width for gaussian integration
!!           or buffer in energy for calculations with tetrahedra (telphint=0)
!!     anaddb_dtset%elph_fermie = input value of Fermi energy
!!           0 means use value from wfk file
!!     anaddb_dtset%enunit = governs the units to be used for the output of
!!           the phonon frequencies and e-ph quantities
!!     anaddb_dtset%gkk2write= flag to write out gkk2 matrix elements to disk
!!     anaddb_dtset%gkk_rptwrite= flag to write out real space gkk_rpt matrix elements to disk
!!     anaddb_dtset%gkqwrite= flag to write out gkq matrix elements to disk
!!     anaddb_dtset%ep_b_min= first band taken into account in FS integration (if telphint==2)
!!     anaddb_dtset%ep_b_max= last band taken into account in FS integration (if telphint==2) 
!!     anaddb_dtset%prtfsurf = integer flag for the output of the Fermi surface (XCrysden file format)
!!     anaddb_dtset%prtnest = integer flag for the calculation of the nesting function
!!     anaddb_dtset%ifcflag = flag for IFC matrices in anaddb calling routine
!!           the IFCs are presumed to be known!
!!     anaddb_dtset%ifltransport= flag for transport properties (no=0: yes=1 )
!!     anaddb_dtset%kptrlatt=kpoint grid generating vectors, as in abinit
!!     anaddb_dtset%kptrlatt_fine=kpoint grid generating vectors, for fine grid used in FS integration
!!     anaddb_dtset%mustar = parameter for Coulombic pseudo-potential in McMillan T_c calculation
!!     anaddb_dtset%ngqpt(3)=integers defining the number of points in the qpt sampling
!!     anaddb_dtset%nqpath=number of vertices in the path in reciprocal space, for band structure
!!           and phonon linewidth output
!!     anaddb_dtset%nqshft= number of shift vectors for defining the sampling of q points
!!     anaddb_dtset%ntemper = number of temperature points to calculate, from tempermin to 
!!           tempermin+ntemper*temperinc
!!     anaddb_dtset%qpath=vertices in the path in reciprocal space, for band structure
!!           and phonon linewidth output
!!     anaddb_dtset%q1shft(3,4) =qpoint shifts considered
!!     anaddb_dtset%telphint = flag for integration over the FS with 0=tetrahedra 1=gaussians
!!     anaddb_dtset%tempermin = minimum temperature at which resistivity etc are calculated (in K)
!!     anaddb_dtset%temperinc = interval temperature grid on which resistivity etc are calculated (in K)
!!     anaddb_dtset%ep_keepbands = flag to keep gamma matrix dependence on electronic bands
!!
!! filnam(7)=character strings giving file names
!!   acell_in(3)= input length scales of cell (bohr)
!!   amu(ntypat)=mass of the atoms (atomic mass unit)
!!   atmfrc  = inter-atomic force constants from anaddb
!!   dielt(3,3) = dielectric tensor
!!   dyewq0(3,3,natom)=atomic self-interaction correction to the dynamical matrix (only when anaddb_dtset%dipdip=1)
!!   gmet(3,3) =metric in reciprocal space (telphint=1)
!!   gprim(3,3) =dimensionless basis vectors of reciprocal space
!!   indsym = mapping of atoms btw themselves under symmetry
!!   mpert =maximum number of ipert
!!   natom=number of atoms in cell
!!   nrpt =number of real space points used to integrate IFC (for interpolation of dynamical matrices)
!!   nsym=number of space group symmetries
!!   ntypat = number of types of atoms
!!   rcan(3,natom) =canonical positions of atoms
!!   rmet(3,3)=metric tensor in real space (bohr^2)
!!   rprim_in(3,3)= input primitive translation vectors
!!   rpt(3,nprt) =canonical positions of R points in the unit cell
!!   symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!   symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!   tnons(3,nsym)=fractional nonsymmorphic translations
!!   trans(3,natom) = Atomic translations : xred = rcan + trans
!!   typat(natom)=type integer for each atom in cell
!!   ucvol=unit cell volume in bohr**3
!!   wghatm(natom,natom,nrpt) =Weight for the pair of atoms and the R vector
!!   xred(3,natom)=fractional dimensionless atomic coordinates
!!   zeff(3,3,natom) =effective charge on each atom, versus electric field and atomic displacement
!!
!! OUTPUT
!!
!! NOTES
!!  inspired to a large extent by epcouple.f from the DecAFT package by J. Kay Dewhurst
!!  most inputs taken from mkifc.f
!!  in anaddb anaddb_dtset%ifcflag must be 1 such that the IFC are calculated in atmfrc prior to calling elphon
!!
!!  brav not taken into account propely in all of the code. (MG?)
!!
!!  could choose to make a full 3 dimensional kpt array (:,:,:). Easier for many operations
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      a2f_dump,a2f_free,a2f_init,bst_init_from_hdr,bstruct_clean
!!      clean_phon_ds,complete_gamma,complete_gamma_tr,copy_kptrank
!!      destroy_crystal,destroy_phondos,eliashberg_1d,elph_ds_clean
!!      elph_ds_nullify,elph_tr_ds_clean,elph_tr_ds_nullify,ep_fs_weights
!!      ep_setupqpt,ftgam,gamma_free,gamma_init,gamma_interp_setup,gamma_linwid
!!      gamma_nullify,get_all_gkk2,get_all_gkq,get_all_gkr,get_fs_bands
!!      get_rank_1kpt,get_veloc_tr,hdr_clean,hdr_io,hdr_nullify,init_crystal
!!      integrate_gamma,integrate_gamma_alt,integrate_gamma_tr,matr3inv,mka2f
!!      mka2f_tr,mka2fqgrid,mkfskgrid,mknesting,mkph_linwid,mkphdos,mkqptequiv
!!      mkrdim,order_fs_kpts,outelph,print_phondos,printbxsf,rchkgsheader
!!      setup_phon_ds,timein,update_occ,wrap2_pmhalf,wrtout
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine elphon(anaddb_dtset,filnam,acell_in,amu,atmfrc,dielt,dyewq0,gmet,&
&  gprim,indsym,mpert,mpi_enreg,natom,nrpt,nsym,ntypat,rcan,rmet,rprim_in,rpt,&
&  symrec,symrel,tnons,trans,typat,ucvol,wghatm,xred,zeff)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon
 use m_kptrank
 use m_errors   

 use m_io_tools,        only : get_unit
 use m_header,          only : hdr_clean, hdr_nullify
 use m_crystal,         only : crystal_structure, init_crystal, destroy_crystal
 !use m_crystal_io,      only : init_crystal_from_hdr
 use m_ebands,          only : bst_init_from_hdr, update_occ, bstruct_clean !bandstructure_type,
 use m_phdos,           only : phonon_dos_type, destroy_phondos, print_phondos, mkphdos !init_phondos, 
 use m_gamma,           only : gamma_t, gamma_nullify, gamma_free, gamma_init, gamma_interp_setup, gamma_linwid, &
&                              a2f_t, a2f_init, a2f_dump, a2f_free
 use m_header,          only : hdr_get_nelect_byocc

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elphon'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_62_occeig
 use interfaces_77_ddb, except_this_one => elphon
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpert,natom
 integer,intent(in) :: nrpt,nsym,ntypat
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(anaddb_dataset_type) :: anaddb_dtset
!arrays
 integer,intent(in) :: indsym(4,nsym,natom)
 integer,intent(in) :: symrec(3,3,nsym),symrel(3,3,nsym),typat(natom)
 real(dp),intent(in) :: acell_in(3),amu(ntypat) 
 real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt) ! inout due to gtdyn9!
 real(dp),intent(inout) :: dyewq0(3,3,natom)              ! inout due to gtdyn9
 real(dp),intent(in) :: dielt(3,3)
 real(dp),intent(in) :: gmet(3,3),gprim(3,3)
 real(dp),intent(in) :: rcan(3,natom),rmet(3,3),rprim_in(3,3)
 real(dp),intent(in) :: rpt(3,nrpt),tnons(3,nsym),trans(3,natom)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt),xred(3,natom),zeff(3,3,natom)
 character(len=fnlen),intent(in) :: filnam(7) 

!Local variables-------------------------------
!scalars
 integer,parameter :: timrev2=2,space_group0=0
 integer :: ikpt_phon,ikpt_fine,ierr,unitgkk
 integer :: iband,ibandp, ieliash,ii,ikpt
 integer :: iqpt,isppol,istat,n1wf,nband
 integer :: neliash,onegkksize,ios
 integer :: timrev,unitfskgrid,qtor,idir
 integer :: iFSkpq, symrankkpt
 integer :: ep_prt_wtk ! eventually to be made into an input variable
 integer :: a2f_nqshift
 integer :: rdwr, fform
 real(dp) :: max_occ,realdp_ex,res,ss
 real(dp) :: tcpu, twall, tcpui, twalli
 logical :: make_gkk2,use_afm,use_tr,use_antiferro
 character(len=500) :: message
 character(len=fnlen) :: fname,elph_base_name,ddkfilename,gkk_fname
 character(len=fnlen) :: nestname
 type(elph_tr_type) :: elph_tr_ds
 type(elph_type) :: elph_ds
 type(hdr_type) :: hdr,hdr1
 type(phon_type) :: phon_ds
 type(crystal_structure) :: Cryst
 type(bandstructure_type) :: Bst
 type(phonon_dos_type) :: PHdos
 type(gamma_t) :: Gam
 type(a2f_t) :: A2f
!arrays
 integer :: qptrlatt(3,3),a2f_qptrlatt(3,3)
 integer,allocatable :: FSfullpqtofull(:,:)
 integer,allocatable :: irredtoGS_phon(:)
 integer,allocatable :: irredtoGS_fine(:)
 integer,allocatable :: qpttoqpt(:,:,:)
 real(dp) :: acell(3),gprimd(3,3),kpt(3),rprim(3,3),rprimd(3,3),shiftk(3)
 real(dp),allocatable :: a2f_1d(:),delta(:,:),dos_phon(:)
 real(dp),allocatable :: eigenGS(:,:,:)
 real(dp),allocatable :: eigenGS_fine(:,:,:)
 real(dp),allocatable :: zz(:,:)
 real(dp),allocatable :: a2f_qshift(:,:)

! *************************************************************************

 write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
& ' Properties based on electron-phonon coupling ',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')


 call timein(tcpui,twalli)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-begin elphon at tcpu',tcpui,'  and twall',twalli,' sec'
 call wrtout(std_out,message,'COLL')

!==================================
!Initialization of some variables
!==================================

 gkk_fname = filnam(5)
 unitgkk   = get_unit()
 open (unit=unitgkk,file=gkk_fname,form='unformatted',status='old',iostat=ios)
 ABI_CHECK(ios==0,"Opening "//TRIM(gkk_fname))

 elph_base_name=trim(filnam(2))//"_ep"
 ddkfilename=trim(filnam(7))

!use time reversal symmetry always when possible for kpoint reduction,
!and suppose it has been used in WF generation
!not used for the moment: values are always taken from input files.
 timrev = 1

 call elph_ds_nullify(elph_ds)

 call elph_tr_ds_nullify(elph_tr_ds)

 elph_ds%mustar       = anaddb_dtset%mustar        ! input mustar
 elph_ds%nbranch      = 3*natom                    ! number of phonon modes = 3 * natom
 elph_ds%ep_keepbands = anaddb_dtset%ep_keepbands  ! flag to sum over bands
 elph_ds%a2fsmear     = anaddb_dtset%a2fsmear      ! smearing for Eliashber functions
 elph_ds%tuniformgrid = 1
 elph_ds%na2f         = 400                        !maximum number of Matsubara frequencies.
!The precise number used depends on the value of Tc:
!they span $w_n = (2n+1) \pi T_c$  where $abs(w_n) < w_{cutoff}$
!ie $|n| < n_{cutoff} = ( \frac{w_{cutoff}}{\pi T_c} ) / 2$

!save gkk data for full kpoints to file on disk

 elph_ds%gkqwrite     = anaddb_dtset%gkqwrite
 elph_ds%gkk_rptwrite = anaddb_dtset%gkk_rptwrite
 elph_ds%gkk2write    = anaddb_dtset%gkk2write

!This should never be turned off: symmetrization of elphon matrix elements in complete_gkk. See get_all_gkq
 elph_ds%symgkq=anaddb_dtset%symgkq

 elph_ds%elph_base_name = trim(elph_base_name)

!normalize input rprim and acell.
 do ii=1,3
   ss = sqrt(rprim_in(1,ii)**2+rprim_in(2,ii)**2+rprim_in(3,ii)**2)
   rprim(:,ii) = rprim_in(:,ii)/ss
   acell(ii) = acell_in(ii) * ss
 end do

!make dimension-ful rprimd and gprimd for transformation of derivatives to cartesian coordinates.
 call mkrdim(acell,rprim,rprimd)
 call matr3inv(rprimd,gprimd)

 call hdr_nullify(hdr)
 call hdr_nullify(hdr1)

!===================
!Check some inputs
!===================

 if (nsym==1) then
   write (message,'(7a)')ch10,&
&   ' elphon: COMMENT- ',ch10,&
&   ' Symmetries are not used! ',ch10,&
&   ' Full matrix elements must be supplied for all perturbations and qpoints!',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   if ( ANY( ABS(tnons(:,1)) > tol10) ) then
     MSG_ERROR('nsym==1 but the symmetry is not the identity')
   end if
 end if

 if (anaddb_dtset%ifcflag/=1) then
   write(message,'(a,i0)')&
&   ' ifcflag should be set to 1 since the IFC matrices are supposed to exist but ifcflag= ',anaddb_dtset%ifcflag
   MSG_ERROR(message)        
 end if

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon begin setup after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

!=================================
!Set up the full grid of qpoints
!=================================
 call ep_setupqpt(anaddb_dtset,elph_ds,gmet,nsym,qptrlatt,rprimd,symrec,symrel,timrev)

!====================================
!Read the GS header of the GKK file
!this will give the phon grid of k
!and the Fermi surface integration weights
!====================================
 call wrtout (std_out,' elphon: reading and checking the GS header of the GKK file','COLL')

 call rchkGSheader(hdr,natom,nband,unitgkk)

!TODO Fix problem with different convention for the treatment of time-reversal symmetry.
!abinit uses timrev==(0,1), Cryst uses (1,2). Here time reversal symmetry is always assumed.

 use_antiferro=(Hdr%nspden==2.and.Hdr%nsppol==1)

 call init_crystal(Cryst,space_group0,natom,hdr%npsp,ntypat,nsym,rprimd,typat,xred,&
& hdr%zionpsp,hdr%znuclpsp,timrev2,use_antiferro,.false.,hdr%title,&
& symrel,tnons,hdr%symafm) ! Optional

!The values of rprimd and gprimd stored in the Hdr slighty differ from the input arguments
!that are read from the DDB, likely due to the use of formatted IO. For the time being
!Ctyst is initialized from the input arguments!
!call init_crystal_from_hdr(Cryst,Hdr,timrev2)
 ABI_CHECK(Cryst%nsym==nsym,"Wrong nsym")
 ABI_CHECK(ALL(Cryst%symrel==symrel),"Wrong symrel")
 ABI_CHECK(ALL(Cryst%symrec==symrec),"Wrong symrel")
 ABI_CHECK(ALL(Cryst%indsym==indsym),"Wrong indsym")

 if (ANY( ABS(Cryst%gprimd - gprimd) > tol16 ) .or.&
& ANY( ABS(Cryst%rprimd - rprimd) > tol16 ) .or.&
& ANY( ABS(Cryst%tnons -  tnons)  > tol16 )     &
& ) then
   write(std_out,*)" DIFFERENCE"
   write(std_out,*)"rprimd ",Cryst%rprimd-rprimd
   write(std_out,*)"gprimd ",Cryst%gprimd-gprimd
   write(std_out,*)"tnons  ",Cryst%tnons-tnons
   
   write(std_out,*)" VALUES "
   write(std_out,*)Cryst%rprimd
   write(std_out,*)rprimd
   write(std_out,*)Cryst%gprimd
   write(std_out,*)gprimd
!  MSG_WARNING("Small difference in rprimd, gprimd")
   MSG_ERROR("Fatal error")
 end if

 elph_ds%nsppol =hdr%nsppol
 elph_ds%nspinor=hdr%nspinor

!in spinor or spin polarized case, orbitals have occupation <= 1 instead of 2
 max_occ = one
 if (hdr%nspinor == 2) max_occ = half ! this accounts for the doubling of the num of bands, even though spin channels are not well defined
 if (elph_ds%nsppol > 1) max_occ = one
 write (std_out,*) ' max_occ factor  ', max_occ

 elph_ds%occ_factor = one
 if (hdr%nspinor == 1 .and. hdr%nsppol == 1) then
   elph_ds%occ_factor = one
 else if (hdr%nspinor == 2) then
   elph_ds%occ_factor = two
 else if (hdr%nsppol == 2) then
   elph_ds%occ_factor = one
 end if

!==================================================
!Read GS eigenvalues for each irreducible kpt and
!number of 1WF files contributing to the GKK file
!==================================================

 ABI_ALLOCATE(eigenGS,(nband,hdr%nkpt,elph_ds%nsppol))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out-of-memory in eigenGS')

 do isppol=1,elph_ds%nsppol
   do ikpt=1,hdr%nkpt
     read(unitgkk) eigenGS(:,ikpt,isppol)
   end do
 end do

!read number of 1WF files contributing to the GKK file
 read(unitgkk) n1wf
 write(message,'(a,i0)')' elphon : number of perturbations in the gkk file = ',n1wf
 call wrtout(std_out,message,'COLL')


!==================================================
!Set elph_ds%fermie: either comes from anaddb input file or from wfk file
!==================================================
 elph_ds%fermie = hdr%fermie
 elph_ds%nelect = hdr_get_nelect_byocc(Hdr)
 if (abs(anaddb_dtset%elph_fermie) > tol10) then
   elph_ds%fermie = anaddb_dtset%elph_fermie
   write(message,'(a,f12.6)')' Fermi level set by the user at :',elph_ds%fermie
   call bst_init_from_hdr(Bst,Hdr,nband,eigenGS)
 else if (abs(anaddb_dtset%ep_extrael) > tol10) then
   write(message,'(a,f12.6)')' Additional electrons per unit cell set by the user at :',&
&   anaddb_dtset%ep_extrael
   elph_ds%nelect = elph_ds%nelect + anaddb_dtset%ep_extrael
   call bst_init_from_hdr(Bst,Hdr,nband,eigenGS,nelect=elph_ds%nelect)
!  call bst_init_from_hdr(Bst,Hdr,nband,eigenGS,elph_ds%nelect)
 else
   call bst_init_from_hdr(Bst,Hdr,nband,eigenGS)
 end if
 call wrtout(std_out,message,'COLL')

!set BSt to use FD occupations:
 BSt%occopt = 3

!will need to be changed as well to get low T reference point
 BSt%tsmear = 0.00001_dp ! is this small enough?
!BSt%tsmear = 0.00285_dp ! 900 K

!Calculate occupation numbers.
 call update_occ(BSt,-99.99_dp)
 write(message,'(a,f12.6)')' Fermi level is now calculated to be :',BSt%fermie
 call wrtout(std_out,message,'COLL')
 
 if (abs(anaddb_dtset%ep_extrael) > tol10) then
   elph_ds%fermie = BSt%fermie
 end if

!====================================================================
!Setup of the phon k-grid :
!1) get bands near Ef
!====================================================================

 call get_fs_bands(eigenGS,hdr,elph_ds%fermie,anaddb_dtset%ep_b_min, anaddb_dtset%ep_b_max,&
& elph_ds%minFSband,elph_ds%maxFSband,elph_ds%k_phon%nkptirr)

 elph_ds%nFSband = elph_ds%maxFSband - elph_ds%minFSband + 1


 if (elph_ds%ep_keepbands == 0) then !we are summing over bands
   elph_ds%ngkkband = 1
 else if (elph_ds%ep_keepbands == 1) then
!  keep the band dependency btw elph_ds%minFSband and elph_ds%maxFSband
   elph_ds%ngkkband = elph_ds%nFSband
 else
   write(message,'(a,i0)')' ep_keepbands must be 0 or 1 while it is: ',elph_ds%ep_keepbands
   MSG_BUG(message)
 end if

 write(message,'(a,i0,2x,i0)')' elphon : minFSband, maxFSband = ',elph_ds%minFSband,elph_ds%maxFSband
 call wrtout(std_out,message,'COLL')

 ABI_ALLOCATE(elph_ds%k_phon%kptirr,(3,elph_ds%k_phon%nkptirr))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"allocating elph_ds%k_phon%kptirr")

 ABI_ALLOCATE(irredtoGS_phon,(elph_ds%k_phon%nkptirr))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"allocating irredtoGS_phon")

!====================================================================
!2) order irred k-points 
!====================================================================
 call order_fs_kpts(elph_ds%k_phon%kptirr,irredtoGS_phon,hdr,elph_ds%k_phon%nkptirr)

!==========================================
!3) reconstruct full kgrid from irred kpoints,
!==========================================
 call mkFSkgrid (elph_ds%k_phon, Cryst%nsym, Cryst%symrec, timrev) 

!====================================================================
!4) setup weights for integration (gaussian or tetrahedron method) 
!====================================================================
 ABI_ALLOCATE(elph_ds%k_phon%wtk,(elph_ds%nFSband,elph_ds%k_phon%nkpt,elph_ds%nsppol))
 istat = ABI_ALLOC_STAT

 call ep_fs_weights(anaddb_dtset%ep_b_min, anaddb_dtset%ep_b_max, eigenGS, anaddb_dtset%elphsmear, &
& elph_ds%fermie, gprimd, irredtoGS_phon, anaddb_dtset%kptrlatt, max_occ, elph_ds%minFSband, nband, elph_ds%nFSband, &
& elph_ds%nsppol, anaddb_dtset%telphint, elph_ds%k_phon)

!=====================================================
!get kpt info from the fine grid part 
!=====================================================
 if (anaddb_dtset%ep_alter_int_gam == 0) then
!  simply copy over _phon variables
   elph_ds%k_fine%nkpt = elph_ds%k_phon%nkpt
   elph_ds%k_fine%nkptirr = elph_ds%k_phon%nkptirr

   ABI_ALLOCATE(elph_ds%k_fine%kptirr,(3,elph_ds%k_fine%nkptirr))
   istat = ABI_ALLOC_STAT
   elph_ds%k_fine%kptirr = elph_ds%k_phon%kptirr
   ABI_ALLOCATE(elph_ds%k_fine%wtkirr,(elph_ds%k_fine%nkptirr))
   elph_ds%k_fine%wtkirr = elph_ds%k_phon%wtkirr
   ABI_ALLOCATE(irredtoGS_fine,(elph_ds%k_fine%nkptirr))
   istat = ABI_ALLOC_STAT
   irredtoGS_fine = irredtoGS_phon

   ABI_ALLOCATE(elph_ds%k_fine%wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol))
   istat = ABI_ALLOC_STAT
   elph_ds%k_fine%wtk = elph_ds%k_phon%wtk
   ABI_ALLOCATE(elph_ds%k_fine%kpt,(3,elph_ds%k_fine%nkpt))
   istat = ABI_ALLOC_STAT
   elph_ds%k_fine%kpt = elph_ds%k_phon%kpt 

   call copy_kptrank(elph_ds%k_phon%kptrank_t, elph_ds%k_fine%kptrank_t)

   ABI_ALLOCATE(elph_ds%k_fine%irr2full,(elph_ds%k_fine%nkptirr))
   istat = ABI_ALLOC_STAT
   elph_ds%k_fine%irr2full = elph_ds%k_phon%irr2full
   ABI_ALLOCATE(elph_ds%k_fine%full2irr,(3,elph_ds%k_fine%nkpt))
   istat = ABI_ALLOC_STAT
   elph_ds%k_fine%full2irr = elph_ds%k_phon%full2irr
   ABI_ALLOCATE(elph_ds%k_fine%full2full,(2,nsym,elph_ds%k_fine%nkpt))
   istat = ABI_ALLOC_STAT
   elph_ds%k_fine%full2full = elph_ds%k_phon%full2full

 else
!  read in the first header for the gkk part
   if (abs(anaddb_dtset%elph_fermie) > tol10) then
     MSG_ERROR("New fermi level + fine grid not coded!")
   end if

   unitfskgrid = get_unit()
   open (unit=unitfskgrid,file='finegrid_GKK',form='unformatted',status='old',iostat=ios)
   ABI_CHECK(ios==0,"opening finegrid_GKK")

   rewind(unitfskgrid)
   rdwr = 5 !read in header of file without rewinding it
   call hdr_io(fform,hdr1,rdwr,unitfskgrid)
   ABI_CHECK(fform/=0,'fine grid GKK header was mis-read. fform == 0')

   ABI_ALLOCATE(eigenGS_fine,(nband,hdr1%nkpt,elph_ds%nsppol))
   istat = ABI_ALLOC_STAT
   if (istat /= 0) stop 'elphon: error in allocating eigenGS_fine'
   do isppol=1,elph_ds%nsppol
     do ikpt=1,hdr1%nkpt
       read(unitfskgrid) eigenGS_fine(:,ikpt,isppol)
     end do
   end do
   close(unitfskgrid)

!  Reinit the structure storing the eigevalues.
!  Be careful. This part has not been tested.
   call bstruct_clean(Bst)
   call bst_init_from_hdr(Bst,Hdr1,nband,eigenGS_fine)

   elph_ds%k_fine%nkptirr = hdr1%nkpt
   ABI_ALLOCATE(elph_ds%k_fine%kptirr,(3,elph_ds%k_fine%nkptirr))
   istat = ABI_ALLOC_STAT
   ABI_ALLOCATE(irredtoGS_fine,(elph_ds%k_fine%nkptirr))
   istat = ABI_ALLOC_STAT

   call order_fs_kpts(elph_ds%k_fine%kptirr,irredtoGS_fine,&
&   hdr1,elph_ds%k_fine%nkptirr)

   call hdr_clean(hdr1)

   call mkFSkgrid (elph_ds%k_fine, Cryst%nsym, Cryst%symrec, timrev)

   ABI_ALLOCATE(elph_ds%k_fine%wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol))
   istat = ABI_ALLOC_STAT
   if (istat /= 0) stop 'elphon: error in allocating elph_ds%k_fine%wtk'

   call ep_fs_weights(anaddb_dtset%ep_b_min, anaddb_dtset%ep_b_max, eigenGS_fine, anaddb_dtset%elphsmear, &
&   elph_ds%fermie, gprimd, irredtoGS_fine, anaddb_dtset%kptrlatt, max_occ, elph_ds%minFSband, nband, elph_ds%nFSband, &
&   elph_ds%nsppol, anaddb_dtset%telphint, elph_ds%k_fine)

   ABI_DEALLOCATE(eigenGS_fine)
   ABI_DEALLOCATE(irredtoGS_fine)

 end if ! alter_int_gam

 ABI_DEALLOCATE(irredtoGS_phon)

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon k and q grids have been setup after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

!====================================================================
!5) calculate DOS at Ef
!====================================================================
!call get_dos(BSt,Kmesh,method,fildos,broad,dosdeltae)
 ABI_ALLOCATE(elph_ds%n0,(elph_ds%nsppol))
 istat = ABI_ALLOC_STAT

!SPPOL sum over spin channels to get total DOS
!channels decoupled => use separate values for DOS_up(Ef) resp down
 do isppol=1,elph_ds%nsppol
   elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
 end do

 if (elph_ds%nsppol == 1) then
   write (std_out,*) ' elphon : the estimated DOS(E_Fermi) = ', elph_ds%n0(1), ' states/Ha/spin '
   write (std_out,*) ' elphon : the total FS weight and # of kpoints = ',sum(elph_ds%k_fine%wtk),elph_ds%k_fine%nkpt
 else if (elph_ds%nsppol == 2) then
   write (std_out,*) ' elphon : the spin up   DOS(E_Fermi) = ', elph_ds%n0(1), ' states/Ha/spin '
   write (std_out,*) ' elphon : the spin down DOS(E_Fermi) = ', elph_ds%n0(2), ' states/Ha/spin '
   write (std_out,*) ' elphon : total DOS(E_Fermi) = ', elph_ds%n0(1)+elph_ds%n0(2), ' states/Ha '
   write (std_out,*) ' elphon : the spin up   FS weight and # of kpoints = ',&
&   sum(elph_ds%k_fine%wtk(:,:,1)),elph_ds%k_fine%nkpt
   write (std_out,*) ' elphon : the spin down FS weight and # of kpoints = ',&
&   sum(elph_ds%k_fine%wtk(:,:,2)),elph_ds%k_fine%nkpt
 else
   write (message,'(a,i0)') 'bad value for nsppol ', elph_ds%nsppol
   MSG_ERROR(message)
 end if

 ABI_ALLOCATE(elph_ds%gkk_intweight,(elph_ds%ngkkband,elph_ds%k_fine%nkpt,elph_ds%nsppol))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"allocating elph_ds%gkk_intweight")

 if (elph_ds%ep_keepbands == 0) then
!  use trivial integration weights  for single band,
!  since average over bands is done in normsq_gkk
   elph_ds%gkk_intweight(1,:,:) = one

 else if (elph_ds%ep_keepbands == 1) then
!  use elph_ds%k_fine%wtk since average over bands is not done in normsq_gkk
   elph_ds%gkk_intweight(:,:,:) = elph_ds%k_fine%wtk(:,:,:)
 else
   write(message,'(a,i0)')' ep_keepbands must be 0 or 1 while it is : ',elph_ds%ep_keepbands
   MSG_ERROR(message)
 end if

 ep_prt_wtk = 0
 if (ep_prt_wtk == 1) then
   do iband=1, elph_ds%ngkkband
     do ikpt_fine=1, elph_ds%k_fine%nkpt
       write (300,*) ikpt_fine, elph_ds%gkk_intweight(iband,ikpt_fine,1)
     end do
   end do
 end if


 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon weights and DOS setup after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

!Output of the Fermi Surface
 if (anaddb_dtset%prtfsurf == 1) then
   fname=trim(elph_ds%elph_base_name) // '_BXSF'

!  FIXME
!  shiftk is defined neither in the anaddb nor in the hdr data type
!  an incorrect FS will be produced in case of a shifted k-grid used during the GS calculation
!  check if we are using a unshifthed kgrid, obviously doesnt work in case
!  of multiple shifts containg a zero translation but in this case prtbxsf should work
   shiftk=one
   do ii=1,hdr%nkpt
     if (all(hdr%kptns(:,ii) == zero)) shiftk=zero
   end do

   use_afm=(hdr%nsppol==1.and.hdr%nspden==2)
!  MG FIXME warning time reversal is always assumed to be present. 
!  the header should report this information.

   use_tr=(timrev==1)

   call printbxsf(eigenGS,zero,elph_ds%fermie,Cryst%gprimd,&
&   anaddb_dtset%kptrlatt,nband,hdr%nkpt,hdr%kptns,&
&   Cryst%nsym,use_afm,Cryst%symrec,Cryst%symafm,use_tr,elph_ds%nsppol,shiftk,1,fname,ierr)

 end if !anaddb_dtset%prtfsurf

 ABI_DEALLOCATE(eigenGS)

!=========================================================    
!Get equivalence between a kpt_phon pair and a qpt in qpt_full
!only works if the qpt grid is complete (identical to
!the kpt one, with a basic shift of (0,0,0)
!=========================================================    

!mapping of k + q onto k' for k and k' in full BZ
 ABI_ALLOCATE(FSfullpqtofull,(elph_ds%k_phon%nkpt,elph_ds%nqpt_full))
 istat = ABI_ALLOC_STAT
 if (istat /= 0) stop 'elphon: error in allocating FSfullpqtofull'

!qpttoqpt(itim,isym,iqpt) = qpoint index which transforms to iqpt under isym and with time reversal itim.
 ABI_ALLOCATE(qpttoqpt,(2,Cryst%nsym,elph_ds%nqpt_full))
 istat = ABI_ALLOC_STAT
 if (istat /= 0) stop 'elphon: error in allocating qpttoqpt'

 call wrtout(std_out,'elphon: calling mkqptequiv to set up the FS qpoint set',"COLL")

 call mkqptequiv (FSfullpqtofull,Cryst,elph_ds%k_phon%kpt,elph_ds%k_phon%nkpt,&
& elph_ds%nqpt_full,qpttoqpt,elph_ds%qpt_full)

!==========================================
!Set up dataset for phonon interpolations
!==========================================

 call setup_phon_ds(phon_ds,Cryst,anaddb_dtset%dipdip,mpert,nrpt,anaddb_dtset%symdynmat,&
& ucvol,acell,amu,atmfrc,dielt,dyewq0,gprim,gmet,&
& xred,zeff,rcan,rmet,rprim,rprimd,rpt,trans,wghatm)

!transfer ifltransport flag to structure
 elph_tr_ds%ifltransport=anaddb_dtset%ifltransport
!transfer name of files file for ddk
 elph_tr_ds%ddkfilename=ddkfilename

!reduce qpt_full to correct zone
 do iqpt=1,elph_ds%nqpt_full
   call wrap2_pmhalf(elph_ds%qpt_full(1,iqpt),kpt(1),res)
   call wrap2_pmhalf(elph_ds%qpt_full(2,iqpt),kpt(2),res)
   call wrap2_pmhalf(elph_ds%qpt_full(3,iqpt),kpt(3),res)
   elph_ds%qpt_full(:,iqpt)=kpt
 end do

!test density of k+q grid: the following should be close to n0 squared
!FIXME: generalize for sppol
 res = zero
 do ikpt_fine = 1, elph_ds%k_fine%nkpt
   do iqpt = 1, elph_ds%nqpt_full
     kpt = elph_ds%k_fine%kpt(:,ikpt_fine) + elph_ds%qpt_full(:,iqpt)
     call get_rank_1kpt (kpt,symrankkpt,elph_ds%k_fine%kptrank_t)
     iFSkpq = elph_ds%k_fine%kptrank_t%invrank(symrankkpt)
     do iband = 1, elph_ds%ngkkband
       do ibandp = 1, elph_ds%ngkkband
         res = res + elph_ds%gkk_intweight(iband,ikpt_fine,1)*elph_ds%gkk_intweight(ibandp,iFSkpq,1)
       end do
     end do 
   end do
 end do
 res = res / elph_ds%k_fine%nkpt/elph_ds%k_fine%nkpt
 write (std_out,*) 'elphon: integrated value of intweight for given k and q grid : ', res, res / elph_ds%n0(1)**2

 res = zero
 do ikpt_fine = 1, elph_ds%k_fine%nkpt
   do iqpt = 1, elph_ds%k_fine%nkpt
     kpt = elph_ds%k_fine%kpt(:,ikpt_fine) + elph_ds%k_fine%kpt(:,iqpt)
     call get_rank_1kpt (kpt,symrankkpt,elph_ds%k_fine%kptrank_t)
     iFSkpq = elph_ds%k_fine%kptrank_t%invrank(symrankkpt)
     do iband = 1, elph_ds%ngkkband
       do ibandp = 1, elph_ds%ngkkband
         res = res + elph_ds%gkk_intweight(iband,ikpt_fine,1)*elph_ds%gkk_intweight(ibandp,iFSkpq,1)
       end do
     end do
   end do
 end do
 res = res / elph_ds%k_fine%nkpt/elph_ds%k_fine%nkpt
 write (std_out,*) 'elphon: integrated value of intweight for double k grid : ', res, res / elph_ds%n0(1)**2

!===================================================
!Allocate all important arrays for FS integrations
!===================================================

!Record sizes for matrices on disk: complex and real versions (for real and recip space resp!)
 onegkksize = 2*elph_ds%nbranch*elph_ds%nbranch*&
& elph_ds%ngkkband*elph_ds%ngkkband*&
& elph_ds%nsppol*kind(realdp_ex)

 elph_tr_ds%onegkksize=onegkksize

 write (message,'(4a)')&
& ' elphon : preliminary setup completed ',ch10,&
& '          calling get_all_gkq to read in all the e-ph matrix elements',ch10
 call wrtout(std_out,message,'COLL')

!flag to do scalar product in gkq before interpolation:
!should also used in interpolate_gkk and mkph_linwid
 elph_ds%ep_scalprod=anaddb_dtset%ep_scalprod
 if (elph_ds%ep_scalprod==0) then
   write (std_out,*) ' elphon: will NOT perform scalar product with phonon'
   write (std_out,*) '  displacement vectors in read_gkk. ep_scalprod==0'
 else if (elph_ds%ep_scalprod==1) then
   write (std_out,*) ' elphon: will perform scalar product with phonon'
   write (std_out,*) '  displacement vectors in read_gkk. ep_scalprod==1'
 else
   MSG_ERROR('illegal value for ep_scalprod')
 end if

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon begin gkq construction after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

 call get_all_gkq (elph_ds,Cryst,Bst,FSfullpqtofull,nband,n1wf,onegkksize,phon_ds,&
& qpttoqpt,anaddb_dtset%ep_prt_yambo,unitgkk)

 close (unitgkk)

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon end gkq construction after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

 if (elph_tr_ds%ifltransport==1 )then
   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-elphon begin gkq_tr construction after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   tcpui = tcpu
   twalli = twall

   call get_veloc_tr(elph_ds,mpi_enreg,nband,elph_tr_ds)

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-elphon end gkq_tr construction after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   tcpui = tcpu
   twalli = twall
 end if

!============================================================================
!Evaluate lambda and omega_log using the weighted sum over the irred q-points
!found in the GKK file. All the data we need are stored in elph_ds%qgrid_data
!============================================================================

 fname=trim(elph_ds%elph_base_name) // '_QPTS'
 call outelph(elph_ds,anaddb_dtset%enunit,fname)

!========================================================
!Get FS averaged gamma matrices and Fourier transform to real space 
!========================================================
 
 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon begin integration of gkq after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

 if (anaddb_dtset%ep_alter_int_gam == 1) then
   call integrate_gamma_alt(elph_ds,elph_tr_ds,Cryst,gprim,anaddb_dtset%kptrlatt,&
&   natom,nrpt,nsym,qpttoqpt,rpt,wghatm)

 else
   call integrate_gamma(elph_ds,FSfullpqtofull,nrpt)

   if (elph_ds%symgkq ==1) then
!    complete the gamma_qpt here instead of the gkk previously
     call complete_gamma(Cryst,elph_ds%nbranch,elph_ds%nsppol,elph_ds%nqptirred,elph_ds%nqpt_full,&
&     elph_ds%ep_scalprod,elph_ds%qirredtofull,qpttoqpt,elph_ds%gamma_qpt)
   end if

!  Now FT to real space too
!  NOTE: gprim (not gprimd) is used for all FT interpolations,
!  to be consistent with the dimensions of the rpt, which come from anaddb.
   qtor = 1 ! q --> r
   do isppol=1,elph_ds%nsppol
     call ftgam(wghatm,elph_ds%gamma_qpt(:,:,isppol,:),elph_ds%gamma_rpt(:,:,isppol,:),gprim,natom,&
&     elph_ds%nqpt_full,nrpt,qtor,rpt,elph_ds%qpt_full)
   end do
 end if

!DEBUG
 if (anaddb_dtset%ep_alter_int_gam == 1) then
   do ii=1, elph_ds%nbranch*elph_ds%nbranch
     write (100,*) ' ibr gamma ', ii, elph_ds%gamma_qpt(1,ii,1,:)
     write (100,*)
   end do
   do ikpt_phon=1, elph_ds%k_phon%nkpt
     write (200,'(I6,3e16.6,I12)') ikpt_phon, elph_ds%k_phon%kpt(:,ikpt_phon),&
&     elph_ds%k_phon%kptrank_t%rank(ikpt_phon)
   end do 
 end if
!ENDDEBUG

!BEGIN NEW_GKK
 if (.FALSE.) then
!  if (.TRUE.) then
   MSG_WARNING("Entering gamma_t part")
   call gamma_nullify(Gam)
   call gamma_init(Gam,gkk_fname,elph_ds,Cryst,Bst,Phon_ds,FSfullpqtofull,gprim,qptrlatt,nrpt,rpt,wghatm)

   call gamma_interp_setup(Gam,Cryst,"INIT")

   call gamma_linwid(Gam,Cryst,anaddb_dtset%ep_alter_int_gam,elph_ds,anaddb_dtset%nqpath,Phon_ds,anaddb_dtset%qpath)

   a2f_qptrlatt = 0
   do ii=1,3
     a2f_qptrlatt(ii,ii) = anaddb_dtset%ng2qpt(ii)
   end do
!  a2f_qptrlatt = RESHAPE( (/4,0,0,0,4,0,0,0,4/),(/3,3/) )

   a2f_nqshift  = 1
   ABI_ALLOCATE(a2f_qshift,(3,a2f_nqshift))
   a2f_qshift(:,1) = anaddb_dtset%q2shft(:)  ! FIXME small inconsistency in the dimension of q1shft

   call a2f_init(A2f,Cryst,Gam,Phon_ds,elph_ds%na2f,elph_ds%a2fsmear,elph_ds%n0,anaddb_dtset%ep_alter_int_gam,&
&   a2f_qptrlatt,a2f_nqshift,a2f_qshift,qptopt=3)

   call a2f_dump(A2f,"NEW_"//TRIM(elph_ds%elph_base_name)//'_A2F_Q3')
   call a2f_free(A2f)

   call a2f_init(A2f,Cryst,Gam,Phon_ds,elph_ds%na2f,elph_ds%a2fsmear,elph_ds%n0,anaddb_dtset%ep_alter_int_gam,&
&   a2f_qptrlatt,a2f_nqshift,a2f_qshift,qptopt=1)

   call a2f_dump(A2f,"NEW_"//TRIM(elph_ds%elph_base_name)//'_A2F_Q1')
   call a2f_free(A2f)

   ABI_DEALLOCATE(a2f_qshift)

   call gamma_interp_setup(Gam,Cryst,"FREE")
   call gamma_free(Gam)

   if (anaddb_dtset%prtdos /=0 ) then
     call mkphdos(PHdos,anaddb_dtset%prtdos,anaddb_dtset%dosdeltae,anaddb_dtset%dossmear,&
&     anaddb_dtset%dipdip,anaddb_dtset%symdynmat,&
&     acell,amu,anaddb_dtset,atmfrc,dielt,dyewq0,gmet,gprim,indsym,&
&     mpert,nsym,natom,nrpt,nsym,ntypat,rmet,rprim,rpt,symrec,symrel,trans,typat,ucvol,wghatm,xred,zeff)

     call print_phondos(PHdos,"NEW_PHDOS")
     call destroy_phondos(PHdos)
   end if
 end if
!END NEW_GKK

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon end integration and completion of gkq after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall


!==========================================================
!calculate transport matrix elements, integrated over FS
!in ep_alter_int_gam == 1 case the gamma_qpt_trout are
!already calculated in integrate_gamma_alt
!==========================================================
 if (elph_tr_ds%ifltransport==1 .and. anaddb_dtset%ep_alter_int_gam == 0 )then
   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-elphon begin integrate gkq_tr after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   tcpui = tcpu
   twalli = twall

   call integrate_gamma_tr(elph_ds,FSfullpqtofull,nrpt,elph_tr_ds)

   call complete_gamma_tr(elph_ds,elph_tr_ds%gamma_qpt_trout, &
&   gprimd,indsym,natom,nsym,qpttoqpt,rprimd, &
&   symrec,symrel)

   call complete_gamma_tr(elph_ds,elph_tr_ds%gamma_qpt_trin, &
&   gprimd,indsym,natom,nsym,qpttoqpt,rprimd, &
&   symrec,symrel)

!  Now FT to real space too
   qtor = 1 ! q --> r
   do isppol=1,elph_ds%nsppol
     do idir=1,9
       call ftgam(wghatm,elph_tr_ds%gamma_qpt_trout(:,idir,:,isppol,:),&
&       elph_tr_ds%gamma_rpt_trout(:,idir,:,isppol,:),gprim,natom,&
&       elph_ds%nqpt_full,nrpt,qtor,rpt,elph_ds%qpt_full)

       call ftgam(wghatm,elph_tr_ds%gamma_qpt_trin(:,idir,:,isppol,:),&
&       elph_tr_ds%gamma_rpt_trin(:,idir,:,isppol,:),gprim,natom,&
&       elph_ds%nqpt_full,nrpt,qtor,rpt,elph_ds%qpt_full)
     end do
   end do

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-elphon end integrate gkq_tr after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   tcpui = tcpu
   twalli = twall
 end if

 ABI_DEALLOCATE(qpttoqpt)
 ABI_DEALLOCATE(FSfullpqtofull)

!==============================================================
!Calculate phonon linewidths, interpolating on chosen qpoints
!==============================================================

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon begin linewidths after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

 call mkph_linwid(Cryst,anaddb_dtset%ep_alter_int_gam, elph_ds,gprim,anaddb_dtset%kptrlatt_fine,&
& nrpt,anaddb_dtset%nqpath,phon_ds,anaddb_dtset%qpath,rpt,wghatm)

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon end linewidths after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

!==============================================================
!the nesting factor calculation
!FIXME: this could go higher up, before the call to get_all_gkq
!you only need the kpt and weight info
!==============================================================
 if (anaddb_dtset%prtnest==1 .or. anaddb_dtset%prtnest==2) then

   nestname = trim(elph_ds%elph_base_name) // "_NEST"
   call mknesting(elph_ds%k_phon%nkpt,elph_ds%k_phon%kpt,anaddb_dtset%kptrlatt,elph_ds%nFSband,&
&   elph_ds%k_phon%wtk,anaddb_dtset%nqpath,anaddb_dtset%qpath,elph_ds%nqpt_full, &
&   elph_ds%qpt_full,nestname,gprimd,gmet,anaddb_dtset%prtnest,qptrlatt)

 end if
 
!======================================================
!Calculate alpha^2 F integrating over fine kpt_phon grid
!======================================================

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon begin mka2f after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

 ABI_ALLOCATE(a2f_1d,(elph_ds%na2f))
 ABI_ALLOCATE(dos_phon,(elph_ds%na2f))
 istat = ABI_ALLOC_STAT
 if (istat /= 0) stop 'elphon: error in allocating a2f_1d,dos_phon'

 call mka2f(Cryst,anaddb_dtset%ep_alter_int_gam,a2f_1d,dos_phon,elph_ds,&
& gprim,anaddb_dtset%kptrlatt_fine,anaddb_dtset%mustar,nrpt,phon_ds,rpt,wghatm)

!MG FIXME  Why do we need to store domega in elph_ds?
!%elph_ds%domega  = (elph_ds%omega_max-elph_ds%omega_min)/(elph_ds%na2f-one) 

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon end mka2f after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

!calculate transport spectral function and coefficients
 if (elph_tr_ds%ifltransport==1 )then
   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-elphon begin mka2f_tr after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   tcpui = tcpu
   twalli = twall

   call mka2f_tr(anaddb_dtset%ep_alter_int_gam, elph_ds,gprim,gprimd,ucvol,natom,nrpt,&
&   anaddb_dtset%ntemper,anaddb_dtset%tempermin,anaddb_dtset%temperinc,&
&   phon_ds,rpt,wghatm,elph_tr_ds)

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-elphon end mka2f_tr after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   tcpui = tcpu
   twalli = twall

 end if

!evaluate a2F only using the input Q-grid (without using interpolated matrices)
!SCOPE: test the validity of the Fourier interpolation
 call wrtout(std_out,' elphon : calling mka2fQgrid',"COLL")

 fname=trim(elph_ds%elph_base_name) // '_A2F_QGRID'
 call mka2fQgrid(elph_ds,fname)

!=============================================
!Eliashberg equation in 1-D (isotropic case)
!=============================================

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon begin eliashberg_1d after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

 call eliashberg_1d(a2f_1d,elph_ds,anaddb_dtset%mustar)

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon end eliashberg_1d after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

 ABI_DEALLOCATE(a2f_1d)
 ABI_DEALLOCATE(dos_phon)

!MJV: 20070805 should exit here. None of the rest is tested or used yet to my knowledge

!========================================================================
!Now gkk contains the matrix elements of dH(1)/dxi i=1,2,3
!for kpoints on the FS but qpoints only in the given grid {Q}.
!
!1.) Need to complete the gkk elements for q and k\prime=k+q not 
!in the set of {k+Q} by Fourier interpolation on the Q.
!
!2.) Need to complete the dynamical matrices and phonon freqs for
!all q between points on the FS.
!
!3.) With the eigenvectors e_ph of the dyn mats, do the scalar product
!e_ph . gkk, which implies the gkk are turned to the eigenbasis of
!the phonons. Before the (non eigen-) modes are ordered
!atom1 xred1 atom1 xred2 atom1 xred3
!atom2 xred1 atom2 xred2 atom2 xred3 ...
!=======================================================================

 make_gkk2=.false.
 
 if (.not. make_gkk2) then
   call wrtout(std_out,' elphon : skipping full g(k,k") interpolation ',"COLL")
 else

!  ==========================================================
!  FT of recip space gkk matrices to real space (gkk_rpt)
!  NOTE: could be made into FFT, couldnt it? If shifts are
!  used with a homogeneous grid
!  ==========================================================
   write (message,'(2a,i0)')ch10,&
&   ' elphon : Fourier transform (q --> r) of the gkk matrices using nrpt = ',nrpt
   call wrtout(std_out,message,'COLL')

   call get_all_gkr(elph_ds,gprim,natom,nrpt,onegkksize,rpt,elph_ds%qpt_full,wghatm)

!  =========================================================
!  complete gkk2 for all qpts between points
!  on full kpt grid (interpolation from real space values)
!  =========================================================

   write(message,'(2a)')ch10,&
&   ' elphon : Calling get_all_gkk2 to calculate gkk2 for q points over the full k grid'
   call wrtout(std_out,message,'COLL')
   
   call get_all_gkk2(elph_ds,elph_ds%k_phon%kptirr,elph_ds%k_phon%kpt,natom,nrpt,phon_ds,rcan,wghatm)
 end if

!=====================================================
!Here should be the anisotropic Eliashberg equations.
!=====================================================

 if (.FALSE.) then
   ABI_ALLOCATE(zz,(elph_ds%na2f,elph_ds%k_phon%nkpt))
   istat = ABI_ALLOC_STAT
   if (istat /= 0) stop 'elphon: error in allocating zz'
   ABI_ALLOCATE(delta,(elph_ds%na2f,elph_ds%k_phon%nkpt))
   istat = ABI_ALLOC_STAT
   if (istat /= 0) stop 'elphon: error in allocating delta'

!  initialize delta function

!  initialize T_c

!  initialize delta

!  iterate for calculation of T_c
   neliash = 10   !maximum number of iterations to converge Tc
   do ieliash=1,neliash

!    ===========================
!    calculate lambda function
!    ===========================

!    ========================================
!    integrate lambda over FS -> Z function
!    ========================================
     
!    ========================================
!    integrate delta*Z over FS -> new delta
!    ========================================
     
!    update T_c

   end do

!  end iterate ieliash

!  output T_c and related quantities
 end if

!clean and deallocate junk
 call destroy_crystal(Cryst)
 call bstruct_clean(Bst)

 call elph_ds_clean(elph_ds)
 call elph_tr_ds_clean(elph_tr_ds)

 call clean_phon_ds(phon_ds)

 call hdr_clean(hdr)

end subroutine elphon
!!***
