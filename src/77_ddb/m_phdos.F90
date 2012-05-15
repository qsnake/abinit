!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_phdos
!! NAME
!! m_phdos
!!
!! FUNCTION
!! Module for the phonon density of states.
!! Container type is defined, and destruction, print subroutines 
!! as well as the central mkphdos 
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MG,MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_phdos

 use m_profiling

 use defs_basis
 use defs_abitypes,  only : anaddb_dataset_type
 use m_errors

 use m_io_tools, only : get_unit

 implicit none

 private
!!***

!!****t* m_phdos/phonon_dos_type
!! NAME
!! phonon_dos_type
!! 
!! FUNCTION
!! Container for phonon DOS and atom projected contributions 
!! 
!! SOURCE

 type,public :: phonon_dos_type

! Integer
  integer :: ntypat
  ! Number of type of atoms.

  integer :: natom
  ! Number of atoms is the unit cell.

  integer :: prtdos
  ! Option of DOS calculation (1 for Gaussian, 2 for tetrahedrons).

  integer :: nomega
  ! Number of frequency points in DOS mesh.

  integer :: nqibz
  ! Number of q-points in the IBZ.

  integer :: ntetra_ibz
  ! Number of tetrahedrons in the IBZ.

! Reals
  real(dp) :: omega_min
  ! Min frequency for DOS calculation.

  real(dp) :: omega_max
  ! Max frequency for DOS calculation.

  real(dp) :: omega_step
  ! Frequency step.

  real(dp) :: dossmear
  ! Gaussian broadening.

! Real pointers
  real(dp), pointer :: omega(:)     SET2NULL
   ! omega(nomega)   
   ! Frequency grid.

  real(dp), pointer :: phdos(:)     SET2NULL
   ! phdos(nomega)   
   ! phonon DOS.

  real(dp), pointer :: phdos_int(:)   SET2NULL
   ! phdos_int(nomega)  
   ! integrated phonon DOS

  real(dp), pointer :: pjdos(:,:,:)   SET2NULL
   ! pjdos(nomega,3,natom)
   ! projected DOS (over atoms)

  real(dp), pointer :: pjdos_int(:,:,:)   SET2NULL
   ! pjdos_int(nomega,3,natom)
   ! Integrate atomic PJDOS

  real(dp), pointer :: pjdos_typ(:,:)    SET2NULL       
   ! pjdos_typ(nomega,ntypat)
   ! phonon DOS contribution arising from a particular atom-type.

  real(dp), pointer :: pjdos_typ_int(:,:)   SET2NULL    
   ! pjdos_typ_int(nomega,ntypat)
   ! Integrate phonon DOS contribution arising from a particular atom-type.

  real(dp), pointer :: pjdos_xyz_typ(:,:,:)   SET2NULL
   ! phdos(nomega,3,ntypat)
   ! phonon DOS contribution arising from a particular atom-type 
   ! decomposed along the three reduced directions.

 end type phonon_dos_type

 public :: print_phondos
 public :: print_phondos_debye
 public :: init_phondos
 public :: nullify_phondos
 public :: destroy_phondos
 public :: mkphdos
!!**

CONTAINS  !=========================================================================================================================
!!***

!!****f* m_phdos/print_phondos
!!
!! NAME
!! print_phondos
!!
!! FUNCTION
!! Print out phonon DOS (and partial DOS etc) in meV units
!!
!! INPUTS
!! PHdos= container object for phonon DOS
!! fname=File name for output
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb,elphon,scphon_interpolate_phonon_and_dos
!!
!! CHILDREN
!!      get_full_kgrid,get_tetra_weight,gtdyn9,init_phondos,init_tetra,matr3inv
!!      mkrdim,phfrq3,simpson_int,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine print_phondos(PHdos,fname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_phondos'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: fname
 type(phonon_dos_type),intent(in) :: PHdos

!Local variables-------------------------------
 integer :: io,itype,ios,unt
 real(dp) :: cfact
 character(len=50) :: frmt
 character(len=500) :: msg
 character(len=3) :: unitname

! *************************************************************************

! === Convert everything into meV ===
 cfact=Ha_eV*1000 ; unitname='meV'
! === Leave everything in Ha      ===
! this should be the abinit default!
 cfact=one        ; unitname='Ha'

! === Open external file and write results ===
! TODO Here I have to rationalize how to write all this stuff!!
! 
 unt=get_unit()
! open(unit=unt,file=fname,form='formatted',status='new')
 open(unit=unt,file=fname,form='formatted',iostat=ios)
 ABI_CHECK(ios==0,'Opening '//TRIM(fname))

 write(msg,'(3a)')'# ',ch10,'# Phonon density of states and projected DOS'
 call wrtout(unt,msg,'COLL')
 write(msg,'(6a)')'# ',ch10,'# Energy in ',unitname,', DOS in states/',unitname
 call wrtout(unt,msg,'COLL')

 select case (PHdos%prtdos)
 case (1)
   write(msg,'(a,es16.8,2a,i4)')&
&   '# Gaussian method with smearing = ',PHdos%dossmear*cfact,unitname,', nqibz =',PHdos%nqibz
 case (2) 
   write(msg,'(a,i5,a,i4)')'# Tetrahedron method, number of irreducible tetrahedrons = ',&
&   PHdos%ntetra_ibz,', nqibz= ',PHdos%nqibz
 case default
   write(msg,'(a,i0)')" Wrong prtdos = ",PHdos%prtdos
   MSG_ERROR(msg)
 end select
 call wrtout(unt,msg,'COLL')

 write(msg,'(5a)')'# ',ch10,'# omega     PHDOS    INT_PHDOS   PJDOS[atom_type=1]  INT_PJDOS[atom_type1] ...  ',ch10,'# '
 call wrtout(unt,msg,'COLL')
 write(frmt,*)'(',PHdos%ntypat,'(2es17.8))'

 do io=1,PHdos%nomega
   write(unt,'(3es17.8)',advance='NO')PHdos%omega(io)*cfact,PHdos%phdos(io)/cfact,PHdos%phdos_int(io)/cfact 
   do itype=1,PHdos%ntypat
     write(unt,frmt,advance='NO')PHdos%pjdos_typ(io,itype)/cfact,PHdos%pjdos_typ_int(io,itype)/cfact
   end do 
   write(unt,*)
 end do

 close(unt)

end subroutine print_phondos
!!***

!----------------------------------------------------------------------
!****f* m_phdos/print_phondos_debye
!!
!! NAME
!! print_phondos_debye
!!
!! FUNCTION
!! Print out global Debye temperature, force constant, etc... from phonon DOS
!!
!! INPUTS
!! phonon_dos= container object for phonon DOS
!! ucvol = unit cell volume
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      get_full_kgrid,get_tetra_weight,gtdyn9,init_phondos,init_tetra,matr3inv
!!      mkrdim,phfrq3,simpson_int,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine print_phondos_debye(PHdos, ucvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_phondos_debye'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

implicit none

!Arguments ------------------------------------
 real(dp), intent(in) :: ucvol
 type(phonon_dos_type),intent(in) :: PHdos

!Local variables-------------------------------
 integer :: io, iomax, iomin
 real(dp) :: avgom2dos, avgspeedofsound
 real(dp) :: debyefreq, meanfreq, meanfreq2

 real(dp), allocatable :: om2dos(:), om1dos(:), intdos(:)
 character(len=500) :: msg


! average speed of sound: coefficient of omega^2 in the DOS is = Volume / 2 pi^2 hbar^3 v_s^3
! first find how far out we can fit with a parabola
 ABI_ALLOCATE(om2dos,(PHdos%nomega))
 ABI_ALLOCATE(om1dos,(PHdos%nomega))
 ABI_ALLOCATE(intdos,(PHdos%nomega))
 avgom2dos = zero
 do io=1,PHdos%nomega
   om1dos(io) = PHdos%omega(io)    * PHdos%phdos(io)
   om2dos(io) = PHdos%omega(io)**2 * PHdos%phdos(io)
 end do

! integrate dos * omega
 intdos = zero
 call simpson_int(PHdos%nomega,PHdos%omega_step,om1dos,intdos)
 meanfreq = intdos(PHdos%nomega)

! integrate dos * omega^2
 intdos = zero
 call simpson_int(PHdos%nomega,PHdos%omega_step,om2dos,intdos)
 meanfreq2 = intdos(PHdos%nomega)

! Debye frequency = sqrt(<omega^2>* 3/2)
 debyefreq = sqrt (meanfreq2 * three * half)
 write (msg,'(a,E20.10,a,E20.10,a)') ' Debye frequency from DOS: ', debyefreq, ' (Ha) = ', &
&    debyefreq*Ha_THz, ' (THz)'
 call wrtout (ab_out,msg,"COLL")
 call wrtout (std_out,msg,"COLL")
! Debye temperature = hbar * Debye frequency / kb
 write (msg,'(a,E20.10,2a)') ' Debye temperature from DOS: ', debyefreq*Ha_K, ' (K)', ch10
 call wrtout (ab_out,msg,"COLL")
 call wrtout (std_out,msg,"COLL")

 iomin = 1
 iomax = PHdos%nomega
 do io = 1, PHdos%nomega
   ! skip eventual negative frequency modes
   if (PHdos%omega(io) <= tol10) then
     iomin = io
     cycle
   end if

   ! accumulate dos * om^2 to make an average
   avgom2dos = avgom2dos + om2dos(io)
   ! first deviation from initial value of more than 10 percent
   if (abs(one-om2dos(1)/om2dos(io)) > 0.1_dp) then
     iomax = io
     exit
   end if
 end do

 avgom2dos = avgom2dos / (iomax-iomin+1)
! this value is also useful for partial atomic DOS, related to kinetic energy and Force constant in Moessbauer

 avgspeedofsound = (ucvol / 2 / pi**2 / avgom2dos)**third
 write (msg,'(a,E20.10,a,F16.4,2a)') '- Average speed of sound: ', avgspeedofsound, ' (at units) = ', &
&    avgspeedofsound * Bohr_Ang * 1.d-10 / Time_Sec, ' (m/s)',ch10
 call wrtout (ab_out,msg,"COLL")
 call wrtout (std_out,msg,"COLL")

! average force constant

! 

 ABI_DEALLOCATE(om2dos)
 ABI_DEALLOCATE(om1dos)
 ABI_DEALLOCATE(intdos)

end subroutine print_phondos_debye
!!***

!----------------------------------------------------------------------

!!****f* m_phdos/nullify_phondos
!!
!! NAME
!! nullify_phondos
!!
!! FUNCTION
!! Set all pointers to null()
!!
!! SIDE EFFECTS
!! PHdos= Pointers in the datastructure set to null()
!!
!! PARENTS
!!      m_phdos
!!
!! CHILDREN
!!      get_full_kgrid,get_tetra_weight,gtdyn9,init_phondos,init_tetra,matr3inv
!!      mkrdim,phfrq3,simpson_int,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine nullify_phondos(PHdos)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_phondos'
!End of the abilint section

 implicit none

!Arguments -------------------------------
 type(phonon_dos_type),intent(inout) :: PHdos

! *************************************************************************

 !@phonon_dos_type
 nullify(PHdos%omega)
 nullify(PHdos%phdos)
 nullify(PHdos%phdos_int)
 nullify(PHdos%pjdos)
 nullify(PHdos%pjdos_int)
 nullify(PHdos%pjdos_typ)
 nullify(PHdos%pjdos_typ_int)
 nullify(PHdos%pjdos_xyz_typ)

end subroutine nullify_phondos
!!***

!----------------------------------------------------------------------

!!****f* m_phdos/init_phondos
!!
!! NAME
!! init_phondos
!!
!! FUNCTION
!! Init function for phonon DOS object
!!
!! INPUTS
!!
!! OUTPUT
!! PHdos= container object for phonon DOS, filled and allocated
!!
!! PARENTS
!!      m_phdos,scphon
!!
!! CHILDREN
!!      get_full_kgrid,get_tetra_weight,gtdyn9,init_phondos,init_tetra,matr3inv
!!      mkrdim,phfrq3,simpson_int,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine init_phondos(PHdos,ntypat,natom,prtdos,nomega,nqibz,ntetra_ibz,omega_max,omega_min,omega_step,dossmear)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_phondos'
!End of the abilint section

 implicit none

!Arguments -------------------------------
 integer, intent(in) :: ntypat,natom,prtdos
 integer, intent(in) :: nomega,nqibz,ntetra_ibz
 real(dp), intent(in) :: omega_step,dossmear
 real(dp), intent(in) :: omega_max,omega_min
 type(phonon_dos_type),intent(out) :: PHdos

! *************************************************************************

 !@phonon_dos_type
 call nullify_phondos(PHdos)

 PHdos%ntypat     = ntypat
 PHdos%natom      = natom
 PHdos%prtdos     = prtdos
 PHdos%nomega     = nomega
 PHdos%nqibz      = nqibz
 PHdos%ntetra_ibz = ntetra_ibz

 PHdos%omega_max  = omega_max
 PHdos%omega_min  = omega_min
 PHdos%omega_step = omega_step
 PHdos%dossmear   = dossmear

 ABI_ALLOCATE(PHdos%omega,(nomega))
 ABI_ALLOCATE(PHdos%phdos,(nomega))
 ABI_ALLOCATE(PHdos%phdos_int,(nomega))
 ABI_ALLOCATE(PHdos%pjdos,(nomega,3,natom))
 ABI_ALLOCATE(PHdos%pjdos_int,(nomega,3,natom))
 ABI_ALLOCATE(PHdos%pjdos_typ,(nomega,ntypat))
 ABI_ALLOCATE(PHdos%pjdos_typ_int,(nomega,ntypat))
 ABI_ALLOCATE(PHdos%pjdos_xyz_typ,(nomega,3,ntypat))

end subroutine init_phondos
!!***

!----------------------------------------------------------------------

!!****f* m_phdos/destroy_phondos
!!
!! NAME
!! destroy_phondos
!!
!! FUNCTION
!! destructor function for phonon DOS object
!!
!! INPUTS
!! PHdos= container object for phonon DOS
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb,elphon
!!
!! CHILDREN
!!      get_full_kgrid,get_tetra_weight,gtdyn9,init_phondos,init_tetra,matr3inv
!!      mkrdim,phfrq3,simpson_int,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine destroy_phondos(PHdos)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_phondos'
!End of the abilint section

 implicit none

!Arguments -------------------------------
 type(phonon_dos_type),intent(inout) ::PHdos

! *************************************************************************

 !@phonon_dos_type       
 if (associated(PHdos%omega))           then
   ABI_DEALLOCATE(PHdos%omega)
 end if
 if (associated(PHdos%phdos))           then
   ABI_DEALLOCATE(PHdos%phdos)
 end if
 if (associated(PHdos%phdos_int))       then
   ABI_DEALLOCATE(PHdos%phdos_int)
 end if
 if (associated(PHdos%pjdos))           then
   ABI_DEALLOCATE(PHdos%pjdos)
 end if
 if (associated(PHdos%pjdos_int))       then
   ABI_DEALLOCATE(PHdos%pjdos_int)
 end if
 if (associated(PHdos%pjdos_typ))       then
   ABI_DEALLOCATE(PHdos%pjdos_typ)
 end if
 if (associated(PHdos%pjdos_typ_int))   then
   ABI_DEALLOCATE(PHdos%pjdos_typ_int)
 end if
 if (associated(PHdos%pjdos_xyz_typ))   then
   ABI_DEALLOCATE(PHdos%pjdos_xyz_typ)
 end if

end subroutine destroy_phondos
!!***

!----------------------------------------------------------------------

!!****f* m_phdos/mkphdos
!!
!! NAME
!! mkphdos
!!
!! FUNCTION
!! Function to calculate the phonon density of states as well as 
!! the contributions associated to the different types of atoms in the unit cell.
!! Two methods are implemented: gaussian method and linear interpolation based on 
!! tetrahedrons.
!!
!! INPUTS
!! prtdos=1 for Gaussian method, 2 for tetrahedra.
!! dosdeltae=Step for the frequency mesh.
!! dossmear=Gaussian broadening used if prtdos==1.
!! dipdip= if 0, no dipole-dipole interaction was subtracted in atmfrc
!!         if 1, atmfrc has been build without dipole-dipole part
!! symdynmat=if 1, (re)symmetrize the dynamical matrix, except if Gamma wavevector with electric field added.
!! acell(3) =length scales by which rprim is to be multiplied
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! anaddb_dtset= (derived datatype) contains all the input variables
!! atmfrc(2,3,natom,3,natom,nrpt) = Interatomic Forces in real space
!! dielt(3,3)=dielectric tensor
!! dyewq0(3,3,natom)=Ewald part of the dynamical matrix, at q=0.
!! gmet(3,3)= metric tensor in reciprocal space.
!! gprim(3,3)=normalized coordinates in reciprocal space
!! indsym(4,nsym,natom)=label given by subroutine symatm, indicating atom
!!  label which gets rotated into given atom by given symmetry
!!  (first three elements are related primitive translation--
!!  see symatm where this is computed)
!! mpert =maximum number of ipert
!! msym = maximum number of symmetries
!! natom=number of atoms in the unit cell
!! nrpt=number of R points in the Big Box
!! nsym=number of symmetries
!! ntypat=number of atom types
!! rmet(3,3)=metric tensor in real space.
!! rprim(3,3)=dimensionless primitive translations in real space
!! rpt(3,nprt)=canonical coordinates of the R points in the unit cell
!!  These coordinates are normalized (=> * acell(3)!!)
!! symrec(3,3,nsym)=symmetry operations
!! symrel(3,3,nsym)=symmetry operations
!! trans(3,natom)=atomic translations : xred = rcan + trans
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume
!! wghatm(natom,natom,nrpt)=weights associated to a pair of atoms and to a R vector
!! xred(3,natom)= relative coords of atoms in unit cell (dimensionless)
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement
!!
!! OUTPUT
!!
!! NOTES
!! 
!! On the use of the q-grids : 
!! Two different q-meshes are used in this subroutine. The first one is the coarse 
!! mesh where the interatomic forces have been calculated during the DFPT run. 
!! This q-grid is used to obtain an initial guess for the max and min frequency 
!! value of the phonon spectrum. These values are, indeed, required to dimension 
!! the array containing the PHDOS. The second (dense) grid is used to perform the 
!! PHDOS calculation. If the Fourier interpolation on the second dense q-grid 
!! generates a phonon frequency outside the initially calculated frequency mesh,
!! the mesh is enlarged and the calculation is restarted.
!!
!! PARENTS
!!      anaddb,elphon,scphon_interpolate_phonon_and_dos
!!
!! CHILDREN
!!      get_full_kgrid,get_tetra_weight,gtdyn9,init_phondos,init_tetra,matr3inv
!!      mkrdim,phfrq3,simpson_int,smpbz,symkpt,wrtout
!!
!! SOURCE

subroutine mkphdos(PHdos,prtdos,dosdeltae,dossmear,dipdip,symdynmat,acell,amu,anaddb_dtset,&
& atmfrc,dielt,dyewq0,gmet,gprim,indsym,&
& mpert,msym,natom,nrpt,nsym,ntypat,rmet,rprim,rpt,symrec,symrel,trans,typat,ucvol,wghatm,xred,zeff)

 use m_tetrahedron
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkphdos'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_56_recipspace
 use interfaces_72_response
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: prtdos,dipdip,symdynmat
 integer,intent(in) :: mpert,msym,natom,nrpt,nsym,ntypat
 real(dp),intent(in) :: dosdeltae,dossmear,ucvol
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
 type(phonon_dos_type),intent(inout) :: PHdos
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym),symrel(3,3,nsym),typat(natom)
 real(dp),intent(in) :: acell(3),amu(ntypat),dielt(3,3),gmet(3,3),gprim(3,3)
 real(dp),intent(in) :: rmet(3,3),rprim(3,3),rpt(3,nrpt)
 real(dp),intent(in) :: trans(3,natom),wghatm(natom,natom,nrpt),xred(3,natom)
 real(dp),intent(in) :: zeff(3,3,natom)
 real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt),dyewq0(3,3,natom)

!Local variables -------------------------
!scalars
 integer,parameter :: brav1=1,chksymbreak0=0
 integer :: facbrv,iat,idir,imesh,imode,io,iq_ibz,istat,itype,mtetra,nkpt_fullbz
 integer :: nmesh,nqbz,nqpt_max,nqshft,option,timrev
 integer :: ierr
 real(dp) :: bzvol,dum,gaussfactor,gaussprefactor
 real(dp) :: gaussval,low_bound,max_occ,pnorm 
 real(dp) :: qphnrm,upr_bound,xx
 logical :: out_of_bounds
 character(len=500) :: msg
 character(len=80) :: errstr
 type(t_tetrahedron) :: tetrahedra_q
!arrays
 integer :: qptrlatt(3,3)
 integer,allocatable :: bz2ibz(:),ibz2bz(:),ngqpt(:,:),tetra_full(:,:,:)
 integer,allocatable :: tetra_mult(:),tetra_wrap(:,:,:)
 real(dp) :: d2cart(2,3,mpert,3,mpert),displ(2*3*natom*3*natom),eigval(3*natom)
 real(dp) :: eigvec(2,3,natom,3*natom),gprimd(3,3),phfrq(3*natom)
 real(dp) :: qlatt(3,3),qphon(3),rlatt(3,3),rprimd(3,3)
 real(dp),allocatable :: dtweightde(:,:),full_eigvec(:,:,:,:,:),full_phfrq(:,:)
 real(dp),allocatable :: kpt_fullbz(:,:),qbz(:,:),qibz(:,:),qshft(:,:),tmp_phfrq(:),tweight(:,:)
 real(dp),allocatable :: wtq(:),wtq_folded(:),wtqibz(:)

! *********************************************************************

 DBG_ENTER("COLL")

 if ( ALL(prtdos/=(/1,2/) ) ) then 
   write(msg,'(a,i0)')&
&   ' The argument prtdos should be 1 or 2, however prtdos = ',prtdos 
   MSG_BUG(msg)
 end if 

 if (dosdeltae<=zero) then 
   write(msg,'(a,es16.8)')&
&   ' The argument dosdeltae should be positive, however dosdeltae = ',dosdeltae 
   MSG_BUG(msg)
 end if 

 if (prtdos==1.and.dossmear<=zero) then 
   write(msg,'(a,es16.8)')&
&   ' The argument anaddb%dossmear should be positive however dossmear = ',dossmear
   MSG_BUG(msg)
 end if 
 !
 call mkrdim(acell,rprim,rprimd)
 call matr3inv(rprimd,gprimd)

 bzvol=ABS ( gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
& -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
& +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))
 !
 ! Initialize container type, but with minimal values
 call init_phondos(PHdos,ntypat,natom,prtdos,1,1,1,smallest_real,greatest_real,dosdeltae,dossmear)
 !
 ! === Parameters defining the gaussian approximant ===
 if (prtdos==1) then 
   gaussprefactor=one/(dossmear*sqrt(two_pi))
   gaussfactor=one/(sqrt2*dossmear)
   write(msg,'(4a,f8.5,2a,f8.5)')ch10,&
&   ' mkphdos : calculating phonon DOS using gaussian method :',ch10,&
&   '  gaussian smearing [meV] = ',dossmear*Ha_meV,ch10,&
&   '  frequency step    [meV] = ',PHdos%omega_step*Ha_meV
 else if (prtdos==2) then 
   write(msg,'(2a)')ch10,' mkphdos : calculating phonon DOS using tetrahedron method '
 end if 
 call wrtout(std_out,msg,'COLL')
 !
 ! Initial lower and upper bound of the phonon spectrum.
 low_bound=greatest_real 
 upr_bound=smallest_real
 !
 ! Save memory during the generation of the q-mesh in the full BZ  
 ! Take into account the type of Bravais lattice
 facbrv=1
 if (brav1==2) facbrv=2
 if (brav1==3) facbrv=4

 nmesh=2
 ABI_ALLOCATE(ngqpt,(3,nmesh))
 do imesh=1,nmesh

   if (imesh==1) then  ! Coarse q-mesh used during RF calculation.
     ngqpt(:,imesh)=anaddb_dtset%ngqpt(1:3)
     nqshft=anaddb_dtset%nqshft 
     ABI_ALLOCATE(qshft,(3,nqshft))
!    TODO this has to be fixed  there is a small inconsistency in the dimension of q1shft
     qshft(:,1:nqshft)=anaddb_dtset%q1shft(:,1:nqshft)
     !
   else ! Dense q-mesh used for the Fourier interpolation. 
     ngqpt(1:3,imesh)=anaddb_dtset%ng2qpt(1:3)
     nqshft=1 !always 1 
     ABI_ALLOCATE(qshft,(3,nqshft))
     qshft(:,1)=anaddb_dtset%q2shft(:)  ! FIXME small inconsistency in the dimension of q1shft
   end if 

   nqpt_max=(ngqpt(1,imesh)*ngqpt(2,imesh)*ngqpt(3,imesh)*nqshft)/facbrv
   ABI_ALLOCATE(qibz,(3,nqpt_max))
   ABI_ALLOCATE(qbz,(3,nqpt_max))

   qptrlatt(:,:)=0
   qptrlatt(1,1)=ngqpt(1,imesh)
   qptrlatt(2,2)=ngqpt(2,imesh)
   qptrlatt(3,3)=ngqpt(3,imesh)
   option=1 
   !
   ! here I noticed a problem in the declaration of q1shft in the anaddb datatype 
   ! FIXME we write on unit std_out just to avoid problem with automatic tests
   call smpbz(brav1,std_out,qptrlatt,nqpt_max,nqbz,nqshft,option,qshft,qbz)
   !  
   !  Reduce the number of such points by symmetrization.
   ABI_ALLOCATE(ibz2bz,(nqbz))
   ABI_ALLOCATE(wtq,(nqbz))
   ABI_ALLOCATE(wtq_folded,(nqbz))
   wtq(:)=one/nqbz         ! Weights sum up to one
   timrev=1; option=1     ! TODO timrev should be input 
   !
   ! This call will set PHdos%nqibz
   call symkpt(chksymbreak0,gmet,ibz2bz,std_out,qbz,nqbz,PHdos%nqibz,nsym,symrec,timrev,wtq,wtq_folded)
   write(std_out,*) 'PHdos%nqibz = ', PHdos%nqibz

   ABI_ALLOCATE(wtqibz,(PHdos%nqibz))
   do iq_ibz=1,PHdos%nqibz
     wtqibz(iq_ibz)=wtq_folded(ibz2bz(iq_ibz))
     qibz(:,iq_ibz)=qbz(:,ibz2bz(iq_ibz))
   end do
   ABI_DEALLOCATE(wtq_folded)
   ABI_DEALLOCATE(qshft)

   if (prtdos==2.and.imesh==2) then
     !    
     ! Second mesh with tetrahedron method
     ! convert kptrlatt to double and invert, qlatt here refer to the shortest qpt vectors
     rlatt(:,:)=qptrlatt(:,:)
     call matr3inv(rlatt,qlatt)

     ABI_ALLOCATE(qshft,(3,nqshft))
     qshft(:,1)=anaddb_dtset%q2shft(:)  ! FIXME small inconsistency in the dimension of q1shft
     nkpt_fullbz=nqbz 
     ABI_ALLOCATE(bz2ibz,(nkpt_fullbz))
     ABI_ALLOCATE(kpt_fullbz,(3,nkpt_fullbz))
     !    
     ! Make full kpoint grid and get equivalence to irred kpoints.
     ! This routines scales badly wrt nkpt_fullbz, should introduce checl on the norm.
     call get_full_kgrid(bz2ibz,qlatt,qibz,kpt_fullbz,qptrlatt,PHdos%nqibz,nkpt_fullbz,nqshft,nsym,qshft,symrel)
     !    
     ! Get tetrahedra, ie indexes of the full q-points at their summits
     ! tetra_full(:,1,i) contains the irred qpt  number
     ! tetra_full(:,2,i) contains the full  qpt number
     ! tetra_wrap(:,:,i) contains a flag to wrap q-points outside the IBZ (+-1) to get the irreducible tetrahedra
     ! the number of equivalent tetrahedra is counted in tetra_mult and the inequivalent few (ntetra < mtetra) are 
     ! packed into the beginning of tetra_full
     mtetra=6*nqbz
     ABI_ALLOCATE(tetra_full,(4,2,mtetra))
     ABI_ALLOCATE(tetra_wrap,(3,4,mtetra))
     ABI_ALLOCATE(tetra_mult,(mtetra))
     

     call init_tetra(bz2ibz, gprimd, qlatt, kpt_fullbz, nqbz, tetrahedra_q, ierr, errstr)
     ABI_CHECK(ierr==0,errstr)
     ABI_DEALLOCATE(bz2ibz)
     !    
     ! Allocate arrays used to store the entire spectrum, Required to calculate tetra weights.
     ABI_ALLOCATE(full_phfrq,(3*natom,PHdos%nqibz))
     ABI_ALLOCATE(full_eigvec,(2,3,natom,3*natom,PHdos%nqibz))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,'out-of-memory in full_eigvec')
   end if  ! prtdos==2.and.imesh==2
   !    
   ! This infinite loop is used to be sure that the frequency mesh is large enough to contain 
   ! the entire phonon spectrum. The mesh is enlarged if, during the Fourier interpolation,
   ! a phonon frequency turns out to be outside the interval [omega_min:omega_max]
   do 
     out_of_bounds=.FALSE.
     if (associated(PHdos%omega))  then
       ABI_DEALLOCATE(PHdos%omega)
     end if
     if (associated(PHdos%phdos))  then
       ABI_DEALLOCATE(PHdos%phdos)
     end if
     if (associated(PHdos%pjdos))  then
       ABI_DEALLOCATE(PHdos%pjdos)
     end if
     !
     ! Frequency mesh.
     PHdos%omega_min=low_bound; if (ABS(PHdos%omega_min)<tol5) PHdos%omega_min=-tol5
     PHdos%omega_max=upr_bound 
     PHdos%nomega=NINT((PHdos%omega_max-PHdos%omega_min)/PHdos%omega_step)+1
     PHdos%nomega=MAX(6,PHdos%nomega) ! Ensure Simpson integration will be ok

     ABI_ALLOCATE(PHdos%omega,(PHdos%nomega))
     do io=1,PHdos%nomega
       PHdos%omega(io)=PHdos%omega_min+PHdos%omega_step*(io-1)
     end do

     if (imesh/=1) then 
       write(std_out,*)&
&         'nomega = ',PHdos%nomega,' omega_min [cm-1] =',PHdos%omega_min*Ha_cmm1,' omega_max [cm-1] =',PHdos%omega_max*Ha_cmm1
     end if 

     ABI_ALLOCATE(PHdos%phdos,(PHdos%nomega))
     PHdos%phdos=zero
     ABI_ALLOCATE(PHdos%pjdos,(PHdos%nomega,3,natom))
     PHdos%pjdos=zero
     !    
     ! === Sum over irreducible q-points ===
     do iq_ibz=1,PHdos%nqibz
       qphon(:)=qibz(:,iq_ibz); qphnrm=one
       !      
       ! Get d2cart using interatomic forces and the long-range coulomb interaction through Ewald summation
       call gtdyn9(acell,atmfrc,dielt,dipdip,dyewq0,d2cart,gmet,gprim,&
&       mpert,natom,nrpt,qphnrm,qphon,rmet,rprim,rpt,trans,ucvol,wghatm,xred,zeff)
       !      
       ! Get eigenvectors and eigenvalues of the dynamical matrix, eigvec are normalized to one
       call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,mpert,msym,natom,nsym,ntypat,&
&       phfrq,qphnrm,qphon,rprimd,symdynmat,symrel,typat,ucvol)
       
       dum=MINVAL(phfrq); PHdos%omega_min=MIN(PHdos%omega_min,dum)
       dum=MAXVAL(phfrq); PHdos%omega_max=MAX(PHdos%omega_max,dum)
       out_of_bounds = (PHdos%omega_min<low_bound .or. PHdos%omega_max>upr_bound) 

       if (imesh>1.and..not.out_of_bounds) then
         select case (prtdos)
         case (1) 
           !
           ! Accumulate PHDOS and PJDOS using gaussian method ===
           do imode=1,3*natom 
             do io=1,PHdos%nomega
               xx=(PHdos%omega(io)-phfrq(imode))*gaussfactor
               gaussval=gaussprefactor*exp(-xx*xx)
               PHdos%phdos(io)=PHdos%phdos(io) + wtqibz(iq_ibz)*gaussval
               do iat=1,natom
                 do idir=1,3
                   pnorm=eigvec(1,idir,iat,imode)**2+eigvec(2,idir,iat,imode)**2
                   PHdos%pjdos(io,idir,iat)=PHdos%pjdos(io,idir,iat)+ pnorm*wtqibz(iq_ibz)*gaussval
                 end do
               end do
             end do 
           end do 
           !
         case (2) 
           !
           ! === Tetrahedrons ===
           !  * Save phonon frequencies and eigenvectors. 
           !  Summation is done after the loops over the two meshes.
           full_phfrq(:,iq_ibz)=phfrq(:)
           full_eigvec(:,:,:,:,iq_ibz)=eigvec(:,:,:,:)
         case default
           write(msg,'(a,i0)')" Wrong value for prtdos= ",prtdos
           MSG_ERROR(msg)
         end select
       end if !Second mesh and not out of boundaries
       !
     end do !irred q-points

     if (out_of_bounds) then 
       upr_bound=PHdos%omega_max+ABS(PHdos%omega_max/ten)
       low_bound=PHdos%omega_min-ABS(PHdos%omega_min/ten)
       write(msg,'(3a)')&
&       ' At least one phonon frequency falls outside the frequency mesh chosen',ch10,&
&       ' restarting the calculation with a larger frequency mesh ' 
       if (imesh>1) then
         MSG_COMMENT(msg)
       end if
     else
       EXIT !infinite loop
     end if 
   end do !infinite loop

   ABI_DEALLOCATE(ibz2bz)
   ABI_DEALLOCATE(qibz)
   ABI_DEALLOCATE(qbz)
   ABI_DEALLOCATE(wtq)
   ABI_DEALLOCATE(wtqibz)
 end do !imesh
 ABI_DEALLOCATE(ngqpt)

 if (associated(PHdos%phdos_int))  then
   ABI_DEALLOCATE(PHdos%phdos_int)
 end if
 if (associated(PHdos%pjdos_int))  then
   ABI_DEALLOCATE(PHdos%pjdos_int)
 end if

 ABI_ALLOCATE(PHdos%phdos_int,(PHdos%nomega))
 PHdos%phdos_int=zero 

 if (prtdos==2) then 
   ! === Integrate using tetrahedrons ===
   !  * All the data are contained in full_phfrq and full_eigvec. 
   !  * low_bound and upr_bound contain the entire spectrum calculated on the dense mesh. 
   ABI_ALLOCATE(tmp_phfrq,(PHdos%nqibz))
   ABI_ALLOCATE(tweight,(PHdos%nqibz,PHdos%nomega))
   ABI_ALLOCATE(dtweightde,(PHdos%nqibz,PHdos%nomega))
   ABI_ALLOCATE(PHdos%pjdos_int,(PHdos%nomega,3,natom))
   PHdos%phdos=zero; PHdos%pjdos=zero; PHdos%pjdos_int=zero
   max_occ=one 

   do imode=1,3*natom 
     tmp_phfrq(:)=full_phfrq(imode,:)
     !    
     ! === Calculate general integration weights at each irred kpoint as in Blochl et al PRB 49 16223 ===
     call get_tetra_weight(tmp_phfrq, low_bound, upr_bound,&
&            max_occ, PHdos%nomega, PHdos%nqibz, tetrahedra_q,&
&            tweight,dtweightde)


     do io=1,PHdos%nomega
       do iq_ibz=1,PHdos%nqibz
         PHdos%phdos(io)=PHdos%phdos(io)+dtweightde(iq_ibz,io)
         PHdos%phdos_int(io)=PHdos%phdos_int(io)+tweight(iq_ibz,io)
         do iat=1,natom
           do idir=1,3
             pnorm=full_eigvec(1,idir,iat,imode,iq_ibz)**2 + full_eigvec(2,idir,iat,imode,iq_ibz)**2
             PHdos%pjdos(io,idir,iat)=PHdos%pjdos(io,idir,iat) + pnorm*dtweightde(iq_ibz,io)
             PHdos%pjdos_int(io,idir,iat)=PHdos%pjdos_int(io,idir,iat) + pnorm*tweight(iq_ibz,io)         
           end do
         end do
       end do
     end do

   end do 
   ABI_DEALLOCATE(tmp_phfrq)
   ABI_DEALLOCATE(tweight)
   ABI_DEALLOCATE(dtweightde)
 end if 
 !
 ! =======================
 ! === calculate IPDOS ===
 ! =======================
 if (associated(PHdos%pjdos_xyz_typ))  then
   ABI_DEALLOCATE(PHdos%pjdos_xyz_typ)
 end if
 if (associated(PHdos%pjdos_typ))  then
   ABI_DEALLOCATE(PHdos%pjdos_typ)
 end if
 if (associated(PHdos%pjdos_typ_int))  then
   ABI_DEALLOCATE(PHdos%pjdos_typ_int)
 end if

 ABI_ALLOCATE(PHdos%pjdos_xyz_typ,(PHdos%nomega,3,ntypat))
 PHdos%pjdos_xyz_typ=zero
 ABI_ALLOCATE(PHdos%pjdos_typ    ,(PHdos%nomega,ntypat))
 PHdos%pjdos_typ    =zero
 ABI_ALLOCATE(PHdos%pjdos_typ_int,(PHdos%nomega,ntypat))
 PHdos%pjdos_typ_int=zero

 do iat=1,natom 
   itype=typat(iat)
   do io=1,PHdos%nomega
     PHdos%pjdos_xyz_typ(io,:,itype)=PHdos%pjdos_xyz_typ(io,:,itype)+PHdos%pjdos(io,:,iat)
     PHdos%pjdos_typ(io,itype)=PHdos%pjdos_typ(io,itype)+sum(PHdos%pjdos(io,:,iat))
   end do
   if (prtdos==2) then 
     do io=1,PHdos%nomega
       PHdos%pjdos_typ_int(io,itype)=PHdos%pjdos_typ_int(io,itype)+SUM(PHdos%pjdos_int(io,:,iat))
     end do
   end if 
 end do
 !
 ! Evaluate IDOS using simple simpson integration
 ! TODO should avoid the simpson rule using derf.F90, just to be consistent
 if (prtdos==1) then 
   call simpson_int(PHdos%nomega,PHdos%omega_step,PHdos%phdos,PHdos%phdos_int)
   do itype=1,ntypat
     call simpson_int(PHdos%nomega,PHdos%omega_step,PHdos%pjdos_typ(:,itype),PHdos%pjdos_typ_int(:,itype))
   end do
 end if 

 if (prtdos==2) then
   ABI_DEALLOCATE(tetra_full)
   ABI_DEALLOCATE(tetra_wrap)
   ABI_DEALLOCATE(tetra_mult)
   ABI_DEALLOCATE(full_phfrq)
   ABI_DEALLOCATE(full_eigvec)
 end if

 DBG_EXIT("COLL")

end subroutine mkphdos
!!***

end module m_phdos
!!***
