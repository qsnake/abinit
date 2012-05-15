!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spectra
!! NAME
!! m_spectra
!!
!! FUNCTION
!! This module contains the definition of the specta_type data type 
!! used to store results related to optical spectra with or without 
!! nonlocal field effects as well as the electron energy loss function
!! for a given q. These quantities are obtained from the dielectric 
!! matrix as calculated in the GW part of ABINIT (screening.F90)
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_spectra

 use m_profiling

 use defs_basis
 use m_errors

 use m_io_tools,      only : get_unit
 use m_fstrings,      only : str_conct

 implicit none

 private
!!***

!!****t* m_spectra/spectra_type
!! NAME
!! spectra_type
!!
!! FUNCTION
!!  Object used to store optical spectra with or without non-local field effects.
!!
!! SOURCE

 type,public :: spectra_type

 !scalars
  integer :: nomega
  
  integer :: nqpts

 !arrays
  real(dp),pointer :: omega(:)  SET2NULL
  ! omega(nomega)
  ! Real frequency mesh for optical spectra.

  real(dp),pointer :: qpts(:,:)  SET2NULL
  ! qpts(3,nqpoints)
  ! The list of q-points used for the spectra

  real(dp),pointer :: eelf(:,:)   SET2NULL
  ! eelf(nomega,nqpoints)
  ! contains the Electron Energy Loss Function i.e. -\Im{ e^{-1}_{G1=0,G2=0}(q-->0,nomega)}

  complex(dpc),pointer :: emacro_lf(:,:)  SET2NULL
  ! emacro_lf(nomega,nqpoints)
  ! contains 1/e^{-1}_{G1=0,G2=0}(q-->0,nomega) (with Local field effects)

  complex(dpc),pointer :: emacro_nlf(:,:)  SET2NULL
  ! emacro_nlf(nomega,nqpoints)
  ! contains e_{G1=0,G2=0}(q-->0,nomega) (without Local field effects)

 end type spectra_type
!!***

 public :: init_spectra       ! Creation method.
 public :: destroy_spectra    ! Destruction method.
 public :: dump_spectra       ! Write results on file.
 public :: repr_dielconst     ! Return info on Macroscopic diel. constant in form of a string.
                                                                                               
 integer,public,parameter :: W_EM_LF  = 1
 integer,public,parameter :: W_EM_NLF = 2
 integer,public,parameter :: W_EELF   = 4

CONTAINS  !========================================================================================
!!***

!!****f* m_spectra/init_spectra
!! NAME
!!  init_spectra
!!
!! FUNCTION
!!  Initialize the object.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_spectra(Spectra,nomega,omega,nqpts,qpts)

 use defs_basis

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_spectra'
!End of the abilint section

 integer,intent(in) :: nomega,nqpts
!arrays
 real(dp),intent(in) :: omega(nomega),qpts(3,nqpts)
 type(spectra_type),intent(out) :: Spectra

! *********************************************************************

 call nullify_spectra_(Spectra)

 Spectra%nomega = nomega
 Spectra%nqpts  = nqpts

 ABI_ALLOCATE(Spectra%qpts,(3,nqpts))
 Spectra%qpts = qpts

 ABI_ALLOCATE(Spectra%omega,(nomega))
 Spectra%omega = omega

 ABI_ALLOCATE(Spectra%emacro_lf,(nomega,nqpts))
 Spectra%emacro_lf = czero

 ABI_ALLOCATE(Spectra%emacro_nlf,(nomega,nqpts))
 Spectra%emacro_nlf = czero

 ABI_ALLOCATE(Spectra%eelf,(nomega,nqpts))
 Spectra%eelf = zero

end subroutine init_spectra
!!***

!----------------------------------------------------------------------

!!****f* m_spectra/nullify_spectra_
!! NAME
!!  nullify_spectra_
!!
!! FUNCTION
!!  Nullify all pointers in the structure {PRIVATE}
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_spectra
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_spectra_(Spectra)

 use defs_basis

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_spectra_'
!End of the abilint section

 type(spectra_type),intent(inout) :: Spectra

! *********************************************************************

 nullify(Spectra%omega     )
 nullify(Spectra%qpts      )
 nullify(Spectra%emacro_lf )
 nullify(Spectra%emacro_nlf)
 nullify(Spectra%eelf      )

end subroutine nullify_spectra_
!!***

!----------------------------------------------------------------------

!!****f* m_spectra/destroy_spectra
!! NAME
!!  destroy_spectra
!!
!! FUNCTION
!!  Deallocate all associated pointers defined in the structure.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_screen,m_screening,screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_spectra(Spectra)

 use defs_basis

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_spectra'
!End of the abilint section

 type(spectra_type),intent(inout) :: Spectra

! *********************************************************************

 if (associated(Spectra%omega     ))   then
   ABI_DEALLOCATE(Spectra%omega)
 end if
 if (associated(Spectra%qpts      ))   then
   ABI_DEALLOCATE(Spectra%qpts)
 end if
 if (associated(Spectra%emacro_lf ))   then
   ABI_DEALLOCATE(Spectra%emacro_lf)
 end if
 if (associated(Spectra%emacro_nlf))   then
   ABI_DEALLOCATE(Spectra%emacro_nlf)
 end if
 if (associated(Spectra%eelf      ))   then
   ABI_DEALLOCATE(Spectra%eelf)
 end if

end subroutine destroy_spectra 
!!***

!----------------------------------------------------------------------

!!****f* m_spectra/dump_spectra
!! NAME
!!  dump_spectra
!!
!! FUNCTION
!!  Write the optical spectra stored in the object on an external formatted file.
!!
!! INPUTS
!!  Spectra=The Object containing the spectra
!!  write_bits=Positive integer defining the quantities to be written (bit representation is used)
!!  fname=Name of the file to be written.
!!
!! OUTPUT
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine dump_spectra(Spectra,write_bits,fname)

 use defs_basis

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dump_spectra'
!End of the abilint section

 integer,intent(in) :: write_bits
 character(len=*),intent(in) :: fname
!arrays
 type(spectra_type),intent(in) :: Spectra

!Local variables-------------------------------
!scalars
 integer :: unt,io,iqpt
 real(dp) :: mino,maxo
 character(len=100) :: fmt

! *********************************************************************

 if (write_bits<0) RETURN

 mino = MINVAL(Spectra%omega)
 maxo = MAXVAL(Spectra%omega)

 unt=get_unit() 
 open(unit=unt,file=fname)

 !write(unt,'(a,i5,2(a,f9.1),a)')'# nomega : ',Spectra%nomega,' from ',mino*Ha_eV,' up to ',maxo*Ha_eV,' [eV] '

 if ( IAND(write_bits,W_EM_NLF ) == W_EM_NLF ) then
   write(unt,'(a)')'#'
   write(unt,'(a)')'# Macroscopic Dielectric Function without local fields'
   call dump_qlist()
   write(unt,'(a)')'# Omega [eV]    Re epsilon_M       IM eps_M '
   write(fmt,'(a,i3,a)') '(1x,f7.3,7x,',Spectra%nqpts,'(2(e12.4),2x))'
   do io=1,Spectra%nomega
     write(unt,fmt) Spectra%omega(io)*Ha_eV, ( Spectra%emacro_nlf(io,iqpt), iqpt=1,Spectra%nqpts )
   end do
 end if
 !
 if ( IAND(write_bits,W_EM_LF ) == W_EM_LF ) then 
   write(unt,'(a)')'#'
   write(unt,'(a)')'# Macroscopic Dielectric Function with local fields included'
   call dump_qlist()
   write(unt,'(a)')'# Omega [eV]    Re epsilon_M       Im eps_M '
   write(fmt,'(a,i3,a)') '(1x,f7.3,7x,',Spectra%nqpts,'(2(e12.4),2x))'
   do io=1,Spectra%nomega
     write(unt,fmt) Spectra%omega(io)*Ha_eV, ( Spectra%emacro_lf(io,iqpt), iqpt=1,Spectra%nqpts )
   end do
 end if
 !
 if ( IAND(write_bits,W_EELF) == W_EELF) then
   write(unt,'(a)')'#'
   write(unt,'(a)')'# Electron Energy Loss Function -Im(1/epsilon_M)'
   call dump_qlist()
   write(unt,'(a)')'# Omega [eV]    -Im(1/epsilon_M)'
   write(fmt,'(a,i3,a)') '(1x,f7.3,7x,',Spectra%nqpts,'(e12.4,2x))'
   do io=1,Spectra%nomega
     write(unt,fmt) Spectra%omega(io)*Ha_eV, ( Spectra%eelf(io,iqpt), iqpt=1,Spectra%nqpts ) ! -AIMAG(chi0(1,1,io))
   end do
 end if

 close(unt)

CONTAINS
!!***

!!****f* m_spectra/dump_Qlist
!! NAME
!!  dump_Qlist
!!
!! FUNCTION
!!  Helper function used to write the list of q-points used for the long-wavelength limit.
!!
!! INPUTS
!!  Spectra=The Object containing the spectra
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_spectra
!!
!! CHILDREN
!!
!! SOURCE

 subroutine dump_Qlist()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dump_Qlist'
!End of the abilint section

  integer :: iqpt
  write(unt,'(a,i3)')'# Q-point list, No. ',Spectra%nqpts
  do iqpt=1,Spectra%nqpts
    write(unt,'(a,i3,a,3f9.6,a)')'# ',iqpt,')  [',Spectra%qpts(:,iqpt),'] r.l.u. '
  end do
 end subroutine dump_Qlist

end subroutine dump_spectra 
!!***

!----------------------------------------------------------------------

!!****f* m_spectra/repr_dielconst
!! NAME
!!  repr_dielconst
!!
!! FUNCTION
!!  Returns a string reporting info on the calculated dielectric constant.
!!
!! INPUTS
!!  Spectra=The Object containing the spectra
!!
!! OUTPUT
!!
!! PARENTS
!!      m_screening,screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine repr_dielconst(Spectra,str)

 use defs_basis

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'repr_dielconst'
!End of the abilint section

 type(spectra_type),intent(in) :: Spectra
 character(len=*),intent(out) :: str

!Local variables-------------------------------
!scalars
 integer :: iqpt
 real(dp) :: epsilon0,epsilon0_nlf
 character(len=500) :: msg

! *********************************************************************

 !istatic = -1
 !do io = Spectra%nomega 
 ! if (ABS(REAL(Spectra%omega(io)))<1.e-3.and.ABS(AIMAG(Spectra%omega(io)))<1.e-3) then
 !  istatic = io 
 !  EXIT
 ! end if
 !end do

 str = ""
 do iqpt=1,Spectra%nqpts
   epsilon0    = REAL(Spectra%emacro_lf (1,iqpt))
   epsilon0_nlf= REAL(Spectra%emacro_nlf(1,iqpt))
   write(msg,'(a,3f9.6,a)')' For q-point: ',Spectra%qpts(:,iqpt),ch10
   str = str_conct(str,msg)
   write(msg,'(1x,a,f8.4,a)')' dielectric constant = ',epsilon0,ch10
   str = str_conct(str,msg)
   write(msg,'(1x,a,f8.4,a)')' dielectric constant without local fields = ',epsilon0_nlf,ch10
   str = str_conct(str,msg)
 end do

end subroutine repr_dielconst

!----------------------------------------------------------------------

END MODULE m_spectra
!!***
