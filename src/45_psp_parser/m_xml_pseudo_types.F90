!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xml_pseudo_types
!! NAME
!! m_xml_pseudo_types
!!
!! FUNCTION
!!  Module with pseudopotential structures for SIESTA/XML format
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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


module m_xml_pseudo_types

 use m_profiling

implicit none
!
! Data structures for a prototype pseudopotential
!
integer, parameter, private    :: MAXN_POTS = 8
integer, parameter, private    :: dp = selected_real_kind(14)
!
public  :: dump_pseudo
!
!-----------------------------------------------------------
type, public :: grid_t
!
!     It should be possible to represent both log and linear
!     grids with a few parameters here.
!
      character(len=20)              :: type
      real(kind=dp)                  :: scale
      real(kind=dp)                  :: step
      integer                        :: npts
end type grid_t
!
type, public :: radfunc_t
      type(grid_t)                            :: grid
      real(kind=dp), dimension(:), pointer    :: data
end type radfunc_t

type, public :: vps_t
      character(len=1)               :: l
      integer                        :: n
      integer                        :: spin
      real(kind=dp)                  :: occupation
      real(kind=dp)                  :: cutoff
      type(radfunc_t)                :: V
end type vps_t

type, public :: pswf_t
      character(len=1)               :: l
      integer                        :: n
      integer                        :: spin
      type(radfunc_t)                :: V
end type pswf_t

type, public :: header_t
        character(len=2)        :: symbol
        real(kind=dp)           :: atomicnumber
        real(kind=dp)           :: zval
        character(len=10)       :: creator
        character(len=10)       :: date
        character(len=40)       :: flavor
        logical                 :: relativistic
        logical                 :: polarized
        character(len=30)       :: xcfunctionaltype
        character(len=30)       :: xcfunctionalparametrization
        character(len=4)        :: core_corrections
end type header_t

type, public :: pseudo_t
      type(header_t)                     :: header
      integer                            :: npots
      integer                            :: npswfs
      integer                            :: npots_down
      integer                            :: npots_up
      type(vps_t), dimension(MAXN_POTS)  :: pot
      type(pswf_t), dimension(MAXN_POTS) :: pswf
      type(radfunc_t)                    :: core_charge
      type(radfunc_t)                    :: valence_charge
end type pseudo_t


CONTAINS !===============================================
!!***

!!****f* m_xml_pseudo_types/dump_pseudo
!! NAME
!! dump_pseudo
!!
!! FUNCTION
!!  Additional routine for dumping output about pseudopotential contents.
!!
!! INPUTS
!!  pseudo = structure containing pseudopotential info to be dumped to stdout.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine dump_pseudo(pseudo)

use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dump_pseudo'
!End of the abilint section

implicit none 

type(pseudo_t), intent(in), target   :: pseudo

integer  :: i
type(vps_t), pointer :: pp
type(pswf_t), pointer :: pw
type(radfunc_t), pointer :: rp

write(std_out,*) "---PSEUDO data:"

write(std_out,*) pseudo%header%xcfunctionaltype
do i = 1, pseudo%npots
      pp =>  pseudo%pot(i)
      rp =>  pseudo%pot(i)%V
      write(std_out,*) "VPS ", i, " angular momentum: ", pp%l
      write(std_out,*) "                 n: ", pp%n
      write(std_out,*) "                 occupation: ", pp%occupation
      write(std_out,*) "                 cutoff: ", pp%cutoff
      write(std_out,*) "                 spin: ", pp%spin
      write(std_out,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
end do
do i = 1, pseudo%npswfs
      pw =>  pseudo%pswf(i)
      rp =>  pseudo%pswf(i)%V
      write(std_out,*) "PSWF ", i, " angular momentum: ", pw%l
      write(std_out,*) "                 n: ", pw%n
      write(std_out,*) "                 spin: ", pw%spin
      write(std_out,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
end do
rp => pseudo%valence_charge
write(std_out,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)
rp => pseudo%core_charge
write(std_out,*) "grid data: ", rp%grid%npts, rp%grid%scale, rp%data(1)

end subroutine dump_pseudo

end module m_xml_pseudo_types
!!***
