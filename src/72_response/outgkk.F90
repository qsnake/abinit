!{\src2tex{textfont=tt}}
!!****f* ABINIT/outgkk
!! NAME
!! outgkk
!!
!! FUNCTION
!! output gkk file for one perturbation (used for elphon calculations in anaddb)
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  bantot0 = total number of bands for all kpoints
!!  bantot1 = total number of matrix elements for 1st order eigenvalues
!!  eigen0 = GS eigenvalues
!!  eigen1 = response function 1st order eigenvalue matrix
!!  hdr0 = GS header
!!  hdr1 = RF header
!!  mpi_enreg=informations about MPI parallelization
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      hdr_io,wrtout,xmaster_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine outgkk(bantot0,bantot1,outfile,eigen0,eigen1,hdr0,hdr1,mpi_enreg,phasecg)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

 use m_io_tools,    only : get_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outgkk'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bantot0,bantot1
 character(len=fnlen),intent(in) :: outfile
 type(MPI_type),intent(inout) :: mpi_enreg
 type(hdr_type),intent(inout) :: hdr0,hdr1
!arrays
 real(dp),intent(in) :: eigen0(bantot0),eigen1(2*bantot1)
 real(dp),intent(in) :: phasecg(2,bantot1)

!Local variables-------------------------------
!scalars
 integer :: fform,iband,ikpt,isppol,master,me,ntot,rdwrout,unitout
 integer :: iband_off, mband
 real(dp), allocatable :: tmpeig(:)
 !character(len=500) :: message

! *************************************************************************

!only master should be writing to disk
!Init me
 call xme_init(mpi_enreg,me)
!Define master
 call xmaster_init(mpi_enreg,master)

 if (master /= me) return

 call wrtout(std_out,' writing gkk file: '//outfile,"COLL")

!initializations
 rdwrout = 6
 unitout = get_unit()
 fform = 42
 ntot = 1

!open gkk file
 open (unit=unitout,file=outfile,form='unformatted',status='unknown')

!output GS header
 call hdr_io(fform,hdr0,rdwrout,unitout)

!output GS eigenvalues
 iband=0
 do isppol=1,hdr0%nsppol
   do ikpt=1,hdr0%nkpt
     write (unitout) eigen0(iband+1:iband+hdr0%nband(ikpt))
     iband=iband+hdr0%nband(ikpt)
   end do
 end do

!output number of gkk in this file (1)
 write (unitout) ntot

!output RF header
 call hdr_io(fform,hdr1,rdwrout,unitout)

!output RF eigenvalues
 mband = maxval(hdr1%nband(:))
 ABI_ALLOCATE(tmpeig,(2*mband**2))
 iband_off = 0
 tmpeig(1) = phasecg(1, 1)
 do isppol = 1, hdr1%nsppol
   do ikpt = 1, hdr1%nkpt
     tmpeig = zero
     do iband = 1, hdr1%nband(ikpt)**2
!      multiply by conjugate of phasecg with minus sign in imag part instead of real
!      tmpeig (2*(iband-1)+1) = eigen1(2*(iband_off+iband-1)+1) * phasecg(1, iband_off+iband) &
!      &                             - eigen1(2*(iband_off+iband-1)+2) * phasecg(2, iband_off+iband)
!      tmpeig (2*(iband-1)+2) = eigen1(2*(iband_off+iband-1)+2) * phasecg(1, iband_off+iband)&
!      &                             + eigen1(2*(iband_off+iband-1)+1) * phasecg(2, iband_off+iband)
       tmpeig (2*(iband-1)+1) = eigen1(2*(iband_off+iband-1)+1)
       tmpeig (2*(iband-1)+2) = eigen1(2*(iband_off+iband-1)+2)
     end do
     write (unitout) tmpeig(1:2*hdr1%nband(ikpt)**2)
     iband_off = iband_off + hdr1%nband(ikpt)**2
   end do
 end do
 ABI_DEALLOCATE(tmpeig)

!close gkk file
 close (unitout)

end subroutine outgkk
!!***
