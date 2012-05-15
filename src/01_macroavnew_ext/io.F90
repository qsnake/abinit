!
! Copyright (C) Alberto Garcia, 1996, 1997, 1998
! This file is distributed under the terms of the
! GNU General Public License, see ~abinit/COPYING
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module implements an interface to the FORTRAN logical unit
! system. Based on code by Richard Maine.
!
!
! Alberto Garcia, December 30, 1996
! Rewritten as a single subroutine 
! with multiple entry points, March 7, 1998
!
! This scheme is actually the closest in spirit to f90 modules, but
! in perfectly legal f77.
!
!---------------------------------------------------------------
!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

      subroutine io
!
!     Logical unit management. Units 0 to min_lun-1 are "reserved",
!     since most of the "typical" files (output, etc) use them.
!
!     Logical units min_lun to min_max are managed by this module.

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'io'
!End of the abilint section

      implicit none 
!
!----------------------------------------------------------------
!     Module variables
!
      integer stdout, stderr
      integer min_lun, max_lun, nunits
      parameter (min_lun=10, max_lun=99, nunits=max_lun-min_lun+1)
      logical lun_is_free(min_lun:max_lun)

      save stdout, stderr, lun_is_free 
!-----------------------------------------------------------------
!
!     Internal and dummy variables
!
      integer i, unit, lun, iostat
      logical used, named, opened
      character filename*50, form*11
!
!-----------------------------------------------------------------
!     Initialization section
!
      data lun_is_free /nunits*.true./
      data stdout, stderr /6,0/ 
!-----------------------------------------------------------------
!
!     Executable routines
!
!     Simple interfaces to modify standard units
!
      entry io_seterr(unit)
      stderr = unit
      return
      entry io_setout(unit)
      stdout = unit
      return

      entry io_geterr(unit)
      unit = stderr
      return
      entry io_getout(unit)
      unit = stdout
      return 
!
!------------------------------------------------------------------     
!
!     Logical unit management
!
      entry io_assign(lun) 
!
!     Looks for a free unit and assigns it to lun
!
      do lun= min_lun, max_lun
         if (lun_is_free(lun)) then
            inquire(unit=lun, opened=used, iostat=iostat)
            if (iostat /= 0) used = .true.
            lun_is_free(lun) = .false.
            if (.not. used) return
         endif
      enddo
      write(stderr,'(a)') 'No luns available in io_assign'
      stop 'LUN' 
!
!===
!
      entry io_reserve(lun) 
!
!     Useful to specify that one needs to use a particular unit number
!
!     For example, assume some legacy code expects to work with unit 15:
!
!     call io_reserve(15)   ! this call at the beginning of the program
!     ...
!     open(15,....)
!
      inquire(unit=lun, opened=used, iostat=iostat)
      if (iostat /= 0) used = .true.
      if (used) then
         write(stderr,'(a,i3,a)')&
     &        'Cannot reserve unit',lun,'. Already connected'
         stop 'LUN'
      endif
      if (lun >= min_lun .and. lun <= max_lun)&
     &                      lun_is_free(lun) = .false.

      return 
!
!===
!
      entry io_close(lun) 
!
!     Use this routine instead of a simple close!!
!
      close(lun)
      if (lun >= min_lun .and. lun <= max_lun)&
     &                     lun_is_free(lun) = .true.
      return 
!
!===
!
      entry io_status 
!
!     Prints a list of the connected logical units and the names of
!     the associated files
!

      write(stdout,'(a)') '******** io_status ********'
      do i = 0, max_lun
         inquire(i,opened=opened,named=named,name=filename,&
     &           form=form,iostat=iostat)
         if (iostat == 0) then
            if (opened) then
               if (named) then
                  write(stdout,'(i4,5x,a,5x,a)') i, form, filename
               else
                  write(stdout,'(i4,5x,a,5x,a)') i, form, 'No name available'
               endif
            endif
         else
            write(stdout,'(i4,5x,a,5x,a)') i, 'Iostat error'
         endif
      enddo
      write(stdout,'(a)') '********           ********'

      return

      end






