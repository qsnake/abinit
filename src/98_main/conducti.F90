!{\src2tex{textfont=tt}}
!!****p* ABINIT/conducti
!! NAME
!! conducti
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula.
!!
!! COPYRIGHT
!! Copyright (C) 2006-2012 ABINIT group (FJ,SMazevet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!!
!! CHILDREN
!!      conducti_nc, conducti_paw
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,conducti_nc,conducti_paw,conducti_paw_core,emispec_paw
!!      linear_optics_paw,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program conducti

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'conducti'
 use interfaces_66_paw
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
 integer :: incpaw,nproc,comm
 character(len=fnlen) :: filnam,filnam_out
 !character(len=500) :: msg
 type(MPI_type) :: mpi_enreg_seq ! needed for calling rwwf

! *********************************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

 call xmpi_init()

 comm = xmpi_world
 nproc = xcomm_size(comm)

 if ( nproc > 1) then
   MSG_ERROR("conducti is not parallelized: run with one processor.")
 end if

!Read data file name
 write(std_out,'(a)')' Please, give the name of the data file ...'
 read(5, '(a)')filnam
 write(std_out,'(a)')' The name of the data file is :',trim(filnam)
 open(15,file=filnam,form='formatted')
 rewind(15)
 read(15,*) incpaw
 close(15)
 write(std_out,'(a)')' Give the name of the output file ...'
 read(5, '(a)') filnam_out
 write(std_out,'(a)')' The name of the output file is :',filnam_out

 if (incpaw==1) then
   call conducti_nc(filnam,filnam_out,mpi_enreg_seq)
 elseif (incpaw==2) then
   call conducti_paw(filnam,filnam_out,mpi_enreg_seq)
 elseif (incpaw==3) then
   call linear_optics_paw(filnam,filnam_out,mpi_enreg_seq)
 elseif (incpaw==4) then
   call conducti_paw(filnam,filnam_out,mpi_enreg_seq)
   call conducti_paw_core(filnam,filnam_out,mpi_enreg_seq)
   call emispec_paw(filnam,filnam_out,mpi_enreg_seq)
 elseif (incpaw==5) then
   call conducti_paw_core(filnam,filnam_out,mpi_enreg_seq)
 elseif (incpaw==6) then
   call emispec_paw(filnam,filnam_out,mpi_enreg_seq)
 end if

 call xmpi_end()

 end program conducti
!!***
