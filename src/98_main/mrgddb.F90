!!****p* ABINIT/mrgddb
!! NAME
!! mrgddb
!!
!! FUNCTION
!! This code merges the derivative databases.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, SP)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt
!.
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!! The heading of the constituted database is read,
!! then the heading of the temporary database to be added is read,
!! the code check their compatibility, and create a new
!! database that mixes the old and the temporary ones.
!! This process can be iterated.
!! The whole database will be stored in
!! central memory. One could introduce a third mode in which
!! only the temporary DDB is in central memory, while the
!! input DDB is read twice : first to make a table of blocks,
!! counting the final number of blocks, and second to merge
!! the two DDBs. This would save memory.
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,destroy_mpi_enreg,herald,init8,inprep8,leave_new
!!      mblktyp1,mblktyp5,mpi_comm_rank,mpi_comm_size,nullify_mpi_enreg,wrtout
!!      xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program mrgddb
 use defs_basis
 use defs_abitypes
 use m_build_info
 use m_xmpi

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mrgddb'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_72_response
 use interfaces_77_ddb
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif
!Arguments -----------------------------------

!Local variables-------------------------------
!no_abirules
!
! Set array dimensions
!  mddb=maximum number of databases (cannot be made dynamic)
 integer,parameter :: mddb=5000,ddbun=2
 integer :: dummy,dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,dummy7
 integer :: vrsddb,iddb,mblktyp,mblktyptmp,nddb,ierr
!  msym=maximum number of symmetry elements in space group
 integer :: msym

 character(len=24) :: codename
 character(len=fnlen) :: dscrpt
 character(len=fnlen) :: filnam(mddb+1)
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg

!******************************************************************
!BEGIN EXECUTABLE SECTION

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

!Initialize MPI : one should write a separate routine -init_mpi_enreg-
!for doing that !!
 call xmpi_init()

!Default for sequential use
 call nullify_mpi_enreg(mpi_enreg)
 mpi_enreg%world_comm=0
 mpi_enreg%world_group=0
 mpi_enreg%me=0
 mpi_enreg%nproc=1
 mpi_enreg%nproc_atom=1
 mpi_enreg%num_group_fft = 0 ! in some cases not initialized but referenced in xdef_comm.F90
 mpi_enreg%paral_compil=0
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%paral_compil_mpio=0
 mpi_enreg%mode_para="n"
 mpi_enreg%flag_ind_kg_mpi_to_seq = 0
!MG080916 If we want to avoid MPI preprocessing options, %proc_distr should be
!always allocated and
!set to mpi_enreg%me. In such a way we can safely test its value inside loops
!parallelized over k-points
!For the time being, do not remove this line since it is needed in outkss.F90.
!nullify(mpi_enreg%proc_distrb)
!nullify(mpi_enreg%bandfft_kpt,mpi_enreg%tab_kpt_distrib)

!Initialize MPI
#if defined HAVE_MPI
 mpi_enreg%world_comm=xmpi_world
 mpi_enreg%world_group=MPI_GROUP_NULL
 call MPI_COMM_RANK(xmpi_world,mpi_enreg%me,ierr)
 call MPI_COMM_SIZE(xmpi_world,mpi_enreg%nproc,ierr)
!write(std_out,*)' abinit : nproc,me=',mpi_enreg%nproc,mpi_enreg%me
 mpi_enreg%paral_compil=1
#endif

!Signal MPI I/O compilation has been activated
#if defined HAVE_MPI_IO
 mpi_enreg%paral_compil_mpio=1
 if(mpi_enreg%paral_compil==0)then
   write(message,'(6a)') ch10,&
&   ' abinit : ERROR -',ch10,&
&   '  In order to use MPI_IO, you must compile with the MPI flag ',ch10,&
&   '  Action : recompile your code with different CPP flags.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
#endif

!Initialize paral_compil_kpt, actually always equal to paral_compil
!(paral_compil_kpt should be suppressed after big cleaning)
 mpi_enreg%paral_compil_kpt=0
 if(mpi_enreg%paral_compil==1) mpi_enreg%paral_compil_kpt=1

!Other values of mpi_enreg are dataset dependent, and should NOT be initialized
!inside mrgddb.F90.


 codename='MRGDDB'//repeat(' ',18)
 call herald(codename,abinit_version,std_out)
!YP: calling dump_config() makes tests fail => commented
!call dump_config()

!Initialise the code : write heading,
!read names of files, operating mode (also check its value),
!and short description of new database.
 call init8(dscrpt,filnam,mddb,nddb)

!Set the ddb version
 vrsddb=100401

!Evaluate the mblktyp of the databases
 mblktyptmp=1
 do iddb=1,nddb
   call inprep8(dummy,filnam(iddb+1),dummy1,dummy2,mblktyp,&
&   msym,dummy3,dummy4,dummy5,dummy6,ddbun,dummy7,vrsddb)
   if(mblktyp > mblktyptmp) then
     mblktyptmp = mblktyp
   end if
 end do

 mblktyp = mblktyptmp

!Debug
!write(std_out,*),'mblktyp',mblktyp
!enddebug

 if(mblktyp==5) then
!  Memory optimized routine
   call mblktyp5(ddbun,dscrpt,filnam,mddb,msym,nddb,vrsddb)
 else
!  Speed optimized routine
   call mblktyp1(ddbun,dscrpt,filnam,mddb,msym,nddb,vrsddb)
 end if

!**********************************************************************

 write(message, '(a)' )'+mrgddb : the run completed successfully '
 call wrtout(std_out,message,'COLL')


 call destroy_mpi_enreg(mpi_enreg)
 call xmpi_end()

 end program mrgddb
!!***

