!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_mover
!! NAME
!! defs_mover
!!
!! FUNCTION
!! This module contains definition of an structured datatypes used
!! by the mover routine.
!! If you want to add one new datatype, please, examine first
!! whether another datatype might meet your need (e.g. adding some
!! records to it).
!! Then, if you are sure your new structured datatype is needed,
!! write it here, and DOCUMENT it properly (not all datastructure
!! here are well documented, it is a shame ...).
!! Do not forget : you will likely be the major winner if you
!! document properly.
!! Proper documentation of a structured datatype means :
!!  (1) Mention it in the list just below
!!  (2) Describe it in the NOTES section
!!  (3) Put it in alphabetical order in the the main section of
!!      this module
!!  (4) Document each of its records, except if they are described
!!      elsewhere
!!      (this exception is typically the case of the dataset
!!      associated with input variables, for which there is a help
!!      file)
!!
!! List of datatypes :
!! * ab_movehistory  : Historical record of previous states of
!!                     during the optimization of geometry
!! * ab_movetype     : Subset of dtset with data for move atoms
!! * ab_xfh_type     : Old history style
!! * mttk_type       : For pred_isothermal
!!
!! List of subroutines :
!! * ab_movehistory_bcast    : Broadcast an ab_movehistory (with MPI)
!! * ab_movehistory_nullify  : Nullify an ab_movehistory
!! * ab_movetype_nullify     : Nullify an ab_movetype
!! * ab_movetype_print       : Print the contents of an ab_movetype
!!
!! COPYRIGHT
!! Copyright (C) 2001-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 module defs_mover

 use m_profiling

 use defs_basis
 use m_delocint

 implicit none

 integer, parameter :: mover_BEFORE=0
 integer, parameter :: mover_AFTER=1

!Structures
!!***

!----------------------------------------------------------------------

!!****t* defs_mover/mttk_type
!! NAME
!! mttk_type
!!
!! FUNCTION
!! For Martyna et al. (TTK) reversible MD integration scheme and related data
!!
!! SOURCE

 type mttk_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Real (double precision) scalars

   real(dp) :: glogv
    !Logarithm of the volume

   real(dp) :: vlogv
    !Derivative of logv

!Real (double precision) arrays

  real(dp) :: gboxg(3,3)
   !Imbalance in pressure (see paper)

  real(dp) :: vboxg(3,3)
   !Velocity of log rprimd (see paper)

  real(dp), pointer :: glogs(:)
   ! glogs(nnos)
   ! Imbalance of kinetic energy

  real(dp), pointer :: vlogs(:)
   ! vlogs(nnos)
   ! Velocities of thermostat variables

  real(dp), pointer :: xlogs(:)
   ! xlogs(nnos)
   ! Positions of thermostat variables

 end type mttk_type
!!***

!----------------------------------------------------------------------

!!****t* defs_mover/ab_movehistory
!! NAME
!! ab_movehistory
!!
!! FUNCTION
!! This function has several vectors, and index scalars to store
!! a proper history of previous evaluations of forces and
!! stresses,velocities,positions and energies
!!
!! NOTES
!! The vectors are not allocated because in some cases
!! not all the vectors are needed, in particular a history
!! of stresses is only needed if optcell/=0, and a history
!! of velocities is needed for ionmov==1
!!
!! Store acell, rprimd and strten even with optcell/=0
!! represent a waste of 12x (dp)[Usually 8 Bytes] per
!! iteration, the reason to store all the records is
!! because some routines (eg bfgs.F90) uses the metric (gmet)
!! for initialize the hessian and we need rprimd for that.
!!
!! SOURCE

type ab_movehistory

! scalars
! Index of the last element on all records
integer :: ihist
! Maximun size of the historical records
integer :: mxhist
! Booleans to know if some arrays are changing
logical :: isVused  ! If velocities are changing
logical :: isARused ! If Acell and Rprimd are changing

! arrays
! Vector of (x,y,z)X(mxhist)
real(dp), pointer :: histA(:,:)
! Vector of (mxhist) values of energy
real(dp), pointer :: histE(:)
! Vector of (mxhist) values of ionic kinetic energy
real(dp), pointer :: histEk(:)
! Vector of (mxhist) values of time (relevant for MD calculations)
real(dp), pointer :: histT(:)
! Vector of (x,y,z)X(x,y,z)X(mxhist)
real(dp), pointer :: histR(:,:,:)
! Vector of (stress [6])X(mxhist)
real(dp), pointer :: histS(:,:)
! Vector of (x,y,z)X(natom)X(mxhist) values of velocity
real(dp), pointer :: histV(:,:,:)
! Vector of (x,y,z)X(natom)X(xcart,xred,fcart,fred)X(mxhist)
real(dp), pointer :: histXF(:,:,:,:)

end type ab_movehistory
!!***

!!****t* defs_mover/ab_forstr
!! NAME
!! ab_forstr
!!
!! FUNCTION
!! Store forces, stress and energy, cartesian and reduced forces
!! one scalar for energy and 6 element array for stress
!!
!! NOTES
!!
!! SOURCE

type ab_forstr

! scalars
real(dp) :: etotal     ! Total energy
! arrays
real(dp),pointer :: fcart(:,:) ! Cartesian forces
real(dp),pointer :: fred(:,:)  ! Reduced forces
real(dp) :: strten(6)  ! Stress tensor (Symmetrical 3x3 matrix)

end type ab_forstr
!!***

!!****t* defs_mover/ab_movetype
!! NAME
!! ab_movetype
!!
!! FUNCTION
!! This datatype has the purpouse of store all the data taked
!! usually from dtset needed for the different predictors
!! to update positions, acell, etc.
!!
!! NOTES
!!  At present 32 variables are present in ab_movetype
!!  if a new variable is added in ab_movetype it should
!!  be added also for nullify in ab_movetype_nullify
!!
!! STATS
!!  integer         => 13
!!  real            =>  7
!!  integer array   =>  4
!!  real array      =>  6
!!  character array =>  1
!!  structures      =>  1
!!  TOTAL              34
!!
!! SOURCE

type ab_movetype

! scalars
! Delay of Permutation (Used by pred_langevin only)
integer,pointer  :: delayperm
! DIIS memory (Used by pred_diisrelax only)
integer,pointer  :: diismemory
! Geometry Optimization Precondition option
integer,pointer  :: goprecon
! include a JELLium SLAB in the cell
integer,pointer  :: jellslab
! Number of ATOMs
integer,pointer  :: natom
! Number of CONstraint EQuations
integer,pointer  :: nconeq
! Use by pred_isothermal only
integer,pointer  :: nnos
! Number of SYMmetry operations
integer,pointer  :: nsym
! Number of Types of atoms
integer,pointer  :: ntypat
! OPTimize the CELL shape and dimensions
integer,pointer  :: optcell
! RESTART Xcart and Fred
integer,pointer  :: restartxf
! Sign of Permutation (Used by pred_langevin only)
integer,pointer  :: signperm
! Ion movement
integer,pointer  :: ionmov

! Use by pred_isothermal only
real(dp),pointer :: bmass
! Delta Time for IONs
real(dp),pointer :: dtion
! Used by pred_langevin only
real(dp),pointer :: friction
! Used by pred_langevin only
real(dp),pointer :: mdwall
! Used by pred_nose only
real(dp),pointer :: noseinert
! STRess PRECONditioner
real(dp),pointer :: strprecon
! VIScosity
real(dp),pointer :: vis

! arrays
! Indices of AToms that are FIXed
integer,pointer  :: iatfix(:,:)         ! iatfix(3,natom)
! SYMmetry in REaL space
integer,pointer  :: symrel(:,:,:)       ! symrel(3,3,nsym)
! TYPe of ATom
integer,pointer  :: typat(:)            ! typat(natom)
! PRTint ATom LIST
integer,pointer  :: prtatlist(:)        ! prtatlist(natom)

! Mass of each atom (NOT IN DTSET)
real(dp),pointer :: amass(:)            ! amass(natom)
! Geometry Optimization Preconditioner PaRaMeters
real(dp),pointer :: goprecprm(:)
! Molecular Dynamics Initial and Final Temperature
real(dp),pointer :: mdtemp(:)           ! mdtemp(2) (initial,final)
! STRess TARGET
real(dp),pointer :: strtarget(:)        ! strtarget(6)
! Use by pred_isothermal only
real(dp),pointer :: qmass(:)
! Z number of each NUCLeus
real(dp),pointer :: znucl(:)            ! znucl(npsp)

! Filename for Hessian matrix
character(len=fnlen), pointer :: fnameabi_hes
! Filename for _HIST file
character(len=fnlen), pointer :: filnam_ds(:)   ! dtfil%filnam_ds(5)

! structure for delocalized internal coordinates
type(ab_delocint) :: deloc

end type ab_movetype
!!***

!!****t* defs_mover/ab_xfh_type
!! NAME
!! ab_xfh_type
!!
!! FUNCTION
!! Datatype with the old structure for storing history
!! used in gstate and brdmin,delocint, and others
!!
!! NOTES
!!
!! This is a transitional structure, to bridge between
!! the old code and the new one base on ab_movehistory
!!
!! SOURCE

type ab_xfh_type

!  mxfh = last dimension of the xfhist array
!  nxfh = actual number of (x,f) history pairs, see xfhist array
integer :: nxfh,nxfhr,mxfh

!  xfhist(3,natom+4,2,mxfh) = (x,f) history array, also including
!   rprim and stress
real(dp),pointer :: xfhist(:,:,:,:)

end type ab_xfh_type
!!***

!!****t* defs_mover/go_bonds
!! NAME
!! go_bonds
!!
!! FUNCTION
!! Datatype all the information relevant to create
!! bonds between atoms inside and outside the
!! cell
!!
!! NOTES
!!
!!
!! SOURCE

type go_bonds

!scalar
real(dp) :: tolerance ! To decide if consider bond the atom or not
                      ! 1.0 means that only consider values lower
                      ! than the sum of covalent radius 

integer  :: nbonds ! Total number of bonds for the system

!arrays

integer,pointer :: nbondi(:)    ! Number of bonds for atom i
integer,pointer :: indexi(:,:)  ! Indices of bonds for atom i
                                ! Positive: Vector from i to j   
                                ! Negative: Vector from j to i

real(dp),pointer :: bond_length(:) ! Bond lengths
real(dp),pointer :: bond_vect(:,:) ! Unitary vectors for bonds

end type go_bonds
!!***

!!****t* defs_mover/go_angles
!! NAME
!! go_angles
!!
!! FUNCTION
!! Datatype all the information relevant to create
!! angles between atoms inside and outside the
!! cell
!!
!! NOTES
!!
!!
!! SOURCE

type go_angles

!scalar

integer  :: nangles ! Total number of bonds for the system

!arrays

integer,pointer  :: angle_vertex(:)  ! Indices of the vertex atom
real(dp),pointer :: angle_value(:)   ! Value of angle in radians
real(dp),pointer :: angle_bonds(:,:) ! Indices of the bonds 
real(dp),pointer :: angle_vect(:,:)  ! Unitary vector perpendicular to the plane

end type go_angles
!!***

!----------------------------------------------------------------------

CONTAINS  !=============================================================
!!***

!!****f* defs_mover/ab_movetype_nullify
!! NAME
!! ab_movetype_nullify
!!
!! FUNCTION
!! Nullify all the pointers in a ab_mover
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ab_mover <type(ab_movetype)> = The ab_mover to nullify
!!
!! PARENTS
!!      defs_mover,mover
!!
!! CHILDREN
!!      xcast_mpi
!!
!! NOTES
!!  At present 32 variables are present in ab_movetype
!!  if a new variable is added in ab_movetype it should
!!  be added also for nullify here
!!
!! STATS
!!  integer         => 13
!!  real            =>  7
!!  integer array   =>  4
!!  real array      =>  6
!!  character array =>  1
!!  structures      =>  1
!!  TOTAL              32
!!
!! SOURCE

subroutine ab_movetype_nullify(ab_mover)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab_movetype_nullify'
!End of the abilint section

implicit none

!Arguments ------------------------------------

 type(ab_movetype),intent(inout) :: ab_mover

! ***************************************************************

! Delay of Permutation (Used by pred_langevin only)
 nullify(ab_mover%delayperm)
! Number of Geometry histories (Use by pred_diisrelax only)
 nullify(ab_mover%diismemory)
! Geometry Optimization precondition option
 nullify(ab_mover%goprecon)
! include a JELLium SLAB in the cell
 nullify(ab_mover%jellslab)
! Number of ATOMs
 nullify(ab_mover%natom)
! Number of CONstraint EQuations
 nullify(ab_mover%nconeq)
! Use by pred_isothermal only
 nullify(ab_mover%nnos)
! Number of SYMmetry operations
 nullify(ab_mover%nsym)
! Number of Types of atoms
 nullify(ab_mover%ntypat)
! OPTimize the CELL shape and dimensions
 nullify(ab_mover%optcell)
! RESTART Xcart and Fred
 nullify(ab_mover%restartxf)
! Sign of Permutation (Used by pred_langevin only)
 nullify(ab_mover%signperm)
! Ion movement
 nullify(ab_mover%ionmov)

! Use by pred_isothermal only
 nullify(ab_mover%bmass)
! Delta Time for IONs
 nullify(ab_mover%dtion)
! Used by pred_langevin only
 nullify(ab_mover%friction)
! Used by pred_langevin only
 nullify(ab_mover%mdwall)
! NOT DOCUMENTED
 nullify(ab_mover%noseinert)
! STRess PRECONditioner
 nullify(ab_mover%strprecon)
! VIScosity
 nullify(ab_mover%vis)

!arrays
! Indices of AToms that are FIXed
 nullify(ab_mover%iatfix)
! SYMmetry in REaL space
 nullify(ab_mover%symrel)
! TYPe of ATom
 nullify(ab_mover%typat)
! TYPe of ATom
 nullify(ab_mover%prtatlist)

! Mass of each atom (NOT IN DTSET)
 nullify(ab_mover%amass)!
! Molecular Dynamics Initial Temperature
 nullify(ab_mover%mdtemp)
! STRess TARGET
 nullify(ab_mover%strtarget)
! Use by pred_isothermal only
 nullify(ab_mover%qmass)
! Z number of each NUCLeus
 nullify(ab_mover%znucl)
! Geometry Optimization Preconditioner PaRaMeters
 nullify(ab_mover%goprecprm)

! Filename for Hessian matrix
 nullify(ab_mover%fnameabi_hes)
! Filename for _HIST file
 nullify(ab_mover%filnam_ds)

 call nullify_delocint(ab_mover%deloc)

end subroutine ab_movetype_nullify
!!***


!!****f* defs_mover/ab_movetype_end
!! NAME
!! ab_movetype_end
!!
!! FUNCTION
!! Deallocate all the pointers in a ab_mover
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ab_mover <type(ab_movetype)> = The ab_mover to destroy
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      xcast_mpi
!!
!! NOTES
!!  At present 32 variables are present in ab_movetype
!!  if a new variable is added in ab_movetype it should
!!  be added also for deallocate here
!!
!! STATS
!!  integer         => 13
!!  real            =>  7
!!  integer array   =>  4
!!  real array      =>  6
!!  character array =>  1
!!  structures      =>  1
!!  TOTAL              32
!!
!! SOURCE

subroutine ab_movetype_end(ab_mover)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab_movetype_end'
!End of the abilint section

implicit none

!Arguments ------------------------------------

 type(ab_movetype),intent(inout) :: ab_mover

! ***************************************************************

 call destroy_delocint(ab_mover%deloc)

! ab_mover is only pointer associated to other data, except possibly amass.
! TODO: check for amass in mover.F90
 call ab_movetype_nullify(ab_mover)

end subroutine ab_movetype_end
!!***

!!****f* defs_mover/ab_movetype_print
!! NAME
!! ab_movetype_print
!!
!! FUNCTION
!! Print all the variables in a ab_mover
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ab_mover <type(ab_movetype)> = The ab_mover to nullify
!!
!! PARENTS
!!
!! CHILDREN
!!      xcast_mpi
!!
!! NOTES
!!  At present 29 variables are present in ab_movetype
!!  if a new variable is added in ab_movetype it should
!!  be added also for print here
!!
!! SOURCE

subroutine ab_movetype_print(ab_mover,iout)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab_movetype_print'
!End of the abilint section

implicit none

!Arguments ------------------------------------
 type(ab_movetype),intent(inout) :: ab_mover
 integer,intent(in) :: iout

!Local variables-------------------------------
!arrays
character(len=1200) :: message
character(len=110)   :: fmt

! ***************************************************************

fmt='(a,e12.5,a,a,I5,a,a,I5,a,a,I5,a,a,I5,a,a,I5,a,a,I5,a,a,e12.5,a,a,e12.5,a,a,e12.5,a,a,e12.5,a,a,e12.5,a)'

write(message,fmt)&
& 'Delta Time for IONs',ab_mover%dtion,ch10, &
& 'include a JELLium SLAB in the cell',ab_mover%jellslab,ch10, &
& 'Number of ATOMs',ab_mover%natom,ch10, &
& 'Number of CONstraint EQuations',ab_mover%nconeq,ch10, &
& 'Number of SYMmetry operations',ab_mover%nsym,ch10, &
& 'OPTimize the CELL shape and dimensions',ab_mover%optcell,ch10, &
& 'RESTART Xcart and Fred',ab_mover%restartxf,ch10, &
& 'Molecular Dynamics Initial Temperature',ab_mover%mdtemp(1),ch10, &
& 'Molecular Dynamics Final Temperature',ab_mover%mdtemp(2),ch10, &
& 'NOT DOCUMENTED',ab_mover%noseinert,ch10, &
& 'STRess PRECONditioner',ab_mover%strprecon,ch10, &
& 'VIScosity',ab_mover%vis,ch10

! ! arrays
! ! Indices of AToms that are FIXed
! integer,  pointer :: iatfix(:,:)
! ! SYMmetry in REaL space
! integer,  pointer :: symrel(:,:,:)
! ! Mass of each atom (NOT IN DTSET)
! real(dp), pointer :: amass(:)
! ! STRess TARGET
! real(dp), pointer :: strtarget(:)
! Filename for Hessian matrix
! character(len=fnlen), pointer :: fnameabi_hes

 write(iout,*) 'CONTENTS of ab_mover'
 write(iout,'(a)') message

end subroutine ab_movetype_print
!!***

!!****f* defs_mover/ab_movehistory_nullify
!! NAME
!! ab_movehistory_nullify
!!
!! FUNCTION
!! Nullify all the pointers in a hist
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  hist <type(ab_movehistory)> = The hist to nullify
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      xcast_mpi
!!
!! NOTES
!!  At present 8 variables are present in ab_movehist
!!  if a new variable is added in ab_movehist it should
!!  be added also for nullify here
!!
!! SOURCE

subroutine ab_movehistory_nullify(hist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab_movehistory_nullify'
!End of the abilint section

implicit none

!Arguments ------------------------------------

 type(ab_movehistory),intent(inout) :: hist

! ***************************************************************

! Vector of (x,y,z)X(mxhist)
 nullify(hist%histA)
! Vector of (mxhist) values of energy
 nullify(hist%histE)
! Vector of (mxhist) values of ionic kinetic energy
 nullify(hist%histEk)
! Vector of (mxhist) values of time (relevant for MD calculations)
 nullify(hist%histT)
! Vector of (x,y,z)X(x,y,z)X(mxhist)
 nullify(hist%histR)
! Vector of (stress [6])X(mxhist)
 nullify(hist%histS)
! Vector of (x,y,z)X(natom)X(mxhist) values of velocity
 nullify(hist%histV)
! Vector of (x,y,z)X(natom)X(xcart,xred,fcart,fred)X(mxhist)
 nullify(hist%histXF)

end subroutine ab_movehistory_nullify
!!***

!----------------------------------------------------------------------

!!****f* defs_mover/ab_movehistory_end
!! NAME
!! ab_movehistory_end
!!
!! FUNCTION
!! Deallocate all the pointers in a hist
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  hist <type(ab_movehistory)> = The hist to deallocate
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      xcast_mpi
!!
!! NOTES
!!  At present 8 variables are present in ab_movehist
!!  if a new variable is added in ab_movehist it should
!!  be added also for deallocate here
!!
!! SOURCE

subroutine ab_movehistory_end(hist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab_movehistory_end'
!End of the abilint section

implicit none

!Arguments ------------------------------------

 type(ab_movehistory),intent(inout) :: hist

! ***************************************************************

! Vector of (x,y,z)X(mxhist)
 if (associated(hist%histA))  then
   ABI_DEALLOCATE(hist%histA)
 end if
! Vector of (mxhist) values of energy
 if (associated(hist%histE))  then
   ABI_DEALLOCATE(hist%histE)
 end if
! Vector of (mxhist) values of ionic kinetic energy
 if (associated(hist%histEk))  then
   ABI_DEALLOCATE(hist%histEk)
 end if
! Vector of (mxhist) values of time (relevant for MD calculations)
 if (associated(hist%histT))  then
   ABI_DEALLOCATE(hist%histT)
 end if
! Vector of (x,y,z)X(x,y,z)X(mxhist)
 if (associated(hist%histR))  then
   ABI_DEALLOCATE(hist%histR)
 end if
! Vector of (stress [6])X(mxhist)
 if (associated(hist%histS))  then
   ABI_DEALLOCATE(hist%histS)
 end if
! Vector of (x,y,z)X(natom)X(mxhist) values of velocity
 if (associated(hist%histV))  then
   ABI_DEALLOCATE(hist%histV)
 end if
! Vector of (x,y,z)X(natom)X(xcart,xred,fcart,fred)X(mxhist)
 if (associated(hist%histXF))  then
   ABI_DEALLOCATE(hist%histXF)
 end if

end subroutine ab_movehistory_end
!!***

!----------------------------------------------------------------------

!!****f* defs_mover/ab_movehistory_bcast
!! NAME
!! ab_movehistory_bcast
!!
!! FUNCTION
!! Broadcast a hist datastructure (from a root process to all others)
!!
!! INPUTS
!!  master=ID of the sending node in spaceComm
!!  spaceComm=MPI Communicator
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  hist <type(ab_movehistory)> = The hist to broadcast
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      xcast_mpi
!!
!! NOTES
!!  At present 8 variables are present in ab_movehist
!!  if a new variable is added in ab_movehist it should
!!  be added also for broadcast here
!!
!! SOURCE

subroutine ab_movehistory_bcast(hist,master,spaceComm)


 use defs_basis
 use defs_datatypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab_movehistory_bcast'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: master,spaceComm
 type(ab_movehistory),intent(inout) :: hist

!Local variables-------------------------------
!scalars
 integer :: bufsize,ierr,indx,nproc,rank
 integer :: sizeA,sizeA1,sizeA2
 integer :: sizeE,sizeEk,sizeT
 integer :: sizeR,sizeR1,sizeR2,sizeR3
 integer :: sizeS,sizeS1,sizeS2
 integer :: sizeV,sizeV1,sizeV2,sizeV3
 integer :: sizeXF,sizeXF1,sizeXF2,sizeXF3,sizeXF4
!arrays
 integer,allocatable :: buffer_i(:)
 real(dp),allocatable :: buffer_r(:)

! ***************************************************************

 ierr=0
 nproc=xcomm_size(spaceComm)
 if (nproc<=1) return

 rank=xcomm_rank(spaceComm)

!=== Broadcast integers and logicals
 ABI_ALLOCATE(buffer_i,(4))
 if (rank==master) then
   buffer_i(1)=hist%ihist
   buffer_i(2)=hist%mxhist
   buffer_i(3)=0;if (hist%isVused)  buffer_i(3)=1
   buffer_i(4)=0;if (hist%isARused) buffer_i(4)=1
 end if
 call xcast_mpi(buffer_i,master,spaceComm,ierr)
 if (rank/=master) then
   hist%ihist=buffer_i(1)
   hist%mxhist=buffer_i(2)
   hist%isVused =(buffer_i(3)==1)
   hist%isARused=(buffer_i(4)==1)
 end if
 ABI_DEALLOCATE(buffer_i)

!If history is empty, return
 if (hist%mxhist==0.or.hist%ihist==0) return

!=== Broadcast sizes of arrays
ABI_ALLOCATE(buffer_i,(17))
 if (rank==master) then
   sizeA1=size(hist%histA,1);sizeA2=size(hist%histA,2)
   sizeE=size(hist%histE,1);sizeEk=size(hist%histEk,1);sizeT=size(hist%histT,1)
   sizeR1=size(hist%histR,1);sizeR2=size(hist%histR,2);sizeR3=size(hist%histR,3)
   sizeS1=size(hist%histS,1);sizeS2=size(hist%histS,2)
   sizeV1=size(hist%histV,1);sizeV2=size(hist%histV,2);sizeV3=size(hist%histV,3)
   sizeXF1=size(hist%histXF,1);sizeXF2=size(hist%histXF,2);sizeXF3=size(hist%histXF,3)
   sizeXF4=size(hist%histXF,4)
   buffer_i(1)=sizeA1  ;buffer_i(2)=sizeA2
   buffer_i(3)=sizeE   ;buffer_i(4)=sizeEk
   buffer_i(5)=sizeT   ;buffer_i(6)=sizeR1
   buffer_i(7)=sizeR2  ;buffer_i(8)=sizeR3
   buffer_i(9)=sizeS1  ;buffer_i(10)=sizeS2
   buffer_i(11)=sizeV1 ;buffer_i(12)=sizeV2
   buffer_i(13)=sizeV3 ;buffer_i(14)=sizeXF1
   buffer_i(15)=sizeXF2;buffer_i(16)=sizeXF3
   buffer_i(17)=sizeXF4
 end if
 call xcast_mpi(buffer_i,master,spaceComm,ierr)
 if (rank/=master) then
   sizeA1 =buffer_i(1) ;sizeA2 =buffer_i(2)
   sizeE  =buffer_i(3) ;sizeEk =buffer_i(4)
   sizeT  =buffer_i(5) ;sizeR1 =buffer_i(6)
   sizeR2 =buffer_i(7) ;sizeR3 =buffer_i(8)
   sizeS1 =buffer_i(9) ;sizeS2 =buffer_i(10)
   sizeV1 =buffer_i(11);sizeV2 =buffer_i(12)
   sizeV3 =buffer_i(13);sizeXF1=buffer_i(14)
   sizeXF2=buffer_i(15);sizeXF3=buffer_i(16)
   sizeXF4=buffer_i(17)
 end if
 ABI_DEALLOCATE(buffer_i)

!=== Broadcast reals
 sizeA=sizeA1*sizeA2;sizeR=sizeR1*sizeR2*sizeR3;sizeS=sizeS1*sizeS2
 sizeV=sizeV1*sizeV2*sizeV3;sizeXF=sizeXF1*sizeXF2*sizeXF3*sizeXF4
 bufsize=sizeA+sizeE+sizeEk+sizeT+sizeR+sizeS+sizeV+sizeXF
 ABI_ALLOCATE(buffer_r,(bufsize))
 if (rank==master) then
   indx=0
   buffer_r(indx+1:indx+sizeA)=reshape(hist%histA(1:sizeA1,1:sizeA2),(/sizeA/))
   indx=indx+sizeA
   buffer_r(indx+1:indx+sizeE)=hist%histE(1:sizeE)
   indx=indx+sizeE
   buffer_r(indx+1:indx+sizeEk)=hist%histEk(1:sizeEk)
   indx=indx+sizeEk
   buffer_r(indx+1:indx+sizeT)=hist%histT(1:sizeT)
   indx=indx+sizeT
   buffer_r(indx+1:indx+sizeR)=reshape(hist%histR(1:sizeR1,1:sizeR2,1:sizeR3),(/sizeR/))
   indx=indx+sizeR
   buffer_r(indx+1:indx+sizeS)=reshape(hist%histS(1:sizeS1,1:sizeS2),(/sizeS/))
   indx=indx+sizeS
   buffer_r(indx+1:indx+sizeV)=reshape(hist%histV(1:sizeV1,1:sizeV2,1:sizeV3),(/sizeV/))
   indx=indx+sizeV
   buffer_r(indx+1:indx+sizeXF)=reshape(hist%histXF(1:sizeXF1,1:sizeXF2,1:sizeXF3,1:sizeXF4),(/sizeXF/))
 else
   if (associated(hist%histA))   then
     ABI_DEALLOCATE(hist%histA)
   end if
   if (associated(hist%histE))   then
     ABI_DEALLOCATE(hist%histE)
   end if
   if (associated(hist%histEk))  then
     ABI_DEALLOCATE(hist%histEk)
   end if
   if (associated(hist%histT))   then
     ABI_DEALLOCATE(hist%histT)
   end if
   if (associated(hist%histR))   then
     ABI_DEALLOCATE(hist%histR)
   end if
   if (associated(hist%histS))   then
     ABI_DEALLOCATE(hist%histS)
   end if
   if (associated(hist%histV))   then
     ABI_DEALLOCATE(hist%histV)
   end if
   if (associated(hist%histXF))  then
     ABI_DEALLOCATE(hist%histXF)
   end if
   ABI_ALLOCATE(hist%histA,(sizeA1,sizeA2))
   ABI_ALLOCATE(hist%histE,(sizeE))
   ABI_ALLOCATE(hist%histEk,(sizeEk))
   ABI_ALLOCATE(hist%histT,(sizeT))
   ABI_ALLOCATE(hist%histR,(sizeR1,sizeR2,sizeR3))
   ABI_ALLOCATE(hist%histS,(sizeS1,sizeS2))
   ABI_ALLOCATE(hist%histV,(sizeV1,sizeV2,sizeV3))
   ABI_ALLOCATE(hist%histXF,(sizeXF1,sizeXF2,sizeXF3,sizeXF4))
 end if
 call xcast_mpi(buffer_r,master,spaceComm,ierr)
 if (rank/=master) then
   indx=0
   hist%histA(1:sizeA1,1:sizeA2)=reshape(buffer_r(indx+1:indx+sizeA),(/sizeA1,sizeA2/))
   indx=indx+sizeA
   hist%histE(1:sizeE)=buffer_r(indx+1:indx+sizeE)
   indx=indx+sizeE
   hist%histEk(1:sizeEk)=buffer_r(indx+1:indx+sizeEk)
   indx=indx+sizeEk
   hist%histT(1:sizeT)=buffer_r(indx+1:indx+sizeT)
   indx=indx+sizeT
   hist%histR(1:sizeR1,1:sizeR2,1:sizeR3)=reshape(buffer_r(indx+1:indx+sizeR),(/sizeR1,sizeR2,sizeR3/))
   indx=indx+sizeR
   hist%histS(1:sizeS1,1:sizeS2)=reshape(buffer_r(indx+1:indx+sizeS),(/sizeS1,sizeS2/))
   indx=indx+sizeS
   hist%histV(1:sizeV1,1:sizeV2,1:sizeV3)=reshape(buffer_r(indx+1:indx+sizeV),(/sizeV1,sizeV2,sizeV3/))
   indx=indx+sizeV
   hist%histXF(1:sizeXF1,1:sizeXF2,1:sizeXF3,1:sizeXF4)=reshape(buffer_r(indx+1:indx+sizeXF), &
&                                                       (/sizeXF1,sizeXF2,sizeXF3,sizeXF4/))
 end if
 ABI_DEALLOCATE(buffer_r)

end subroutine ab_movehistory_bcast

!----------------------------------------------------------------------

end module defs_mover
!!***
