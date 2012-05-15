!{\src2tex{textfont=tt}}
!!****f* ABINIT/initwf
!!
!! NAME
!! initwf
!!
!! FUNCTION
!! Initialization of wavefunctions, fform=2 .
!! If formeig==1, and partially filled case, I am not sure that the eig_k
!! are initialized properly ...
!! formeig option (format of the eigenvalues and eigenvector) :
!!   0 => ground-state format (initialisation of
!!        eigenvectors with random numbers, vector of eigenvalues)
!!   1 => respfn format (initialisation of
!!        eigenvectors with 0 s, hermitian matrix of eigenvalues)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! formeig=see above
!! headform=header format (might be needed to read the block of wfs)
!! icg=shift to be given to the location of the data in the array cg
!! ikpt= number of the k point of which the wf is initialised
!! isppol=spin index
!! mcg=dimension of the cg array
!! mpi_enreg=informations about MPI parallelization
!! nband_k=number of bands at this particular k point
!! nkpt=number of k points
!! npw=number of plane waves
!! nspinor=number of spinorial components of the wavefunctions (on current proc)
!! wff1=structure info for file containing wavefunctions (when needed)
!!
!! OUTPUT
!! cg(2,mcg)=complex wf array
!! if ground state format (formeig=0):
!!   eig_k(nband_k)=list of eigenvalues (input or init to large number), hartree
!! if respfn format (formeig=1):
!!   eig_k(2*nband_k*nband_k)=
!!             matrix of eigenvalues (input or init to large number), hartree
!!
!! SIDE EFFECTS
!! Input/output:
!! occ_k(nband_k)=list of occupations (input or left to their initial value)
!! ikptsp_old=number of the previous spin-k point, or 0 if first call of present file
!!
!! PARENTS
!!      wfsinp
!!
!! CHILDREN
!!      rwwf,timab,wffreadskipk,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initwf(cg,eig_k,formeig,headform,icg,ikpt,ikptsp_old,&
&  isppol,mcg,mpi_enreg,&
&  nband_k,nkpt,npw,nspinor,occ_k,wff1)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initwf'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot, except_this_one => initwf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!  The size of eigen depends on formeig (ground state or respfn)
!scalars
 integer,intent(in) :: formeig,headform,icg,ikpt,isppol,mcg,nband_k,nkpt,npw
 integer,intent(in) :: nspinor
 integer,intent(inout) :: ikptsp_old
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wffile_type),intent(inout) :: wff1
!arrays
 real(dp),intent(inout) :: occ_k(nband_k)
 real(dp),intent(out) :: cg(2,mcg),eig_k((2*nband_k)**formeig*nband_k)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt_max=50
 integer :: ikpt0,nband_disk,tim_rwwf
 character(len=500) :: message
!arrays
 integer,allocatable :: kg_dum(:,:)
 real(dp) :: tsec(2)

! *************************************************************************

!DEBUG
!write(std_out,*)' initwf : enter, ikptsp_old,ikpt,isppol,nkpt= ',&
!& ikptsp_old,ikpt,isppol,nkpt
!stop
!ENDDEBUG

 call timab(770,1,tsec)
 call timab(771,1,tsec)

 ABI_ALLOCATE(kg_dum,(3,0))

!Skip wavefunctions for k-points not treated by this proc.
!(from ikptsp_old+1 to ikpt+(isppol-1)*nkpt-1)
 if(ikptsp_old<ikpt+(isppol-1)*nkpt-1)then

!  DEBUG
!  write(std_out,*)' initwf : skip some k point'
!  ENDDEBUG

   do ikpt0=ikptsp_old+1,ikpt+(isppol-1)*nkpt-1
     call WffReadSkipK(formeig,headform,ikpt0,isppol,mpi_enreg,wff1)
   end do
 end if

!DEBUG
!write(std_out,*)' initwf : before rwwf'
!write(std_out,*)' formeig,icg,ikpt,isppol=',formeig,icg,ikpt,isppol
!write(std_out,*)' nband_k,nband_disk,npw,nspinor=',nband_k,nband_disk,npw,nspinor
!write(std_out,*)' unwff1=',unwff1
!stop
!ENDDEBUG

 call timab(771,2,tsec)

 if(mpi_enreg%paralbd==0)tim_rwwf=2
 if(mpi_enreg%paralbd==1)tim_rwwf=20

 call rwwf(cg,eig_k,formeig,headform,icg,ikpt,isppol,kg_dum,nband_k,mcg,mpi_enreg,nband_k,nband_disk,&
& npw,nspinor,occ_k,1,0,tim_rwwf,wff1)

 call timab(772,1,tsec)

!DEBUG
!write(std_out,*)' initwf : after rwwf'
!stop
!ENDDEBUG

 if(ikpt<=nkpt_max)then
   write(message, '(a,i6,a,i6,a,i5)' ) &
&   ' initwf : disk file gives npw=',npw,&
&   ' nband=',nband_disk,' for kpt number=',ikpt
   call wrtout(std_out,  message,'PERS')
 else if(ikpt==nkpt_max+1)then
   write(message, '(a)' )&
&   ' initwf : the number of similar message is sufficient... stop printing them'
   call wrtout(std_out,message,'PERS')
 end if

!Check the number of bands on disk file against desired number
!(These are not required to agree)
 if (nband_disk/=nband_k) then
   write(message, '(a,a,a,a,i4,a,i6,a,a,a,i6,a,a,a,a,a)' ) ch10,&
&   ' initwf: COMMENT -',ch10,&
&   '  For kpt number',ikpt,' disk file has',nband_disk,' bands',ch10,&
&   '  but input file gave nband=',nband_k,'.',ch10,&
&   '  This is not fatal.',ch10, &
&   '  Bands are skipped or filled with random numbers.'
   call wrtout(std_out,message,'PERS')
 end if

 if(ikpt<=nkpt_max)then
   write(message, '(a,i6,a)' ) &
&   ' initwf :',nband_disk,' bands have been initialized from disk'
   call wrtout(std_out,  message,'PERS')
 end if

 ikptsp_old=ikpt+(isppol-1)*nkpt

 ABI_DEALLOCATE(kg_dum)

 call timab(772,2,tsec)
 call timab(770,2,tsec)

!DEBUG
!write(std_out,*)' initwf : exit '
!stop
!ENDDEBUG

end subroutine initwf
!!***
