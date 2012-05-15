!{\src2tex{textfont=tt}}
!!****f* ABINIT/orthonormalize
!! NAME
!! orthonormalize
!!
!! FUNCTION
!! This routine computes the overlap of two wavefunctions (for a given number of bands)
!! and orthonormalizes it:
!!      - Computes the products of two rectangular matrices
!!         containing the wavefunctions psi and S.psi (where S is the
!!         overlap (with the PAW terms if necessary)).
!!      - Does a Cholesky decomposition of this overlap
!!      - rotates the initial matrix blockvectorx by the triangular matrix to
!!         have an orthonormal set of wavefunctions
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (GZ,AR,MT)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blockvectorbx = matrix of dimension (blocksize,vectsize)
!!                  (e.g. block of overlap*wavefunction)
!!  blocksize     = dimension of matrices (e.g number of bands)
!!  mpi_enreg     = informations about MPI parallelization
!!  vectsize      = dimension of matrices (e.g number of G vector)
!!
!! OUTPUT
!!  sqgram        = Choleski decomposition of transpose(blockvector)*blockvectorx
!!
!! SIDE EFFECTS
!!  blockvectorx  = on input, matrix of dimension (vectsize,blocksize)
!!                  (e.g block of wavefunction)
!!  blockvectorx  = on output, orthonormalized wavefunction.
!!
!! PARENTS
!!      lobpcgIIwf,m_lobpcg,m_lobpcgIIIwf,pw_orthon
!!
!! CHILDREN
!!      dgemm,dpotrf,dtrsm,wrtout,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine orthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,sqgram,vectsize)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'orthonormalize'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,vectsize
 type(mpi_type) :: mpi_enreg
!arrays
 real(dp) :: blockvectorbx(vectsize,blocksize),blockvectorx(vectsize,blocksize)
 real(dp) :: sqgram(blocksize,blocksize)

!Local variables-------------------------------
!scalars
 integer :: ierr,info,old_paral_level,spaceComm
 character(len=500) :: message
!arrays

! *********************************************************************

 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorx,&
& vectsize,blockvectorbx,vectsize,zero,sqgram,blocksize)
 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_3d)
 call xsum_mpi(sqgram,spaceComm,ierr)
 mpi_enreg%paral_level= old_paral_level

!Cholesky factorization of sqgram (ouside upper Triangular of sqgram)
 call dpotrf('u',blocksize,sqgram,blocksize,info)

 if (info /= 0 )  then
   write(message,'(a,i3)') 'WARNING in dpotrf, info=',info
   call wrtout(std_out,message,'COLL')
 end if

!Find X  X*sqgram=blockvectorx
 call dtrsm('r','u','n','n',vectsize,blocksize,one,sqgram,blocksize,&
& blockvectorx,vectsize)

end subroutine orthonormalize
!!***
