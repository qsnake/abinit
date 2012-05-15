!{\src2tex{textfont=tt}}
!!****f* ABINIT/zorthonormalize
!! NAME
!! zorthonormalize
!!
!! FUNCTION
!! This routine computes the overlap of two complex wavefunctions (for a given number of bands)
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
!!
!! PARENTS
!!      lobpcgccIIIwf,lobpcgccIIwf,m_lobpcg,pw_orthon
!!
!! CHILDREN
!!      wrtout,xcomm_init,xsum_mpi,zgemm,zpotrf,ztrsm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine zorthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,sqgram,vectsize)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'zorthonormalize'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,vectsize
 type(mpi_type) :: mpi_enreg
!arrays
 complex(dpc) :: blockvectorbx(vectsize,blocksize)
 complex(dpc) :: blockvectorx(vectsize,blocksize),sqgram(blocksize,blocksize)

!Local variables-------------------------------
!scalars
 integer :: ierr,info,old_paral_level,spaceComm
 character(len=500) :: message
!arrays

! *********************************************************************

 call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
& vectsize,blockvectorbx,vectsize,czero,sqgram,blocksize)

 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_3d)
 call xsum_mpi(sqgram,spaceComm,ierr)
 mpi_enreg%paral_level= old_paral_level

 call zpotrf('u',blocksize,sqgram,blocksize,info)

 if (info /= 0 )  then
   write(message,'(a,i3)') 'WARNING in zpotrf, info=',info
   call wrtout(std_out,message,'COLL')
 end if

 call ztrsm('r','u','n','n',vectsize,blocksize,cone,sqgram,blocksize,blockvectorx,vectsize)


end subroutine zorthonormalize
!!***
