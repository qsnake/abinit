!{\src2tex{textfont=tt}}
!!****f* ABINIT/set_twind
!! NAME
!! set_twind
!!
!! FUNCTION
!! set index tables for mappings used in magnetization
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtefield <type(efield_type)> = variables related to Berry phase
!! 
!! NOTES
!!
!! PARENTS
!!      initberry
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine set_twind(dtefield)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'set_twind'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 type(efield_type), intent(inout) :: dtefield

!Local

! *********************************************************************
!fill in indexing information
!for some quantities only six combinations of bra and ket shift vectors are stored
!(1) -k1 -k2
!(2) +k1 -k2
!(3) -k2 -k3
!(4) +k2 -k3
!(5) -k3 -k1
!(6) +k3 -k1

 dtefield%twind(:,:) = 0
 dtefield%indhk(:,:) = 0

#define PIND(vdir,vsig) 2*vdir-1+(vsig+1)/2

 dtefield%twind(PIND( 1,-1),PIND( 2, 1)) =  1  ! -k1 - k2
 dtefield%twind(PIND( 1, 1),PIND( 2, 1)) =  2  ! +k1 - k2
 dtefield%twind(PIND( 1,-1),PIND( 2,-1)) = -2  ! -k1 + k2
 dtefield%twind(PIND( 1, 1),PIND( 2,-1)) = -1  ! +k1 + k2
 dtefield%twind(PIND( 1,-1),PIND( 3, 1)) =  5  ! -k1 - k3
 dtefield%twind(PIND( 1, 1),PIND( 3, 1)) = -6  ! +k1 - k3
 dtefield%twind(PIND( 1,-1),PIND( 3,-1)) =  6  ! -k1 + k3
 dtefield%twind(PIND( 1, 1),PIND( 3,-1)) = -5  ! +k1 + k3
 
 dtefield%twind(PIND( 2,-1),PIND( 1, 1)) =  1  ! -k2 - k1
 dtefield%twind(PIND( 2, 1),PIND( 1, 1)) = -2  ! +k2 - k1
 dtefield%twind(PIND( 2,-1),PIND( 1,-1)) =  2  ! -k2 + k1
 dtefield%twind(PIND( 2, 1),PIND( 1,-1)) = -1  ! +k2 + k1
 dtefield%twind(PIND( 2,-1),PIND( 3, 1)) =  3  ! -k2 - k3
 dtefield%twind(PIND( 2, 1),PIND( 3, 1)) =  4  ! +k2 - k3
 dtefield%twind(PIND( 2,-1),PIND( 3,-1)) = -4  ! -k2 + k3
 dtefield%twind(PIND( 2, 1),PIND( 3,-1)) = -3  ! +k2 + k3

 dtefield%twind(PIND( 3,-1),PIND( 1, 1)) =  5  ! -k3 - k1
 dtefield%twind(PIND( 3, 1),PIND( 1, 1)) =  6  ! +k3 - k1
 dtefield%twind(PIND( 3,-1),PIND( 1,-1)) = -6  ! -k3 + k1
 dtefield%twind(PIND( 3, 1),PIND( 1,-1)) = -5  ! +k3 + k1
 dtefield%twind(PIND( 3,-1),PIND( 2, 1)) =  3  ! -k3 - k2
 dtefield%twind(PIND( 3, 1),PIND( 2, 1)) = -4  ! +k3 - k2
 dtefield%twind(PIND( 3,-1),PIND( 2,-1)) =  4  ! -k3 + k2
 dtefield%twind(PIND( 3, 1),PIND( 2,-1)) = -3  ! +k3 + k2

 dtefield%indhk(PIND( 1,-1),PIND( 2, 1)) =  1  ! -k1 - k2
 dtefield%indhk(PIND( 1, 1),PIND( 2, 1)) =  2  ! +k1 - k2
 dtefield%indhk(PIND( 1,-1),PIND( 2,-1)) =  3  ! -k1 + k2
 dtefield%indhk(PIND( 1, 1),PIND( 2,-1)) =  4  ! +k1 + k2
 dtefield%indhk(PIND( 1,-1),PIND( 3, 1)) =  5  ! -k1 - k3
 dtefield%indhk(PIND( 1, 1),PIND( 3, 1)) =  6  ! +k1 - k3
 dtefield%indhk(PIND( 1,-1),PIND( 3,-1)) =  7  ! -k1 + k3
 dtefield%indhk(PIND( 1, 1),PIND( 3,-1)) =  8  ! +k1 + k3
 
 dtefield%indhk(PIND( 2,-1),PIND( 1, 1)) =  9  ! -k2 - k1
 dtefield%indhk(PIND( 2, 1),PIND( 1, 1)) = 10  ! +k2 - k1
 dtefield%indhk(PIND( 2,-1),PIND( 1,-1)) = 11  ! -k2 + k1
 dtefield%indhk(PIND( 2, 1),PIND( 1,-1)) = 12  ! +k2 + k1
 dtefield%indhk(PIND( 2,-1),PIND( 3, 1)) = 13  ! -k2 - k3
 dtefield%indhk(PIND( 2, 1),PIND( 3, 1)) = 14  ! +k2 - k3
 dtefield%indhk(PIND( 2,-1),PIND( 3,-1)) = 15  ! -k2 + k3
 dtefield%indhk(PIND( 2, 1),PIND( 3,-1)) = 16  ! +k2 + k3

 dtefield%indhk(PIND( 3,-1),PIND( 1, 1)) = 17  ! -k3 - k1
 dtefield%indhk(PIND( 3, 1),PIND( 1, 1)) = 18  ! +k3 - k1
 dtefield%indhk(PIND( 3,-1),PIND( 1,-1)) = 19  ! -k3 + k1
 dtefield%indhk(PIND( 3, 1),PIND( 1,-1)) = 20  ! +k3 + k1
 dtefield%indhk(PIND( 3,-1),PIND( 2, 1)) = 21  ! -k3 - k2
 dtefield%indhk(PIND( 3, 1),PIND( 2, 1)) = 22  ! +k3 - k2
 dtefield%indhk(PIND( 3,-1),PIND( 2,-1)) = 23  ! -k3 + k2
 dtefield%indhk(PIND( 3, 1),PIND( 2,-1)) = 24  ! +k3 + k2

end subroutine set_twind

!!***
