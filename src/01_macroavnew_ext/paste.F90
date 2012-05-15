! Copyright (C) 1999-2003 (P. Ordejon, J. Junquera)
! Modified by TD (2008)
! This file is distributed under the terms of the
! GNU General Public License, see ~abinit/COPYING
! or http://www.gnu.org/copyleft/gpl.txt .

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

      CHARACTER(LEN=26) FUNCTION PASTE( STR1, STR2 )


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'PASTE'
!End of the abilint section

      IMPLICIT NONE
 
! CONCATENATES THE STRINGS STR1 AND STR2 REMOVING BLANKS IN BETWEEN

      CHARACTER(LEN=*) :: STR1, STR2
      INTEGER :: L
      DO L = LEN( STR1 ), 1, -1
         IF (STR1(L:L) /= ' ') EXIT
      END DO
      PASTE = STR1(1:L)//STR2
      END FUNCTION PASTE


      CHARACTER(LEN=26) FUNCTION PASTEB( STR1, STR2 )
 
! CONCATENATES THE STRINGS STR1 AND STR2 LEAVING ONLY ONE BLANK IN BETWEEN

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'PASTEB'
 use interfaces_01_macroavnew_ext, except_this_one => PASTEB
!End of the abilint section

      IMPLICIT NONE

      CHARACTER(LEN=*) :: STR1, STR2 
      CHARACTER(LEN=1), PARAMETER :: BLANK=' '
      INTEGER :: L
      DO L = LEN( STR1 ), 1, -1
         IF (STR1(L:L) /= ' ') EXIT
      END DO
      PASTEB = STR1(1:L)//BLANK
      PASTEB = PASTEB(1:L+1)//STR2
      END FUNCTION PASTEB

