/** -*- C source file -*-
 *
 * Copyright (C) 2009-2012 ABINIT Group (Damien Caliste)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#include "ab6_base.h"

/*
  AB6_NO_ERROR,
  AB6_ERROR_OBJ,
  AB6_ERROR_ARG,
  AB6_ERROR_INVARS_ATT,
  AB6_ERROR_INVARS_ID,
  AB6_ERROR_INVARS_SIZE,
  AB6_ERROR_SYM_NOT_PRIMITIVE,
  AB6_ERROR_SYM_BRAVAIS_XRED
*/
static char* messages[8] = {
  "No error.",
  "Wrong pointer to object.",
  "Wrong value for one argument of the calling routine.",
  "Unknown attribute or wrong attribute type from Dtset structure.",
  "Out of bounds dtset number.",
  "Wrong input size for array.",
  "The cell is not primitive.",
  "The bravais lattice has more symmetries than the lattice"
  " system found from the atom coordinates."};

char* ab6_error_string_from_id(Ab6Error errno)
{
  if (errno < 0 || errno >= 8)
    return "Unknown error id.";
  else
    return messages[errno];
}
