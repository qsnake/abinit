/** -*- C header file -*-
 *
 * Copyright (C) 2009-2012 ABINIT Group (Damien Caliste)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef AB6_BASE
#define AB6_BASE

/**
 * Ab6Error:
 * @AB6_NO_ERROR: no error.
 * @AB6_ERROR_INVARS_OBJ: wrong dataset object.
 * @AB6_ERROR_INVARS_ATT: wrong attribute in dataset.
 * @AB6_ERROR_INVARS_ID: wrong dataset index.
 * @AB6_ERROR_INVARS_SIZE: wrong size when accessing arrays.
 *
 * An error code.
 */
typedef enum
  {
    AB6_NO_ERROR,
    AB6_ERROR_OBJ,
    AB6_ERROR_ARG,
    AB6_ERROR_INVARS_ATT,
    AB6_ERROR_INVARS_ID,
    AB6_ERROR_INVARS_SIZE,
    AB6_ERROR_SYM_NOT_PRIMITIVE,
    AB6_ERROR_SYM_BRAVAIS_XRED
  } Ab6Error;

char* ab6_error_string_from_id(Ab6Error errno);

#ifndef GLIB_MAJOR_VERSION
typedef int gboolean;
#define g_malloc(A) malloc(A)
#define g_free(A) free(A)
#endif

#endif
