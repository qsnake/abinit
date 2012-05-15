/** - C source file -
 *
 * Copyright (C) 2009-2012 ABINIT Group (Damien Caliste)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#include "ab6_invars.h"
#include <string.h>
#include <stdlib.h>
#include <config.h>

#define VARS_NAME(a,A)     ABI_FC_MOD(m_ab6_invars,M_AB6_INVARS,	\
                                      ab6_invars_ ## a,AB6_INVARS_ ## A)
#define VARS_CALL(a,A,...) VARS_NAME(a,A)(__VA_ARGS__)

/* Fortran interface. */
void VARS_NAME(new_from_file,NEW_FROM_FILE)
     (int *dt, const char *filename, int *len, const char *pseudo, int *npsp);
void VARS_NAME(new_from_string,NEW_FROM_STRING) 
  (int *dt, const char *string, int *len);
void VARS_NAME(free,FREE)            
  (int *dt);
void VARS_NAME(get_ndtset,GET_NDTSET)   
  (int *dt, int *ndtset, unsigned int *ab6_errno);
void VARS_NAME(get_integer,GET_INTEGER) 
  (int *dt, int *value, int *att, int *idtset, unsigned int *ab6_errno);
void VARS_NAME(get_real,GET_REAL)
  (int *dt, double *value, int *att, int *idtset, unsigned int *ab6_errno);
void VARS_NAME(get_shape,GET_SHAPE) 
  (int *dt, int *dims, int *size, int *att, int *idtset, unsigned int *ab6_errno);
void VARS_NAME(get_integer_array,GET_INTEGER_ARRAY)
  (int *dt, int *values, int *size, int *att, int *idtset, unsigned int *ab6_errno);
void VARS_NAME(get_real_array,GET_REAL_ARRAY)
  (int *dt, double *values, int *size, int *att, int *idtset, unsigned int *ab6_errno);

#include "ab6_invars_c.h"

Ab6InvarsTypes ab6_invars_get_type_from_id(Ab6InvarsIds id)
{
  return ((id >= 0 && id < AB6_INVARS_N_IDS)?ab6_invars_types[id]:_OTHER);
}

Ab6Invars* ab6_invars_new_from_file(const char *filename)
{
  int n, npsp;
  Ab6Invars *dt, tmpDt;
  char buf[8] = "        ";
  
  n = strlen(filename);
  npsp = 0;

  VARS_CALL(new_from_file, NEW_FROM_FILE, &tmpDt, filename, &n, buf, &npsp);
  if (tmpDt > 0)
    {
      dt = malloc(sizeof(Ab6Invars));
      *dt = tmpDt;
    }
  else
    dt = (Ab6Invars*)0;
  return dt;
}
Ab6Invars* ab6_invars_new_from_file_with_pseudo(const char *filename, const char **pspfiles)
{
  int n, npsp;
  Ab6Invars *dt, tmpDt;
  char *psp;
  size_t i, len;
  
  n = strlen(filename);
  npsp = 0;
  if (pspfiles)
    for (npsp = 0; pspfiles[npsp]; npsp++);

  if (npsp)
    {
      len = sizeof(char) * npsp * 264;
      psp = malloc(len);
      for (npsp = 0; pspfiles[npsp]; npsp++)
        strncpy(psp + 264 * npsp, pspfiles[npsp], 264);
      for (i = 0; i < len; i++)
        if (psp[i] == '\0')
          psp[i] = ' ';
    }
  
  VARS_CALL(new_from_file, NEW_FROM_FILE, &tmpDt, filename, &n, psp, &npsp);
  if (tmpDt > 0)
    {
      dt = malloc(sizeof(Ab6Invars));
      *dt = tmpDt;
    }
  else
    dt = (Ab6Invars*)0;
  return dt;
}
Ab6Invars* ab6_invars_new_from_string(const char *string)
{
  int n;
  Ab6Invars *dt, tmpDt;
  
  n = strlen(string);
  
  VARS_CALL(new_from_string, NEW_FROM_STRING, &tmpDt, string, &n);
  if (tmpDt > 0)
    {
      dt = malloc(sizeof(Ab6Invars));
      *dt = tmpDt;
    }
  else
    dt = (Ab6Invars*)0;
  return dt;
}

void ab6_invars_free(Ab6Invars *dt)
{
  VARS_CALL(free, FREE, dt);
  free(dt);
}

Ab6Error ab6_invars_get_ndtset(Ab6Invars *dt, int *ndtset)
{
  Ab6Error ab6_errno;

  VARS_CALL(get_ndtset, GET_NDTSET, dt, ndtset, &ab6_errno);
  return ab6_errno;
}

Ab6Error ab6_invars_get_integer(Ab6Invars *dt, Ab6InvarsIds id, int idtset, int *value)
{
  Ab6Error ab6_errno;

  if (AB6_INVARS_TYPE(id) != _INT_SCALAR) return AB6_ERROR_INVARS_ATT;
  VARS_CALL(get_integer, GET_INTEGER, dt, value, (int*)&id, &idtset, &ab6_errno);
  return ab6_errno;
}

Ab6Error ab6_invars_get_real(Ab6Invars *dt, Ab6InvarsIds id, int idtset, double *value)
{
  Ab6Error ab6_errno;

  if (AB6_INVARS_TYPE(id) != _DOUBLE_SCALAR) return AB6_ERROR_INVARS_ATT;
  VARS_CALL(get_real, GET_REAL, dt, value, (int*)&id, &idtset, &ab6_errno);
  return ab6_errno;
}

Ab6Error ab6_invars_get_shape(Ab6Invars *dt, int *n, int dims[7], Ab6InvarsIds id, int idtset)
{
  Ab6Error ab6_errno;

  VARS_CALL(get_shape, GET_SHAPE, dt, dims, n, (int*)&id, &idtset, &ab6_errno);
  return ab6_errno;
}

Ab6Error ab6_invars_get_integer_array(Ab6Invars *dt, int *values, size_t n, Ab6InvarsIds id, int idtset)
{
  Ab6Error ab6_errno;

  if (AB6_INVARS_TYPE(id) != _INT_ARRAY) return AB6_ERROR_INVARS_ATT;
  VARS_CALL(get_integer_array, GET_INTEGER_ARRAY,
	    dt, values, (int*)&n, (int*)&id, &idtset, &ab6_errno);
  return ab6_errno;
}

Ab6Error ab6_invars_get_real_array(Ab6Invars *dt, double *values, size_t n, Ab6InvarsIds id, int idtset)
{
  Ab6Error ab6_errno;

  if (AB6_INVARS_TYPE(id) != _DOUBLE_ARRAY) return AB6_ERROR_INVARS_ATT;
  VARS_CALL(get_real_array, GET_REAL_ARRAY,
	    dt, values, (int*)&n, (int*)&id, &idtset, &ab6_errno);
  return ab6_errno;
}
