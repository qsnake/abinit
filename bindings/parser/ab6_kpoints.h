/** C header file
 *
 * Copyright (C) 2008-2012 ABINIT Group (Damien Caliste)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef AB6_KPOINTS
#define AB6_KPOINTS

#include "ab6_base.h"
#include "ab6_symmetry.h"

Ab6Error ab6_kpoints_get_irreductible_zone(Ab6Symmetry *sym,
					   int *irrzon, double *phnons,
					   int n1, int n2, int n3,
					   int nsppol, int nspden);
Ab6Error ab6_kpoints_get_mp_k_grid   (Ab6Symmetry *sym, int *nkpt, double **kpt,
				      double **wkpt, const int ngkpt[3],
				      const int nshiftk, const double *shiftk);
Ab6Error ab6_kpoints_get_auto_k_grid (Ab6Symmetry *sym, int *nkpt, double **kpt,
				      double **wkpt, const double kptrlen);
#endif
