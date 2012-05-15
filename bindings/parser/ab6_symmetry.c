/** - C source file -
 *
 * Copyright (C) 2009-2012 ABINIT Group (Damien Caliste)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#include <stdlib.h>

#include <config.h>
#include <ab6_symmetry.fortran.h>
#include "ab6_symmetry.h"

Ab6Symmetry* ab6_symmetry_new()
{
  Ab6Symmetry *sym, tmpSym;
  
  SYM_CALL(new,NEW, &tmpSym);
  if (tmpSym > 0)
    {
      sym = g_malloc(sizeof(Ab6Symmetry));
      *sym = tmpSym;
    }
  else
    sym = (Ab6Symmetry*)0;
  return sym;
}
void ab6_symmetry_free(Ab6Symmetry *sym)
{
  SYM_CALL(free,FREE, sym);
  g_free(sym);
}
Ab6Error ab6_symmetry_set_tolerance(Ab6Symmetry *sym, double tolsym)
{
  Ab6Error ab6_errno;

  SYM_CALL(set_tolerance,SET_TOLERANCE, sym, &tolsym, &ab6_errno);
  return ab6_errno;
}
Ab6Error ab6_symmetry_set_lattice(Ab6Symmetry *sym, double rprimd[3][3])
{
  Ab6Error ab6_errno;
  double rprimd_[9];

  rprimd_[0] = rprimd[0][0];
  rprimd_[1] = rprimd[1][0];
  rprimd_[2] = rprimd[2][0];
  rprimd_[3] = rprimd[0][1];
  rprimd_[4] = rprimd[1][1];
  rprimd_[5] = rprimd[2][1];
  rprimd_[6] = rprimd[0][2];
  rprimd_[7] = rprimd[1][2];
  rprimd_[8] = rprimd[2][2];
  SYM_CALL(set_lattice,SET_LATTICE, sym, rprimd_, &ab6_errno);
  return ab6_errno;
}
Ab6Error ab6_symmetry_set_structure(Ab6Symmetry *sym, int natoms,
				    int *typeAt, double *xRed)
{
  Ab6Error ab6_errno;

  SYM_CALL(set_structure,SET_STRUCTURE, sym, &natoms, typeAt, xRed, &ab6_errno);
  return ab6_errno;
}
Ab6Error ab6_symmetry_set_spin(Ab6Symmetry *sym, double *spinAt)
{
  Ab6Error ab6_errno;
  int nAtoms;

  SYM_CALL(get_n_atoms,GET_N_ATOMS, sym, &nAtoms, &ab6_errno);
  if (ab6_errno != AB6_NO_ERROR)
    return ab6_errno;

  SYM_CALL(set_spin,SET_SPIN, sym, &nAtoms, spinAt, &ab6_errno);
  return ab6_errno;
}
Ab6Error ab6_symmetry_set_collinear_spin(Ab6Symmetry *sym, int *spinAt)
{
  Ab6Error ab6_errno;
  int nAtoms;

  SYM_CALL(get_n_atoms,GET_N_ATOMS, sym, &nAtoms, &ab6_errno);
  if (ab6_errno != AB6_NO_ERROR)
    return ab6_errno;

  SYM_CALL(set_collinear_spin,SET_COLLINEAR_SPIN, sym, &nAtoms, spinAt, &ab6_errno);
  return ab6_errno;
}
Ab6Error ab6_symmetry_set_spin_orbit(Ab6Symmetry *sym, gboolean status)
{
  Ab6Error ab6_errno;

  SYM_CALL(set_spin_orbit,SET_SPIN_ORBIT, sym, &status, &ab6_errno);
  return ab6_errno;
}
Ab6Error ab6_symmetry_set_field(Ab6Symmetry *sym, double field[3])
{
  Ab6Error ab6_errno;

  SYM_CALL(set_field,SET_FIELD, sym, field, &ab6_errno);
  return ab6_errno;
}
Ab6Error ab6_symmetry_set_jellium(Ab6Symmetry *sym, gboolean jellium)
{
  Ab6Error ab6_errno;

  SYM_CALL(set_jellium,SET_JELLIUM, sym, &jellium, &ab6_errno);
  return ab6_errno;
}
Ab6Error ab6_symmetry_set_periodicity(Ab6Symmetry *sym, gboolean periodic[3])
{
  Ab6Error ab6_errno;

  SYM_CALL(set_periodicity,SET_PERIODICITY, sym, periodic, &ab6_errno);
  return ab6_errno;
}

Ab6Error ab6_symmetry_get_n_atoms(Ab6Symmetry *sym, int *nAtoms)
{
  Ab6Error ab6_errno;

  SYM_CALL(get_n_atoms,GET_N_ATOMS, sym, nAtoms, &ab6_errno);
  return ab6_errno;
}
Ab6Error ab6_symmetry_get_n_sym(Ab6Symmetry *sym, int *nSym)
{
  Ab6Error ab6_errno;

  SYM_CALL(get_n_sym,GET_N_SYM, sym, nSym, &ab6_errno);
  return ab6_errno;
}
Ab6Error ab6_symmetry_get_multiplicity(Ab6Symmetry *sym, int *multiplicity)
{
  Ab6Error ab6_errno;

  SYM_CALL(get_multiplicity,GET_MULTIPLICITY, sym, multiplicity, &ab6_errno);
  return ab6_errno;
}
Ab6Error ab6_symmetry_get_bravais(Ab6Symmetry *sym, int bravais[3][3],
				  int *holohedry, int *center, int *nBravSym,
				  Ab6SymmetryMat **bravSym)
{
  Ab6Error ab6_errno;
  int bravais_[9];
  int syms[3 * 3 * AB6_MAX_SYMMETRIES];
  int i;

  SYM_CALL(get_bravais,GET_BRAVAIS, sym, bravais_, holohedry, center, nBravSym, syms, &ab6_errno);
  if (ab6_errno != AB6_NO_ERROR) return ab6_errno;

  *bravSym = g_malloc(sizeof(Ab6SymmetryMat) * (*nBravSym));
  for (i = 0; i < *nBravSym; i++)
    {
      (*bravSym)[i].mat[0][0] = syms[9 * i + 0];
      (*bravSym)[i].mat[0][1] = syms[9 * i + 1];
      (*bravSym)[i].mat[0][2] = syms[9 * i + 2];
      (*bravSym)[i].mat[1][0] = syms[9 * i + 3];
      (*bravSym)[i].mat[1][1] = syms[9 * i + 4];
      (*bravSym)[i].mat[1][2] = syms[9 * i + 5];
      (*bravSym)[i].mat[2][0] = syms[9 * i + 6];
      (*bravSym)[i].mat[2][1] = syms[9 * i + 7];
      (*bravSym)[i].mat[2][2] = syms[9 * i + 8];
    }
  bravais[0][0] = bravais_[0];
  bravais[0][1] = bravais_[1];
  bravais[0][2] = bravais_[2];
  bravais[1][0] = bravais_[3];
  bravais[1][1] = bravais_[4];
  bravais[1][2] = bravais_[5];
  bravais[2][0] = bravais_[6];
  bravais[2][1] = bravais_[7];
  bravais[2][2] = bravais_[8];

  return ab6_errno;
}
Ab6Error ab6_symmetry_get_matrices(Ab6Symmetry *sym, int *nSym,
				   Ab6SymmetryMat **syms,
				   Ab6SymmetryTrans **transNon,
				   int **symAfm)
{
  Ab6Error ab6_errno;
  int syms_[3 * 3 * AB6_MAX_SYMMETRIES];
  int symAfm_[AB6_MAX_SYMMETRIES];
  double trans[3 * AB6_MAX_SYMMETRIES];
  int i;

  SYM_CALL(get_matrices,GET_MATRICES, sym, nSym, syms_, trans, symAfm_, &ab6_errno);
  if (ab6_errno != AB6_NO_ERROR && ab6_errno != AB6_ERROR_SYM_BRAVAIS_XRED) return ab6_errno;

  if (syms)
    *syms = g_malloc(sizeof(Ab6SymmetryMat) * (*nSym));
  if (transNon)
    *transNon = g_malloc(sizeof(Ab6SymmetryTrans) * (*nSym));
  if (symAfm)
    *symAfm = g_malloc(sizeof(int) * (*nSym));
  for (i = 0; i < *nSym; i++)
    {
      if (syms)
	{
	  (*syms)[i].mat[0][0] = syms_[9 * i + 0];
	  (*syms)[i].mat[0][1] = syms_[9 * i + 1];
	  (*syms)[i].mat[0][2] = syms_[9 * i + 2];
	  (*syms)[i].mat[1][0] = syms_[9 * i + 3];
	  (*syms)[i].mat[1][1] = syms_[9 * i + 4];
	  (*syms)[i].mat[1][2] = syms_[9 * i + 5];
	  (*syms)[i].mat[2][0] = syms_[9 * i + 6];
	  (*syms)[i].mat[2][1] = syms_[9 * i + 7];
	  (*syms)[i].mat[2][2] = syms_[9 * i + 8];
	}
      if (transNon)
	{
	  (*transNon)[i].vect[0] = trans[0];
	  (*transNon)[i].vect[1] = trans[1];
	  (*transNon)[i].vect[2] = trans[2];
	}
      if (symAfm)
	(*symAfm)[i] = symAfm_[i];
    }

  return ab6_errno;
}
Ab6Error ab6_symmetry_get_group(Ab6Symmetry *sym,
				char **spaceGroup, int *spaceGroupId,
				int *pointGroupMagn, double genAfm[3])
{
  Ab6Error ab6_errno;
  char name[15];
  int i;

  SYM_CALL(get_group,GET_GROUP, sym, name, spaceGroupId,
	      pointGroupMagn, genAfm, &ab6_errno);
  if (ab6_errno != AB6_NO_ERROR && ab6_errno != AB6_ERROR_SYM_BRAVAIS_XRED) return ab6_errno;

  *spaceGroup = g_malloc(sizeof(char) * 16);
  for (i = 0; i < 15; i++)
    (*spaceGroup)[i] = name[i];
  for (i = 14; i >= 0; i--)
    if ((*spaceGroup)[i] == ' ')
      (*spaceGroup)[i] = '\0';
    else
      break;
  (*spaceGroup)[15] = '\0';

  return ab6_errno;
}
Ab6Error ab6_symmetry_get_equivalent_atom(Ab6Symmetry *sym, int **equiv,
					  int *nSym, int iAtom)
{
  Ab6Error ab6_errno;
  int equiv_[4 * AB6_MAX_SYMMETRIES];
  int i;

  SYM_CALL(get_equivalent_atom,GET_EQUIVALENT_ATOM, sym, equiv_, &iAtom, &ab6_errno);
  if (ab6_errno != AB6_NO_ERROR && ab6_errno != AB6_ERROR_SYM_BRAVAIS_XRED) return ab6_errno;

  ab6_errno = ab6_symmetry_get_n_sym(sym, nSym);
  if (ab6_errno != AB6_NO_ERROR) return ab6_errno;

  *equiv = g_malloc(sizeof(int) * 4 * (*nSym));
  for (i = 0; i < *nSym; i++)
    {
      (*equiv)[i * 4 + 0] = equiv_[i * 4 + 0];
      (*equiv)[i * 4 + 1] = equiv_[i * 4 + 1];
      (*equiv)[i * 4 + 2] = equiv_[i * 4 + 2];
      (*equiv)[i * 4 + 3] = equiv_[i * 4 + 3];
    }

  return ab6_errno;
}
Ab6Error ab6_symmetry_get_type(Ab6Symmetry *sym, int *type, char **label, int iSym)
{
  Ab6Error ab6_errno;
  char label_[128];
  int i;

  SYM_CALL(get_type, GET_TYPE, sym, &iSym, label_, type, &ab6_errno);
  if (ab6_errno != AB6_NO_ERROR && ab6_errno != AB6_ERROR_SYM_BRAVAIS_XRED)
    return ab6_errno;

  *label = g_malloc(sizeof(char) * 129);
  for (i = 0; i < 128; i++)
    (*label)[i] = label_[i];
  for (i = 127; i >= 0; i--)
    if ((*label)[i] == ' ')
      (*label)[i] = '\0';
    else
      break;
  (*label)[128] = '\0';

  return ab6_errno;
}
