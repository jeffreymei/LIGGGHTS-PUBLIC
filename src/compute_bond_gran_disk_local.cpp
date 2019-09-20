/* ----------------------------------------------------------------------

   This file is part of the sea ice toolbox for LIGGGHTS.
   It was created based on the analogous code of LAMMPS and LIGGGHTS.
   See the documentation of the toolbox for details.

   Author: Agnieszka Herman, University of Gdansk, Poland
   http://herman.ocean.ug.edu.pl/LIGGGHTSseaice.html
   agnieszka.herman@ug.edu.pl

------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "compute_bond_gran_disk_local.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "bond.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

enum{THICKNESS,LENGTH,FORCEN,FORCET,TORQUEN,TORQUETH,TORQUETZ};

/* ---------------------------------------------------------------------- */

ComputeBondGranDiskLocal::ComputeBondGranDiskLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute bond/gran/disk/local command");

  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Compute bond/gran/disk/local used when bonds are not allowed");

  local_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  int numbondhist = atom->n_bondhist;

  bstyle = new int[nvalues];
  nvalues = 0;
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"thickness") == 0) bstyle[nvalues++] = THICKNESS;
    else if (strcmp(arg[iarg],"length") == 0) bstyle[nvalues++] = LENGTH;
    else if (strcmp(arg[iarg],"forcen") == 0) bstyle[nvalues++] = FORCEN;
    else if (strcmp(arg[iarg],"forcet") == 0) bstyle[nvalues++] = FORCET;    
    else if (strcmp(arg[iarg],"torquen") == 0) {
      if (numbondhist==6)  // 2D
        error->all(FLERR,"Keyword TORQUEN in compute bond/gran/disk/local command can only be used with quasi-3D atom style");
      else {
        bstyle[nvalues++] = TORQUEN;
      }
    }
    else if (strcmp(arg[iarg],"torqueth") == 0) {
      if (numbondhist==6)  // 2D
        error->all(FLERR,"Keyword TORQUETH in compute bond/gran/disk/local command can only be used with quasi-3D atom style");
      else bstyle[nvalues++] = TORQUETH;
    }
    else if (strcmp(arg[iarg],"torquetz") == 0) bstyle[nvalues++] = TORQUETZ;
    else error->all(FLERR,"Invalid keyword in compute bond/gran/disk/local command");
  }

  // set singleflag if need to call bond->single()
  singleflag = 0;
  for (int i = 0; i < nvalues; i++)
    if (bstyle[i] == LENGTH) singleflag = 1;

  nmax = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBondGranDiskLocal::~ComputeBondGranDiskLocal()
{
  memory->destroy(vector);
  memory->destroy(array);
  delete [] bstyle;
}

/* ---------------------------------------------------------------------- */

void ComputeBondGranDiskLocal::init()
{
  if (force->bond == NULL)
    error->all(FLERR,"No bond style is defined for compute bond/gran/disk/local");

  // do initial memory allocation so that memory_usage() is correct

  ncount = compute_bonds(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeBondGranDiskLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute bond info

  ncount = compute_bonds(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_bonds(1);
}

/* ----------------------------------------------------------------------
   count bonds and compute bond info on this proc
   only count bond once if newton_bond is off
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if bond is deleted (type = 0), do not count
   if bond is turned off (type < 0), still count
   if flag is set, compute requested info about bond
   if bond is turned off (type < 0), energy = 0.0
------------------------------------------------------------------------- */

int ComputeBondGranDiskLocal::compute_bonds(int flag)
{
  int i,m,n,atom1,atom2;
  double rsum;
  double *ptr; 

  double *radius = atom->radius;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  double ***bondhist = atom->bond_hist;
  int numbondhist = atom->n_bondhist;

  Bond *bond = force->bond;
  double thbond,lbond,fbondn,fbondt,tbondn,tbondth,tbondtz,tmp;

  m = n = 0;
  for (atom1 = 0; atom1 < nlocal; atom1++) {
    if (!(mask[atom1] & groupbit)) continue;
    for (i = 0; i < num_bond[atom1]; i++) {
      atom2 = atom->map(bond_atom[atom1][i]);
      if (atom2 < 0 || !(mask[atom2] & groupbit)) continue;
      if (newton_bond == 0 && tag[atom1] > tag[atom2]) continue;
      if (bond_type[atom1][i] == 0) continue;

      if (flag) {
        rsum = radius[atom1]+radius[atom2];

        if (bond_type[atom1][i] > 0) { 
          if (singleflag) lbond = bond->single(bond_type[atom1][i],rsum,atom1,atom2,tmp);
//          fprintf(screen,"atom1 = %d, i = %d, bondhist[5] = %f\n",atom1,i,bondhist[atom1][i][5]);
          fbondn = sqrt(bondhist[atom1][i][0]*bondhist[atom1][i][0]+bondhist[atom1][i][1]*bondhist[atom1][i][1]);
          if (numbondhist==6) { // 2D
            thbond = bondhist[atom1][i][5];
            fbondt = sqrt(bondhist[atom1][i][3]*bondhist[atom1][i][3]+bondhist[atom1][i][4]*bondhist[atom1][i][4]);
            tbondtz = bondhist[atom1][i][4];
          } else {              // quasi-3D
            thbond = bondhist[atom1][i][10];
            fbondt = sqrt(bondhist[atom1][i][3]*bondhist[atom1][i][3]+bondhist[atom1][i][4]*bondhist[atom1][i][4]
                                                                     +bondhist[atom1][i][5]*bondhist[atom1][i][5]);
            tbondn  = sqrt(bondhist[atom1][i][5]*bondhist[atom1][i][5]+bondhist[atom1][i][6]*bondhist[atom1][i][6]);
            tbondth = sqrt(bondhist[atom1][i][7]*bondhist[atom1][i][7]+bondhist[atom1][i][8]*bondhist[atom1][i][8]);
            tbondtz = bondhist[atom1][i][9];
            //fprintf(screen,"thbond = %f, fbondt = %f, tbondn = %f\n",thbond,fbondt,tbondn);
          }
        }
        else lbond = thbond = fbondn = fbondt = tbondn = tbondth = tbondtz = 0.0;

        if (nvalues == 1) ptr = &vector[m];
        else ptr = array[m];

        for (n = 0; n < nvalues; n++) {
          switch (bstyle[n]) {
          case THICKNESS:
            ptr[n] = thbond;
            break;
          case LENGTH:
            ptr[n] = lbond;
            break;
          case FORCEN:
            ptr[n] = fbondn;
            break;
          case FORCET:
            ptr[n] = fbondt;
            break;
          case TORQUEN:
            ptr[n] = tbondn;
            break;
          case TORQUETH:
            ptr[n] = tbondth;
            break;
          case TORQUETZ:
            ptr[n] = tbondtz;
            break;
          }
        }
      }

      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeBondGranDiskLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vector);
    memory->create(vector,nmax,"bond/gran/disk/local:vector");
    vector_local = vector;
  } else {
    memory->destroy(array);
    memory->create(array,nmax,nvalues,"bond/gran/disk/local:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeBondGranDiskLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}
