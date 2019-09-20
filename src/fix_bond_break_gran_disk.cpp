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
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_bond_break_gran_disk.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBondBreakGranDisk::FixBondBreakGranDisk(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix bond/break/gran/disk command");

  MPI_Comm_rank(world,&me);

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix bond/break/gran/disk command");

  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;

  // error check

  if (atom->disk_flag == 0)
    error->all(FLERR,"fix bond/break/gran/disk can be used only with disk-shaped atoms");

  // set comm sizes needed by this fix

  comm_forward = 2;
  comm_reverse = 2;

  // allocate arrays local to this fix

  nmax = 0;
  partner = NULL;

  // zero out stats

  breakcount = 0;
  breakcounttotal = 0;
}

/* ---------------------------------------------------------------------- */

FixBondBreakGranDisk::~FixBondBreakGranDisk()
{
  // delete locally stored arrays
  memory->destroy(partner);
}

/* ---------------------------------------------------------------------- */

int FixBondBreakGranDisk::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondBreakGranDisk::init()
{
  // require special bonds = 0,1,1
/*
  int flag = 0;
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0) flag = 1;
  if (force->special_coul[1] != 0.0 || force->special_coul[2] != 1.0 ||
      force->special_coul[3] != 1.0) flag = 1;
  if (flag) error->all(FLERR,"Fix bond/break/gran/disk requires special_bonds = 0,1,1");
*/
  // warn if angles, dihedrals, impropers are being used

  if (force->angle || force->dihedral || force->improper) {
    if (me == 0)
      error->warning(FLERR,"Broken bonds will not alter angles, "
                     "dihedrals, or impropers");
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixBondBreakGranDisk::post_integrate()
{
  int i,j,k,m,n,i1,i2,n1,n3,type;
  int *slist;

  if (update->ntimestep % nevery) return;

  // need updated ghost atom positions

  comm->forward_comm();

  // resize bond partner list and initialize it
  // probability array overlays distsq array
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(partner);
    nmax = atom->nmax;
    memory->create(partner,nmax,"bond/break/gran/disk:partner");
  }

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    partner[i] = 0;
  }

  // loop over bond list
  // setup partner list of bonds to break

  double **x = atom->x;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    if (!(mask[i1] & groupbit)) continue;
    if (!(mask[i2] & groupbit)) continue;
    if (type != 0) continue;  // bonds of type==0 are broken

    partner[i1] = tag[i2];
    partner[i2] = tag[i1];
  }

  // reverse comm of partner info

  if (force->newton_bond) comm->reverse_comm_fix(this);

  // each atom now knows its winning partner

  comm->forward_comm_fix(this);

  // break bonds
  // if both atoms list each other as winning bond partner
  // and probability constraint is satisfied

  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  int **nspecial = atom->nspecial;
  int **special = atom->special;

  int nbreak = 0;
  for (i = 0; i < nlocal; i++) {
    if (partner[i] == 0) continue;
    j = atom->map(partner[i]);
    if (partner[j] != tag[i]) continue;

    // delete bond from atom I if I stores it
    // atom J will also do this

    for (m = 0; m < num_bond[i]; m++) {
      if (bond_atom[i][m] == partner[i]) {
        for (k = m; k < num_bond[i]-1; k++) {
          bond_atom[i][k] = bond_atom[i][k+1];
          bond_type[i][k] = bond_type[i][k+1];
        }
        num_bond[i]--;
        break;
      }
    }

    // remove J from special bond list for atom I
    // atom J will also do this

    slist = special[i];
    n1 = nspecial[i][0];
    n3 = nspecial[i][2];
    for (m = 0; m < n1; m++)
      if (slist[m] == partner[i]) break;
    for (; m < n3-1; m++) slist[m] = slist[m+1];
    nspecial[i][0]--;
    nspecial[i][1]--;
    nspecial[i][2]--;

    // count the broken bond once

    if (tag[i] < tag[j]) nbreak++;
  }

  // tally stats

  MPI_Allreduce(&nbreak,&breakcount,1,MPI_INT,MPI_SUM,world);
  breakcounttotal += breakcount;
  atom->nbonds -= breakcount;

  // trigger reneighboring if any bonds were formed

  if (breakcount) next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixBondBreakGranDisk::post_integrate_respa(int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */

int FixBondBreakGranDisk::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = partner[j];
  }
  return 2;
}

/* ---------------------------------------------------------------------- */

void FixBondBreakGranDisk::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    partner[i] = static_cast<int> (buf[m++]);
  }
}

/* ---------------------------------------------------------------------- */

int FixBondBreakGranDisk::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = partner[i];
  }
  return 2;
}

/* ---------------------------------------------------------------------- */

void FixBondBreakGranDisk::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
//    if (buf[m+1] > distsq[j]) {
//      partner[j] = static_cast<int> (buf[m++]);
//      distsq[j] = buf[m++];
//    } else m += 2;
      partner[j] = static_cast<int> (buf[m++]);
  }
}

/* ---------------------------------------------------------------------- */

double FixBondBreakGranDisk::compute_vector(int n)
{
  if (n == 1) return (double) breakcount;
  return (double) breakcounttotal;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixBondBreakGranDisk::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  return bytes;
}
