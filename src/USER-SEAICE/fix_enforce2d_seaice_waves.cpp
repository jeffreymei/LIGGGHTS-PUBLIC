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

#include "string.h"
#include "fix_enforce2d_seaice_waves.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEnforce2DSeaiceWaves::FixEnforce2DSeaiceWaves(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal fix enforce2d/seaice/waves command");
  if (!atom->tilt_flag) error->all(FLERR,"Fix enforce2d/seaice/waves requires atom style disk/waves");
}

/* ---------------------------------------------------------------------- */

int FixEnforce2DSeaiceWaves::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEnforce2DSeaiceWaves::init()
{
  if (domain->dimension == 3)
    error->all(FLERR,"Cannot use fix enforce2d/seaice/waves with 3d simulation");
}

/* ---------------------------------------------------------------------- */

void FixEnforce2DSeaiceWaves::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    int nlevels_respa = ((Respa *) update->integrate)->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixEnforce2DSeaiceWaves::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEnforce2DSeaiceWaves::post_force(int vflag)
{
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][2] = 0.0;
      f[i][2] = 0.0;
    }
/*
AH: omega/angmom/torque can be nonzero due to the presence of waves

  // for systems with omega/angmom/torque, zero x and y components

  if (atom->omega_flag) {
    double **omega = atom->omega;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        omega[i][0] = 0.0;
        omega[i][1] = 0.0;
      }
  }

  if (atom->angmom_flag) {
    double **angmom = atom->angmom;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        angmom[i][0] = 0.0;
        angmom[i][1] = 0.0;
      }
  }

  if (atom->torque_flag) {
    double **torque = atom->torque;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        torque[i][0] = 0.0;
        torque[i][1] = 0.0;
      }
  }
*/
}
/* ---------------------------------------------------------------------- */

void FixEnforce2DSeaiceWaves::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEnforce2DSeaiceWaves::min_post_force(int vflag)
{
  post_force(vflag);
}
