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
#include "stdio.h"
#include "string.h"
#include "fix_nve_disk.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "domain.h" 

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIAZ  0.5                 // moment of inertia prefactor for disk (around z-axis)
#define INERTIAXY 0.08333333333333333 // moment of inertia prefactor for disk (around x- and y-axes)

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

FixNVEDisk::FixNVEDisk(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nve/disk command");

  time_integrate = 1;

  // process extra keywords

  extra = NONE;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"update") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix nve/disk command");
      if (strcmp(arg[iarg+1],"dipole") == 0) extra = DIPOLE;
      else error->all(FLERR,"Illegal fix nve/disk command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix nve/disk command");
  }

  // error checks

  if (!atom->disk_flag)
    error->all(FLERR,"Fix nve/disk requires atom style disk");
  if (extra == DIPOLE && !atom->mu_flag)
    error->all(FLERR,"Fix nve/disk requires atom attribute mu");
}

/* ---------------------------------------------------------------------- */

void FixNVEDisk::init()
{
  FixNVE::init();

  // check that all particles are finite-size disks
  // no point particles allowed

  double *radius = atom->radius;
  double *thickness = atom->thickness;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0 || thickness[i] == 0.0)
        error->one(FLERR,"Fix nve/disk requires extended particles");
}

/* ---------------------------------------------------------------------- */

void FixNVEDisk::initial_integrate(int vflag)
{
  double dtfm,dtirotate,msq,scale;
  double g[3];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *thickness = atom->thickness;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotatez,dtfrotatexy; 
  dtfrotatez  = dtf / INERTIAZ;
  dtfrotatexy = dtf / INERTIAXY;

  // update v,x,omega for all particles
  // d_omega/dt = torque / inertia

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      dtirotate = dtfrotatexy / ((3.0*radius[i]*radius[i]+thickness[i]*thickness[i])*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      dtirotate = dtfrotatez / (radius[i]*radius[i]*rmass[i]);
      omega[i][2] += dtirotate * torque[i][2];
    }
  }

  if (atom->tilt_flag) {
    double **tilt = atom->tilt;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        tilt[i][0] += dtv * omega[i][0];
        tilt[i][1] += dtv * omega[i][1];
      }
    }
  }
/*
  // update mu for dipoles
  // d_mu/dt = omega cross mu
  // renormalize mu to dipole length

  if (extra == DIPOLE) {
    double **mu = atom->mu;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        if (mu[i][3] > 0.0) {
          g[0] = mu[i][0] + dtv * (omega[i][1]*mu[i][2]-omega[i][2]*mu[i][1]);
          g[1] = mu[i][1] + dtv * (omega[i][2]*mu[i][0]-omega[i][0]*mu[i][2]);
          g[2] = mu[i][2] + dtv * (omega[i][0]*mu[i][1]-omega[i][1]*mu[i][0]);
          msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
          scale = mu[i][3]/sqrt(msq);
          mu[i][0] = g[0]*scale;
          mu[i][1] = g[1]*scale;
          mu[i][2] = g[2]*scale;
        }
  }
*/
}

/* ---------------------------------------------------------------------- */

void FixNVEDisk::final_integrate()
{
  double dtfm,dtirotate;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  double *thickness = atom->thickness;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotatez,dtfrotatexy; 
  dtfrotatez  = dtf / INERTIAZ;
  dtfrotatexy = dtf / INERTIAXY;

  // update v,omega for all particles
  // d_omega/dt = torque / inertia

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      dtirotate = dtfrotatexy / ((3.0*radius[i]*radius[i]+thickness[i]*thickness[i])*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      dtirotate = dtfrotatez / (radius[i]*radius[i]*rmass[i]);
      omega[i][2] += dtirotate * torque[i][2];
    }

  if (atom->tilt_flag) {
    double **tilt = atom->tilt;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        tilt[i][0] += dtv * omega[i][0];
        tilt[i][1] += dtv * omega[i][1];
      }
    }
  }    
}
