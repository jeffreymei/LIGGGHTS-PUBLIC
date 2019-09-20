/* ----------------------------------------------------------------------

   This file is part of the sea ice toolbox for LIGGGHTS.
   It was created based on the analogous code of LAMMPS and LIGGGHTS.
   See the documentation of the toolbox for details.

   Author: Agnieszka Herman, University of Gdansk, Poland
   http://herman.ocean.ug.edu.pl/LIGGGHTSseaice.html
   agnieszka.herman@ug.edu.pl

------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "compute_erotate_disk.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h" 
#include "fix_multisphere.h" 
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INERTIAZ  0.5               // moment of inertia prefactor for disk (around z-axis)
#define INERTIAXY 0.083333333333333 // moment of inertia prefactor for disk (around x-and y-axes)

/* ---------------------------------------------------------------------- */

ComputeERotateDisk::ComputeERotateDisk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute erotate/disk command");

  scalar_flag = 1;
  extscalar = 1;

  // error check

  if (!atom->disk_flag)
    error->all(FLERR,"Compute erotate/disk requires atom style disk");

  fix_ms = 0; 
}

/* ---------------------------------------------------------------------- */

void ComputeERotateDisk::init()
{
  pfactor = 0.5 * force->mvv2e;  // mvv2e = 1 for SI

  fix_ms =  static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0)); 
}

/* ---------------------------------------------------------------------- */

double ComputeERotateDisk::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **omega = atom->omega;
  double *radius = atom->radius;
  double *thickness = atom->thickness;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // sum rotational energy for each particle
  // point particles will not contribute, due to radius = 0.0

  double erotate = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && (!fix_ms || fix_ms->belongs_to(i) < 0)) 
      erotate += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1])
                 *INERTIAXY*rmass[i]*(3.0*radius[i]*radius[i]+thickness[i]*thickness[i]) +
                  omega[i][2]*omega[i][2]*INERTIAZ*radius[i]*radius[i]*rmass[i];

  MPI_Allreduce(&erotate,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}
