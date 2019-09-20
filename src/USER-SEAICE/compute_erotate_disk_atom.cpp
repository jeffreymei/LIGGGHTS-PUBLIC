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
#include "string.h"
#include "compute_erotate_disk_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "fix_multisphere.h" 

using namespace LAMMPS_NS;

#define INERTIAZ  0.5                 // moment of inertia prefactor for disk (around z-axis)
#define INERTIAXY 0.08333333333333333 // moment of inertia prefactor for disk (around x- and y-axes)

/* ---------------------------------------------------------------------- */

ComputeErotateDiskAtom::
ComputeErotateDiskAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3)
    error->all(FLERR,"Illegal compute erotate/disk/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  // error check

  if (!atom->disk_flag)
    error->all(FLERR,"Compute erotate/disk/atom requires atom style disk");

  nmax = 0;
  erot = NULL;

  fix_ms = 0; 
}

/* ---------------------------------------------------------------------- */

ComputeErotateDiskAtom::~ComputeErotateDiskAtom()
{
  memory->destroy(erot);
}

/* ---------------------------------------------------------------------- */

void ComputeErotateDiskAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"erotate/disk/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute erotate/disk/atom");

  pfactor = 0.5 * force->mvv2e; // mvv2e = 1 for SI

  fix_ms =  static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0)); 
}

/* ---------------------------------------------------------------------- */

void ComputeErotateDiskAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow erot array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(erot);
    nmax = atom->nmax;
    memory->create(erot,nmax,"erotate/disk/atom:erot");
    vector_atom = erot;
  }

  // compute rotational kinetic energy for each atom in group
  // point particles will have erot = 0.0, due to radius = 0.0

  double **omega = atom->omega;
  double *radius = atom->radius;
  double *thickness = atom->thickness;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && (!fix_ms || fix_ms->belongs_to(i) < 0)) 
    {
      erot[i] = (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1])*
                INERTIAXY*rmass[i]*(3.0*radius[i]*radius[i]+thickness[i]*thickness[i]) +
                INERTIAZ*omega[i][2]*omega[i][2]*radius[i]*radius[i]*rmass[i];
      erot[i] *= pfactor;
    } else erot[i] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeErotateDiskAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
