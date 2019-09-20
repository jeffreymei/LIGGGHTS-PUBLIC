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

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "atom_vec_disk.h"
#include "domain_wedge.h"

#ifndef DOMAIN_WEDGE_REAL_H
#define DOMAIN_WEDGE_REAL_H

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

int AtomVecDisk::pack_border_vel_wedge(int n, int *list, double *buf,
                                     int pbc_flag, int *pbc)
{
    return 0;
}

/* ---------------------------------------------------------------------- */

int AtomVecDisk::pack_comm_vel_wedge(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
    return 0;
}

#endif
