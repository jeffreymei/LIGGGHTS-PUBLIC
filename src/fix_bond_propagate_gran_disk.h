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

#ifdef FIX_CLASS

FixStyle(bond/propagate/gran/disk,FixBondPropagateGranDisk)

#else

#ifndef LMP_FIX_BOND_PROPAGATE_GRAN_DISK_H
#define LMP_FIX_BOND_PROPAGATE_GRAN_DISK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondPropagateGranDisk : public Fix {
 public:
  FixBondPropagateGranDisk(class LAMMPS *, int, char **);
  ~FixBondPropagateGranDisk();
  int setmask();
  void post_integrate();
  void post_integrate_respa(int,int);  
  void pre_exchange();
  void write_restart(FILE *);
  void restart(char *);

 private:
  void remove_bond(int ilocal,int ibond, int bondnumber);
  int nlevels_respa;
};

}

#endif
#endif
