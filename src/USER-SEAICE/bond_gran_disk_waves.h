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

#ifdef BOND_CLASS

BondStyle(gran/disk/waves,BondGranDiskWaves)

#else

#ifndef LMP_BOND_GRAN_DISK_WAVES_H
#define LMP_BOND_GRAN_DISK_WAVES_H

#include "stdio.h"
#include "bond.h"

namespace LAMMPS_NS {

class BondGranDiskWaves : public Bond {
 public:
  BondGranDiskWaves(class LAMMPS *);
  ~BondGranDiskWaves();
  void init_style();
  void compute(int, int);
  void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, double, int, int, double &);

 protected:
  int breakmode;
  double *rb;        // half bond width, relative to the radius of the smaller grain in pair (number between 0 and 1)
  double *lb;        // bond length, relative to the sum of the grain radii (number between 0 and 1)
  double *thb;       // mean bond thickness
  double *thstdb;    // std.dev. of bond thickness  
  double *Ec,*kn2kt; // Young modulus of bonds and their kn/kt ratio
  double *sigman_break_c,*sigman_break_t,*tau_break;
  void allocate();

 private:
  int j0, j1, j3, j4, j5, j6, j7, j9, j10, j11, j13;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/
