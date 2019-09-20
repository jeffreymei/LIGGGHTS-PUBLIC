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

#ifdef FIX_CLASS

FixStyle(bond/break/gran/disk,FixBondBreakGranDisk)

#else

#ifndef LMP_FIX_BOND_BREAK_GRAN_DISK_H
#define LMP_FIX_BOND_BREAK_GRAN_DISK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondBreakGranDisk : public Fix {
 public:
  FixBondBreakGranDisk(class LAMMPS *, int, char **);
  ~FixBondBreakGranDisk();
  int setmask();
  void init();
  void post_integrate();
  void post_integrate_respa(int,int);

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double compute_vector(int);
  double memory_usage();

 private:
  int me;
  //int btype,seed;
  //double cutsq,fraction;

  int breakcount,breakcounttotal;
  int nmax;
  int *partner;
  //double *distsq,*probability;

  //class RanMars *random;
  int nlevels_respa;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid bond type in fix bond/break command

Self-explanatory.

E: Cannot use fix bond/break with non-molecular systems

Self-explanatory.

E: Fix bond/break requires special_bonds = 0,1,1

This is a restriction of the current fix bond/break implementation.

W: Broken bonds will not alter angles, dihedrals, or impropers

See the doc page for fix bond/break for more info on this
restriction.

*/
