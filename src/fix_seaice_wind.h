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

FixStyle(seaice/wind,FixSeaiceWind)

#else

#ifndef LMP_FIX_SEAICE_WIND_H
#define LMP_FIX_SEAICE_WIND_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSeaiceWind : public Fix {
 public:
  FixSeaiceWind(class LAMMPS *, int, char **);
  ~FixSeaiceWind();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  double Uxvalue,Uyvalue,Csvalue,Cfvalue;
  int varflag,iregion;
  char *Uxstr,*Uystr,*Csstr,*Cfstr;
  char *idregion;
  int Uxvar,Uyvar,Csvar,Cfvar,Uxstyle,Uystyle,Csstyle,Cfstyle;
  double foriginal[4],foriginal_all[4];
  int force_flag;
  int nlevels_respa;

  int maxatom;
  double **sforce;
  
  double wspd,sstrx,sstry,fstrx,fstry; 
  double rhoair,rhowater;
  double pirhoa,freeboardrhoa2;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix addforce does not exist

Self-explanatory.

E: Variable name for fix addforce does not exist

Self-explanatory.

E: Variable for fix addforce is invalid style

Self-explanatory.

E: Cannot use variable energy with constant force in fix addforce

This is because for constant force, LAMMPS can compute the change
in energy directly.

E: Must use variable energy with fix addforce

Must define an energy vartiable when applyting a dynamic
force during minimization.

*/
