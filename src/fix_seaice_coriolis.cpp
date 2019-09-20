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
#include "stdlib.h"
#include "fix_seaice_coriolis.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "math.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

#define CORCONST 0.0001458 // 2 Omega_Z
#define DEGTORAD 0.01745329251994329576923690768489 // pi/180

/* ---------------------------------------------------------------------- */

FixSeaiceCoriolis::FixSeaiceCoriolis(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix seaice/coriolis command: not enough arguments");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3; //AH ??????
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  phistr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    phistr = new char[n];
    strcpy(phistr,&arg[3][2]);
  } else {
    phivalue = force->numeric(FLERR,arg[3]);
    phistyle = CONSTANT;
  }

  // optional args

  iregion = -1;
  idregion = NULL;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix seaice/coriolis command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix seaice/coriolis does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix seaice/coriolis command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  maxatom = 0;
  sforce = NULL;
}

/* ---------------------------------------------------------------------- */

FixSeaiceCoriolis::~FixSeaiceCoriolis()
{
  delete [] phistr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixSeaiceCoriolis::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSeaiceCoriolis::init()
{
  // check variables
  if (phistr) {
    phivar = input->variable->find(phistr);
    if (phivar < 0)
      error->all(FLERR,"Variable name for fix seaice/coriolis does not exist");
    if (input->variable->equalstyle(phivar)) phistyle = EQUAL;
    else if (input->variable->atomstyle(phivar)) phistyle = ATOM;
    else error->all(FLERR,"Variable for fix seaice/coriolis is invalid style");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix seaice/coriolis does not exist");
  }

  if (phistyle == ATOM)
    varflag = ATOM;
  else if (phistyle == EQUAL)
    varflag = EQUAL;
  else { 
    varflag = CONSTANT;
    fcoeff = CORCONST*sin(phivalue*DEGTORAD); 
  }
//AH??  if ((varflag == EQUAL || varflag == ATOM) &&
//      update->whichflag == 2 && estyle == NONE)
//    error->all(FLERR,"Must use variable energy with fix coriolis");

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // error checks on coarsegraining
  if(force->cg_active())
    error->cg(FLERR,this->style);
}

/* ---------------------------------------------------------------------- */

void FixSeaiceCoriolis::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) 
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixSeaiceCoriolis::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSeaiceCoriolis::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;
  double *rmass = atom->rmass;
  double **v = atom->v;

  // reallocate sforce array if necessary

  if (varflag == ATOM && nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,1,"seaice/coriolis:sforce");  
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // constant force

  if (varflag == CONSTANT) {
    double unwrap[3];
    double tmp;
    
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (iregion >= 0 &&
            !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;

        tmp = fcoeff*rmass[i];
        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];
        f[i][0] += tmp*v[i][1]; // +fv
        f[i][1] -= tmp*v[i][0]; // -fu
      }

  } else {

    modify->clearstep_compute();

    if (phistyle == EQUAL) phivalue = input->variable->compute_equal(phivar);
    else if (phistyle == ATOM && sforce)
      input->variable->compute_atom(phivar,igroup,&sforce[0][0],1,0);

    modify->addstep_compute(update->ntimestep + 1);

    double locphi;
    double tmp;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (iregion >= 0 &&
            !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;

        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];
        if (phistyle == ATOM) locphi = sforce[i][0];
        else if (phistyle) locphi = phivalue;
        tmp = CORCONST*rmass[i]*sin(phivalue*DEGTORAD);
        f[i][0] += tmp*v[i][1]; // +fv
        f[i][1] -= tmp*v[i][0]; // -fu
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSeaiceCoriolis::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSeaiceCoriolis::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixSeaiceCoriolis::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixSeaiceCoriolis::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSeaiceCoriolis::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*4 * sizeof(double);
  return bytes;
}
