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
#include "fix_seaice_wind.h"
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

#define PI 3.1415926536 //pi

/* ---------------------------------------------------------------------- */

FixSeaiceWind::FixSeaiceWind(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix seaice/wind command: not enough arguments");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3; //AH ??????
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  Uxstr = Uystr = Csstr = Cfstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    Uxstr = new char[n];
    strcpy(Uxstr,&arg[3][2]);
  } else {
    Uxvalue = force->numeric(FLERR,arg[3]);
    Uxstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    Uystr = new char[n];
    strcpy(Uystr,&arg[4][2]);
  } else {
    Uyvalue = force->numeric(FLERR,arg[4]);
    Uystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    Csstr = new char[n];
    strcpy(Csstr,&arg[5][2]);
  } else {
    Csvalue = force->numeric(FLERR,arg[5]);
    Csstyle = CONSTANT;
  }
  if (strstr(arg[6],"v_") == arg[6]) {
    int n = strlen(&arg[6][2]) + 1;
    Cfstr = new char[n];
    strcpy(Cfstr,&arg[6][2]);
  } else {
    Cfvalue = force->numeric(FLERR,arg[6]);
    Cfstyle = CONSTANT;
  }

  // optional args

  iregion = -1;
  idregion = NULL;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix seaice/wind command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix seaice/wind does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix seaice/wind command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  maxatom = 0;
  sforce = NULL;
}

/* ---------------------------------------------------------------------- */

FixSeaiceWind::~FixSeaiceWind()
{
  delete [] Uxstr;
  delete [] Uystr;
  delete [] Csstr;
  delete [] Cfstr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixSeaiceWind::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSeaiceWind::init()
{
  // check variables

  if (Uxstr) {
    Uxvar = input->variable->find(Uxstr);
    if (Uxvar < 0)
      error->all(FLERR,"Variable name for fix seaice/wind does not exist");
    if (input->variable->equalstyle(Uxvar)) Uxstyle = EQUAL;
    else if (input->variable->atomstyle(Uxvar)) Uxstyle = ATOM;
    else error->all(FLERR,"Variable for fix seaice/wind is invalid style");
  }
  if (Uystr) {
    Uyvar = input->variable->find(Uystr);
    if (Uyvar < 0)
      error->all(FLERR,"Variable name for fix seaice/wind does not exist");
    if (input->variable->equalstyle(Uyvar)) Uystyle = EQUAL;
    else if (input->variable->atomstyle(Uyvar)) Uystyle = ATOM;
    else error->all(FLERR,"Variable for fix seaice/wind is invalid style");
  }
  if (Csstr) {
    Csvar = input->variable->find(Csstr);
    if (Csvar < 0)
      error->all(FLERR,"Variable name for fix seaice/wind does not exist");
    if (input->variable->equalstyle(Csvar)) Csstyle = EQUAL;
    else if (input->variable->atomstyle(Csvar)) Csstyle = ATOM;
    else error->all(FLERR,"Variable for fix seaice/wind is invalid style");
  }
  if (Cfstr) {
    Cfvar = input->variable->find(Cfstr);
    if (Cfvar < 0)
      error->all(FLERR,"Variable name for fix seaice/wind does not exist");
    if (input->variable->equalstyle(Cfvar)) Cfstyle = EQUAL;
    else if (input->variable->atomstyle(Cfvar)) Cfstyle = ATOM;
    else error->all(FLERR,"Variable for fix seaice/wind is invalid style");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix seaice/wind does not exist");
  }

  rhoair=static_cast<FixPropertyGlobal*>(modify->find_fix_property("rhoAir","property/global","scalar",0,0,style))->compute_scalar();
  if (rhoair<=0) error->all(FLERR,"Fix seaice/wind requires a global property rhoAir > 0");
  
  rhowater=static_cast<FixPropertyGlobal*>(modify->find_fix_property("rhoWater","property/global","scalar",0,0,style))->compute_scalar();
  if (rhowater<=0) error->all(FLERR,"Fix seaice/wind requires a global property rhoWater > 0");  
  
  pirhoa = PI*rhoair;
  freeboardrhoa2 = 2.0*rhoair/rhowater;
  
  if (Uxstyle == ATOM || Uystyle == ATOM || Csstyle == ATOM || Cfstyle == ATOM)
    varflag = ATOM;
  else if (Uxstyle == EQUAL || Uystyle == EQUAL || Csstyle == EQUAL || Cfstyle == EQUAL)
    varflag = EQUAL;
  else { 
    varflag = CONSTANT;
    wspd = sqrt(Uxvalue*Uxvalue+Uyvalue*Uyvalue);
    sstrx = pirhoa*Csvalue*wspd*Uxvalue;
    sstry = pirhoa*Csvalue*wspd*Uyvalue;
    fstrx = freeboardrhoa2*Cfvalue*wspd*Uxvalue;
    fstry = freeboardrhoa2*Cfvalue*wspd*Uyvalue;
  }

//AH??  if ((varflag == EQUAL || varflag == ATOM) &&
//      update->whichflag == 2 && estyle == NONE)
//    error->all(FLERR,"Must use variable energy with fix wind");

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // error checks on coarsegraining
  if(force->cg_active())
    error->cg(FLERR,this->style);
}

/* ---------------------------------------------------------------------- */

void FixSeaiceWind::setup(int vflag)
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

void FixSeaiceWind::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSeaiceWind::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;
  double *r = atom->radius;
  double *thickness = atom->thickness;
  double *density = atom->density;

  // reallocate sforce array if necessary

  if (varflag == ATOM && nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,4,"seaice/wind:sforce");  
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  if (varflag == CONSTANT) {
    double unwrap[3];
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (iregion >= 0 &&
            !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;

        //AH domain->unmap(x[i],image[i],unwrap);
        //AH foriginal[0] -= tmpx*unwrap[0] + tmpy*unwrap[1];
        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];
        f[i][0] += r[i]*(sstrx*r[i] + fstrx*thickness[i]*(rhowater-density[i]));
        f[i][1] += r[i]*(sstry*r[i] + fstry*thickness[i]*(rhowater-density[i]));
      }

  // variable force, wrap with clear/add
  // potential energy = evar if defined, else 0.0
  // wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (Uxstyle == EQUAL) Uxvalue = input->variable->compute_equal(Uxvar);
    else if (Uxstyle == ATOM && sforce)
      input->variable->compute_atom(Uxvar,igroup,&sforce[0][0],4,0);
    if (Uystyle == EQUAL) Uyvalue = input->variable->compute_equal(Uyvar);
    else if (Uystyle == ATOM && sforce)
      input->variable->compute_atom(Uyvar,igroup,&sforce[0][1],4,0);
    if (Csstyle == EQUAL) Csvalue = input->variable->compute_equal(Csvar);
    else if (Csstyle == ATOM && sforce)
      input->variable->compute_atom(Csvar,igroup,&sforce[0][2],4,0);
    if (Cfstyle == ATOM && sforce)
      input->variable->compute_atom(Cfvar,igroup,&sforce[0][3],4,0);

    modify->addstep_compute(update->ntimestep + 1);

    double locUx,locUy,locCs,locCf;
    double tmp;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (iregion >= 0 &&
            !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;

        //AH foriginal[0] += sforce[i][3]; //AH WRONG!!!
        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];
        if (Uxstyle == ATOM) locUx = sforce[i][0];
        else if (Uxstyle) locUx = Uxvalue;
        if (Uystyle == ATOM) locUy = sforce[i][1];
        else if (Uystyle) locUy = Uyvalue;
        if (Csstyle == ATOM) locCs = sforce[i][2];
        else if (Csstyle) locCs = Csvalue;
        if (Cfstyle == ATOM) locCf = sforce[i][3];
        else if (Cfstyle) locCf = Cfvalue;
        wspd = sqrt(locUx*locUx+locUy*locUy);
        tmp = r[i]*wspd*(r[i]*pirhoa*locCs + freeboardrhoa2*locCf*thickness[i]*(rhowater-density[i]));
        f[i][0] += tmp*locUx;
        f[i][1] += tmp*locUy;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSeaiceWind::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSeaiceWind::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixSeaiceWind::compute_scalar()
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

double FixSeaiceWind::compute_vector(int n)
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

double FixSeaiceWind::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*4 * sizeof(double);
  return bytes;
}
