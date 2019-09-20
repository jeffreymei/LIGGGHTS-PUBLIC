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
#include "fix_seaice_current.h"
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
#include "fix_seaice_waves_extra.h"
#include "random_park.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

#define deg2rad 0.0174532925 	// pi/180
#define coefk 4.024303537457	// 4pi^2/g

/* ---------------------------------------------------------------------- */

FixSeaiceCurrent::FixSeaiceCurrent(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix seaice/current command: not enough arguments");

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

  wavecur = wavever = 0;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix seaice/current command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix seaice/current does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"wavecur") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix seaice/current command");
      if (strcmp(arg[iarg+1],"jonswap") == 0) {
        wavecur = 1;
        wavever = 1;
      } else if (strcmp(arg[iarg+1],"regWaves") == 0) {
        wavecur = 1;
        wavever = 2;
      } else error->all(FLERR,"Option wavecur in fix seaice/current requires style jonswap or regWaves");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix seaice/current command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  maxatom = 0;
  sforce = NULL;

}

/* ---------------------------------------------------------------------- */

FixSeaiceCurrent::~FixSeaiceCurrent()
{
  delete [] Uxstr;
  delete [] Uystr;
  delete [] Csstr;
  delete [] Cfstr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixSeaiceCurrent::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSeaiceCurrent::init()
{
  // check variables

  if (Uxstr) {
    Uxvar = input->variable->find(Uxstr);
    if (Uxvar < 0)
      error->all(FLERR,"Variable name for fix seaice/current does not exist");
    if (input->variable->equalstyle(Uxvar)) Uxstyle = EQUAL;
    else if (input->variable->atomstyle(Uxvar)) Uxstyle = ATOM;
    else error->all(FLERR,"Variable for fix seaice/current is invalid style");
  }
  if (Uystr) {
    Uyvar = input->variable->find(Uystr);
    if (Uyvar < 0)
      error->all(FLERR,"Variable name for fix seaice/current does not exist");
    if (input->variable->equalstyle(Uyvar)) Uystyle = EQUAL;
    else if (input->variable->atomstyle(Uyvar)) Uystyle = ATOM;
    else error->all(FLERR,"Variable for fix seaice/current is invalid style");
  }
  if (Csstr) {
    Csvar = input->variable->find(Csstr);
    if (Csvar < 0)
      error->all(FLERR,"Variable name for fix seaice/current does not exist");
    if (input->variable->equalstyle(Csvar)) Csstyle = EQUAL;
    else if (input->variable->atomstyle(Csvar)) Csstyle = ATOM;
    else error->all(FLERR,"Variable for fix seaice/current is invalid style");
  }
  if (Cfstr) {
    Cfvar = input->variable->find(Cfstr);
    if (Cfvar < 0)
      error->all(FLERR,"Variable name for fix seaice/current does not exist");
    if (input->variable->equalstyle(Cfvar)) Cfstyle = EQUAL;
    else if (input->variable->atomstyle(Cfvar)) Cfstyle = ATOM;
    else error->all(FLERR,"Variable for fix seaice/current is invalid style");
  }
  
  // wave-induced current
  
  if (wavecur==1) {
    int i;
    if (wavever==1) {
      jonspar = static_cast<FixPropertyGlobal*>(modify->find_fix_property("jonswap","property/global","vector",11,0,style))->values;
      if (jonspar[0]<=0.0)   error->all(FLERR,"Fix seaice/waves with option jonswap requires Hs > 0");  
      if (jonspar[1]<=0.0)   error->all(FLERR,"Fix seaice/waves with option jonswap requires Tp > 0");  
      if (jonspar[2]<-180.0) error->all(FLERR,"Fix seaice/waves with option jonswap requires theta_m >= -180");  
      if (jonspar[2]>180.0)  error->all(FLERR,"Fix seaice/waves with option jonswap requires theta_m <= 180");  
      if (jonspar[3]<=0.0)   error->all(FLERR,"Fix seaice/waves with option jonswap requires gamma > 0");  
      if (jonspar[4]<=0.0)   error->all(FLERR,"Fix seaice/waves with option jonswap requires dir_spread > 0"); 
      if (jonspar[5]<=0.0)   error->all(FLERR,"Fix seaice/waves with option jonswap requires Nfreqs > 0");       
      nw = (int)jonspar[5]*(int)jonspar[8]; // number elementary of waves = Nf*Nd
      if (jonspar[6]<=0.0)   error->all(FLERR,"Fix seaice/waves with option jonswap requires fmin > 0");  
      if (jonspar[7]<=0.0)   error->all(FLERR,"Fix seaice/waves with option jonswap requires fmax > 0");  
      if (jonspar[7]<=jonspar[6]) error->all(FLERR,"Fix seaice/waves with option jonswap requires fmax > fmin");  
      if (jonspar[8]<=0.0)   error->all(FLERR,"Fix seaice/waves with option jonswap requires Ndirs > 0");      
      if (jonspar[9]<-180.0) error->all(FLERR,"Fix seaice/waves with option jonswap requires thetamin >= -180");  
      if (jonspar[10]>180.0) error->all(FLERR,"Fix seaice/waves with option jonswap requires thetamax <= 180");  
      if (jonspar[10]<=jonspar[9]) error->all(FLERR,"Fix seaice/waves with option jonswap requires thetamax > thetamin");      
      if (jonspar[2]<jonspar[9]) error->all(FLERR,"Fix seaice/waves with option jonswap: theta_m < thetamin");      
      if (jonspar[2]>jonspar[10]) error->all(FLERR,"Fix seaice/waves with option jonswap: theta_m > thetamax");              
      memory->create(regwave,nw,4,"seaice/waves:regwave");  
      RanPark *random = new RanPark(lmp,100+comm->me);
      SeaiceWavesExtra::jonswapToRegWave(nw,jonspar,regwave,random);
    } else if (wavever==2) {
      nw = static_cast<FixPropertyGlobal*>(modify->find_fix_property("regWaves","property/global","matrix",0,0,style))->size_array_cols;
      memory->create(regwave,nw,4,"seaice/waves:regwave");  
      double** tmp = static_cast<FixPropertyGlobal*>(modify->find_fix_property("regWaves","property/global","matrix",nw,4,style))->array;
      for (i = 0; i < nw; i++) {    
        regwave[i][0] = tmp[0][i];
        regwave[i][1] = tmp[1][i];
        regwave[i][2] = tmp[2][i];
        regwave[i][3] = tmp[3][i];
        if (regwave[i][0]<=0.0)   error->all(FLERR,"Fix seaice/waves with option regWaves requires a > 0");  
        if (regwave[i][1]<=0.0)   error->all(FLERR,"Fix seaice/waves with option regWaves requires T > 0");  
        if (regwave[i][2]<-180.0) error->all(FLERR,"Fix seaice/waves with option regWaves requires theta >= -180");  
        if (regwave[i][2]>180.0)  error->all(FLERR,"Fix seaice/waves with option regWaves requires theta <= 180");  
        if (regwave[i][3]<-180.0) error->all(FLERR,"Fix seaice/waves with option regWaves requires phi >= -180");  
        if (regwave[i][3]>180.0)  error->all(FLERR,"Fix seaice/waves with option regWaves requires phi <= 180");  
      }  
    }
    memory->create(wavenum,nw,2,"seaice/waves:wavenum");
    memory->create(wccoefx,nw,"seaice/waves:wccoefx");
    memory->create(wccoefy,nw,"seaice/waves:wccoefy");
    for (i = 0; i < nw; i++) {    
      regwave[i][2] *= deg2rad;
      regwave[i][3] *= deg2rad;
      wavenum[i][0] = coefk/regwave[i][1]/regwave[i][1]*cos(regwave[i][2]); // kx
      wavenum[i][1] = coefk/regwave[i][1]/regwave[i][1]*sin(regwave[i][2]); // ky
      wccoefx[i] = MY_2PI*regwave[i][0]/regwave[i][1]*cos(regwave[i][2]);
      wccoefy[i] = MY_2PI*regwave[i][0]/regwave[i][1]*sin(regwave[i][2]);
      fprintf(screen,"i = %d, a = %f, T = %f, kx = %f\n",i,regwave[i][0],regwave[i][1],wavenum[i][0]);
    }
    
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix seaice/current does not exist");
  }

  rhowater=static_cast<FixPropertyGlobal*>(modify->find_fix_property("rhoWater","property/global","scalar",0,0,style))->compute_scalar();
  if (rhowater<=0) error->all(FLERR,"Fix seaice/current requires a global property rhoWater > 0");  
  
  pirhow = MY_PI*rhowater;

  if (Uxstyle == ATOM || Uystyle == ATOM || Csstyle == ATOM || Cfstyle == ATOM)
    varflag = ATOM;
  else if (Uxstyle == EQUAL || Uystyle == EQUAL || Csstyle == EQUAL || Cfstyle == EQUAL)
    varflag = EQUAL;
  else { 
    varflag = CONSTANT;
    sstrx = pirhow*Csvalue;
    sstry = pirhow*Csvalue;
    fstrx = 2.0*Cfvalue;
    fstry = 2.0*Cfvalue;
    sstrang = pirhow*Csvalue/2.0;
    fstrang = MY_PI*Cfvalue/2.0;
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

void FixSeaiceCurrent::setup(int vflag)
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

void FixSeaiceCurrent::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSeaiceCurrent::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  double **omega = atom->omega;
  double **torque = atom->torque;
  int *mask = atom->mask;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;
  double *r = atom->radius;
  double *thickness = atom->thickness;
  double *density = atom->density;
  double relspd,relspdx,relspdy;

  // reallocate sforce array if necessary

  if (varflag == ATOM && nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,4,"seaice/current:sforce");  
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  if (varflag == CONSTANT) {
    double unwrap[3];
    // without wave-induced current:
    if (wavecur==0) {
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
          relspdx = Uxvalue - v[i][0];
          relspdy = Uyvalue - v[i][1];
          relspd = sqrt(relspdx*relspdx+relspdy*relspdy);
          f[i][0] += r[i]*relspd*relspdx*(sstrx*r[i] + fstrx*thickness[i]*density[i]);
          f[i][1] += r[i]*relspd*relspdy*(sstry*r[i] + fstry*thickness[i]*density[i]);
          torque[i][2] -= (sstrang*r[i]+fstrang*thickness[i]*density[i])*r[i]*r[i]*r[i]*omega[i][2];
        }
    } else { 
    // with wave-induced current:
      double locUx,locUy,krcoef,cosphi;
      double tnow = (update->ntimestep-update->beginstep)*update->dt;
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
          locUx = Uxvalue;
          locUy = Uyvalue;
          for (int k = 0; k < nw; k++) {
            krcoef = 1.0;
            if (wavenum[k][0]>0.0) krcoef *= sin(wavenum[k][0]*r[i])/(wavenum[k][0]*r[i]);
            if (wavenum[k][1]>0.0) krcoef *= sin(wavenum[k][1]*r[i])/(wavenum[k][1]*r[i]);
            cosphi = cos(wavenum[k][0]*x[i][0]+wavenum[k][1]*x[i][1]-MY_2PI/regwave[k][1]*tnow+regwave[k][3]);
            locUx += wccoefx[k]*krcoef*cosphi;
            locUy += wccoefy[k]*krcoef*cosphi;
          }
          relspdx = locUx - v[i][0];
          relspdy = locUy - v[i][1];
          relspd = sqrt(relspdx*relspdx+relspdy*relspdy);
          f[i][0] += r[i]*relspd*relspdx*(sstrx*r[i] + fstrx*thickness[i]*density[i]);
          f[i][1] += r[i]*relspd*relspdy*(sstry*r[i] + fstry*thickness[i]*density[i]);
          torque[i][2] -= (sstrang*r[i]+fstrang*thickness[i]*density[i])*r[i]*r[i]*r[i]*omega[i][2];
        }
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

    // without wave-induced current:
    if (wavecur==0) {
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
          relspdx = locUx - v[i][0];
          relspdy = locUy - v[i][1];
          relspd = sqrt(relspdx*relspdx+relspdy*relspdy);
          tmp = r[i]*relspd*(r[i]*pirhow*locCs + 2.0*thickness[i]*density[i]*locCf);
          f[i][0] += tmp*locUx;
          f[i][1] += tmp*locUy;
          torque[i][2] -= MY_PI/2.0*(rhowater*locCs*r[i]+density[i]*thickness[i]*locCf)*r[i]*r[i]*r[i]*omega[i][2];
        }
    } else {
    // with wave-induced current:
      double krcoef,cosphi,tnow;    
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
          tnow = (update->ntimestep-update->beginstep)*update->dt;
          for (int k = 0; k < nw; k++) {
            krcoef = sin(wavenum[k][0]*r[i])*sin(wavenum[k][1]*r[i])/(wavenum[k][0]*wavenum[k][1]*r[i]*r[i]);
            cosphi = cos(wavenum[k][0]*x[i][0]+wavenum[k][1]*x[i][1]-MY_2PI/regwave[k][1]*tnow+regwave[k][3]);
            locUx += wccoefx[k]*krcoef*cosphi;
            locUy += wccoefy[k]*krcoef*cosphi;
          }
          relspdx = locUx - v[i][0];
          relspdy = locUy - v[i][1];
          relspd = sqrt(relspdx*relspdx+relspdy*relspdy);
          tmp = r[i]*relspd*(r[i]*pirhow*locCs + 2.0*thickness[i]*density[i]*locCf);
          f[i][0] += tmp*locUx;
          f[i][1] += tmp*locUy;
          torque[i][2] -= MY_PI/2.0*(rhowater*locCs*r[i]+density[i]*thickness[i]*locCf)*r[i]*r[i]*r[i]*omega[i][2];
        }      
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSeaiceCurrent::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSeaiceCurrent::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixSeaiceCurrent::compute_scalar()
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

double FixSeaiceCurrent::compute_vector(int n)
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

double FixSeaiceCurrent::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*4 * sizeof(double);
  return bytes;
}
