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
#include "fix_seaice_waves.h"
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
#include "fix_seaice_waves_extra.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

#define PIG4 7.70475598292 		// pi*g/4
#define deg2rad 0.01745329251994330 	// pi/180
#define coefk 4.024303537457	// 4pi^2/g
#define twopi 6.28318530718     // 2pi

/* ---------------------------------------------------------------------- */

FixSeaiceWaves::FixSeaiceWaves(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix seaice/waves command: not enough arguments");
  
  if (!atom->disk_flag) error->all(FLERR,"Fix seaice/waves requires atom style disk");
  if (!atom->tilt_flag) error->all(FLERR,"Fix seaice/waves requires atom style disk/waves");  

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3; //AH ??????
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  wavever = 0;
  if (strcmp(arg[3],"jonswap") == 0) {
    wavever = 1;
  } else if (strcmp(arg[3],"regWaves") == 0) {
    wavever = 2;
  } else {
    error->all(FLERR,"Fix seaice/waves requires style jonswap or regWaves");
  }
  
  // optional args

  iregion = -1;
  idregion = NULL;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix seaice/waves command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix seaice/waves does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix seaice/waves command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

}

/* ---------------------------------------------------------------------- */

FixSeaiceWaves::~FixSeaiceWaves()
{
  delete [] idregion;
}

/* ---------------------------------------------------------------------- */

int FixSeaiceWaves::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSeaiceWaves::init()
{
  int i;
  // check variables

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
    RanPark *random = new RanPark(lmp,200+comm->me);
    SeaiceWavesExtra::jonswapToRegWave(nw,jonspar,regwave,random);
    delete random;
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
  for (i = 0; i < nw; i++) {    
    printf("%f %f %f %f\n",regwave[i][0],regwave[i][1],regwave[i][2],regwave[i][3]);
    regwave[i][2] *= deg2rad;
    regwave[i][3] *= deg2rad;
        fprintf(screen," regwave = %e\n",regwave[i][2]);            
    wavenum[i][0] = coefk/regwave[i][1]/regwave[i][1]*cos(regwave[i][2]); // kx
    wavenum[i][1] = coefk/regwave[i][1]/regwave[i][1]*sin(regwave[i][2]); // ky
        fprintf(screen," wavenumx = %e, wavenumy = %f\n",wavenum[i][0],wavenum[i][1]);    
      // UWAGA!!! Dzielimy T przez n dla otrzymania szybciej/wolniej propagujących się fal!!!
    //  regwave[i][1] /= 0.5;
  }
  
  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix seaice/waves does not exist");
  }

  rhowater=static_cast<FixPropertyGlobal*>(modify->find_fix_property("rhoWater","property/global","scalar",0,0,style))->compute_scalar();
  if (rhowater<=0) error->all(FLERR,"Fix seaice/waves requires a global property rhoWater > 0");  
  
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // error checks on coarsegraining
  if(force->cg_active())
    error->cg(FLERR,this->style);

}

/* ---------------------------------------------------------------------- */

void FixSeaiceWaves::setup(int vflag)
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

void FixSeaiceWaves::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSeaiceWaves::post_force(int vflag)
{
  double **x = atom->x;
  int *mask = atom->mask;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;
  double *r = atom->radius;
  double **tilt = atom->tilt;
  double **torque = atom->torque;
  double *density = atom->density;

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  modify->clearstep_compute();
  modify->addstep_compute(update->ntimestep + 1);

    double locSLPX,locSLPY;
    double tnow,tmp;
    
    tnow = (update->ntimestep - update->beginstep) * update->dt;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (iregion >= 0 &&
            !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;

        foriginal[1] += torque[i][0];
        foriginal[2] += torque[i][1];
        foriginal[3] += torque[i][2];
        locSLPX = 0.0;
        locSLPY = 0.0;
        for (int k = 0; k < nw; k++) {
//          fprintf(screen,"a = %f,kx = %f\n",regwave[k][0],wavenum[k][0]);
          locSLPX += regwave[k][0]*cos(wavenum[k][0]*(x[i][0]+r[i])+wavenum[k][1]*x[i][1]-twopi/regwave[k][1]*tnow+regwave[k][3]);
          locSLPX -= regwave[k][0]*cos(wavenum[k][0]*(x[i][0]-r[i])+wavenum[k][1]*x[i][1]-twopi/regwave[k][1]*tnow+regwave[k][3]);
          locSLPY += regwave[k][0]*cos(wavenum[k][0]*x[i][0]+wavenum[k][1]*(x[i][1]+r[i])-twopi/regwave[k][1]*tnow+regwave[k][3]);
          locSLPY -= regwave[k][0]*cos(wavenum[k][0]*x[i][0]+wavenum[k][1]*(x[i][1]-r[i])-twopi/regwave[k][1]*tnow+regwave[k][3]);
        }
        locSLPX /= 2.0*r[i];
        locSLPY /= 2.0*r[i];
        tmp = PIG4*(rhowater+density[i])/2.0*r[i]*r[i]*r[i]*r[i];
        torque[i][0] += tmp*tan(atan(locSLPY)-tilt[i][0]);
        torque[i][1] += tmp*tan(atan(locSLPX)-tilt[i][1]);
//        fprintf(screen,"i = %d, slpx = %f, tiltx = %f, torquex = %f\n",i,locSLPY,tilt[i][0],torque[i][0]);
      }
}

/* ---------------------------------------------------------------------- */

void FixSeaiceWaves::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSeaiceWaves::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixSeaiceWaves::compute_scalar()
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

double FixSeaiceWaves::compute_vector(int n)
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

double FixSeaiceWaves::memory_usage()
{
  double bytes = atom->nmax*2 * sizeof(double);
  return bytes;
}
