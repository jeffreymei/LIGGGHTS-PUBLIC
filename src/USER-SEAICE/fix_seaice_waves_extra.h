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

/* ----------------------------------------------------------------------
   Contributing author: Agnieszka Herman, IOUG, Poland
   
   Functions in this file are based on the code of the SWAN wave model
   (version 41.01, file swanser.ftn: subroutine SSHAPE, functions GAMMA 
   and GAMMALN)
------------------------------------------------------------------------- */


#ifndef LMP_FIX_SEAICE_WAVES_EXTRA_H
#define LMP_FIX_SEAICE_WAVES_EXTRA_H

#include "math.h"
#include "random_park.h"  
#include "math_const.h"
#include "pointers.h"

#define SIGMA1 0.07
#define SIGMA2 0.09
#define BETA  -1.25 // 5/4

using namespace MathConst;

namespace SeaiceWavesExtra {

  inline void jonswapToRegWave(int nw, double jonspar[], double **regwave, LAMMPS_NS::RanPark *random);
  inline double gamma(double xx);
  inline double gammaln(double xx);

}

/* ----------------------------------------------------------------------
express a 2D jonswap spectrum as a superposition of regular waves
------------------------------------------------------------------------- */

void SeaiceWavesExtra::jonswapToRegWave(int nw, double jonspar[], double **regwave, LAMMPS_NS::RanPark *random)
{
  int nf = (int)jonspar[5]; // no. of frequencies
  int nd = (int)jonspar[8]; // no. of directions
  int i;

  // logarithmically spaced frequency vector:
  double *f = new double[nf];
  double *df = new double[nf];
  double dlogf = (log10(jonspar[7])-log10(jonspar[6]))/((double)nf-1.0);
  double tmp;
  f[0] = jonspar[6];
  for (i = 1; i < (nf-1); i++) {
    tmp = log10(jonspar[6]) + (double)i*dlogf;
    f[i] = pow(10.0,tmp);
  }
  f[nf-1] = jonspar[7];
  df[0] = f[1]-f[0];
  for (i = 1; i < (nf-1); i++) 
    df[i] = (f[i+1]-f[i-1])/2.;
  df[nf-1] = f[nf-1]-f[nf-2];

  // uniformly spaced direction vector:
  double *theta = new double[nd];
  double dth = (jonspar[10]-jonspar[9])/(jonspar[8]-1);
  theta[0] = jonspar[9];
  for (i = 1; i < (nd-1); i++) {
    theta[i] = theta[i-1] + dth;
  }
  theta[nd-1] = jonspar[10];

  // 1D frequency spectrum, E1D:
  double fp = 1.0/jonspar[1]; // peak frequency
  double *E1D = new double[nf];
  double fp4 = pow(fp,4);
  double salpha = jonspar[0]*jonspar[0] * fp4 / ((0.06533*pow(jonspar[3],0.8015)+0.13467)*16.);
  double cpshap = BETA*fp4;
  double sigma,a,syf,ra;
  for (i = 0; i < nf; i++) {
    if (f[i]<=fp) sigma = SIGMA1;
    else sigma = SIGMA2;
    tmp = (f[i]/fp-1.0)/sigma;
    a = tmp*tmp/2.0;
    syf = 1.0;
    if (a<=10.0) syf = pow(jonspar[3],exp(-a));
    tmp = pow(f[i],4);
    ra = 0.0;
    if ((cpshap/tmp)<=10.0)
      ra = salpha/(f[i]*tmp)*exp(cpshap/tmp/f[i]);
    E1D[i] = ra*syf/(MY_2PI*f[i]);
    printf("%f %f\n",f[i],E1D[i]);
  }

  // 2D frequency-direction spectrum, E2D:
  double *E2D = new double[nd*nf];
  double adir = MY_PI*jonspar[2]/180.0;
  double dspr,ms,ctot,acos,cdir;
  int j,ij = 0;
// this code would be valid for directional spreading expressed in terms of directional standard deviation:
//  if (dshapl==1) {
//    dspr = MY_PI * jonspar[4]/180.0;
//    ms = MAX(pow(dspr,-2)-2.0,1.0);
//  } else
// this code is valid for directional spreading expressed with the power m in cos^m(theta-theta_peak):
  ms = jonspar[4];
  if (ms<12.0) {
    ctot = pow(2.0,ms) * pow(SeaiceWavesExtra::gamma(0.5*ms+1.0),2) / (MY_PI*SeaiceWavesExtra::gamma(ms+1.0));
  } else
    ctot = sqrt(ms/MY_2PI) / (1.0-0.25/ms);
  for (i = 0; i < nd; i++) {
    acos = cos(MY_PI*theta[i]/180.-adir);
    if (acos>0.0) {
      cdir = ctot * MAX(pow(acos,ms),1.E-10);
    } else 
      cdir = 0.0;
    for (j = 0; j < nf; j++) {
      E2D[ij] = cdir * E1D[j];
      ij++;
    }
    printf("%f %f\n",theta[i],cdir);
  }

  // find the characteristics of the regular waves:
  ij = 0;
  for (i = 0; i < nd; i++)
    for (j = 0; j < nf; j++) {
      regwave[ij][0] = sqrt(2./9.81/1025*E2D[ij]*dth*MY_PI/180.*MY_2PI*df[j]); //?????????
      regwave[ij][1] = 1.0/f[j];
      regwave[ij][2] = theta[i]; 
      regwave[ij][3] = random->uniform()*360.0 - 180.0;
      ij++;
    }
}

/* ----------------------------------------------------------------------
Compute the transcendental function Gamma
------------------------------------------------------------------------- */

double SeaiceWavesExtra::gamma(double xx) 
{ 
  double abig = 30.0;
  double yy = gammaln(xx);
  if (yy>abig)  yy = abig;
  if (yy<-abig) yy = -abig;
  return exp(yy);
}

double SeaiceWavesExtra::gammaln(double xx) 
{
  const double cof[6] = {76.18009173,-86.50532033,24.01409822,-1.231739516,.120858003e-2,-.536382e-5};
  double x,tmp,ser,res;
  int j;
  x = xx - 1.0;
  tmp = x + 5.5;
  tmp = (x+0.5)*log(tmp)-tmp;
  ser = 1.0;
  for (j = 0; j<6; j++) {
    x = x + 1.0;
    ser = ser + cof[j]/x;
  }
  res = tmp + log(2.50662827465*ser);
  return res;
}

#endif

