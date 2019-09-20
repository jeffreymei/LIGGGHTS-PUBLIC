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

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom_vec_disk_waves.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "domain_wedge.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecDiskWaves::AtomVecDiskWaves(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;

  comm_x_only = 1;
  comm_f_only = 1;
  size_forward = 3;  // see pack_comm
  size_reverse = 3;  // see pack_reverse
  size_border = 8;   // see pack_border
  size_velocity = 3; // see pack_comm_vel
  size_data_atom = 5;// number of entries for atoms in input data files
  size_data_vel = 4; // number of entries for vels in input data files
  xcol_data = 3;     // column where x is located in input data files

  atom->sphere_flag = 0;
  atom->disk_flag = 1;
  atom->tilt_flag = 1;
}

/* ---------------------------------------------------------------------- */

void AtomVecDiskWaves::init()
{
/*
  AtomVec::init();

  // set radvary if particle diameters are time-varying due to fix adapt

  radvary = 0;
  comm_x_only = 1;
  size_forward = 3;

  // set radvary if particle diameters are time-varying due to some fix
  for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->rad_mass_vary_flag) {
        radvary = 1;
        size_forward = 7; // see pack_comm
        comm_x_only = 1;
      }

  if(radvary) atom->radvary_flag = 1;
*/  
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecDiskWaves::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  tilt = memory->grow(atom->tilt,nmax,2,"atom:tilt");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecDiskWaves::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  tilt = atom->tilt;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecDiskWaves::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  tilt[j][0] = tilt[i][0];
  tilt[j][1] = tilt[i][1];  

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecDiskWaves::pack_comm(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0];
        buf[m++] = x[j][1];
        buf[m++] = x[j][2];
      }
    } else {
      if (domain->triclinic == 0) {
        dx = pbc[0]*domain->xprd;
        dy = pbc[1]*domain->yprd;
        dz = pbc[2]*domain->zprd;
      } else {
        dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
        dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
        dz = pbc[2]*domain->zprd;
      }
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
      }
    }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDiskWaves::pack_comm_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
//  if(dynamic_cast<DomainWedge*>(domain))
//    return pack_comm_vel_wedge(n,list,buf,pbc_flag,pbc);

  int i,j,m;

  double dx,dy,dz,dvx,dvy,dvz;

    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0];
        buf[m++] = x[j][1];
        buf[m++] = x[j][2];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
    } else {
      if (domain->triclinic == 0) {
        dx = pbc[0]*domain->xprd;
        dy = pbc[1]*domain->yprd;
        dz = pbc[2]*domain->zprd;
      } else {
        dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
        dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
        dz = pbc[2]*domain->zprd;
      }
      if (!deform_vremap) {
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = x[j][0] + dx;
          buf[m++] = x[j][1] + dy;
          buf[m++] = x[j][2] + dz;
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      } else {
        dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
        dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
        dvz = pbc[2]*h_rate[2];
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = x[j][0] + dx;
          buf[m++] = x[j][1] + dy;
          buf[m++] = x[j][2] + dz;
          if (mask[i] & deform_groupbit) {
            buf[m++] = v[j][0] + dvx;
            buf[m++] = v[j][1] + dvy;
            buf[m++] = v[j][2] + dvz;
          } else {
            buf[m++] = v[j][0];
            buf[m++] = v[j][1];
            buf[m++] = v[j][2];
          }
        }
      }
    }
  return m;
}

/* ---------------------------------------------------------------------- */
/*
int AtomVecDiskWaves::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = tilt[j][0]; 	//AH: ????
    buf[m++] = tilt[j][1]; 	//AH: ????   
  }
  return m;
}
*/
/* ---------------------------------------------------------------------- */

void AtomVecDiskWaves::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
    }
}

/* ---------------------------------------------------------------------- */

void AtomVecDiskWaves::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
      v[i][0] = buf[m++];
      v[i][1] = buf[m++];
      v[i][2] = buf[m++];
    }
}

/* ---------------------------------------------------------------------- */
/*
int AtomVecDiskWaves::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    tilt[i][0] = buf[m++];
    tilt[i][1] = buf[m++];
  }
  return m;
}
*/
/* ---------------------------------------------------------------------- */

int AtomVecDiskWaves::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;
  error->warning(FLERR,"in pack_reverse");
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDiskWaves::pack_reverse_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = tilt[i][0];  //AH: ????
    buf[m++] = tilt[i][1];  //AH: ????
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecDiskWaves::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecDiskWaves::unpack_reverse_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    tilt[j][0] += buf[m++];
    tilt[j][1] += buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDiskWaves::pack_border(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = tilt[j][0];
      buf[m++] = tilt[j][1];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = tilt[j][0];
      buf[m++] = tilt[j][1];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDiskWaves::pack_border_vel(int n, int *list, double *buf,
                                     int pbc_flag, int *pbc)
{
//  if(dynamic_cast<DomainWedge*>(domain))
//    return pack_border_vel_wedge(n,list,buf,pbc_flag,pbc);

  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = tilt[j][0];
      buf[m++] = tilt[j][1];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = tilt[j][0];
        buf[m++] = tilt[j][1];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = tilt[j][0];
        buf[m++] = tilt[j][1];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDiskWaves::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = tilt[j][0];
    buf[m++] = tilt[j][1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecDiskWaves::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (int) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    tilt[i][0] = buf[m++];
    tilt[i][1] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecDiskWaves::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (int) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    tilt[i][0] = buf[m++];
    tilt[i][1] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecDiskWaves::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    tilt[i][0] = buf[m++];
    tilt[i][1] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecDiskWaves::pack_exchange(int i, double *buf)
{

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  buf[m++] = tilt[i][0];
  buf[m++] = tilt[i][1];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDiskWaves::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (int) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (tagint) ubuf(buf[m++]).i;

  tilt[nlocal][0] = buf[m++];
  tilt[nlocal][1] = buf[m++];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecDiskWaves::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 18 * nlocal; 

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecDiskWaves::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = tilt[i][0];
  buf[m++] = tilt[i][1];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecDiskWaves::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (int) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (tagint) ubuf(buf[m++]).i;
  
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  tilt[nlocal][0] = buf[m++];
  tilt[nlocal][1] = buf[m++];

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecDiskWaves::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((tagint) IMGMAX << IMG2BITS) |
    ((tagint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  tilt[nlocal][0] = 0.0;
  tilt[nlocal][1] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecDiskWaves::data_atom(double *coord, tagint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = atoi(values[0]);
  if (tag[nlocal] <= 0)
    error->one(FLERR,"Invalid atom ID in Atoms section of data file");

  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  tilt[nlocal][0] = 0.0;
  tilt[nlocal][1] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecDiskWaves::pack_data(double **buf)
{
/*
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = 2.0*radius[i];
    buf[i][3] = thickness[i];
    if (radius[i] == 0.0) buf[i][4] = rmass[i];
    else
      buf[i][4] = rmass[i] / (MY_PI*radius[i]*radius[i]*thickness[i]);
    buf[i][5] = x[i][0];
    buf[i][6] = x[i][1];
    buf[i][7] = x[i][2];
    buf[i][8] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][10] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
 */
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecDiskWaves::pack_data(double **buf,int tag_offset)
{
/*
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]+tag_offset).d; 
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = 2.0*radius[i];
    buf[i][3] = thickness[i];
    if (radius[i] == 0.0) buf[i][4] = rmass[i];
    else
      buf[i][4] = rmass[i] / (MY_PI*radius[i]*radius[i]*thickness[i]);
    buf[i][5] = x[i][0];
    buf[i][6] = x[i][1];
    buf[i][7] = x[i][2];
    buf[i][8] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][10] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
*/  
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecDiskWaves::pack_data_hybrid(int i, double *buf)
{
  buf[0] = tilt[i][0];
  buf[1] = tilt[i][1];
  return 2;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecDiskWaves::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (int) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],buf[i][4],
            buf[i][5],buf[i][6],buf[i][7],
            (int) ubuf(buf[i][8]).i,(int) ubuf(buf[i][9]).i,
            (int) ubuf(buf[i][10]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecDiskWaves::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 2;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */
/*
void AtomVecDiskWaves::pack_vel(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
  }
}
*/
/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */
/*
void AtomVecDiskWaves::pack_vel(double **buf,int tag_offset) 
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]+tag_offset).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
  }
}
*/
/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */
/*
int AtomVecDiskWaves::pack_vel_hybrid(int i, double *buf)
{
  return 3;
}
*/
/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */
/*
void AtomVecDiskWaves::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (int) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6]);
}
*/
/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */
/*
int AtomVecDiskWaves::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}
*/
/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecDiskWaves::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("tilt")) bytes += memory->usage(tilt,nmax,2);

  return bytes;
}
