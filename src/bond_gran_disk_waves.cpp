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

#include "math.h"
#include "stdlib.h"
#include "bond_gran_disk_waves.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "fix_property_atom.h"
#include "error.h"
#include "update.h"
#include "vector_liggghts.h"
#include "random_park.h"  //AH

using namespace LAMMPS_NS;

/*NP
large TODO list for granular bonds:  (could be a diploma thesis?)

+ need a better dissipative formulation than the hardcoded
  'dissipate' value which produces plastic deformation
  need some vel-dependant damping
+ need to carefully debug and validate this bond style
  valiation against fix rigid
+ check whether run this bond style w/ or w/o gran pair style active,
  (neigh_modify command)
+ need to store bond radii per particle, not per type
+ parallel debugging and testing not yet done
+ need evtally implemetation
*/

enum{
     BREAKSTYLE_STRESS,
    };

/* ---------------------------------------------------------------------- */

BondGranDiskWaves::BondGranDiskWaves(LAMMPS *lmp) : Bond(lmp)
{
    // 
    n_granhistory(11); 
    j0 = 0;
    j1 = 1;
    j3 = 2;
    j4 = 3;
    j5 = 4;
    j6 = 5;
    j7 = 6;
    j9 = 7;
    j10 = 8;
    j11 = 9;
    j13 = 10;  //AH: additional place for variable bond thickness
    /*	NP
    /* number of entries in bondhistlist. bondhistlist[number of bond][number of value (from 0 to number given here)]
    /* so with this number you can modify how many pieces of information you savae with every bond
    /* following dependencies and processes for saving,copying,growing the bondhistlist: */
     
    /* NP
    /* gibt groesse der gespeicherten werte  pro bond wieder 
    /* neighbor.cpp:       memory->create(bondhistlist,maxbond,atom->n_bondhist,"neigh:bondhistlist");
    /* neigh_bond.cpp:     memory->grow(bondhistlist,maxbond,atom->n_bondhist,"neighbor:bondhistlist");
    /* bond.cpp: void Bond::n_granhistory(int nhist) {ngranhistory = nhist;     atom->n_bondhist = ngranhistory; if(){FLERR}}
    /* atom_vec_bond_gran.cpp:  memory->grow(atom->bond_hist,nmax,atom->bond_per_atom,atom->n_bondhist,"atom:bond_hist");

    /* 
     */
    if(!atom->style_match("bond/gran/disk"))
      error->all(FLERR,"A granular bond/disk style can only be used together with atom style bond/gran");
}

/* ---------------------------------------------------------------------- */

BondGranDiskWaves::~BondGranDiskWaves()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(rb);
    memory->destroy(lb);
    memory->destroy(thb);
    memory->destroy(thstdb);
    memory->destroy(Ec);
    memory->destroy(kn2kt);
    memory->destroy(sigman_break_c);
    memory->destroy(sigman_break_t);
    memory->destroy(tau_break);
  }
}

/* ---------------------------------------------------------------------- */

void  BondGranDiskWaves::init_style()
{
//    if(breakmode == BREAKSTYLE_STRESS_TEMP)
//       fix_Temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",1,0,"bond gran waves"));
}

/* ---------------------------------------------------------------------- */

void BondGranDiskWaves::compute(int eflag, int vflag)
{

  double rsq,r,rinv,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,vtr1,vtr2,vtr3,vrel,tor1,tor2,tor3;
  double wnnr,wn1,wn2,wn3,wt1,wt2,wt3;
  double fs1,fs2,fs3;

  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double dnforce[2],dtforce[3];
  double dntorque[2],dttorque[3];
  double rot;
  double A,Ir,Ih,J;          //AH: Ir, Ih is for bending aroundvert and horiz axis, repsectively, J for twisting
  double thbloc,rbloc,lbloc; //AH: "local" bond thickness (vertical), half-width (horiz.) and length
  double knloc,ktloc;        //AH: kn=Ec/(r1+r2), kt=kn/kn2kt

  ebond = fbond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *radius = atom->radius;
  double *thickness = atom->thickness;
  double **torque = atom->torque;
  int *tag = atom->tag; // tag of atom is their ID number
  double **omega = atom->omega;
  int **bondlist = neighbor->bondlist;
  double **bondhistlist = neighbor->bondhistlist;
  //AH: new code, transferred from fix_bond_propagate_gran_disk:
  double ***bond_hist = atom->bond_hist; 
  int n_bondhist = atom->n_bondhist;
  bool found;
  //AH  

  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  double dt = update->dt;
  
  int seed = 100;
  RanPark *random = new RanPark(lmp,seed); //AH

  for (n = 0; n < nbondlist; n++) {
    if (bondlist[n][3]) continue;
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    /*NL*/ //fprintf(screen,"ts %d: handling id %d and %d\n",update->ntimestep,tag[i1],tag[i2]);
    type = bondlist[n][2];
    
    if(!bondlist[n][3] & bondhistlist[n][j13]==0) {
      random->reset(seed,x[i1]);
      bondhistlist[n][j13] = thb[type] - thstdb[type]/2.0 + random->uniform()*thstdb[type];
    }

    thbloc = bondhistlist[n][j13];		 // thbloc = MIN(thickness[i1],thickness[i2]);    //AH
    lbloc = lb[type]*(radius[i1]+radius[i2]);    //AH
    rbloc = rb[type]*MIN(radius[i1],radius[i2]); //AH
    A = 2.0 * rbloc * thbloc;                    //AH (orig: A = M_PI * rb[type] * rb[type];)
    Ir = A * rbloc * rbloc / 3.0;                //AH
    Ih = A * thbloc * thbloc / 12.0;             //AH
    J = Ir + Ih;                                 //AH (orig: J = A * 0.5 * rb[type] * rb[type];)

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = 0.0;
    domain->minimum_image(delx,dely,delz);

    rsq = delx*delx + dely*dely;
    rsqinv = 1./rsq;
    r = sqrt(rsq);
    rinv = 1./r;

    // relative translational velocity

        vr1 = v[i1][0] - v[i2][0];
        vr2 = v[i1][1] - v[i2][1];

        // normal component

        vnnr = vr1*delx + vr2*dely;
        vn1 = delx*vnnr * rsqinv;
        vn2 = dely*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;

        // relative rotational velocity for shear

	wr1 = (radius[i1]*omega[i1][0] + radius[i2]*omega[i2][0]) * rinv;
	wr2 = (radius[i1]*omega[i1][1] + radius[i2]*omega[i2][1]) * rinv;
        wr3 = (radius[i1]*omega[i1][2] + radius[i2]*omega[i2][2]) * rinv;

        // relative velocities for shear

        vtr1 = vt1 + dely*wr3;
        vtr2 = vt2 - delx*wr3;
        vtr3 = -(dely*wr1 - delx*wr2);

        // relative rotational velocity for torsion and bending

        wr1 = (radius[i1]*omega[i1][0] - radius[i2]*omega[i2][0]) * rinv;
        wr2 = (radius[i1]*omega[i1][1] - radius[i2]*omega[i2][1]) * rinv;
        wr3 = (radius[i1]*omega[i1][2] - radius[i2]*omega[i2][2]) * rinv;

        // normal component

        wnnr = wr1*delx + wr2*dely;
        wn1 = delx*wnnr * rsqinv;
        wn2 = dely*wnnr * rsqinv;

        // tangential component

        wt1 = wr1 - wn1;
        wt2 = wr2 - wn2;
        wt3 = wr3;
    
    // calc normal and shear stiffness: //AH
    knloc = Ec[type] / lbloc;	//AH
    ktloc = knloc/kn2kt[type];	//AH

    // calc change in normal forces
    dnforce[0] = - vn1 * knloc * A * dt;
    dnforce[1] = - vn2 * knloc * A * dt;

    // calc change in shear forces
    dtforce[0] = - vtr1 * ktloc * A * dt;
    dtforce[1] = - vtr2 * ktloc * A * dt;
    dtforce[2] = - vtr3 * ktloc * A * dt;

    // calc change in normal torque
    dntorque[0] = - wn1 * ktloc * J * dt;
    dntorque[1] = - wn2 * ktloc * J * dt;

    // calc change in tang torque
    dttorque[0] = - wt1 * knloc * Ih * dt;
    dttorque[1] = - wt2 * knloc * Ih * dt;
    dttorque[2] = - wt3 * knloc * Ir * dt;

    // rotate forces

    //rotate normal force
    rot = bondhistlist[n][j0]*delx + bondhistlist[n][j1]*dely;
    rot *= rsqinv;
    bondhistlist[n][j0] = rot*delx;
    bondhistlist[n][j1] = rot*dely;

    //rotate tangential force
    rot = bondhistlist[n][j3]*delx + bondhistlist[n][j4]*dely;
    rot *= rsqinv;
    bondhistlist[n][j3] -= rot*delx;
    bondhistlist[n][j4] -= rot*dely;

    //rotate normal torque
    rot = bondhistlist[n][j6]*delx + bondhistlist[n][j7]*dely;
    rot *= rsqinv;
    bondhistlist[n][j6] = rot*delx;
    bondhistlist[n][j7] = rot*dely;

    //rotate tangential torque
    rot = bondhistlist[n][j9]*delx + bondhistlist[n][j10]*dely;
    rot *= rsqinv;
    bondhistlist[n][ j9] -= rot*delx;
    bondhistlist[n][j10] -= rot*dely;

    //increment normal and tangential force and torque
    double dissipate = 0.995;
    bondhistlist[n][j0] = dissipate * bondhistlist[n][j0] + dnforce[0];
    bondhistlist[n][j1] = dissipate * bondhistlist[n][j1] + dnforce[1];
    bondhistlist[n][j3] = dissipate * bondhistlist[n][j3] + dtforce[0];
    bondhistlist[n][j4] = dissipate * bondhistlist[n][j4] + dtforce[1];
    bondhistlist[n][j5] = dissipate * bondhistlist[n][j5] + dtforce[2];
    bondhistlist[n][j6] = dissipate * bondhistlist[n][j6] + dntorque[0];
    bondhistlist[n][j7] = dissipate * bondhistlist[n][j7] + dntorque[1];
    bondhistlist[n][j9] = dissipate * bondhistlist[n][j9] + dttorque[0];
    bondhistlist[n][j10] = dissipate * bondhistlist[n][j10] + dttorque[1];
    bondhistlist[n][j11] = dissipate * bondhistlist[n][j11] + dttorque[2];

    tor1 = - rinv * dely*bondhistlist[n][j5];
    tor2 =   rinv * delx*bondhistlist[n][j5];
    tor3 = - rinv * (delx*bondhistlist[n][j4] - dely*bondhistlist[n][j3]);

//    if(breakmode == BREAKSTYLE_STRESS)
//    {
        double nforce = sqrt(bondhistlist[n][j0]*bondhistlist[n][j0]+bondhistlist[n][j1]*bondhistlist[n][j1]);
        rot = bondhistlist[n][j0]*delx + bondhistlist[n][j1]*dely; // rot<0 => tension; rot>0 => compression
        if (rot < 0.0) nforce = -1.0*nforce; // negative force means tension
        double tforce_mag = vectorMag3D(&bondhistlist[n][j3]);
        double ntorque_mag = sqrt(bondhistlist[n][j6]*bondhistlist[n][j6]+bondhistlist[n][j7]*bondhistlist[n][j7]);
        double ttorque_mag01 = sqrt(bondhistlist[n][j9]*bondhistlist[n][j9]+bondhistlist[n][j10]*bondhistlist[n][j10])/Ih*thbloc;
        double ttorque_mag2 = abs(bondhistlist[n][j11])/Ir*rbloc;

	// this version treats stress acting on the horizontal and vertical "walls" of the bond separately:
        //bool ntensstress = (sigman_break_t[type] < (-nforce/A + ttorque_mag01)) || (sigman_break_t[type] < (-nforce/A + ttorque_mag2)); 
        //bool ncompstress = (sigman_break_c[type] < ( nforce/A + ttorque_mag01)) || (sigman_break_c[type] < ( nforce/A + ttorque_mag2)); 
        // this version assumes that the maximum stress occurs along the edges ("corners") of the bond:
        bool ntensstress = sigman_break_t[type] < (-nforce/A + ttorque_mag01 + ttorque_mag2); 
        bool ncompstress = sigman_break_c[type] < ( nforce/A + ttorque_mag01 + ttorque_mag2); 
        bool tstress     = tau_break[type]      < (tforce_mag/A + ntorque_mag/J*rbloc); //AH

        if(ntensstress || ncompstress || tstress)
        {
//            bondhistlist[n][j13] *= 0.99;
//            fprintf(screen,"n = %d, thickness = %f\n",n,bondhistlist[n][j13]);
//            if (bondhistlist[n][j13]<=0.1)
//            {
              bondlist[n][3] = 1;
              int m;
              int histmax = atom->n_bondhist;
              for (m = 0; m < histmax; m++)		//AH
                bondhistlist[n][m] = 0.0;		//AH
              for (m = 0; m < atom->num_bond[i1]; m++)
                if (atom->bond_atom[i1][m] == atom->tag[i2])
                  atom->bond_type[i1][m] = 0;
              if (i2 < atom->nlocal)
                for (m = 0; m < atom->num_bond[i2]; m++)
                  if (atom->bond_atom[i2][m] == atom->tag[i1])
                    atom->bond_type[i2][m] = 0;
//            }
        }
//    }

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += (bondhistlist[n][j0] + bondhistlist[n][j3]);
      f[i1][1] += (bondhistlist[n][j1] + bondhistlist[n][j4]);
      torque[i1][0] += radius[i1]*tor1 + (bondhistlist[n][j6] + bondhistlist[n][j9]);
      torque[i1][1] += radius[i1]*tor2 + (bondhistlist[n][j7] + bondhistlist[n][j10]);
      torque[i1][2] += radius[i1]*tor3 + bondhistlist[n][j11];
      //AH: new code, transferred from fix_bond_propagate_gran_disk:
      found = false;
      for(int k = 0; k < atom->num_bond[i1]; k++) {
        if (atom->bond_atom[i1][k] == atom->tag[i2])
          {
            found = true;
            for (int q = 0; q < n_bondhist; q++) bond_hist[i1][k][q] = bondhistlist[n][q];
            break;
          }
      }
      if (!found) error->all(FLERR,"Failed to operate on bond history during copy");      
      //      
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= (bondhistlist[n][j0] + bondhistlist[n][j3]);
      f[i2][1] -= (bondhistlist[n][j1] + bondhistlist[n][j4]);
      torque[i2][0] += radius[i2]*tor1 - (bondhistlist[n][j6] + bondhistlist[n][ j9]);
      torque[i2][1] += radius[i2]*tor2 - (bondhistlist[n][j7] + bondhistlist[n][j10]);
      torque[i2][2] += radius[i2]*tor3 - bondhistlist[n][j11];
      //AH: new code, transferred from fix_bond_propagate_gran_disk:
      found = false;
      for(int k = 0; k < atom->num_bond[i2]; k++) {
        if (atom->bond_atom[i2][k] == atom->tag[i1])
          {
            found = true;
            for (int q = 0; q < n_bondhist; q++) bond_hist[i2][k][q] = bondhistlist[n][q];
            break;
          }
      }
      if (!found) error->all(FLERR,"Failed to operate on bond history during copy");      
      //            
    }

    fbond = nforce + tforce_mag;
    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
  delete random;
}

/* ---------------------------------------------------------------------- */

void BondGranDiskWaves::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(rb,n+1,"bond:rb");
  memory->create(thb,n+1,"bond:thb");
  memory->create(thstdb,n+1,"bond:thstdb");  
  memory->create(lb,n+1,"bond:lb");
  memory->create(Ec,n+1,"bond:Ec");
  memory->create(kn2kt,n+1,"bond:kn2kt");
  memory->create(sigman_break_c,(n+1),"bond:sigman_break_c");
  memory->create(sigman_break_t,(n+1),"bond:sigman_break_t");
  memory->create(tau_break,(n+1),"bond:tau_break");

  memory->create(setflag,(n+1),"bond:setflag");
  for (int i = 1; i <= n; i++)
    setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondGranDiskWaves::coeff(int narg, char **arg)
{
  if(narg < 7)  error->all(FLERR,"Incorrect args for bond coefficients");

  double rb_one = force->numeric(FLERR,arg[1]);
  double thb_one = force->numeric(FLERR,arg[2]);
  double thstdb_one = force->numeric(FLERR,arg[3]);  
  double lb_one = force->numeric(FLERR,arg[4]);
  double Ec_one = force->numeric(FLERR,arg[5]);
  double kn2kt_one = force->numeric(FLERR,arg[6]);

  if(force->numeric(FLERR,arg[7]) == 1. )
  {
      breakmode = BREAKSTYLE_STRESS;
      if (narg != 11) error->all(FLERR,"Incorrect args for bond coefficients");
  }
  else  error->all(FLERR,"Only breakmode = 1 is allowed in this version");

  if (!allocated) allocate();

  double sigman_break_c_one,sigman_break_t_one,tau_break_one;

  sigman_break_c_one = force->numeric(FLERR,arg[8]);
  sigman_break_t_one = force->numeric(FLERR,arg[9]);
  tau_break_one = force->numeric(FLERR,arg[10]);

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rb[i]    = rb_one;
    thb[i]   = thb_one;
    thstdb[i]= thstdb_one;
    lb[i]    = lb_one;
    Ec[i]    = Ec_one;
    kn2kt[i] = kn2kt_one;
    sigman_break_c[i] = sigman_break_c_one;
    sigman_break_t[i] = sigman_break_t_one;
    tau_break[i] = tau_break_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients - or the bonds are not initialized in create_atoms");

}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondGranDiskWaves::equilibrium_distance(int i)
{
  //NP ATTENTION: this is _not_ correct - and has implications on fix shake, pair_lj_cut_coul_long and pppm
  //NP it is not possible to define a general equilibrium distance for this bond model
  //NP as rotational degree of freedom is present
  return 0.;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondGranDiskWaves::write_restart(FILE *fp)
{
  fwrite(&rb[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&thb[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&thstdb[1],sizeof(double),atom->nbondtypes,fp);  
  fwrite(&lb[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Ec[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&kn2kt[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondGranDiskWaves::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&rb[1],sizeof(double),atom->nbondtypes,fp);
    fread(&thb[1],sizeof(double),atom->nbondtypes,fp);
    fread(&thstdb[1],sizeof(double),atom->nbondtypes,fp);    
    fread(&lb[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Ec[1],sizeof(double),atom->nbondtypes,fp);
    fread(&kn2kt[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&rb[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&thb[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&thstdb[1],atom->nbondtypes,MPI_DOUBLE,0,world);  
  MPI_Bcast(&lb[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Ec[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kn2kt[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double BondGranDiskWaves::single(int type, double rsum, int i, int j, double &fval)
{ 
    return rsum*lb[type];
}