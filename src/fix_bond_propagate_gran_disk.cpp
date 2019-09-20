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
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_bond_propagate_gran_disk.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include <list>

using namespace LAMMPS_NS;
using namespace FixConst;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixBondPropagateGranDisk::FixBondPropagateGranDisk(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    restart_global = 1;
    global_freq = 1;

    force_reneighbor = 1;
    next_reneighbor = -1;
}

/* ---------------------------------------------------------------------- */

FixBondPropagateGranDisk::~FixBondPropagateGranDisk()
{
}

/* ---------------------------------------------------------------------- */

int FixBondPropagateGranDisk::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= POST_INTEGRATE;  
  mask |= POST_INTEGRATE_RESPA;    
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondPropagateGranDisk::post_integrate()
{
  int n;
  int broken = 0;
  int broken_global = 0;
  int **bondlist = neighbor->bondlist;  
  int nbondlist = neighbor->nbondlist;


  for (n = 0; n < nbondlist; n++) {
    if(bondlist[n][3]) {
      broken = 1;
      break;
    }
  }
  MPI_Allreduce(&broken,&broken_global,1,MPI_INT,MPI_LOR,world);
  if (broken_global) next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixBondPropagateGranDisk::post_integrate_respa(int ilevel, int iloop)
{
 if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */
//NP  copy partner info from neighbor lists to atom arrays
//NP  so can be exchanged with atoms
//NP  IMPORTANT: newtonflag = 1 assumed for all history values
void FixBondPropagateGranDisk::pre_exchange()
{
  int i1,i2,n;
  bool found;
  int **bondlist = neighbor->bondlist;
  double **bondhistlist = neighbor->bondhistlist;
  int nbondlist = neighbor->nbondlist;

  int **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  double ***bond_hist = atom->bond_hist;// schickt eig. in. neighlist zu per atom arrays
  int n_bondhist = atom->n_bondhist;
  int nlocal = atom->nlocal;
  int *tag = atom->tag;

  int newton_bond = force->newton_bond;

/*
  //NP task 1
  //NP propagate bond contact history
  if(!n_bondhist) return;

  //fprintf(screen,"nbondlist_propagate %d\n",nbondlist);

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    int broken = bondlist[n][3];
    if(broken) continue;
    if (newton_bond || i1 < nlocal)
    {
        found = false;
        for(int k = 0; k < num_bond[i1]; k++) {
           if(bond_atom[i1][k] == tag[i2])
           {
	      found = true;
	      for(int q = 0; q < n_bondhist; q++) bond_hist[i1][k][q] = bondhistlist[n][q];
              break;
           }
        }
        if(!found) error->all(FLERR,"Failed to operate on granular bond history during copy");
    }

    if (newton_bond || i2 < nlocal)
    {
        found = false;
        for(int k = 0; k < num_bond[i2]; k++)
           if(bond_atom[i2][k] == tag[i1])
           {
	      found = true;
	      for(int q = 0; q < n_bondhist; q++) bond_hist[i2][k][q] = bondhistlist[n][q];
              break;
           }
        if(!found) error->all(FLERR,"Failed to operate on granular bond history during copy");
      
    }
  }
*/
  //NP task 2
  //NP remove broken bonds
  //NP should be done equally on all processors
  
  //NP P.F. create list for all BOND-IDs,erase all the broken bonds compress bondlist with new value list
  std::list<unsigned int> list_bond_id;
  std::list<unsigned int>::iterator it1;
  
//AH: There was a major bug here!!!
//AH  for (n=0; n<nbondlist; n++) if(!bondlist[n][3])list_bond_id.push_back(n);// load all unbroken bond ids into the list
  for (n=0; n<nbondlist; n++) if(bondlist[n][3])list_bond_id.push_back(n);// load all BROKEN  bond ids into the list  
  for (it1=list_bond_id.begin(); it1!=list_bond_id.end(); ++it1) 
  { 
    n = *it1;
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];

    
    int broken = bondlist[n][3];
    if(!broken) continue;
    
    //NP if the bond is broken, we remove it from
    //NP both atom data

    // delete bond from atom I if I stores it
    // atom J will also do this

    if (newton_bond || i1 < nlocal)
    {
        found = false;
        for(int k = 0; k < num_bond[i1]; k++)
           if(bond_atom[i1][k] == tag[i2])
           {
	      found = true;
	      remove_bond(i1,k,n);
	      break;
           }
        if(!found) error->all(FLERR,"Failed to operate on granular bond history during deletion1");
    }
    if (newton_bond || i2 < nlocal)
    {
        found = false;
        for(int k = 0; k < num_bond[i2]; k++)
           if(bond_atom[i2][k] == tag[i1])
           {
	     found = true;
	       remove_bond(i2,k,n);
              int nbondlist = neighbor->nbondlist;
		
		if(n<nbondlist)
		{
		  neighbor->nbondlist = (nbondlist-1); // delete one bond -> reduce nbondlist by one
		  for(int i = 0; i <= 3; i++) 		
		  neighbor->bondlist[n][i] = neighbor->bondlist[nbondlist-1][i]; // NP P.F. added also change in neighbor bondlist
		  break;
		}
		else if(n == nbondlist) 
		{
		  neighbor->nbondlist = (nbondlist-1);
		  break;
		}break;
           }
        if(!found) error->all(FLERR,"Failed to operate on granular bond history during deletion2");
    }
  }
}

inline void FixBondPropagateGranDisk::remove_bond(int ilocal,int ibond, int bondnumber) //NP P.F. added bondnumber
{
    /*NL*///fprintf(screen,"removing bond\n");
    /*NL*///error->one(FLERR,"romoving bond");
    int nbond = atom->num_bond[ilocal];
    atom->bond_atom[ilocal][ibond] = atom->bond_atom[ilocal][nbond-1];
    atom->bond_type[ilocal][ibond] = atom->bond_type[ilocal][nbond-1];
    for(int k = 0; k < atom->n_bondhist; k++) atom->bond_hist[ilocal][ibond][k] = atom->bond_hist[ilocal][nbond-1][k];
    atom->num_bond[ilocal]--;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixBondPropagateGranDisk::write_restart(FILE *fp)
{
  error->warning(FLERR,"Restart functionality not yet tested for granular bonds...");

  //NP write a dummy value
  int n = 0;
  double list[1];
  list[n++] = 1.;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }

  //NP write data to atom arrays where it can then be stored
  //NP can be done this way bc modify writes before avec
  pre_exchange();
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixBondPropagateGranDisk::restart(char *buf)
{
//  int n = 0;
  //double *list = (double *) buf;

//  double dummy = static_cast<int> (list[n++]);

  error->warning(FLERR,"Restart functionality not yet tested for granular bonds...");
}
