/* ----------------------------------------------------------------------

   This file is part of the sea ice toolbox for LIGGGHTS.
   It was created based on the analogous code of LAMMPS and LIGGGHTS.
   See the documentation of the toolbox for details.

   Author: Agnieszka Herman, University of Gdansk, Poland
   http://herman.ocean.ug.edu.pl/LIGGGHTSseaice.html
   agnieszka.herman@ug.edu.pl

------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */
#ifdef NORMAL_MODEL
NORMAL_MODEL(HERTZ_DISK,hertz/disk,5)
#else
#ifndef NORMAL_MODEL_HERTZ_DISK_H_
#define NORMAL_MODEL_HERTZ_DISK_H_
#include "global_properties.h"
#include "math.h"
#include "math_const.h"

using namespace MathConst;

namespace LIGGGHTS {

namespace ContactModels
{
  template<>
  class NormalModel<HERTZ_DISK> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

    NormalModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp),
      Yeff(NULL),
      Geff(NULL),
      betaeff(NULL),
      limitForce(false),
      displayedSettings(false)
    {
      
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
    }

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff,"model hertz/disk");
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff,"model hertz/disk");
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff,"model hertz/disk");

      registry.connect("Yeff", Yeff,"model hertz/disk");
      registry.connect("Geff", Geff,"model hertz/disk");
      registry.connect("betaeff", betaeff,"model hertz/disk");
    }

    // effective exponent for stress-strain relationship
    
    inline double stressStrainExponent()
    {
      return 1.5;
    }

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;
      double ri = cdata.radi;
      double rj = cdata.radj;
      double reff=cdata.is_wall ? cdata.radi : (ri*rj/(ri+rj));
      double meff=cdata.meff;
      
      double hmin = cdata.hmin;
      const double p1 = 0.9117;
      const double p2 = 0.2722;
      const double p3 = 0.003324;
      const double q1 = 1.524;
      const double q2 = 0.03159;
      
      double x = cdata.deltan*reff/(2.*hmin*hmin);
      double fx = (p1*x*x+p2*x+p3)/(x*x+q1*x+q2);
      double kn = MY_PI*Yeff[itype][jtype]*hmin*fx;
      double kt = 6.*Geff[itype][jtype]/Yeff[itype][jtype]*kn;
      double tmp = sqrt(5.*kn*meff);
      double gamman=-betaeff[itype][jtype]*tmp;
      double gammat=-2.*betaeff[itype][jtype]*tmp*sqrt(Geff[itype][jtype]/Yeff[itype][jtype]);

      if (!tangential_damping) gammat = 0.0;

      if(!displayedSettings)
      {
        displayedSettings = true;

        /*
        if(limitForce)
            if(0 == comm->me) fprintf(screen," NormalModel<HERTZ_DISK>: will limit normal force.\n");
        */
      }
      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double Fn_damping = -gamman*cdata.vn;
      const double Fn_contact = kn*(cdata.radsum-cdata.r);
      double Fn = Fn_damping + Fn_contact;

      //limit force to avoid the artefact of negative repulsion force
      if(limitForce && (Fn<0.0) )
      {
          Fn = 0.0;
      }

      cdata.Fn = Fn;
      cdata.kn = kn;
      cdata.kt = kt;
      cdata.gamman = gamman;
      cdata.gammat = gammat;

      // apply normal force
      if(cdata.is_wall) {
        const double Fn_ = Fn * cdata.area_ratio;
        i_forces.delta_F[0] = Fn_ * cdata.en[0];
        i_forces.delta_F[1] = Fn_ * cdata.en[1];
        i_forces.delta_F[2] = Fn_ * cdata.en[2];
      } else {
        i_forces.delta_F[0] = cdata.Fn * cdata.en[0];
        i_forces.delta_F[1] = cdata.Fn * cdata.en[1];
        i_forces.delta_F[2] = cdata.Fn * cdata.en[2];

        j_forces.delta_F[0] = -i_forces.delta_F[0];
        j_forces.delta_F[1] = -i_forces.delta_F[1];
        j_forces.delta_F[2] = -i_forces.delta_F[2];
      }
    }

    void noCollision(ContactData&, ForceData&, ForceData&){}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
    double ** Yeff;
    double ** Geff;
    double ** betaeff;

    bool tangential_damping;
    bool limitForce;
    bool displayedSettings;
  };

}

}
#endif
#endif
