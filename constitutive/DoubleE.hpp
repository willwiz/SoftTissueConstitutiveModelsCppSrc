#pragma once

#include "../kinematics/kinematics.hpp"
#include "../interfaces.hpp"

namespace constitutive_models {
  class DoubleE2D: public MatLaw<4>
  {
    public:
      double b1, b2;
      double k1, k2, k3;
      double mxm[4], nxn[4], mxn[4];

      DoubleE2D();
      DoubleE2D(double b1, double b2, double k1, double k2, double k3, double theta);
      ~DoubleE2D();

      void set_pars(double b1, double b2, double k1, double k2, double k3, double theta);
      double stress(const kinematics::kinematics<4> &kin, double stress[4]);
      void stress(double args[], double stress[4]);
  };

}

