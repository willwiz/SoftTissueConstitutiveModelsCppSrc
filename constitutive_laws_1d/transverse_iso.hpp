#pragma once

#include "../kinematics/kinematics.hpp"
#include "../interfaces.hpp"

namespace constitutive_models
{

  class TransverseIso1D : public MatLaw<1>
  {
  public:
    double k1, k2;
    double fiber, iso;

    TransverseIso1D();
    ~TransverseIso1D();
    TransverseIso1D(double k1, double k2, double kappa);
    TransverseIso1D(double k1, double k2, double kappa, double Cmax);

    void set_pars(double kappa);
    double stress(const kinematics::kinematics<1> &kin, double stress[]);
    double stress(double args);
  };

}