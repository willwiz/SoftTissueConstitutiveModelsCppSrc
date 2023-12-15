#pragma once

#include "../kinematics/kinematics.hpp"
#include "../interfaces.hpp"

namespace constitutive_models {

  class TransverseHog1D: public MatLaw<1>
  {
  public:
    double k1, k2;
    double fiber, iso;
    double E1;
    double E2;

    TransverseHog1D();
    ~TransverseHog1D();
    TransverseHog1D(double k1, double k2, double kappa);
    TransverseHog1D(double k1, double k2, double kappa, double Cmax);

    void set_pars(double kappa);
    void set_pars(double kappa, double Cmax);
    double stress(const kinematics::kinematics<1> &kin, double stress[]);
    double stress(double args);
  };

}