#pragma once

#include "../kinematics/kinematics.hpp"
#include "../interfaces.hpp"

namespace constitutive_models
{
  class LinearFiber: public MatLaw<4>
  {

  public:
    double k;
    double m4[4];
    double m6[4];
    LinearFiber();
    ~LinearFiber();
    LinearFiber(double mu, double theta, double alpha, double beta);

    void set_pars(double mu, double theta, double alpha, double beta);
    double stress(const kinematics::kinematics<4> &kin, double stress[4]);
    void stress(double args[4], double stress[4]);
  };

} // namespace constitutive_models
