#pragma once

#include "../../interfaces.hpp"
#include "../../kinematics/kinematics.hpp"

namespace constitutive_models {

class NeoHookean1D : public MatLaw<1> {
  public:
    double mu;

    NeoHookean1D() {};
    NeoHookean1D(double mu);
    ~NeoHookean1D() {};
    void set_pars(double mu);
    double stress(const kinematics::kinematics<1> &kin, double stress[]);
    double stress(double args);
};

} // namespace constitutive_models
