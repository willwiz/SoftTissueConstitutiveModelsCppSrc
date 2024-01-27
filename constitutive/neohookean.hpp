#pragma once

#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

class NeoHookean : public MatLaw<4> {
  public:
    double mu;

    NeoHookean(){};
    NeoHookean(double mu);
    ~NeoHookean(){};
    void set_pars(double mu);
    double stress(const kinematics::kinematics<4> &kin, double stress[]);
    void stress(double args[], double stress[]);
};

} // namespace constitutive_models
