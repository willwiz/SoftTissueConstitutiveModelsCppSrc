#pragma once

#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

class NeoHookean : public MatLaw3D<9> {
  public:
    double mu;

    NeoHookean(){};
    NeoHookean(double mu);
    ~NeoHookean(){};
    void set_pars(double mu);
    void stress(const kinematics::kinematics<9> &kin, double stress[]);
};

} // namespace constitutive_models
