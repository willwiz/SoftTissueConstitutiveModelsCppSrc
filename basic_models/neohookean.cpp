#include "../kinematics/tensor_algebra.hpp"
#include "neohookean.hpp"
#include <cmath>

namespace constitutive_models {

/*----------------------------------------------------------------------
 |  Standard Neohookean, no much to say
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

NeoHookean::NeoHookean(double mu) {
    this->set_pars(mu);
}

void NeoHookean::set_pars(double mu) {

    this->mu = mu;
}

void NeoHookean::stress(const kinematics::kinematics<9> &kin, double stress[]) {

    for (int i = 0; i < 9; i++) {
        stress[i] = mu * id3d[i];
    }
}

} // namespace constitutive_models
