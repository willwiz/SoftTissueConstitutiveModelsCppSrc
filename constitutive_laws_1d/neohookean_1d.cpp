#include <cmath>
#include "../kinematics/tensor_algebra.hpp"
#include "../kinematics/kinematics.hpp"
#include "neohookean_1d.hpp"

namespace constitutive_models {

/*----------------------------------------------------------------------
 |  Standard Neohookean, no much to say
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

    NeoHookean1D::NeoHookean1D(double mu)
    {
        this->set_pars(mu);
    }

    void NeoHookean1D::set_pars(double mu) {

        this -> mu = mu;
    }

    double NeoHookean1D::stress(const kinematics::kinematics<1> &kin, double stress[]){
        stress[0] = mu;
        return mu;
    }

    double NeoHookean1D::stress(double args){
        double stress[1];
        kinematics::deformation1D kin(args);
        double p = this->stress(kin, stress);
        return stress[0] - p * kin.I_n * kin.Cinv[0];
    }

}
