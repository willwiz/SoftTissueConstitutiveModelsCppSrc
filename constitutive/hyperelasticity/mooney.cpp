#include "mooney.hpp"
#include "../kinematics/tensor_algebra.hpp"
#include <cmath>

namespace constitutive_models {

/*----------------------------------------------------------------------
 |  Standard Neohookean, no much to say
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

// Fung form
mooney::mooney(double mu, double b) : m_mu(mu), m_b(b) {};

void mooney::set_pars(double mu, double b) {
    m_mu = mu;
    m_b = b;
}

double mooney::get_scaled_modulus() {
    return m_mu;
}

double mooney::stress(const kinematics::kinematics<4> &kin, double stress[]) {

    double exponent = m_mu * pow(kin.I_1 - 3, m_b);

    for (int i = 0; i < 4; i++) {
        stress[i] = exponent * id2d[i];
    }
    return exponent;
}

void mooney::stress(const double args[], double stress[]) {

    kinematics::deformation2D kin(args);
    double p = this->stress(kin, stress);
    for (int i = 0; i < 4; i++) {
        stress[i] = stress[i] - p * kin.I_n * kin.Cinv[i];
    }
}
} // namespace constitutive_models
