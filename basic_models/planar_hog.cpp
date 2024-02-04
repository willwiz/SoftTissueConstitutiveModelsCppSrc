#include "planar_hog.hpp"
#include "../kinematics/tensor_algebra.hpp"
#include <algorithm>
#include <cmath>

/*----------------------------------------------------------------------
 |  This file provides the holzapfel class of models in 2D. Single and
 |  dual family models are available
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

namespace constitutive_models {

/*******************************************************************************
 * Standard holzapfel ogden with one fiber family
 *
 * COMMENTS:
 *
 *******************************************************************************/

const double id_planar[9] = {1, 0, 0, 0, 1, 0, 0, 0, 0};

PlanarHog3D::PlanarHog3D() : E1{}, E2{} {};
PlanarHog3D::~PlanarHog3D(){};

PlanarHog3D::PlanarHog3D(double k1, double k2, double theta, double kappa)
    : k1(k1), k2(k2), E1{}, E2{}, m{} {
    this->set_pars(theta, kappa);
};

void PlanarHog3D::set_pars(double theta, double kappa) {
    double c = cos(theta);
    double s = sin(theta);
    m[0] = c * c;
    m[1] = c * s;
    m[3] = m[1];
    m[4] = s * s;
    for (size_t i = 0; i < 9; i++) {
        m[i] = (1.0 - 2.0 * kappa) * m[i] + kappa * id_planar[i];
    }
}

// Stress functions
void PlanarHog3D::stress(const kinematics::kinematics<9> &kin, double stress[4]) {

    double I_4 = ddot(m, kin.C, 9) - 1.0;
    double dWd4 = k1 * I_4 * exp(k2 * (I_4 * I_4));

    for (int i = 0; i < 9; i++) {
        stress[i] = dWd4 * m[i];
    }
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/

} // namespace constitutive_models