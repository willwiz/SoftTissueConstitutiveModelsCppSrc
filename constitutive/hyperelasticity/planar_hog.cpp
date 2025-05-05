#include "planar_hog.hpp"
#include "../../kinematics/tensor_algebra.hpp"
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
PlanarHog2D::PlanarHog2D() : E1{}, E2{} {};
PlanarHog2D::~PlanarHog2D() {};

PlanarHog2D::PlanarHog2D(double k1, double k2, double theta, double kappa)
    : k1(k1), k2(k2), E1{}, E2{} {
    this->set_pars(theta, kappa);
};

PlanarHog2D::PlanarHog2D(double k1, double k2, double theta, double kappa, const double Cmax[])
    : k1(k1), k2(k2) {
    this->set_pars(theta, kappa, Cmax);
};

void PlanarHog2D::set_pars(double theta, double kappa) {
    double c = cos(theta);
    double s = sin(theta);
    m[0] = c * c;
    m[1] = c * s;
    m[2] = m[1];
    m[3] = s * s;
    for (size_t i = 0; i < 4; i++) {
        m[i] = (1.0 - 2.0 * kappa) * m[i] + kappa * id2d[i];
    }
}

// Scaled version
void PlanarHog2D::set_pars(double theta, double kappa, const double Cmax[]) {
    this->set_pars(theta, kappa);
    E1 = ddot2D(m, Cmax) - 1;
    k1 = k1 / E1;
    E2 = E1 * E1;
}

double PlanarHog2D::get_scaled_modulus() {
    return k1 * exp(-k2 * E2);
}

// Stress functions
double PlanarHog2D::stress(const kinematics::kinematics<4> &kin, double stress[4]) {

    double I_4 = ddot2D(m, kin.C) - 1;
    double dWd4 = k1 * I_4 * exp(k2 * (I_4 * I_4 - E2));

    for (int i = 0; i < 4; i++) {
        stress[i] = dWd4 * m[i];
    }
    return 0.0;
}

void PlanarHog2D::stress(const double args[4], double stress[4]) {

    kinematics::deformation2D kin(args);
    (void)this->stress(kin, stress);
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/

} // namespace constitutive_models