#include "struc_hog_2D.hpp"
#include "../kinematics/kinematics.hpp"
#include "../kinematics/tensor_algebra.hpp"
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
 * Structural tensor forms
 *
 * COMMENTS:
 *
 *******************************************************************************/
StrucHOG2D::StrucHOG2D() : E1{}, E2{} {};
StrucHOG2D::~StrucHOG2D(){};

StrucHOG2D::StrucHOG2D(
    double k1, double k2, double theta, double alpha, double beta, double kip, double kop
)
    : E1{}, E2{} {
    this->set_pars(k1, k2, theta, alpha, beta, kip, kop);
};

StrucHOG2D::StrucHOG2D(
    double k1, double k2, double theta, double alpha, double beta, double kip, double kop,
    double Cmax[]
) {
    this->set_pars(k1, k2, theta, alpha, beta, kip, kop, Cmax);
};

void StrucHOG2D::set_pars(
    double k1, double k2, double theta, double alpha, double beta, double kip, double kop
) {

    this->k1 = k1;
    this->k2 = k2;
    this->A = 2.0 * kop * kip;
    this->B = 2.0 * kop * (1.0 - 2.0 * kip);
    this->C = 1.0 - 3.0 * A - B;

    double ca4 = cos(theta + alpha);
    double sa4 = sin(theta + alpha);
    double ca6 = cos(theta + beta);
    double sa6 = sin(theta + beta);

    this->m4[0] = ca4 * ca4;
    this->m4[1] = ca4 * sa4;
    this->m4[2] = m4[1];
    this->m4[3] = sa4 * sa4;

    this->m6[0] = ca6 * ca6;
    this->m6[1] = ca6 * sa6;
    this->m6[2] = m6[1];
    this->m6[3] = sa6 * sa6;

    for (int i = 0; i < 4; i++) {
        this->H4[i] = A * id2d[i] + B * m4[i];
        this->H6[i] = A * id2d[i] + B * m6[i];
    }
}

// Scaled version
void StrucHOG2D::set_pars(
    double k1, double k2, double theta, double alpha, double beta, double kip, double kop,
    double Cmax[]
) {
    this->set_pars(k1, k2, theta, alpha, beta, kip, kop);

    double det = Cmax[0] * Cmax[3] - Cmax[1] * Cmax[1];
    // double I_n = 1.0 / det;
    // double I_1 = Cmax[0] + Cmax[3] + I_n;
    double I_4 = ddot(H4, Cmax);
    double I_6 = ddot(H6, Cmax);
    // double E6 = A * I_1 + C * I_n - 1.0;
    // double E4 = E6 + B * I_4;
    // E6 = E6 + B * I_6;
    E1 = 0.5 * (I_4 + I_6);
    this->k1 = k1 / E1;
    this->E2 = E1 * E1;
}

double StrucHOG2D::get_scaled_modulus() {
    return k1 * exp(-k2 * E2);
}

// Stress functions
double StrucHOG2D::stress(const kinematics::kinematics<4> &kin, double stress[4]) {
    double I_n = (A + C) * kin.I_n;
    double I_4 = ddot(H4, kin.C) + I_n - 1.0;
    double I_6 = ddot(H6, kin.C) + I_n - 1.0;
    // double E6 = A * kin.I_1 + C * kin.I_n - 1.0;
    // double E4 = E6 + B * I_4;
    // E6 = E6 + B * I_6;
    double dWd4 = k1 * I_4 * exp(k2 * (I_4 * I_4 - E2));
    double dWd6 = k1 * I_6 * exp(k2 * (I_6 * I_6 - E2));

    for (int i = 0; i < 4; i++) {
        stress[i] = dWd4 * H4[i] + dWd6 * H6[i];
    }
    return (A + C) * (dWd4 + dWd6);
}

void StrucHOG2D::stress(double args[4], double stress[4]) {

    kinematics::deformation2D kin(args);
    double p = this->stress(kin, stress);
    for (int i = 0; i < 4; i++) {
        stress[i] = stress[i] - p * kin.I_n * kin.Cinv[i];
    }
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/

} // namespace constitutive_models