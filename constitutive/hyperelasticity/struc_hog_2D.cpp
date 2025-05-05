#include "struc_hog_2D.hpp"
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
 * Structural tensor forms
 *
 * COMMENTS:
 *
 *******************************************************************************/
StrucHOG2D::StrucHOG2D() : z{1.0} {};
StrucHOG2D::~StrucHOG2D() {};

StrucHOG2D::StrucHOG2D(
    double k, double b, double theta, double alpha, double beta, double kip, double kop
)
    : z{1.0} {
    this->set_pars(k, b, theta, alpha, beta, kip, kop);
};

StrucHOG2D::StrucHOG2D(
    double k, double b, double theta, double alpha, double beta, double kip, double kop,
    const double Cmax[]
) {
    this->set_pars(k, b, theta, alpha, beta, kip, kop, Cmax);
};

void StrucHOG2D::set_pars(
    double k, double b, double theta, double alpha, double beta, double kip, double kop
) {
    double m4[4], m6[4];
    double ca, sa;

    this->k = k;
    this->b = b;
    double A = 2.0 * kop * kip;
    double B = 2.0 * kop * (1.0 - 2.0 * kip);
    double C = 1.0 - 2.0 * A - B; // Because this is 2D we add back the A, i.e. C = A+C

    ca = cos(theta + alpha);
    sa = sin(theta + alpha);

    m4[0] = ca * ca;
    m4[1] = ca * sa;
    m4[2] = m4[1];
    m4[3] = sa * sa;

    ca = cos(theta + beta);
    sa = sin(theta + beta);
    m6[0] = ca * ca;
    m6[1] = ca * sa;
    m6[2] = m6[1];
    m6[3] = sa * sa;

    for (int i = 0; i < 4; i++) {
        this->H4[i] = this->H6[i] = A * id2d[i];
        this->H4[i] = this->H4[i] + B * m4[i];
        this->H6[i] = this->H6[i] + B * m6[i];
    }
    this->H_33 = A + C;
    this->low_slope = k * exp(-b);
    this->high_slope = (1.0 + 2.0 * b) * k;
}

// Scaled version
void StrucHOG2D::set_pars(
    double k, double b, double theta, double alpha, double beta, double kip, double kop,
    const double Cmax[]
) {
    this->set_pars(k, b, theta, alpha, beta, kip, kop);

    double det = Cmax[0] * Cmax[3] - Cmax[1] * Cmax[1];
    double I_n = H_33 / det;

    double I_4 = ddot2D(H4, Cmax) + I_n - 1.0;
    double I_6 = ddot2D(H6, Cmax) + I_n - 1.0;
    z = 1.0 / std::max(I_4, I_6);

    // this->k1 = k1 / z;
    // this->E2 = z * z;
}

double StrucHOG2D::get_scaled_modulus() {
    return k * z * exp(-b);
}

// Stress functions
double StrucHOG2D::stress(const kinematics::kinematics<4> &kin, double stress[4]) {
    double I_n = H_33 * kin.I_n;
    double I_4 = (ddot2D(H4, kin.C) + I_n - 1.0) * z;
    double I_6 = (ddot2D(H6, kin.C) + I_n - 1.0) * z;
    double dWd4, dWd6;
    if (I_4 < 0.0) {
        dWd4 = low_slope * I_4;
    } else if (I_4 < 1.0) {
        dWd4 = k * I_4 * exp(b * (I_4 * I_4 - 1.0));
    } else {
        dWd4 = high_slope * (I_4 - 1.0) + k;
    }
    if (I_6 < 0.0) {
        dWd6 = low_slope * I_6;
    } else if (I_6 < 1.0) {
        dWd6 = k * I_6 * exp(b * (I_6 * I_6 - 1.0));
    } else {
        dWd6 = high_slope * (I_6 - 1.0) + k;
    }
    for (int i = 0; i < 4; i++) {
        stress[i] = dWd4 * H4[i] + dWd6 * H6[i];
    }
    return H_33 * (dWd4 + dWd6);
}

void StrucHOG2D::stress(const double args[4], double stress[4]) {

    kinematics::deformation2D kin(args);
    double p = this->stress(kin, stress);
    for (int i = 0; i < 4; i++) {
        stress[i] = stress[i] - (p * kin.I_n) * kin.Cinv[i];
    }
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/

} // namespace constitutive_models