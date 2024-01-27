#include "../kinematics/tensor_algebra.hpp"
#include "struc_hog.hpp"
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
StrucHOG2D::StrucHOG2D() : H4{}, H6{} {};
StrucHOG2D::~StrucHOG2D(){};

StrucHOG2D::StrucHOG2D(
    double k1, double k2, double theta, double alpha, double beta, double kip, double kop
)
    : H4{}, H6{} {
    this->set_pars(k1, k2, theta, alpha, beta, kip, kop);
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

    this->H4[0] = ca4 * ca4;
    this->H4[1] = ca4 * sa4;
    this->H4[3] = H4[1];
    this->H4[4] = sa4 * sa4;

    this->H6[0] = ca6 * ca6;
    this->H6[1] = ca6 * sa6;
    this->H6[3] = H6[1];
    this->H6[4] = sa6 * sa6;

    for (int i = 0; i < 9; i++) {
        this->H4[i] = A * id3d[i] + B * H4[i];
        this->H6[i] = A * id3d[i] + B * H6[i];
    }
    H4[8] = H4[8] + C;
    H6[8] = H6[8] + C;
}

// Stress functions
void StrucHOG2D::stress(const kinematics::kinematics<9> &kin, double stress[9]) {
    double I_4 = std::max(ddot(H4, kin.C, 9) - 1.0, 0.0);
    double I_6 = std::max(ddot(H6, kin.C, 9) - 1.0, 0.0);
    // double E6 = A * kin.I_1 + C * kin.I_n - 1.0;
    // double E4 = E6 + B * I_4;
    // E6 = E6 + B * I_6;
    double dWd4 = k1 * I_4 * exp(k2 * (I_4 * I_4));
    double dWd6 = k1 * I_6 * exp(k2 * (I_6 * I_6));

    for (int i = 0; i < 9; i++) {
        stress[i] = dWd4 * H4[i] + dWd6 * H6[i];
    }
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/

} // namespace constitutive_models