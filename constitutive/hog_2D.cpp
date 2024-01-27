#include "../kinematics/kinematics.hpp"
#include "../kinematics/tensor_algebra.hpp"
#include "hog_2D.hpp"
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
Hog2D::Hog2D() : E1{}, E2{} {};
Hog2D::~Hog2D(){};

Hog2D::Hog2D(double k1, double k2, double theta) : E1{}, E2{} {
    this->set_pars(k1, k2, theta);
};

Hog2D::Hog2D(double k1, double k2, double theta, double Cmax[]) {
    this->set_pars(k1, k2, theta, Cmax);
};

void Hog2D::set_pars(double k1, double k2, double theta) {
    this->k1 = k1;
    this->k2 = k2;
    double c = cos(theta);
    double s = sin(theta);
    this->m[0] = c * c;
    this->m[1] = c * s;
    this->m[2] = this->m[1];
    this->m[3] = s * s;
}

// Scaled version
void Hog2D::set_pars(double k1, double k2, double theta, double Cmax[]) {
    this->set_pars(k1, k2, theta);
    E1 = ddot2D(m, Cmax) - 1;
    this->k1 = k1 / E1;
    this->E2 = E1 * E1;
}

double Hog2D::get_scaled_modulus() {
    return k1 * exp(-k2 * E2);
}

// Stress functions
double Hog2D::stress(const kinematics::kinematics<4> &kin, double stress[4]) {

    double I_4 = ddot2D(m, kin.C) - 1;
    double dWd4 = k1 * I_4 * exp(k2 * (I_4 * I_4 - E2));

    for (int i = 0; i < 4; i++) {
        stress[i] = dWd4 * m[i];
    }
    return 0.0;
}

void Hog2D::stress(double args[4], double stress[4]) {

    kinematics::deformation2D kin(args);
    (void)this->stress(kin, stress);
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/

} // namespace constitutive_models