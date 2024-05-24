#include "hog_2D.hpp"
#include "../kinematics/kinematics.hpp"
#include "../kinematics/tensor_algebra.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

/*----------------------------------------------------------------------
 |  This file provides the holzapfel class of models in 2D. Single and
 |  dual family models are available
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

namespace constitutive_models {

constexpr double std_strain[4] = {2.0, 0.0, 0.0, 2.0};
/*******************************************************************************
 * Standard holzapfel ogden with one fiber family
 *
 * COMMENTS:
 *
 *******************************************************************************/
Hog2D::Hog2D() : z{1.0} {};
Hog2D::~Hog2D(){};

Hog2D::Hog2D(const double k1, const double k2, const double theta) : z{1.0} {
    this->set_pars(k1, k2, theta);
};

Hog2D::Hog2D(const double k1, const double k2, const double theta, const double Cmax[]) {
    this->set_pars(k1, k2, theta, Cmax);
};

void Hog2D::set_pars(const double k1, const double k2, const double theta) {
    this->k = k1;
    this->b = k2;
    double c = cos(theta);
    double s = sin(theta);
    this->m[0] = c * c;
    this->m[1] = c * s;
    this->m[2] = this->m[1];
    this->m[3] = s * s;
    low_slope = k1 * exp(-k2);
    high_slope = (1.0 + 2.0 * k2) * k1;
}

// Scaled version
void Hog2D::set_pars(const double k1, const double k2, const double theta, const double Cmax[]) {
    this->k = k1;
    this->b = k2;
    double c = cos(theta);
    double s = sin(theta);
    this->m[0] = c * c;
    this->m[1] = c * s;
    this->m[2] = this->m[1];
    this->m[3] = s * s;
    low_slope = k1 * exp(-k2);
    high_slope = (1.0 + 2.0 * k2) * k1;
    z = 1.0 / (ddot2D(m, Cmax) - 1.0);
}

double Hog2D::get_scaled_modulus() {
    return k * z * exp(-b);
}

// Stress functions
double Hog2D::stress(const kinematics::kinematics<4> &kin, double stress[4]) {
    double dWd4;
    double x = (ddot2D(m, kin.C) - 1.0) * z;
    // std::cout << "z = " << z << ", x = " << x << std::endl;
    if (x < 0.0) {
        dWd4 = low_slope * x;
    } else if (x < 1.0) {
        dWd4 = k * x * exp(b * (x * x - 1.0));
    } else {
        dWd4 = high_slope * (x - 1.0) + k;
    }
    for (int i = 0; i < 4; i++) {
        stress[i] = dWd4 * m[i];
    }

    return 0.0;
}

void Hog2D::stress(const double args[4], double stress[4]) {

    kinematics::deformation2D kin(args);
    (void)this->stress(kin, stress);
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/

} // namespace constitutive_models