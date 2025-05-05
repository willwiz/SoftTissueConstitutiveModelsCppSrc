#include "transverse_iso.hpp"
#include "../../kinematics/tensor_algebra.hpp"
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

/*******************************************************************************
 * Standard holzapfel ogden with one fiber family
 *
 * COMMENTS:
 *
 *******************************************************************************/
TransverseIso1D::TransverseIso1D() {};
TransverseIso1D::~TransverseIso1D() {};

TransverseIso1D::TransverseIso1D(double k1, double k2, double kappa)
    : k1(k1), k2(k2), fiber{1 - kappa + kappa / 3.0}, iso{kappa / 3.0} {};

TransverseIso1D::TransverseIso1D(double k1, double k2, double kappa, double Cmax)
    : k1{k1}, k2{k2}, fiber{1 - kappa + kappa / 3.0}, iso{kappa / 3.0} {};

inline double calc_W1(double b, double lambda2, double lambda) {
    double exponent = b * (lambda2 + 2.0 / lambda - 3.0);
    return exp(exponent);
}

void TransverseIso1D::set_pars(double kappa) {
    fiber = (1 - kappa);
    iso = kappa / 3.0;
}

// Stress functions
double TransverseIso1D::stress(const kinematics::kinematics<1> &kin, double stress[]) {

    double lambda = std::sqrt(kin.C[0]);
    double W1 = calc_W1(k2, kin.C[0], lambda);

    stress[0] = fiber * (W1 * kin.C[0] - 1.0);
    // std::cout << kin.C[0] << " " << kin.I_n << "\n";
    // std::cout << dWd4 << " " << I_4 << " " << E2 << " " << E1 << "\n";

    return iso * (W1 * lambda - 1.0);
}

double TransverseIso1D::stress(double args) {
    double stress[1];
    kinematics::deformation1D kin(args);
    double p = this->stress(kin, stress);
    return stress[0] - p * kin.I_n * kin.Cinv[0];
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/

} // namespace constitutive_models