#include "../kinematics/kinematics.hpp"
#include "../kinematics/tensor_algebra.hpp"
#include "hog_double_2D.hpp"
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
 * Standard holzapfel ogden with two fiber families
 *
 * COMMENTS:
 *
 *******************************************************************************/
HOGDouble2D::HOGDouble2D(double k1, double k2, double theta, double alpha) {
    this->set_pars(k1, k2, theta, alpha);
};

void HOGDouble2D::set_pars(double k1, double k2, double theta, double alpha) {

    this->k1 = k1;
    this->k2 = k2;

    double ca4 = cos(theta + alpha);
    double sa4 = sin(theta + alpha);
    double ca6 = cos(theta - alpha);
    double sa6 = sin(theta - alpha);

    this->m4[0] = ca4 * ca4;
    this->m4[1] = ca4 * sa4;
    this->m4[2] = m4[1];
    this->m4[3] = sa4 * sa4;

    this->m6[0] = ca6 * ca6;
    this->m6[1] = ca6 * sa6;
    this->m6[2] = m6[1];
    this->m6[3] = sa6 * sa6;
}

// Stress functions
double HOGDouble2D::stress(const kinematics::kinematics<4> &kin, double stress[4]) {

    double I_4 = ddot2D(m4, kin.C) - 1;
    double I_6 = ddot2D(m6, kin.C) - 1;
    double dWd4 = k1 * I_4 * exp(k2 * I_4 * I_4);
    double dWd6 = k1 * I_6 * exp(k2 * I_6 * I_6);
    for (int i = 0; i < 4; i++) {
        stress[i] = dWd4 * m4[i] + dWd6 * m6[i];
    }
    return 0.0;
}

void HOGDouble2D::stress(double args[4], double stress[4]) {

    kinematics::deformation2D kin(args);
    (void)this->stress(kin, stress);
    // for (int i = 0; i < 4; i++)
    // {
    //   stress[i] = stress[i] - p*kin.I_n*kin.Cinv[i];
    // }
}
/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/

} // namespace constitutive_models