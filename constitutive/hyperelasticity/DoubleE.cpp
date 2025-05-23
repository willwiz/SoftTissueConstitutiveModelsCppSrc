#define _USE_MATH_DEFINES
#include "DoubleE.hpp"
#include "../../kinematics/tensor_algebra.hpp"
#include <cmath>

namespace constitutive_models {

DoubleE2D::DoubleE2D() {};
DoubleE2D::~DoubleE2D() {};
DoubleE2D::DoubleE2D(double b1, double b2, double k1, double k2, double k3, double theta) {
    this->set_pars(b1, b2, k1, k2, k3, theta);
};

void DoubleE2D::set_pars(double b1, double b2, double k1, double k2, double k3, double theta) {

    this->b1 = b1;
    this->b2 = b2;
    this->k1 = k1;
    this->k2 = k2;
    this->k3 = k3;

    double ca4 = cos(theta);
    double sa4 = sin(theta);
    double ca6 = cos(theta + M_PI_2);
    double sa6 = sin(theta + M_PI_2);

    this->mxm[0] = ca4 * ca4;
    this->mxm[1] = ca4 * sa4;
    this->mxm[2] = mxm[1];
    this->mxm[3] = sa4 * sa4;

    this->nxn[0] = ca6 * ca6;
    this->nxn[1] = ca6 * sa6;
    this->nxn[2] = nxn[1];
    this->nxn[3] = sa6 * sa6;

    this->mxn[0] = ca4 * ca6;
    this->mxn[1] = 0.5 * (ca4 * sa6 + sa4 * ca6);
    this->mxn[2] = mxn[1];
    this->mxn[3] = sa4 * sa6;
}

double DoubleE2D::stress(const kinematics::kinematics<4> &kin, double stress[4]) {
    double I_ff = ddot2D(mxm, kin.C);
    double I_ss = ddot2D(nxn, kin.C);
    double I_fs = ddot2D(mxn, kin.C);

    double W1 = exp(b1 * kin.I_1m3);
    double W2 = exp(b2 * I_fs * I_fs);

    double W1ff = k1 * (W1 * I_ff - 1.0);
    double W1ss = k2 * (W1 * I_ss - 1.0);
    double W2fs = k3 * (W2 * I_fs);
    for (int i = 0; i < 4; i++) {
        stress[i] = W1ff * mxm[i] + W1ss * nxn[i] + W2fs * mxn[i];
    }
    return 0.0;
}

void DoubleE2D::stress(const double args[], double stress[4]) {

    kinematics::deformation2D kin(args);
    double p = this->stress(kin, stress);
    for (int i = 0; i < 4; i++) {
        stress[i] = stress[i] - p * kin.I_n * kin.Cinv[i];
    }
}

} // namespace constitutive_models