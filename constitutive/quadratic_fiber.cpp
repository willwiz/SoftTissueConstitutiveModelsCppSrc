#include "../kinematics/kinematics.hpp"
#include "../kinematics/tensor_algebra.hpp"
#include "quadratic_fiber.hpp"
#include <cmath>

namespace constitutive_models {

QuadraticFiber::QuadraticFiber(/* args */) : k{}, m4{}, m6{} {
}
QuadraticFiber::~QuadraticFiber(){};
QuadraticFiber::QuadraticFiber(double mu, double theta, double alpha, double beta) : k{mu} {
    this->set_pars(mu, theta, alpha, beta);
}

void QuadraticFiber::set_pars(double k, double theta, double alpha, double beta) {
    this->k = k;
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
}

double QuadraticFiber::stress(const kinematics::kinematics<4> &kin, double stress[4]) {
    double I_4 = ddot2D(m4, kin.C) - 1;
    double I_6 = ddot2D(m6, kin.C) - 1;
    double dWd4 = k * I_4 * I_4;
    double dWd6 = k * I_6 * I_6;
    for (int i = 0; i < 4; i++) {
        stress[i] = dWd4 * m4[i] + dWd6 * m6[i];
    }
    return 0.0;
}

void QuadraticFiber::stress(double args[4], double stress[4]) {

    kinematics::deformation2D kin(args);
    (void)this->stress(kin, stress);
}

} // namespace constitutive_models
