#define _USE_MATH_DEFINES

#include "femoral_model.hpp"
#include "constants.hpp"
#include <cmath>

/*----------------------------------------------------------------------
 |  This file provides the definitions of the different model forms
 |  which combines the primative constitive model from the other cpp files
 |  in the constitutive_models namespace
 |
 |  These models are mainly used for the collaboration with Alexey Kamenskiy
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

using namespace constitutive_models;

namespace femoral {

/*----------------------------------------------------------------------
 |  This provides the main models in the full constitutive model
 -----------------------------------------------------------------------*/
FemoralBase::FemoralBase(const double pars[], const double fiber[])
    : m_matrix(pars[0]), m_elastin(pars[1], 0.0, fiber[0]),
      m_muscle(pars[2], pars[3], fiber[0] + M_PI_2),
      m_collagen(pars[4], pars[5], fiber[0], fiber[1], -fiber[1], M_kip, M_kop) {
}

FemoralBase::FemoralBase(const double pars[], const double fiber[], const double Cmax[])
    : m_matrix(pars[0]), m_elastin(pars[1], 0.0, fiber[0]),
      m_muscle(pars[2], pars[3], fiber[0] + M_PI_2, Cmax),
      m_collagen(pars[4], pars[5], fiber[0], fiber[1], -fiber[1], M_kip, M_kop, Cmax) {
}

FemoralBase::~FemoralBase() {
}

void FemoralBase::get_scaled_pars(double pars[]) {
    pars[0] = m_matrix.mu;
    pars[1] = m_elastin.k1;
    pars[2] = m_muscle.get_scaled_modulus();
    pars[3] = m_muscle.k2;
    pars[4] = m_collagen.get_scaled_modulus();
    pars[5] = m_collagen.k2;
}

void FemoralBase::stress(const kinematics::kinematics<4> &kin, double dt, double stress[]) {
    double p = 0.0;
    double iso[4], el[4], smc[4], col[4];
    p = m_matrix.stress(kin, iso);
    (void)m_elastin.stress(kin, el);
    (void)m_muscle.stress(kin, smc);
    p = p + m_collagen.stress(kin, col);
    for (int j = 0; j < 4; j++) {
        stress[j] = iso[j] + col[j] + el[j] + smc[j] - p * kin.I_n * kin.Cinv[j];
    }
}

FemoralVEBase::FemoralVEBase(
    const double pars[], const double fiber[], const double visco[], double Tf
)
    : FemoralBase(pars, fiber), muscle(m_muscle, visco[0], Tf), collagen(m_collagen, visco[1], Tf) {
}

FemoralVEBase::FemoralVEBase(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[]
)
    : FemoralBase(pars, fiber, Cmax), muscle(m_muscle, visco[0], Tf),
      collagen(m_collagen, visco[1], Tf) {
}
FemoralVEBase::~FemoralVEBase() {
}

void FemoralVEBase::stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]) {
    double p = 0.0;
    double iso[4], el[4], smc[4], col[4];
    p = m_matrix.stress(kin, iso);
    (void)m_elastin.stress(kin, el);
    p = p + muscle.stress(kin, dt, smc);
    p = p + collagen.stress(kin, dt, col);
    for (int j = 0; j < 4; j++) {
        stress[j] = iso[j] + col[j] + el[j] + smc[j] - p * kin.I_n * kin.Cinv[j];
    }
}
} // namespace femoral