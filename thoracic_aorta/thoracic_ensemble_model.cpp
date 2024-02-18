#define _USE_MATH_DEFINES
#include "thoracic_ensemble_model.hpp"
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

namespace thoracic {

/*----------------------------------------------------------------------
 |  This provides the main models in the full constitutive model
 -----------------------------------------------------------------------*/
ThoracicFullEnsembleBase::ThoracicFullEnsembleBase(double pars[])
    : m_matrix(0), m_elastin(pars[1], 0.0, 0.0, 0.5), m_muscle(pars[2], pars[3], 0.0, 0.5),
      m_collagen(pars[4], pars[5], 0.0, M_ideal_alpha, -M_ideal_alpha, M_kip, M_kop) {
}

ThoracicFullEnsembleBase::ThoracicFullEnsembleBase(double pars[], double Cmax[])
    : m_matrix(0), m_elastin(pars[1], 0.0, 0.0, 0.5), m_muscle(pars[2], pars[3], 0.0, 0.5, Cmax),
      m_collagen(pars[4], pars[5], 0.0, M_ideal_alpha, -M_ideal_alpha, M_kip, M_kop, Cmax) {
}

void ThoracicFullEnsembleBase::get_scaled_pars(double pars[]) {
    pars[0] = m_matrix.mu;
    pars[1] = m_collagen.get_scaled_modulus();
    pars[2] = m_collagen.k2;
    pars[3] = m_elastin.get_scaled_modulus();
}

void ThoracicFullEnsembleBase::stress(
    const kinematics::kinematics<4> &kin, const double dt, double stress[]
) {
    double p = 0.0;
    double mat[4], el[4], smc[4], col[4];
    p = m_matrix.stress(kin, mat);
    (void)m_elastin.stress(kin, el);
    (void)m_muscle.stress(kin, smc);
    p = p + m_collagen.stress(kin, col);
    for (int j = 0; j < 4; j++) {
        stress[j] = mat[j] + el[j] + smc[j] + col[j] - p * kin.I_n * kin.Cinv[j];
    }
}

void ThoracicFullEnsembleVE::stress(
    const kinematics::kinematics<4> &kin, const double dt, double stress[]
) {
    double p = 0.0;
    double mat[4], el[4], smc[4], col[4];
    p = m_matrix.stress(kin, mat);
    (void)m_elastin.stress(kin, el);
    (void)muscle.stress(kin, dt, smc);
    p = p + collagen.stress(kin, dt, col);
    for (int j = 0; j < 4; j++) {
        stress[j] = mat[j] + el[j] + smc[j] + col[j] - p * kin.I_n * kin.Cinv[j];
    }
}

} // namespace thoracic