#define _USE_MATH_DEFINES

#include "thoracic_circ_model.hpp"
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

namespace sim {

/*----------------------------------------------------------------------
 |  This provides the main models in the full constitutive model
 -----------------------------------------------------------------------*/
void ThoracicCircBase::get_scaled_pars(double pars[]) {
    pars[0] = m_matrix.mu;
    pars[1] = m_elastin.k1;
    pars[2] = m_muscle.get_scaled_modulus();
    pars[3] = m_muscle.k2;
}

void ThoracicCircBase::stress(
    const kinematics::kinematics<4> &kin, const double dt, double stress[]
) {
    double p;
    double mat[4], el[4], smc[4];
    p = m_matrix.stress(kin, mat);
    (void)m_elastin.stress(kin, el);
    (void)m_muscle.stress(kin, smc);
    for (int j = 0; j < ctv::prob_dim; j++) {
        stress[j] = mat[j] + el[j] + smc[j] - p * kin.I_n * kin.Cinv[j];
    }
}

void ThoracicCircVEBase::stress(
    const kinematics::kinematics<4> &kin, const double dt, double stress[]
) {
    double p;
    double mat[4], el[4], smc[4];
    p = m_matrix.stress(kin, mat);
    (void)m_elastin.stress(kin, el);
    (void)muscle.stress(kin, dt, smc);
    for (int j = 0; j < ctv::prob_dim; j++) {
        stress[j] = mat[j] + el[j] + smc[j] - p * kin.I_n * kin.Cinv[j];
    }
}

} // namespace sim