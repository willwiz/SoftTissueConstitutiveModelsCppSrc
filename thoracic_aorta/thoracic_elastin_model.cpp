#include "thoracic_elastin_model.hpp"

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

void PlanarElastinMatrix::get_scaled_pars(double pars[]) {
    pars[0] = m_matrix.mu;
    pars[1] = m_elastin.get_scaled_modulus();
}

void PlanarElastinMatrix::stress(const kinematics::kinematics<4> &kin, double dt, double stress[]) {
    double p = 0.0;
    double mat[4], el[4];
    p = m_matrix.stress(kin, mat);
    (void)m_elastin.stress(kin, el);
    for (int j = 0; j < 4; j++) {
        stress[j] = mat[j] + el[j] - p * kin.I_n * kin.Cinv[j];
    }
}

} // namespace thoracic