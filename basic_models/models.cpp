#define _USE_MATH_DEFINES

#include "models.hpp"
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

void ModelMatrix::stress(const kinematics::kinematics<9> &kin, const double dt, double stress[]) {
    m_model.stress(kin, stress);
}

void ModelElastin::stress(const kinematics::kinematics<9> &kin, const double dt, double stress[]) {
    m_model.stress(kin, stress);
}

void ModelSMC::stress(const kinematics::kinematics<9> &kin, const double dt, double stress[]) {
    m_model.stress(kin, stress);
}

void ModelCollagen::stress(const kinematics::kinematics<9> &kin, const double dt, double stress[]) {
    m_model.stress(kin, stress);
}

void ModelSMCVE::stress(const kinematics::kinematics<9> &kin, const double dt, double stress[]) {
    ve_model.stress(kin, dt, stress);
}

void ModelCollagenVE::stress(
    const kinematics::kinematics<9> &kin, const double dt, double stress[]
) {
    ve_model.stress(kin, dt, stress);
}

} // namespace sim