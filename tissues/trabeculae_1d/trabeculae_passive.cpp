#define _USE_MATH_DEFINES
#include "trabeculae_passive.hpp"
#include <cmath>

namespace optimization_1d {

void TrabeculaePassive::stress(
    const kinematics::kinematics<1> &kin, const double dt, double stress[]
) {
    double p = 0.0;
    double iso[1], smc[1];
    p = m_matrix.stress(kin, iso);
    p = p + muscle.stress(kin, dt, smc);
    stress[0] = iso[0] + smc[0] - p * kin.I_n * kin.Cinv[0];
}

} // namespace optimization_1d
