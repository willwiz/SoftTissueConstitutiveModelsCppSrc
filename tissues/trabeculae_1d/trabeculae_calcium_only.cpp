#define _USE_MATH_DEFINES
#include "trabeculae_calcium_only.hpp"
#include <cmath>

namespace optimization_1d {

void TrabeculaeCalciumOnly::stress(
    const kinematics::kinematics<1> &kin, const double dt, double stress[]
) {
    double p = 0.0;
    double iso[1], smc[1], tit[1], cal[1];
    p = matrix.stress(kin, iso);
    p = p + muscle.stress(kin, dt, smc);
    p = p + titin.stress(kin, dt, tit);
    p = p + calcium.stress(kin, dt, cal);
    stress[0] = iso[0] + smc[0] + tit[0] + cal[0] - p * kin.I_n * kin.Cinv[0];
}

} // namespace optimization_1d
