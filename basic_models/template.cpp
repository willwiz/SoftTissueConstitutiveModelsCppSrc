#include "template.hpp"
#include "../CTvalues_optimization.hpp"
#include "../kinematics/kinematics.hpp"
#include <cmath>

namespace templates {

/* ------------------------------------------------------------------------------
 |  This file provides the definitions for calculating the residuals for
 |  the different model forms
 |
 |  These models are mainly used for the collaboration with Alexey Kamenskiy
 |
 |  Author: Will Zhang
 |  Dependencies: None
 ----------------------------------------------------------------------------- */

template <class matlaw>
void simulate(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
) {
    int strd_i;

    double vals[9];

    kinematics::deformation3D kin;
    matlaw law(pars, fiber, caputo, Tf);

    for (int i = 0; i < n; i++) {
        strd_i = 9 * i;
        // Calculate deformation
        kin.precompute(&args[strd_i]);
        // Compute Final Stress
        law.stress(kin, dt[i], vals);
        for (int j = 0; j < 9; j++) {
            stress[strd_i + j] = vals[j];
        }
    }
}
/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/
} // namespace templates
