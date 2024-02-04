// #include <stdio.h>
// #include <iostream>
// #include <stdlib.h>
// #include <algorithm>
#include "objective.hpp"
#include "models.hpp"
#include "templates.hpp"

namespace residuals {

using namespace sim;
/* ------------------------------------------------------------------------------
 |  This file provides the definitions for calculating the residuals for
 |  the different model forms
 |
 |  These models are mainly used for the collaboration with Alexey Kamenskiy
 |
 |  Author: Will Zhang
 |  Dependencies: None
 ----------------------------------------------------------------------------- */

/*******************************************************************************
 * Calculating the residual for the hyperelastic functions
 *
 * COMMENTS:
 * The main form is predetermined, the forms differs by whether the model is
 * scaled.
 *******************************************************************************/

double planar_elastin_matrix_residual(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double args[],
    double stress[], double dt[], double weights[], double deltaCG[], double hysteresis[],
    double alphas[], int index[], int select[], int n, int dim, int nprot, int skip
) {
    return calc_residual_general<
        PlanarElastinMatrix, residual_body_below, hysteresis_body_null, penalty_body_kappa,
        calculate_viscopart_body_null>(
        pars, fiber, visco, Tf, Cmax, args, stress, dt, weights, deltaCG, hysteresis, alphas, index,
        select, n, dim, nprot, skip
    );
}

double thoracic_default_residual_VE_scaled_hyst_relax(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double args[],
    double stress[], double dt[], double weights[], double deltaCG[], double hysteresis[],
    double alphas[], int index[], int select[], int n, int dim, int nprot, int skip
) {
    return calc_residual_general<
        ThoracicDefaultVEScaled, residual_body, hysteresis_body, penalty_body_kappa_3,
        calculate_viscopart_body_null>(
        pars, fiber, visco, Tf, Cmax, args, stress, dt, weights, deltaCG, hysteresis, alphas, index,
        select, n, dim, nprot, skip
    );
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/
} // namespace residuals
