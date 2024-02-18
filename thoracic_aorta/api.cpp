// #include <stdio.h>
// #include <iostream>
// #include <stdlib.h>
// #include <algorithm>
#include "api.hpp"
#include "../simulate/templates.hpp"
#include "constants.hpp"
#include "planar_elastin_matrix_model.hpp"
#include "thoracic_default_model.hpp"

namespace thoracic {

/* ------------------------------------------------------------------------------
 |  This file provides the definitions for calculating the residuals for
 |  the different model forms
 |
 |  These models are mainly used for the collaboration with Alexey Kamenskiy
 |
 |  Author: Will Zhang
 |  Dependencies: None
 ----------------------------------------------------------------------------- */

double penalty_body_null(
    const double pars[], const double fiber[], const double visco[], const double data[]
) {
    return 0.0;
}
/*******************************************************************************
 * Calculating the residual for the hyperelastic functions
 *
 * COMMENTS:
 * The main form is predetermined, the forms differs by whether the model is
 * scaled.
 *******************************************************************************/

void planar_elastin_matrix_get_model_parameters_scaled(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    double pars_out[2]
) {
    PlanarElastinMatrix psi(pars, fiber, visco, Tf, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
}

void planar_elastin_matrix_simulate(
    const double pars[], const double fiber[], const double caputo[], double Tf,
    const double Cmax[], const double args[], const double dt[], double stress[], int n
) {
    simulate::simulate<PlanarElastinMatrix, kinematics::deformation2D, 2>(
        pars, fiber, caputo, Tf, Cmax, args, dt, stress, n
    );
}

double planar_elastin_matrix_residual(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    const double args[], const double stress[], const double dt[], const double weights[],
    const double deltaCG[], const double hysteresis[], const double data[], const int index[],
    const int select[], int n, int nprot, int skip
) {
    return simulate::calc_residual<
        PlanarElastinMatrix, kinematics::deformation2D, simulate::quart_quad_residual,
        simulate::hysteresis_body_null, penalty_body_null, 2>(
        pars, fiber, visco, Tf, Cmax, args, stress, dt, weights, deltaCG, hysteresis, data, index,
        select, n, nprot, skip, M_w_hyst
    );
}

/*******************************************************************************
 * Calculating the residual for the hyperelastic functions
 *
 * COMMENTS:
 * The main form is predetermined, the forms differs by whether the model is
 * scaled.
 *******************************************************************************/

double thoracic_default_residual_VE_scaled_hyst_relax(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    const double args[], const double stress[], const double dt[], const double weights[],
    const double deltaCG[], const double hysteresis[], const double alphas[], const int index[],
    const int select[], int n, int nprot, int skip
) {
    return simulate::calc_residual<
        ThoracicDefaultVEScaled, kinematics::deformation2D, simulate::quart_quad_residual,
        simulate::hysteresis_body<2>, penalty_body_null, 2>(
        pars, fiber, visco, Tf, Cmax, args, stress, dt, weights, deltaCG, hysteresis, alphas, index,
        select, n, nprot, skip, M_w_hyst
    );
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/
} // namespace thoracic
