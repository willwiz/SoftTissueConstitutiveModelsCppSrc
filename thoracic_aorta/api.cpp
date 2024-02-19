// #include <stdio.h>
// #include <iostream>
// #include <stdlib.h>
// #include <algorithm>
#include "api.hpp"
#include "../simulate/templates.hpp"
#include "constants.hpp"
#include "ensemble_template.hpp"
#include "thoracic_elastin_model.hpp"
#include "thoracic_ensemble_model.hpp"
#include "thoracic_ve_model.hpp"

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
inline double quart_quad_residual(double sim, double data, double strain) {
    double difference = sim - data;
    double r2 = difference * difference;
    double mult = (difference > 0) ? r2 : 1.0;
    return r2 * mult;
}

inline double quadratic_residual(double sim, double data, double strain) {
    double difference = sim - data;
    return difference * difference;
}

double penalty_body_null(
    const double pars[], const double fiber[], const double visco[], const double data[]
) {
    return 0.0;
}

double penalty_ensemble_null(const double pars[], const double visco[], const double data[]) {
    return 0.0;
};

double penalty_ensemble3(const double pars[], const double visco[], const double data[]) {
    double b_s = pars[3];
    // double b_c = pars[5] - 2*pars[13];
    double v_c = visco[0] - data[7];
    double v_s = visco[1] - data[8];
    return 1.0 + 0.1 * b_s * b_s + 10.0 * v_c * v_c + 10.0 * v_s * v_s;
}

double hysteresis_body_null(const double sims[], const double deltaCG[], int n, double hysteresis) {
    return 0.0;
}

// Hysteresis Calculation
double hysteresis_body(const double sims[], const double deltaCG[], int n, double hysteresis) {

    double hyst = 0;
    for (int i = 0; i < n; i++) {
        hyst = hyst + sims[i] * deltaCG[i];
    }
    double ds = hyst - hysteresis;

    return ds * ds;
}

/*******************************************************************************
 * Calculating the residual for the hyperelastic functions
 *
 * COMMENTS:
 * The main form is predetermined, the forms differs by whether the model is
 * scaled.
 *******************************************************************************/

void thoracic_elastin_get_parameters(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    double pars_out[2]
) {
    PlanarElastinMatrix psi(pars, fiber, visco, Tf, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
}

void thoracic_elastin_simulate(
    const double pars[], const double fiber[], const double caputo[], double Tf,
    const double Cmax[], const double args[], const double dt[], double stress[], int n
) {
    simulate::simulate<PlanarElastinMatrix, kinematics::deformation2D, 2>(
        pars, fiber, caputo, Tf, Cmax, args, dt, stress, n
    );
}

double thoracic_elastin_residual(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    const double args[], const double stress[], const double dt[], const double weights[],
    const double deltaCG[], const double hysteresis[], const double data[], const int index[],
    const int select[], int n, int nprot, int skip
) {
    return simulate::calc_residual<
        PlanarElastinMatrix, kinematics::deformation2D, simulate::quart_quad_residual,
        simulate::hysteresis_body_null<2>, penalty_body_null, 2>(
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
void thoracic_ensemble_get_parameters(
    const double pars[], const double visco[], double Tf, const double Cmax[], double pars_out[6]
) {
    ThoracicEnsembleVE psi(pars, visco, Tf, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
}

void thoracic_ensemble_simulate(
    const double pars[], const double visco[], const double Tf, const double Cmax[],
    const double strain[], const double dt[], double stress[], int n
) {
    ensemble_simulate<ThoracicEnsembleVE>(pars, visco, Tf, Cmax, strain, dt, &stress[0], n);
}

double thoracic_ensemble_residual(
    const double pars[], const double visco[], double Tf, const double Cmax[],
    const double strain[], const double stress[], const double dt[], const double deltaCG[],
    double hysteresis, const double data[], int n, int skip
) {
    return calc_ensemble_objective<
        ThoracicEnsembleVE, quart_quad_residual, hysteresis_body, penalty_ensemble3>(
        pars, visco, Tf, Cmax, strain, stress, dt, deltaCG, hysteresis, data, n, skip
    );
}

/*******************************************************************************
 * Calculating the residual for the hyperelastic functions
 *
 * COMMENTS:
 * The main form is predetermined, the forms differs by whether the model is
 * scaled.
 *******************************************************************************/

void thoracic_ve_get_parameters(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    double pars_out[10]
) {
    ThoracicVEScaled psi(pars, fiber, visco, Tf, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[6] = fiber[0];
    pars_out[7] = fiber[1];
    pars_out[8] = visco[0];
    pars_out[9] = visco[1];
}

void thoracic_ve_simulate(
    const double pars[], const double fiber[], const double caputo[], double Tf,
    const double Cmax[], const double args[], const double dt[], double stress[], int n
) {
    simulate::simulate<ThoracicVEScaled, kinematics::deformation2D, 2>(
        pars, fiber, caputo, Tf, Cmax, args, dt, &stress[0], n
    );
}

double thoracic_ve_residual(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    const double args[], const double stress[], const double dt[], const double weights[],
    const double deltaCG[], const double hysteresis[], const double data[], const int index[],
    const int select[], int n, int nprot, int skip
) {
    return simulate::calc_residual<
        ThoracicVEScaled, kinematics::deformation2D, simulate::quart_quad_residual,
        simulate::hysteresis_body<2>, penalty_body_null, 2>(
        pars, fiber, visco, Tf, Cmax, args, stress, dt, weights, deltaCG, hysteresis, data, index,
        select, n, nprot, skip, M_w_hyst
    );
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/
} // namespace thoracic
