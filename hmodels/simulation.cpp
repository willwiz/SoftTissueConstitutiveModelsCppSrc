#include "simulation.hpp"
#include "models.hpp"
#include "objective_template.hpp"

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

void femoral_get_model_parameters(
    double pars[], double fiber[], double visco[], double Tf, double pars_out[11]
) {
    FemoralBase psi(pars, fiber);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[7] = fiber[0];
    pars_out[8] = fiber[1];
    pars_out[9] = visco[0];
    pars_out[10] = visco[1];
}

void femoral_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[11]
) {
    FemoralBase psi(pars, fiber, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[7] = fiber[0];
    pars_out[8] = fiber[1];
    pars_out[9] = visco[0];
    pars_out[10] = visco[1];
}

void thoracic_iso_get_model_parameters(
    double pars[], double fiber[], double visco[], double Tf, double pars_out[9]
) {
    ThoracicIsoBase psi(pars, fiber);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[5] = fiber[0];
    pars_out[6] = fiber[1];
    pars_out[7] = visco[0];
    pars_out[8] = visco[1];
}

void thoracic_default_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[]
) {
    ThoracicDefaultBase psi(pars, fiber, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[4] = fiber[0];
    pars_out[5] = fiber[1];
    pars_out[6] = visco[0];
}

void thoracic_defaultfung_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[6]
) {
    ThoracicDefaultFungBase psi(pars, fiber, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[3] = fiber[0];
    pars_out[4] = fiber[1];
    pars_out[5] = visco[0];
}

void thoracic_iso_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[9]
) {
    ThoracicIsoBase psi(pars, fiber, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[4] = fiber[0];
    pars_out[5] = fiber[1];
    pars_out[6] = visco[0];
    pars_out[7] = visco[1];
}

void thoracic_lin_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[8]
) {
    ThoracicLinBase psi(pars, fiber, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[4] = fiber[0];
    pars_out[5] = fiber[1];
    pars_out[6] = visco[0];
    pars_out[7] = visco[1];
}

void thoracic_quad_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[8]
) {
    ThoracicQuadBase psi(pars, fiber, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[4] = fiber[0];
    pars_out[5] = fiber[1];
    pars_out[6] = visco[0];
    pars_out[7] = visco[1];
}

void thoracic_circ_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[]
) {
    ThoracicCircBase psi(pars, fiber, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[5] = fiber[0];
    pars_out[6] = fiber[1];
    pars_out[7] = fiber[2];
    pars_out[8] = visco[0];
    pars_out[9] = visco[1];
}

/*----------------------------------------------------------------------
|  Some basic models
-----------------------------------------------------------------------*/

void planar_elastin_matrix_simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<PlanarElastinMatrix>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}

void thoracic_circ_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<ThoracicCircVEScaled>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}

void thoracic_default_ve_simulate_standard(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<ThoracicDefaultVEStandard>(
        pars, fiber, caputo, Tf, Cmax, args, dt, stress, n
    );
}

void thoracic_default_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<ThoracicDefaultVEScaled>(
        pars, fiber, caputo, Tf, Cmax, args, dt, stress, n
    );
}

/*----------------------------------------------------------------------
|  The femoral artery models
-----------------------------------------------------------------------*/

void femoral_he_simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<FemoralHE>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}

void femoral_he_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<FemoralHEScaled>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}

void femoral_ve_simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<FemoralVE>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}

void femoral_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<FemoralVEScaled>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}
/*----------------------------------------------------------------------
|  The thoracic models
-----------------------------------------------------------------------*/

void thoracic_iso_he_simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<ThoracicIsoHE>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}

void thoracic_iso_he_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<ThoracicIsoHEScaled>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}

void thoracic_iso_ve_simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<ThoracicIsoVE>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}

void thoracic_defaultfung_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<ThoracicDefaultFungVEScaled>(
        pars, fiber, caputo, Tf, Cmax, args, dt, stress, n
    );
}

void thoracic_iso_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<ThoracicIsoVEScaled>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}

void thoracic_lin_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<ThoracicLinVEScaled>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}

void thoracic_quad_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    residuals::simulate<ThoracicQuadVEScaled>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
}

} // namespace sim