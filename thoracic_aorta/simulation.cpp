#include "simulation.hpp"
#include "../simulate/templates.hpp"
#include "thoracic_default_model.hpp"

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

namespace thoracic {

/*----------------------------------------------------------------------
 |  This provides the main models in the full constitutive model
 -----------------------------------------------------------------------*/

void thoracic_default_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[]
) {
    ThoracicDefaultBase psi(pars, fiber, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[6] = fiber[0];
    pars_out[7] = fiber[1];
    pars_out[8] = visco[0];
    pars_out[9] = visco[1];
}

/*----------------------------------------------------------------------
|  Some basic models
-----------------------------------------------------------------------*/

// void planar_elastin_matrix_simulate(
//     double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
//     double dt[], double stress[], int n
// ) {
//     simulate::simulate<PlanarElastinMatrix>(pars, fiber, caputo, Tf, Cmax, args, dt, stress, n);
// }

// void thoracic_default_ve_simulate_standard(
//     const double pars[], const double fiber[], const double caputo[], double Tf,
//     const double Cmax[], const double args[], const double dt[], double stress[], int n
// ) {
//     residuals::simulate<ThoracicDefaultVEStandard>(
//         pars, fiber, caputo, Tf, Cmax, args, dt, stress, n
//     );
// }

// void thoracic_default_ve_simulate_scaled(
//     const double pars[], const double fiber[], const double caputo[], double Tf,
//     const double Cmax[], const double args[], const double dt[], double stress[], int n
// ) {
//     residuals::simulate<ThoracicDefaultVEScaled>(
//         pars, fiber, caputo, Tf, Cmax, args, dt, stress, n
//     );
// }

} // namespace thoracic