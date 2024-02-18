#include "api.hpp"
#include "../simulate/templates.hpp"
#include "femoral_model.hpp"

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

namespace femoral {

/*----------------------------------------------------------------------
|  Some basic models
-----------------------------------------------------------------------*/
double penalty_body_4(
    const double pars[], const double fiber[], const double visco[], const double data[]
) {
    double d_alpha = fiber[1] - M_ideal_alpha;
    double delta_alpha_L = visco[0] - data[0];
    double delta_alpha_C = 0.5 * (visco[0] + visco[1]) - data[1];
    double frac_alpha_weight = delta_alpha_L * delta_alpha_L + delta_alpha_C * delta_alpha_C;
    return (
        M_p_fiber * fiber[0] * fiber[0] + M_p_alpha * d_alpha * d_alpha +
        M_p_elastin * pars[4] * pars[4] + M_w_visco * frac_alpha_weight
    );
}

/*----------------------------------------------------------------------
 |  This provides the main models in the full constitutive model
 -----------------------------------------------------------------------*/

void femoral_get_model_parameters_scaled(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    double pars_out[11]
) {
    FemoralBase psi(pars, fiber, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[6] = fiber[0];
    pars_out[7] = fiber[1];
    pars_out[8] = visco[0];
    pars_out[9] = visco[1];
}

void femoral_ve_simulate_scaled(
    const double pars[], const double fiber[], const double caputo[], double Tf,
    const double Cmax[], const double args[], const double dt[], double stress[], int n
) {
    simulate::simulate<FemoralVEScaled, kinematics::deformation2D, 2>(
        pars, fiber, caputo, Tf, Cmax, args, dt, stress, n
    );
}

double femoral_residual_VE_scaled_hyst_relax(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    const double args[], const double stress[], const double dt[], const double weights[],
    const double deltaCG[], const double hysteresis[], const double data[], const int index[],
    const int select[], int n, int nprot, int skip
) {
    return simulate::calc_residual<
        FemoralVEScaled, kinematics::deformation2D, simulate::quadratic_residual,
        simulate::hysteresis_body<2>, penalty_body_4, 2>(
        pars, fiber, visco, Tf, Cmax, args, stress, dt, weights, deltaCG, hysteresis, data, index,
        select, n, nprot, skip, M_w_hyst
    );
}

/*----------------------------------------------------------------------
|  The femoral artery models
-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
|  The thoracic models
-----------------------------------------------------------------------*/

} // namespace femoral