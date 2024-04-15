#pragma once

namespace femoral {

void femoral_get_model_parameters_scaled(
    const double pars[10], const double Cmax[], double pars_out[10]
);

void femoral_ve_simulate_scaled(
    const double pars[], const double fiber[], const double caputo[], double Tf,
    const double Cmax[], const double args[], const double dt[], double stress[], int n
);

double femoral_residual_VE_scaled_hyst_relax(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    const double args[], const double stress[], const double dt[], const double weights[],
    const double deltaCG[], const double hysteresis[], const double data[], const int index[],
    const int select[], int n, int nprot, int skip
);

} // namespace femoral