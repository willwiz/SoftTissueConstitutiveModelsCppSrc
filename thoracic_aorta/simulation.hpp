#pragma once

namespace sim {

void thoracic_default_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[]
);

void planar_elastin_matrix_simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void thoracic_default_ve_simulate_standard(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void thoracic_default_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

} // namespace sim