#pragma once

namespace thoracic {

void thoracic_default_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[]
);

void planar_elastin_matrix_simulate(
    const double pars[], const double fiber[], const double caputo[], double Tf,
    const double Cmax[], const double args[], const double dt[], double stress[], int n
);

void thoracic_default_ve_simulate_standard(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void thoracic_default_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

} // namespace thoracic