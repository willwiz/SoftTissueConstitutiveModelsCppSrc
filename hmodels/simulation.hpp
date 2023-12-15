#pragma once

namespace sim {

void femoral_get_model_parameters(
    double pars[], double fiber[], double visco[], double Tf, double pars_out[11]
);

void femoral_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[11]
);

void thoracic_iso_get_model_parameters(
    double pars[], double fiber[], double visco[], double Tf, double pars_out[9]
);

void thoracic_default_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[]
);

void thoracic_defaultfung_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[9]
);

void thoracic_iso_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[9]
);

void thoracic_lin_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[8]
);

void thoracic_quad_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[8]
);

void thoracic_circ_get_model_parameters_scaled(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double pars_out[]
);

void planar_elastin_matrix_simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void femoral_he_simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void femoral_he_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void femoral_ve_simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void femoral_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void thoracic_iso_he_simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void thoracic_iso_he_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void thoracic_iso_ve_simulate(
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

void thoracic_defaultfung_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void thoracic_iso_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void thoracic_lin_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void thoracic_quad_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

void thoracic_circ_ve_simulate_scaled(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

} // namespace sim