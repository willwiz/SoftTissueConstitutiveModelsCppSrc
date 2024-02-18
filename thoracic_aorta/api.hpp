#pragma once

namespace thoracic {

// Planar Elastin Model
void planar_elastin_matrix_get_model_parameters_scaled(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    double pars_out[2]
);

void planar_elastin_matrix_simulate(
    const double pars[], const double fiber[], const double caputo[], double Tf,
    const double Cmax[], const double args[], const double dt[], double stress[], int n
);

double planar_elastin_matrix_residual(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    const double args[], const double stress[], const double dt[], const double weights[],
    const double deltaCG[], const double hysteresis[], const double alphas[], const int index[],
    const int select[], int n, int nprot, int skip
);

// Ensemble Model

void thoracic_ensemble_get_model_parameters_scaled(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    double pars_out[2]
);

void thoracic_ensemble_simulate(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double dt[],
    double out_stress[], int n
);

double thoracic_ensemble_residual(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double stress[],
    double dt[], double deltaCG[], double hysteresis, int n, int skip
);

} // namespace thoracic
