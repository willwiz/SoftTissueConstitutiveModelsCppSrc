#pragma once

namespace thoracic {

// Planar Elastin Model
void thoracic_elastin_get_parameters(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    double pars_out[2]
);

void thoracic_elastin_simulate(
    const double pars[], const double fiber[], const double caputo[], double Tf,
    const double Cmax[], const double args[], const double dt[], double stress[], int n
);

double thoracic_elastin_residual(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    const double args[], const double stress[], const double dt[], const double weights[],
    const double deltaCG[], const double hysteresis[], const double alphas[], const int index[],
    const int select[], int n, int nprot, int skip
);

// Ensemble Model

void thoracic_ensemble_get_parameters(
    const double pars[], const double visco[], double Tf, const double Cmax[], double pars_out[6]
);

void thoracic_ensemble_simulate(
    const double pars[], const double visco[], double Tf, const double Cmax[],
    const double strain[], const double dt[], double out_stress[], int n
);

double thoracic_ensemble_residual(
    const double pars[], const double visco[], double Tf, const double Cmax[],
    const double strain[], const double stress[], const double dt[], const double deltaCG[],
    double hysteresis, const double data[], int n, int skip
);

// Full Model

void thoracic_ve_get_parameters(const double pars[10], const double Cmax[], double pars_out[10]);

void thoracic_ve_simulate(
    const double pars[], const double fiber[], const double caputo[], double Tf,
    const double Cmax[], const double args[], const double dt[], double stress[], int n
);

double thoracic_ve_residual(
    const double pars[], const double fiber[], const double visco[], double Tf, const double Cmax[],
    const double args[], const double stress[], const double dt[], const double weights[],
    const double deltaCG[], const double hysteresis[], const double alphas[], const int index[],
    const int select[], int n, int nprot, int skip
);

} // namespace thoracic
