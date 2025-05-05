#pragma once
#include "../../interfaces.hpp"

namespace ensemble {

typedef double (*ResidualFunction)(double *, double *, int, int, double *, double);
typedef double (*HysteresisFunction)(double *, double *, int, double);
typedef double (*PenaltyFunction)(double *, double *);

// template<class matlaw>
// void ensemble_simulate(
//     double pars[], double visco[], double Tf, double Cmax[],
//     double strain[], double dt[], double stress[], int n);

template <
    class matlaw, ResidualFunction resfunc, HysteresisFunction hystfunc, PenaltyFunction penfunc>
double calc_ensemble_objective(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double stress[],
    double dt[], double deltaCG[], double hysteresis, int n, int skip
);

double residual_body_below(
    double args[], double stress[], int n, int skip, double sims[], double baseline
);

double penalty_function_null(double pars[], double visco[]);

double penalty_body_ensemble2(double pars[], double visco[]);

void thoracic_smc_ensemble_simulate(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double dt[],
    double out_stress[], int n
);

void thoracic_full_ensemble_simulate(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double dt[],
    double out_stress[], int n
);

double thoracic_smc_ensemble_residual(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double stress[],
    double dt[], double deltaCG[], double hysteresis, int n, int skip
);

double thoracic_full_ensemble_residual(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double stress[],
    double dt[], double deltaCG[], double hysteresis, int n, int skip
);

} // namespace ensemble