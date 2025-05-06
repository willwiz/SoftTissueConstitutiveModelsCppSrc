#pragma once
#include "../../interfaces.hpp"

namespace optimization_1d {

inline double quadratic_residual(double sim, double data, double strain);

typedef double (*ResidualNorm)(double, double, double);
typedef double (*PenaltyFunction)(double *, double *, double *);

double penalty_function_null(double pars[], double fiber[], double visco[]);

template <class matlaw>
void simulate_general(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double dt[], double stress_out[], int n
);

template <ResidualNorm norm_func>
double residual_body_general(
    double strain[], double stress[], double weight[], int nset, int index[], int select[],
    double sims[]
);

template <class matlaw, ResidualNorm norm_func, PenaltyFunction pen_func>
double calc_objective_general(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double stress[], double dt[], double weight[], int n, int nset, int index[],
    int select[]
);

} // namespace optimization_1d