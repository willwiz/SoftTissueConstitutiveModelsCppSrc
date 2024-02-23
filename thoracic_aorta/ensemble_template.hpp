#pragma once

#include "../interfaces.hpp"

namespace thoracic {

typedef double (*NormFunction)(double, double, double);
typedef double (*HysteresisFunction)(const double *, const double *, int, double);
typedef double (*PenaltyFunction)(const double *, const double *, const double *);

template <class matlaw>
void ensemble_simulate(
    const double pars[], const double visco[], double Tf, const double Cmax[],
    const double strain[], const double dt[], double stress[], int n
);

template <class matlaw, NormFunction nfunc, HysteresisFunction hfunc, PenaltyFunction pfunc>
double calc_ensemble_objective(
    const double pars[], const double visco[], double Tf, const double Cmax[],
    const double strain[], const double stress[], const double dt[], const double deltaCG[],
    double hysteresis, const double data[], int n, int skip
);

template <NormFunction nfunc>
double ensemble_residual(
    const double strain[], const double stress[], int n, int skip, const double sims[],
    double baseline
);

inline double quart_quad_residual(double sim, double data, double strain);
double hysteresis_body(const double sims[], const double deltaCG[], int n, double hysteresis);
double penalty_ensemble3(const double pars[], const double visco[], const double data[]);

} // namespace thoracic