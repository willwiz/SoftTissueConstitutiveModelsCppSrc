#pragma once

#include "../interfaces.hpp"

namespace thoracic {

typedef double (*ResidualFunction)(
    const double *, const double *, int, int, const double *, double
);
typedef double (*HysteresisFunction)(const double *, const double *, int, double);
typedef double (*PenaltyFunction)(const double *, const double *);

template <class matlaw>
void ensemble_simulate(
    const double pars[], const double visco[], double Tf, const double Cmax[],
    const double strain[], const double dt[], double stress[], int n
);

template <
    class matlaw, ResidualFunction resfunc, HysteresisFunction hystfunc, PenaltyFunction penfunc>
double calc_ensemble_objective(
    const double pars[], const double visco[], double Tf, const double Cmax[],
    const double strain[], const double stress[], const double dt[], const double deltaCG[],
    double hysteresis, int n, int skip
);

} // namespace thoracic