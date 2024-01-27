#pragma once

#include "../interfaces.hpp"

namespace templates {

typedef double (*PenaltyBodyCalculation)(double *, double *, double *);
typedef double (*ViscoBodyCalculation)(double *, double *, double);

template <class matlaw>
void simulate(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
);
} // namespace templates