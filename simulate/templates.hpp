#pragma once

#include "../interfaces.hpp"

namespace simulate {

typedef double (*ResidualNorm)(const double, const double, const double);
typedef double (*HysteresisFunc)(
    const double *, const double *, const double *, const double *, const int *, const int *,
    const int, const int
);
typedef double (*Penalty)(const double *, const double *, const double *, const double *);
typedef double (*ViscoCalculation)(const double *, const double *, const double);

// main loops
template <class matlaw, class kine, int dim>
void simulate(
    const double pars[], const double fiber[], const double caputo[], const double Tf,
    const double Cmax[], const double args[], const double dt[], double stress[], const int n
);

template <
    class matlaw, class kine, ResidualNorm normfunc, HysteresisFunc hystfunc, Penalty penfunc,
    int dim>
double calc_residual(
    const double pars[], const double fiber[], const double visco[], const double Tf,
    const double Cmax[], const double args[], const double stress[], const double dt[],
    const double weights[], const double deltaCG[], const double hysteresis[], const double data[],
    const int index[], const int select[], const int n, const int nprot, const int skip,
    const double w_hyst
);

template <ResidualNorm norm_func, int dim>
double residual_general(
    const double strain[], const double stress[], const double weights[], const int index[],
    const int select[], const int nprot, const int skip, const double sims[]
);

inline double quart_quad_residual(const double sim, const double data, const double strain);

inline double quadratic_residual(const double sim, const double data, const double strain);

double hysteresis_body_null(
    const double sims[], const double deltaCG[], const double hysteresis[], const double weights[],
    const int index[], const int select[], const int nprot, const int skip
);

template <int dim>
double hysteresis_body(
    const double sims[], const double deltaCG[], const double hysteresis[], const double weights[],
    const int index[], const int select[], const int nprot, const int skip
);

} // namespace simulate