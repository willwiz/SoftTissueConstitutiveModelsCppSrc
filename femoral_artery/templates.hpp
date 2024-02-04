#pragma once

#include "../interfaces.hpp"

namespace residuals {

typedef double (*ResidualNorm)(double, double, double);
typedef double (*ResidualCalculation)(double *, double *, double *, int *, int *, int, int, int, double *);
typedef double (*HysteresisCalculation)(
    double *, double *, double *, double *, int *, int *, int, int, int
);
typedef double (*PenaltyCalculation)(double *, double *, double *);
typedef double (*ViscoCalculation)(double *, double *, double);

template <class matlaw>
void simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
);

template <
    class matlaw, ResidualCalculation resfunc, HysteresisCalculation hystfunc,
    PenaltyCalculation penfunc, ViscoCalculation viscofunc>
double calc_residual_general(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double args[],
    double stress[], double dt[], double weights[], double deltaCG[], double hysteresis[],
    double alphas[], int index[], int select[], int n, int dim, int nprot, int skip
);

// main loops

inline double top_sided_residual(double sim, double data, double strain);

inline double quadratic_residual(double sim, double data, double strain);

double residual_body(
    double strain[], double stress[], double weights[], int index[], int select[], int dim,
    int nprot, int skip, double sims[]
);

double residual_body_below(
    double strain[], double stress[], double weights[], int index[], int select[], int dim,
    int nprot, int skip, double sims[]
);

double residual_body_below2(
    double strain[], double stress[], double weights[], int index[], int select[], int dim,
    int nprot, int skip, double sims[]
);

double hysteresis_body_null(
    double sims[], double deltaCG[], double hysteresis[], double weights[], int index[],
    int select[], int dim, int nprot, int skip
);
double hysteresis_body(
    double sims[], double deltaCG[], double hysteresis[], double weights[], int index[],
    int select[], int dim, int nprot, int skip
);

double penalty_body_null(double pars[], double fiber[], double visco[]);
double penalty_body_default(double pars[], double fiber[], double visco[]);
double penalty_body_1(double pars[], double fiber[], double visco[]);
double penalty_body_4(double pars[], double fiber[], double visco[]);
double penalty_body_1_2(double pars[], double fiber[], double visco[]);
double penalty_body_1_3(double pars[], double fiber[], double visco[]);
double penalty_body_2_4(double pars[], double fiber[], double visco[]);
double penalty_body_kappa(double pars[], double fiber[], double visco[]);
double penalty_body_kappa_2(double pars[], double fiber[], double visco[]);
double penalty_body_kappa_3(double pars[], double fiber[], double visco[]);
double penalty_body_kappa_1(double pars[], double fiber[], double visco[]);

double calculate_viscopart_body_null(double visco[], double datas[], double alpha);
double calculate_viscopart_body(double visco[], double datas[], double alpha);
double calculate_viscopart_body1(double visco[], double datas[], double alpha);
double calculate_viscopart_body2(double visco[], double datas[], double alpha);
double calculate_viscopart_body3(double visco[], double datas[], double alpha);

// Residual templates

template <class matlaw, ResidualNorm norm_func>
double calc_residual(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double args[],
    double stress[], double dt[], double weights[], double deltaCG[], double hysteresis[],
    double alphas[], int index[], int select[], int n, int dim, int nprot, int skip
);

template <class matlaw, ResidualNorm norm_func, PenaltyCalculation penalty>
double calc_residual_hyst(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double args[],
    double stress[], double dt[], double weights[], double deltaCG[], double hysteresis[],
    double alphas[], int index[], int select[], int n, int dim, int nprot, int skip
);

template <class matlaw, PenaltyCalculation penalty, ViscoCalculation viscoform>
double calc_residual_hyst_relax(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double args[],
    double stress[], double dt[], double weights[], double deltaCG[], double hysteresis[],
    double alphas[], int index[], int select[], int n, int dim, int nprot, int skip
);

} // namespace residuals