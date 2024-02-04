#include "objective_template.hpp"
#include "../CTvalues_optimization.hpp"
#include <algorithm>
#include <cmath>

namespace residuals {

/* ------------------------------------------------------------------------------
 |  This file provides the definitions for calculating the residuals for
 |  the different model forms
 |
 |  These models are mainly used for the collaboration with Alexey Kamenskiy
 |
 |  Author: Will Zhang
 |  Dependencies: None
 ----------------------------------------------------------------------------- */

template <class matlaw>
void simulate(
    double pars[], double fiber[], double caputo[], double Tf, double Cmax[], double args[],
    double dt[], double stress[], int n
) {
    int strd_i;

    double vals[ctv::prob_dim];

    kinematics::deformation2D kin;
    matlaw law(pars, fiber, caputo, Tf, Cmax);

    for (int i = 1; i < n; i++) {
        strd_i = ctv::prob_dim * i;
        // Calculate deformation
        kin.precompute(&args[strd_i]);
        // Compute Final Stress
        law.stress(kin, dt[i], vals);
        for (int j = 0; j < ctv::prob_dim; j++) {
            stress[strd_i + j] = vals[j];
        }
    }
}

template <
    class matlaw, ResidualCalculation resfunc, HysteresisCalculation hystfunc,
    PenaltyCalculation penfunc, ViscoCalculation viscofunc>
double calc_residual_general(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double args[],
    double stress[], double dt[], double weights[], double deltaCG[], double hysteresis[],
    double alphas[], int index[], int select[], int n, int dim, int nprot, int skip
) {
    double *sims = new double[n * dim]();

    simulate<matlaw>(pars, fiber, visco, Tf, Cmax, args, dt, &sims[0], n);

    double res = resfunc(args, stress, weights, index, select, dim, nprot, skip, sims);

    double hyst = hystfunc(sims, deltaCG, hysteresis, weights, index, select, dim, nprot, skip);

    double alp = viscofunc(visco, alphas, fiber[0]);

    double scale =
        1.0 / std::max(weights[dim * select[nprot - 1]], weights[dim * select[nprot - 1] + 3]);

    delete[] sims;

    return (res + ctv::w_hyst * scale * hyst) * (penfunc(pars, fiber, visco) + ctv::w_visco * alp);
}

/* ******************************************************************************
 * Basic functions, e.g. calculating the residual, penalty, hysteresis etc.
 *
 * COMMENTS:
 * This is the same for all of the models, so they share these same codes
 ****************************************************************************** */
// Residual Calculation

inline double top_sided_residual(double sim, double data, double strain) {
    double difference = sim - data;
    difference = difference * difference;
    double weight = (difference > 0) ? difference : 1.0;
    return weight * difference;
}

inline double quadratic_residual(double sim, double data, double strain) {
    double difference = sim - data;
    return difference * difference;
}

inline double cubic_residual(double sim, double data, double strain) {
    double difference = sim - data;
    return std::abs(difference) * difference * difference;
}

inline double quartic_residual(double sim, double data, double strain) {
    double difference = sim - data;
    return difference * difference * difference * difference;
}

template <ResidualNorm norm_func>
double residual_general(
    double strain[], double stress[], double weights[], int index[], int select[], int dim,
    int nprot, int skip, double sims[]
) {
    int kid, strd_i, strd_k;
    // double ds;
    double res, eps;

    res = 0;
    for (int k = 0; k < nprot; k++) {
        kid = select[k];
        strd_k = dim * kid;
        for (int j = 0; j < dim; j += 3) {
            eps = 0;
            for (int i = index[kid]; i < index[kid + 1] + 1; i += skip) {
                strd_i = i * dim + j;
                eps = eps + norm_func(sims[strd_i], stress[strd_i], strain[strd_i]);
            }
            res = res + weights[strd_k + j] * eps;
        }
    }
    return res;
}

double residual_body(
    double strain[], double stress[], double weights[], int index[], int select[], int dim,
    int nprot, int skip, double sims[]
) {

    int strd_i, strd_k;
    double ds, res, eps;
    int kid, start, stop;

    res = 0;
    for (int k = 0; k < nprot; k++) {
        kid = select[k];
        strd_k = dim * kid;
        start = index[kid];
        stop = index[kid + 1] + 1;
        for (int j = 0; j < dim; j += 3) {
            eps = 0;

            for (int i = start; i < stop; i += skip) {
                strd_i = i * dim + j;
                ds = sims[strd_i] - stress[strd_i];
                eps = eps + ds * ds;
            }
            res = res + weights[strd_k + j] * eps;
            // res = res + eps;
        }
    }
    return res;
}

// Residual Calculation
double residual_body_below(
    double strain[], double stress[], double weights[], int index[], int select[], int dim,
    int nprot, int skip, double sims[]
) {

    int strd_i;
    // int strd_k;
    double ds, res, eps;
    int kid, start, stop;

    res = 0;
    for (int k = 0; k < nprot; k++) {
        kid = select[k];
        // strd_k = dim * kid;
        start = index[kid];
        stop = index[kid + 1];
        // stop   = std::min(index[kid] + 150, index[kid + 1] + 1);
        for (int j = 0; j < dim; j += 3) {
            eps = 0;

            for (int i = start; i < stop; i += skip) {
                strd_i = i * dim + j;
                ds = sims[strd_i] - stress[strd_i];
                if (ds > 0.0) {
                    ds = ds * ds * ds * ds;
                } else {
                    ds = ds * ds;
                }
                // eps = eps + weights[strd_k + j] * ds / std::pow(strain[k], 4);
                eps = eps + ds / std::pow(strain[k], 4);
            }

            res = res + eps;
        }
    }
    return res;
}

double residual_body_below2(
    double strain[], double stress[], double weights[], int index[], int select[], int dim,
    int nprot, int skip, double sims[]
) {

    int strd_i;
    double ds, res, eps;
    int kid, start, stop;

    res = 0;
    for (int k = 0; k < nprot; k++) {
        kid = select[k];
        start = index[kid];
        stop = index[kid + 1];
        // stop   = std::min(index[kid] + 150, index[kid + 1] + 1);
        for (int j = 0; j < dim; j += 3) {
            eps = 0;
            for (int i = start; i < stop; i += skip) {
                strd_i = i * dim + j;
                ds = sims[strd_i] - stress[strd_i];
                if (ds > 0.0) {
                    ds = 10.0 * ds * ds * ds * ds;
                } else {
                    ds = ds * ds;
                }
                eps = eps + ds / std::pow(strain[k], 4);
            }

            res = res + eps;
        }
    }
    return res;
}

double hysteresis_body_null(
    double sims[], double deltaCG[], double hysteresis[], double weights[], int index[],
    int select[], int dim, int nprot, int skip
) {
    return 0.0;
}

// Hysteresis Calculation
double hysteresis_body(
    double sims[], double deltaCG[], double hysteresis[], double weights[], int index[],
    int select[], int dim, int nprot, int skip
) {

    int strd_i, strd_k;
    double ds, res, hyst;
    int kid, start, stop;

    res = 0;
    for (int k = 0; k < nprot; k++) {
        kid = select[k];
        strd_k = dim * kid;
        start = index[kid];
        stop = index[kid + 1] + 1;
        for (int j = 0; j < dim; j += 3) {
            hyst = 0;

            for (int i = start; i < stop; i += skip) {
                strd_i = i * dim + j;
                hyst = hyst + sims[strd_i] * deltaCG[strd_i];
            }
            ds = hyst - hysteresis[strd_k + j];
            // res = res + weights[strd_k + j] * ds*ds;
            res = res + ds * ds;
        }
    }

    return res;
}
/* ******************************************************************************
 * Basic functions, e.g. calculating the residual, penalty, hysteresis etc.
 *
 * COMMENTS:
 * This is the same for all of the models, so they share these same codes
 ****************************************************************************** */

double penalty_body_null(double pars[], double fiber[], double visco[]) {
    return 1.0;
}

double penalty_body_default(double pars[], double fiber[], double visco[]) {
    double d_alpha = fiber[1] - ctv::ideal_alpha;

    return 1.0 + ctv::p_fiber * fiber[0] * fiber[0] + ctv::p_alpha * d_alpha * d_alpha;
}

double penalty_body_4(double pars[], double fiber[], double visco[]) {
    double d_alpha = fiber[1] - ctv::ideal_alpha;

    return 1.0 + ctv::p_fiber * fiber[0] * fiber[0] + ctv::p_alpha * d_alpha * d_alpha +
           ctv::p_elastin * pars[4] * pars[4];
}

double penalty_body_1(double pars[], double fiber[], double visco[]) {
    double d_alpha = fiber[1] - 0.5;

    return 1.0 + ctv::p_fiber * fiber[0] * fiber[0] + ctv::p_alpha * d_alpha * d_alpha;
}

double penalty_body_kappa(double pars[], double fiber[], double visco[]) {
    double d_kappa = fiber[1] - 0.5;

    return 1.0 + 0.1 * fiber[0] * fiber[0] + 0.1 * d_kappa * d_kappa;
}

double penalty_body_kappa_2(double pars[], double fiber[], double visco[]) {

    double d_kappa2 = fiber[2] - 0.5;
    double d_elastin = pars[1] - fiber[3];
    double d_fiber = fiber[0] - fiber[4];
    double d_kappa = fiber[1] - fiber[5];

    return 1.0 + 1.0 * d_fiber * d_fiber + 0.1 * d_elastin * d_elastin + 40.0 * d_kappa * d_kappa +
           10.0 * (d_kappa2 * d_kappa2);
}

// double penalty_body_kappa_3(double pars[], double fiber[], double visco[])
// {

//   double d_kappa = fiber[1] - 0.5;
//   double d_kappa2 = fiber[2] - 0.5;
//   double d_alpha = fiber[3] - ctv::ideal_alpha;

//   return 1.0 + 0.01 * fiber[0] * fiber[0] + 0.1 *d_kappa * d_kappa + 1.0 * (d_kappa2*d_kappa2)
//   + 1.0 * d_alpha * d_alpha;
// }

double penalty_body_kappa_3(double pars[], double fiber[], double visco[]) {

    // double d_bs = pars[1] - fiber[1];
    // double d_bc = pars[2] - fiber[2];
    // double d_kappa = fiber[1] - 0.5;
    double d_kappa2 = fiber[0] - 0.5;
    double d_alpha = fiber[1] - ctv::ideal_alpha;

    return 1.0 + 1.0 * (d_kappa2 * d_kappa2) + 1.0 * d_alpha * d_alpha;
}

double penalty_body_kappa_1(double pars[], double fiber[], double visco[]) {
    double d_alpha = fiber[1] - ctv::ideal_alpha;
    double d_kappa = fiber[2] - 0.5;

    return 1.0 + ctv::p_fiber * fiber[0] * fiber[0] + 0.001 * d_kappa * d_kappa +
           ctv::p_alpha * d_alpha * d_alpha + ctv::p_collagen * pars[1] * pars[1];
}

double penalty_body_1_3(double pars[], double fiber[], double visco[]) {
    double d_alpha = fiber[1] - ctv::ideal_alpha;
    double d_weight = pars[1] - pars[3];

    return 1.0 + ctv::p_fiber * fiber[0] * fiber[0] + ctv::p_alpha * d_alpha * d_alpha +
           ctv::p_elastin * d_weight * d_weight;
}

double penalty_body_2_4(double pars[], double fiber[], double visco[]) {
    double d_kappa = fiber[1] - 0.5;
    double d_alpha = fiber[2] - ctv::ideal_alpha;
    double d_weight = pars[2] - pars[4];
    double d_elastin = pars[1] - fiber[3];

    return 1.0 + ctv::p_fiber * fiber[0] * fiber[0] + 0.1 * d_alpha * d_alpha +
           2.0 * d_kappa * d_kappa + 0.001 * d_elastin * d_elastin + 0.001 * d_weight * d_weight;
}

double penalty_body_1_2(double pars[], double fiber[], double visco[]) {

    double d_alpha = fiber[1] - ctv::ideal_alpha;
    double d_weight = pars[1] - pars[2];

    return 1.0 + ctv::p_fiber * fiber[0] * fiber[0] + ctv::p_alpha * d_alpha * d_alpha +
           ctv::p_elastin * d_weight * d_weight;
}

double calculate_viscopart_body_null(double visco[], double datas[], double alpha) {
    return 0.0;
}

double calculate_viscopart_body1(double visco[], double datas[], double alpha) {
    double cosa = cos(alpha);
    double sina = sin(alpha);
    double delta_alpha_L = 0.5 * (cosa * visco[0]) - datas[0];
    double delta_alpha_C = 0.5 * (sina * visco[0]) - datas[1];

    return delta_alpha_C * delta_alpha_C + delta_alpha_L * delta_alpha_L;
}

// Relaxation penalties
double calculate_viscopart_body(double visco[], double datas[], double alpha) {

    double delta_alpha_L = visco[0] - datas[0];
    double delta_alpha_C = 0.5 * (visco[0] + visco[1]) - datas[1];

    return delta_alpha_L * delta_alpha_L + delta_alpha_C * delta_alpha_C;
}

double calculate_viscopart_body2(double visco[], double datas[], double alpha) {
    double cosa = cos(alpha);
    double sina = sin(alpha);

    double delta_alpha_L = 0.5 * (cosa * visco[0] + visco[1]) - datas[0];
    double delta_alpha_C = 0.5 * (sina * visco[0] + visco[1]) - datas[1];

    return delta_alpha_C * delta_alpha_C + delta_alpha_L * delta_alpha_L;
}

double calculate_viscopart_body3(double visco[], double datas[], double alpha) {
    // double cosa = cos(alpha);
    // double sina = sin(alpha);

    // double delta_alpha_L = 0.5 * cosa * (visco[0] + visco[1]) - datas[0];
    // double delta_alpha_C = 0.5 * sina * (visco[0] + visco[1]) - datas[1];

    double delta = visco[0] - visco[1];

    return delta * delta;
}

double lowerbound(double var, double value) {

    double rat = var / value;

    return 1.0 / (exp(rat) - exp(-rat));
}

/*******************************************************************************
 * Template residual functions
 *
 * COMMENTS:
 * The main form is predetermined, the forms differs by whether the model is
 * scaled.
 *******************************************************************************/

template <class matlaw, ResidualNorm norm_func>
double calc_residual(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double args[],
    double stress[], double dt[], double weights[], double deltaCG[], double hysteresis[],
    double alphas[], int index[], int select[], int n, int dim, int nprot, int skip
) {
    double *sims = new double[n * dim]();

    simulate<matlaw>(pars, fiber, visco, Tf, Cmax, args, dt, &sims[0], n);

    double res =
        residual_general<norm_func>(args, stress, weights, index, select, dim, nprot, skip, sims);

    delete[] sims;

    return res * (1.0 + ctv::p_fiber * fiber[0] * fiber[0]);
}

template <class matlaw, ResidualNorm norm_func, PenaltyCalculation penalty>
double calc_residual_hyst(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double args[],
    double stress[], double dt[], double weights[], double deltaCG[], double hysteresis[],
    double alphas[], int index[], int select[], int n, int dim, int nprot, int skip
) {
    double *sims = new double[n * dim]();

    simulate<matlaw>(pars, fiber, visco, Tf, Cmax, args, dt, &sims[0], n);

    double res =
        residual_general<norm_func>(args, stress, weights, index, select, dim, nprot, skip, sims);

    double hyst =
        hysteresis_body(sims, deltaCG, hysteresis, weights, index, select, dim, nprot, skip);

    double scale =
        1.0 / std::max(weights[dim * select[nprot - 1]], weights[dim * select[nprot - 1] + 3]);

    delete[] sims;

    return (res + ctv::w_hyst * scale * hyst) * penalty(pars, fiber, visco);
}

template <class matlaw, PenaltyCalculation penalty, ViscoCalculation viscoform>
double calc_residual_hyst_relax(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double args[],
    double stress[], double dt[], double weights[], double deltaCG[], double hysteresis[],
    double alphas[], int index[], int select[], int n, int dim, int nprot, int skip
) {
    double *sims = new double[n * dim]();

    simulate<matlaw>(pars, fiber, visco, Tf, Cmax, args, dt, &sims[0], n);

    double res = residual_body(args, stress, weights, index, select, dim, nprot, skip, sims);

    double hyst =
        hysteresis_body(sims, deltaCG, hysteresis, weights, index, select, dim, nprot, skip);

    double alp = viscoform(visco, alphas, fiber[1]);

    double scale =
        1.0 / std::max(weights[dim * select[nprot - 1]], weights[dim * select[nprot - 1] + 3]);

    delete[] sims;

    return (res + ctv::w_hyst * scale * hyst) * (penalty(pars, fiber, visco) + ctv::w_visco * alp);
}

// template<class matlaw, PenaltyCalculation penalty, ViscoCalculation viscoform>
// double calc_residual_hyst_relax(double pars[], double fiber[],
//   double visco[], double Tf, double Cmax[],
//   double args[], double stress[], double dt[], double weights[],
//   double deltaCG[], double hysteresis[], double alphas[],
//   int index[], int select[], int n, int dim, int nprot, int skip)
// {
//   double * sims = new double[n*dim]();

//   simulate<matlaw>(pars, fiber, visco, Tf, Cmax, args, dt, &sims[0], n);

//   double res = residual_body(args, stress, weights, index, select, dim, nprot, skip, sims);

//   double hyst = hysteresis_body(sims, deltaCG, hysteresis, weights, index,
//     select, dim, nprot, skip);

//   double alp  = viscoform(visco, alphas, fiber[1]);

//   double scale = 1.0 /std::max(weights[dim*select[nprot - 1]], weights[dim*select[nprot - 1] +
//   3]);

//   delete [] sims;

//   return (res + ctv::w_hyst * scale * hyst) * (penalty(pars, fiber, visco) + ctv::w_visco * alp);
// }

// template<class matlaw>
// double calc_residual_hyst_relax2(double pars[], double fiber[],
//   double visco[], double Tf, double Cmax[],
//   double args[], double stress[], double dt[], double weights[],
//   double deltaCG[], double hysteresis[], double alphas[],
//   int index[], int select[], int n, int dim, int nprot, int skip)
// {
//   double * sims = new double[n*dim]();

//   simulate<matlaw>(pars, fiber, visco, Tf, Cmax, args, dt, &sims[0], n);

//   double res = residual_body(stress, weights, index, select, dim, nprot, skip, sims);

//   double hyst = hysteresis_body(sims, deltaCG, hysteresis, weights, index,
//     select, dim, nprot, skip);

//   double alp  = calculate_viscopart_body2(visco, alphas, fiber[1]);

//   double scale = 1.0 /std::max(weights[dim*select[nprot - 1]], weights[dim*select[nprot - 1] +
//   3]);

//   delete [] sims;

//   return (res + ctv::w_hyst * scale * hyst) * (penalty_body(pars, fiber, visco) + ctv::w_visco *
//   alp);
// }

// template<class matlaw>
// double calc_residual_hyst_relax3(double pars[], double fiber[],
//   double visco[], double Tf, double Cmax[],
//   double args[], double stress[], double dt[], double weights[],
//   double deltaCG[], double hysteresis[], double alphas[],
//   int index[], int select[], int n, int dim, int nprot, int skip)
// {
//   double * sims = new double[n*dim]();

//   simulate<matlaw>(pars, fiber, visco, Tf, Cmax, args, dt, &sims[0], n);

//   double res = residual_body(stress, weights, index, select, dim, nprot, skip, sims);

//   double hyst = hysteresis_body(sims, deltaCG, hysteresis, weights, index,
//     select, dim, nprot, skip);

//   double alp  = calculate_viscopart_body3(visco, alphas, fiber[1]);

//   double scale = 1.0 /std::max(weights[dim*select[nprot - 1]], weights[dim*select[nprot - 1] +
//   3]);

//   delete [] sims;

//   return (res + ctv::w_hyst * scale * hyst) * (penalty_body(pars, fiber, visco) + ctv::w_visco *
//   alp);
// }
/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/
} // namespace residuals
