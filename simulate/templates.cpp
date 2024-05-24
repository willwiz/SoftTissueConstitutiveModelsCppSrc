#include "templates.hpp"
#include <algorithm>
#include <cmath>

namespace simulate {

template <class matlaw, class kine, int dim>
void simulate(
    const double pars[], const double fiber[], const double caputo[], double Tf,
    const double Cmax[], const double args[], const double dt[], double stress[], int n
) {
    constexpr int size = dim * dim;
    int strd_i;
    double vals[size];

    kine kin;
    matlaw law(pars, fiber, caputo, Tf, Cmax);

    for (int i = 1; i < n; i++) {
        strd_i = size * i;
        // Calculate deformation
        kin.precompute(&args[strd_i]);
        // Compute Final Stress
        law.stress(kin, dt[i], vals);
        for (int j = 0; j < size; j++) {
            stress[strd_i + j] = vals[j];
        }
    }
}

/* ------------------------------------------------------------------------------
 |  This file provides the definitions for calculating the residuals for
 |  the different model forms
 |
 |  These models are mainly used for the collaboration with Alexey Kamenskiy
 |
 |  Author: Will Zhang
 |  Dependencies: None
 ----------------------------------------------------------------------------- */

template <
    class matlaw, class kine, ResidualNorm normfunc, HysteresisFunc hystfunc, Penalty penfunc,
    int dim>
double calc_residual(
    const double pars[], const double fiber[], const double visco[], const double Tf,
    const double Cmax[], const double args[], const double stress[], const double dt[],
    const double weights[], const double deltaCG[], const double hysteresis[], const double data[],
    const int index[], const int select[], const int n, const int nprot, const int skip,
    const double w_hyst
) {
    double *sims = new double[n * dim * dim]();
    double *diff = new double[n * dim * dim]();

    simulate<matlaw, kine, dim>(pars, fiber, visco, Tf, Cmax, args, dt, &sims[0], n);
    compute_residual_difference<dim>(sims, stress, n, diff);
    for (int k = 0; k < nprot; k++) {
        normalize_residuals<dim>(diff, index[2 * select[k]], index[2 * select[k] + 1]);
    }

    double res =
        residual_term_general<normfunc, dim>(args, diff, weights, index, select, nprot, skip);

    double hyst = hystfunc(sims, deltaCG, hysteresis, weights, index, select, nprot, skip);

    delete[] sims;
    delete[] diff;

    return (res + w_hyst * hyst) * (1.0 + penfunc(pars, fiber, visco, data));
}

/* ******************************************************************************
 * Basic functions, e.g. calculating the residual, penalty, hysteresis etc.
 *
 * COMMENTS:
 * This is the same for all of the models, so they share these same codes
 ****************************************************************************** */
// Residual Calculation
template <int dim>
void compute_residual_difference(
    const double sim[], const double stress[], const int n, double res[]
) {
    constexpr int size = dim * dim;
    for (int i = 0; i < size * n; i++) {
        res[i] = sim[i] - stress[i];
    }
}

template <int dim> void normalize_residuals(double res[], const int start, const int end) {
    constexpr int size = dim * dim;
    int strd_i;
    double n = (double)end - start;
    double mean[size] = {0.0};

    for (int i = start; i < end; i++) {
        for (int j = 0; j < size; j++) {
            strd_i = size * i;
            mean[j] = mean[j] + res[strd_i + j];
        }
    }
    for (int j = 0; j < size; j++) {
        mean[j] = mean[j] / n;
    }
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < size; j++) {
            strd_i = size * i;
            res[strd_i + j] = res[strd_i + j] - mean[j];
        }
    }
}

template <ResidualNorm norm_func, int dim>
double residual_term_general(
    const double strain[], const double residual[], const double weights[], const int index[],
    const int select[], int nprot, int skip
) {
    int kid, strd_i, strd_k;
    // double ds;
    double res, eps;
    int size = dim * dim;

    res = 0;
    for (int k = 0; k < nprot; k++) {
        kid = 2 * select[k];
        strd_k = size * select[k];
        for (int j = 0; j < size; j += (dim + 1)) {
            eps = 0;
            for (int i = index[kid]; i < index[kid + 1] + 1; i += skip) {
                strd_i = i * size + j;
                eps = eps + norm_func(residual[strd_i], strain[strd_i]);
            }
            res = res + weights[strd_k + j] * eps;
        }
    }
    return res;
}

inline double quart_quad_residual(const double residual, const double strain) {
    double r2 = residual * residual;
    double mult = (residual > 0) ? r2 * r2 : 1.0;
    return r2 * mult;
}

inline double quadratic_residual(const double residual, const double strain) {
    return residual * residual;
}

template <int dim>
double hysteresis_body_null(
    const double sims[], const double deltaCG[], const double hysteresis[], const double weights[],
    const int index[], const int select[], int nprot, int skip
) {
    return 0.0;
}

// Hysteresis Calculation
template <int dim>
double hysteresis_body(
    const double sims[], const double deltaCG[], const double hysteresis[], const double weights[],
    const int index[], const int select[], int nprot, int skip
) {
    constexpr int size = dim * dim;
    int strd_i, strd_k;
    double ds, res, hyst;
    int kid;

    res = 0;
    for (int k = 0; k < nprot; k++) {
        kid = 2 * select[k];
        strd_k = size * select[k];
        for (int j = 0; j < size; j += (dim + 1)) {
            hyst = 0;

            for (int i = index[kid]; i < index[kid + 1]; i += skip) {
                strd_i = i * size + j;
                hyst = hyst + sims[strd_i] * deltaCG[strd_i];
            }
            ds = hyst - hysteresis[strd_k + j];
            res = res + weights[strd_k + j] * ds * ds;
            // res = res + ds * ds;
        }
    }

    return res;
}

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/
} // namespace simulate
