#define _USE_MATH_DEFINES

#include "caputo.hpp"
#include "caputo_precomp.hpp"
#include <cmath>
#include <iostream>
#include <stdio.h>

/*----------------------------------------------------------------------
 |  Calculations of the caputo derivative
 |  Contains mainly of a scalar and vector version
 |  To Accomadate the dimensions of space, a template the used
 |  The base class contains the precompute for the caputo parameters, which
 |    is the same for all classes
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

namespace caputo {

/*----------------------------------------------------------------------
 |  Caputo Derivative Base Class
 -----------------------------------------------------------------------*/
template <int n_prony>
caputo_init<n_prony>::caputo_init(double alpha, double Tf, double delta) : dt{} {
    set_pars(alpha, Tf, delta);
}

template <int n_prony> void caputo_init<n_prony>::set_pars(double alpha, double Tf, double delta) {
    double freq;
    double val;

    freq = 2.0 * M_PI / Tf;

    this->alpha = alpha;
    this->Tf = Tf;
    this->delta = delta;

    val = interpolate_caputo_parameter_beta(alpha, caputo_pars[n_prony].betam[n_prony]);

    this->beta0 = val * pow(freq, alpha - 1.0);

    for (int i = 0; i < n_prony; i++) {
        this->betas[i] = interpolate_caputo_parameter_beta(alpha, caputo_pars[n_prony].betam[i]) *
                         pow(freq, alpha);
        this->taus[i] =
            interpolate_caputo_parameter_taus(alpha, caputo_pars[n_prony].taum[i]) / freq;
    }
}
// Set parameters for dt, if dt is the same then there is no need to change
template <int n_prony> void caputo_init<n_prony>::update_dt(double dt) {
    if (this->dt == dt) return;
    if (dt < 2.0e-16) return;

    double ek;
    this->dt = dt;
    // C0 = beta0 / dt;
    C0 = 0.0;
    for (int k = 0; k < n_prony; k++) {
        ek = exp(-0.5 * dt / taus[k]);
        e2[k] = ek * ek;
        bek[k] = betas[k] * ek;
    }

    if (delta == 0.0) return;
    K0 = C0;
    for (int k = 0; k < n_prony; k++) {
        K0 = K0 + bek[k];
    }
    K0 = delta * K0;
    K1 = K0 + 1;
}

template <int n_prony> void caputo_init<n_prony>::update_dt_lin(double dt) {
    if (this->dt == dt) return;
    if (dt < 2.0e-16) return;

    this->dt = dt;
    // C0 = beta0 / dt;
    C0 = 0.0;
    for (int k = 0; k < n_prony; k++) {
        e2[k] = taus[k] / (taus[k] + dt);
        bek[k] = betas[k] * e2[k];
    }

    if (delta == 0.0) return;
    K0 = C0;
    for (int k = 0; k < n_prony; k++) {
        K0 = K0 + bek[k];
    }
    K0 = delta * K0;
    K1 = K0 + 1.0;
}

/*----------------------------------------------------------------------
 |  Implementation for the derivatives for scalar valued functions
 -----------------------------------------------------------------------*/
// The caputo derivative
template <int n_prony> double caputo_init_scl<n_prony>::caputo_iter(double fn, double dt) {

    this->caputo_init<n_prony>::update_dt_lin(dt);

    df = fn - f_prev;

    // double v = C0 * df;
    double v = 0.0;

    for (int k = 0; k < n_prony; k++) {
        Q[k] = this->e2[k] * Q[k] + this->bek[k] * df;
        v = v + Q[k];
    }

    f_prev = fn;

    return v;
}

// The solution to the fractional differential equation
template <int n_prony> double caputo_init_scl<n_prony>::diffeq_iter(double fn, double dt) {

    this->caputo_init<n_prony>::update_dt_lin(dt);
    // Stress calculations
    double v = fn + this->K0 * f_prev;

    for (int k = 0; k < n_prony; k++) {
        v = v - this->delta * this->e2[k] * Q[k];
    }
    v = v / this->K1;

    // Updates
    df = v - f_prev;
    f_prev = v;
    for (int k = 0; k < n_prony; k++) {
        Q[k] = this->e2[k] * Q[k] + this->bek[k] * df;
    }

    return v;
}

/*----------------------------------------------------------------------
 |  Implementation for the derivatives for vector valued functions
 |  Template provides the size for the dimensions. Typical values are 1, 4, 9
 -----------------------------------------------------------------------*/

// Fractional Derivative

template <int dim, int n_prony>
void caputo_init_vec<dim, n_prony>::caputo_iter(const double fn[], const double dt, double v[]) {

    int krow;

    this->caputo_init<n_prony>::update_dt_lin(dt);

    for (int i = 0; i < dim; i++) {
        df[i] = fn[i] - f_prev[i];
        // v[i] = C0 * df[i];
        v[i] = 0.0;
        f_prev[i] = fn[i];
    }

    krow = 0;
    for (int k = 0; k < n_prony; k++) {
        for (int i = 0; i < dim; i++) {
            Q[krow + i] = this->e2[k] * Q[krow + i] + this->bek[k] * df[i];
            v[i] = v[i] + Q[krow + i];
        }
        krow += dim;
    }
}

// The fractional differential equation
template <int dim, int n_prony>
void caputo_init_vec<dim, n_prony>::diffeq_iter(const double fn[], const double dt, double v[]) {

    int krow;

    this->caputo_init<n_prony>::update_dt_lin(dt);

    // Stress calculations
    for (int i = 0; i < dim; i++) {
        v[i] = fn[i] + this->K0 * f_prev[i];
    }

    krow = 0;
    for (int k = 0; k < n_prony; k++) {

        for (int i = 0; i < dim; i++) {
            v[i] = v[i] - this->delta * this->e2[k] * Q[krow + i];
        }
        krow += dim;
    }
    for (int i = 0; i < dim; i++) {
        v[i] = v[i] / this->K1;
    }

    // Updates
    for (int i = 0; i < dim; i++) {
        df[i] = v[i] - f_prev[i];
        f_prev[i] = v[i];
    }

    krow = 0;
    for (int k = 0; k < n_prony; k++) {

        for (int i = 0; i < dim; i++) {
            Q[krow + i] = this->e2[k] * Q[krow + i] + this->bek[k] * df[i];
        }
        krow += dim;
    }
}

/*----------------------------------------------------------------------
 |  Tools for interpolation of the precomputed betas and taus
 -----------------------------------------------------------------------*/

/***********************************************************************
 * AUTHOR COMMENTS:
 *
 * These values are given for alpha in the range of [0.01, 1.0] in the
 *    steps of 0.01. There are no values at alpha = 0.0
 *  For betas, beta(0) = 0.0
 *  For taus, tau(0) is unknown, thus extrapolation is used
 ************************************************************************/

double interpolate1D_newton_linear(double p1, double p2, double t) {
    return t * (p2 - p1) + p1;
}

double extrapolate1D_newton_linear(double p1, double p2, double t) {
    return p1 - (1.0 - t) * (p2 - p1);
}

double interpolate_caputo_parameter_arr(double alpha, const double arr[100]) {
    double as, temp;
    int an;

    temp = alpha * 100.0;

    an = (int)temp;
    as = temp - an;
    if (alpha <= 0.01) {
        return arr[0];
    } else if (alpha >= 1.0) {
        return arr[99];
    } else {
        return interpolate1D_newton_linear(arr[an - 1], arr[an], as);
    }
}

double interpolate_caputo_parameter_beta(double alpha, const double arr[100]) {
    double as, temp;
    int an;

    temp = alpha * 100.0;

    an = (int)temp;
    as = temp - an;
    if (alpha <= 0.01) {
        return interpolate1D_newton_linear(0.0, arr[0], as);
    } else if (alpha >= 1.0) {
        return arr[99];
    } else {
        return interpolate1D_newton_linear(arr[an - 1], arr[an], as);
    }
}

double interpolate_caputo_parameter_taus(double alpha, const double arr[100]) {
    double as, temp;
    int an;

    temp = alpha * 100.0;

    an = (int)temp;
    as = temp - an;
    if (alpha <= 0.01) {
        return extrapolate1D_newton_linear(arr[0], arr[1], as);
    } else if (alpha >= 1.0) {
        return arr[99];
    } else {
        return interpolate1D_newton_linear(arr[an - 1], arr[an], as);
    }
}

/*----------------------------------------------------------------------
 |  The precomputed beta and tau values are given here.
 -----------------------------------------------------------------------*/

} // namespace caputo

/*----------------------------------------------------------------------
 |  The precomputed beta and tau values are given here. For now, only
 |  N_p = 9 Prony terms are included (for 500 Fourier series approximations)
 |
 |  More should be imported from the matlab data files at a later time
 -----------------------------------------------------------------------*/