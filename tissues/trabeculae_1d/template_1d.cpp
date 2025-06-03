#include "template_1d.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "../../CTvalues_optimization.hpp"

namespace optimization_1d {

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

double penalty_function_null(double pars[], double fiber[], double visco[]) {
    return 1.0;
};

template <class matlaw>
void simulate_general(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double dt[], double stress_out[], int n
) {
    kinematics::deformation1D kin;
    matlaw law(pars, fiber, visco, Tf, Cmax);
    for (int i = 1; i < n; i++) {
        kin.precompute(strain[i]);
        law.stress(kin, dt[i], stress_out[i]);
    }
}

template <class matlaw>
void simulate_drift(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double dt[], double stress_out[], int n
) {
    double vals[1];
    kinematics::deformation1D kin;
    matlaw law(pars, fiber, visco, Tf, Cmax);
    double d_drift = drift[0];
    for (int i = 1; i < n; i++) {
        kin.precompute(strain[i]);
        d_drift = d_drift + dt[i] * drift[1];
        law.stress(kin, dt[i], vals);
        stress_out[i] = vals[0] + d_drift * (1 / sqrt(strain[i]));
    }
}

template <ResidualNorm norm_func>
double residual_body_general(
    double strain[], double stress[], double weight[], int nset, int index[], int select[],
    double sims[]
) {
    double res = 0.0;
    double res_protocol;
    for (int k = 0; k < nset; k++) {
        res_protocol = 0.0;
        for (int i = index[select[k]]; i <= index[select[k] + 1]; i++) {
            res_protocol = res_protocol + norm_func(sims[i], stress[i], strain[i]);
        }
        res = res + weight[k] * res_protocol;
    }
    return res;
}

double residual_body_drift(
    double strain[], double stress[], double dt[], double weight[], int nset, int index[],
    int select[], double sims[]
) {
    // Get the maximum length of the time vector
    int n = index[select[nset] + 1] + 1;
    // Compute the time vectors for the residual projection matrix
    double *t_ptr = new double[n]();
    double sum_t = 0.0, sum_t2 = 0.0;
    for (int k = 0; k < n; k++) {
        t_ptr[k + 1] = t_ptr[k] + dt[k];
    }
    for (int k = 0; k < nset; k++) {
        for (int m = index[select[k]]; m <= index[select[k] + 1]; m++) {
            sum_t = sum_t + t_ptr[m + 1];
            sum_t2 = sum_t2 + t_ptr[m + 1] * t_ptr[m + 1];
        }
    }
    // Compute the epsilon for sim data difference
    double *eps_ptr = new double[n]();
    for (int k = 0; k < n; k++) {
        eps_ptr[k] = sims[k] - stress[k];
    }
    // Compute the projected residuals
    double *res_ptr = new double[n]();
    double w_projection = 1.0 / (sum_t2 - sum_t * sum_t);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res_ptr[i] += (static_cast<double>(i == j) + (t_ptr[i] + t_ptr[j]) * sum_t -
                           t_ptr[i] * t_ptr[j] - sum_t2) *
                          eps_ptr[j];
        }
        res_ptr[i] *= w_projection;
    }
    // Compute norms
    double Linf_norm = 0.0;
    double L2_norm = 0.0;
    double L2_norm_p;
    for (int k = 0; k < nset; k++) {
        L2_norm_p = 0.0;
        for (int m = index[select[k]]; m <= index[select[k] + 1]; m++) {
            L2_norm_p += eps_ptr[m] * res_ptr[m];
            Linf_norm = std::max(Linf_norm, std::abs(eps_ptr[m]));
        }
        L2_norm += weight[k] * L2_norm_p;
    }
    // Clean up
    delete[] t_ptr;
    delete[] eps_ptr;
    delete[] res_ptr;
    return L2_norm + Linf_norm;
}

template <class matlaw, ResidualNorm norm_func, PenaltyFunction pen_func>
double calc_objective_drift(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double stress[], double dt[], double weight[], int n, int nset, int index[],
    int select[]
) {
    double *sims = new double[n]();
    simulate_drift<matlaw>(pars, fiber, visco, Tf, drift, Cmax, strain, dt, &sims[0], n);
    double res =
        residual_body_general<norm_func>(strain, stress, weight, nset, index, select, sims);
    delete[] sims;
    return res * (pen_func(pars, fiber, visco));
}

template <class matlaw, PenaltyFunction pen_func>
double calc_objective_lgres(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double stress[], double dt[], double weight[], int n, int nset, int index[], int select[]
) {
    double *sims = new double[n]();
    simulate_general<matlaw>(pars, fiber, visco, Tf, Cmax, strain, dt, &sims[0], n);
    double res = residual_body_drift(strain, stress, dt, weight, nset, index, select, sims);
    delete[] sims;
    return res * (pen_func(pars, fiber, visco));
}

} // namespace optimization_1d