#include "template_1d.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "../../CTvalues_optimization.hpp"
#include "../../residual/frechet/onesided.hpp"

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

double penalty_function_13(double pars[], double fiber[], double visco[]) {
    return 1.0 + (pars[12] - 0.2) * (pars[12] - 0.2);
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
        law.stress(kin, dt[i], &stress_out[i]);
    }
}

template <class matlaw>
void simulate_drift(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double dt[], double stress_out[], int n
) {
    // std::cout << "Simulating with parameters: ";
    // for (int i = 0; i < 6; i++) {
    //     std::cout << pars[i] << " ";
    // }
    // std::cout << std::endl;
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
        for (int i = index[select[k]]; i < index[select[k] + 1]; i++) {
            res_protocol = res_protocol + norm_func(sims[i], stress[i], strain[i]);
        }
        res = res + weight[k] * res_protocol;
    }
    return res;
}

double residual_body_zeroed(
    double strain[], double stress[], double dt[], double weight[], int nset, int index[],
    int select[], double sims[]
) {
    double l2norm = 0.0;
    double linfnorm = 0.0;
    double eps;
    double res_protocol;
    double zeroed_stress = stress[index[select[0]]];
    for (int k = 0; k < nset; k++) {
        res_protocol = 0.0;
        for (int i = index[select[k]]; i < index[select[k] + 1]; i++) {
            eps = (sims[i] - stress[i] + zeroed_stress);
            linfnorm = std::max(linfnorm, std::abs(eps));
            res_protocol = res_protocol + eps * eps;
        }
        l2norm = l2norm + weight[k] * res_protocol;
    }
    return l2norm + 0.01 * linfnorm;
}

double residual_body_drift(
    const double strain[], const double stress[], const double dt[], const double weight[],
    int nset, const int index[], const int select[], const double sims[]
) {
    // Get the maximum length of the time vector
    int n = index[select[nset - 1] + 1];
    double *eps_ptr = new double[n + 1]();
    double zeroed_stress = stress[index[select[0]]];
    double zeroed_sim = sims[index[select[0]]];
    for (int k = 0; k < nset; k++) {
        for (int i = index[select[k]]; i < index[select[k] + 1]; i++) {
            // Compute the residual for sim data difference
            eps_ptr[i] = stress[i] - sims[i] - zeroed_stress + zeroed_sim;
        }
    }
    // Compute the normalized time vectors for the residual projection matrix
    double *t_ptr = new double[n + 1]();
    for (int k = 0; k < n; k++) {
        t_ptr[k + 1] = t_ptr[k] + dt[k] / (dt[2] * n);
    }
    // Compute the projection inner matrix
    double sum_t = 0.0, sum_t2 = 0.0, sum_ones = 0.0;
    for (int k = 0; k < nset; k++) {
        for (int i = index[select[k]]; i < index[select[k] + 1]; i++) {
            // sum_ones = sum_ones + 1.0;
            // sum_t = sum_t + t_ptr[i + 1];
            sum_t2 = sum_t2 + t_ptr[i + 1] * t_ptr[i + 1];
        }
    }
    // double w_projection = 1.0 / (sum_t2 * sum_ones - sum_t * sum_t);
    // double w_projection = 1.0 / (sum_t2);
    // sum_t2 = sum_t2 * w_projection;
    // sum_t = sum_t * w_projection;
    // sum_ones = sum_ones * w_projection;
    // Compute the projected regression in 2 steps
    double sum_ty = 0.0;
    // double sum_y = 0.0;
    for (int k = 0; k < nset; k++) {
        for (int i = index[select[k]]; i < index[select[k] + 1]; i++) {
            // sum_y = sum_y + eps_ptr[i];
            sum_ty = sum_ty + t_ptr[i] * eps_ptr[i];
        }
    }
    // double slope = (sum_ones * sum_ty - sum_t * sum_y);
    // double intercept = (-sum_t * sum_ty + sum_t2 * sum_y);
    double slope = sum_ty / sum_t2;
    double *fix_data = new double[n + 1]();
    for (int k = 0; k < nset; k++) {
        for (int i = index[select[k]]; i < index[select[k] + 1]; i++) {
            // Subtract the linear regression from the residual
            fix_data[i] = stress[i] - (slope * t_ptr[i]) - zeroed_stress + zeroed_sim;
        }
    }
    //
    for (int k = 0; k < nset; k++) {
        res::window_residual_onesided(
            &fix_data[0], &sims[0], index[select[k]], index[select[k] + 1], &eps_ptr[0]
        );
    }
    // Compute norms
    double Linf_norm = 0.0;
    double L2_norm = 0.0;
    double L2_norm_p;
    for (int k = 0; k < nset; k++) {
        // if (k % 2 == 0) {
        //     continue; // Skip even numbers
        // }
        L2_norm_p = 0.0;
        for (int i = index[select[k]]; i < index[select[k] + 1]; i++) {
            Linf_norm = std::max(Linf_norm, std::abs(eps_ptr[i]));
            L2_norm_p += eps_ptr[i] * eps_ptr[i];
        }
        L2_norm += weight[k] * L2_norm_p;
    }
    // Clean up
    delete[] fix_data;
    delete[] eps_ptr;
    delete[] t_ptr;
    return L2_norm;
}

double residual_body_window(
    const double strain[], const double stress[], const double dt[], const double weight[],
    int nset, const int index[], const int select[], const double sims[]
) {
    int n = index[select[nset - 1] + 1];
    double *eps_ptr = new double[n + 1]();
    for (int k = 0; k < nset; k++) {
        res::window_residual_onesided(
            &stress[0], &sims[0], index[select[k]], index[select[k] + 1], &eps_ptr[0]
        );
    }
    // Compute norms
    double L2_norm = 0.0;
    double L2_norm_p;
    for (int k = 0; k < nset; k++) {
        // if (k > 5) {
        //     continue; // Skip even numbers
        // }
        L2_norm_p = 0.0;
        for (int i = index[select[k]]; i < index[select[k] + 1]; i++) {
            L2_norm_p += eps_ptr[i] * eps_ptr[i];
        }
        L2_norm += weight[k] * L2_norm_p;
    }
    // Clean up
    delete[] eps_ptr;
    return L2_norm;
}

double residual_protocol_l1(
    const double stress[], const double sims[], int start, int end, int shift
) {
    // Compute norms
    double L1_norm = 0.0;
    for (int i = start; i < end; i++) {
        L1_norm += (stress[i] - sims[i + shift]) * (stress[i] - sims[i + shift]);
    }
    return L1_norm;
}

double residual_body_l2(
    const double strain[], const double stress[], const double dt[], const double weight[],
    int nset, const int index[], const int select[], const double sims[]
) {
    // Compute norms
    double L2_norm = 0.0;
    for (int k = 0; k < nset; k++) {
        // if (k % 2 == 0) {
        //     continue; // Skip even numbers
        // }
        // L2_norm +=
        //     weight[k] * std::min(
        //                     residual_protocol_l1(
        //                         &stress[0], &sims[0], index[select[k]], index[select[k] + 1], 0
        //                     ),
        //                     residual_protocol_l1(
        //                         &stress[0], &sims[0], index[select[k]], index[select[k] + 1], 1
        //                     )
        //                 );
        L2_norm +=
            weight[k] *
            residual_protocol_l1(&stress[0], &sims[0], index[select[k]], index[select[k] + 1], 0);
    }
    // Clean up
    return L2_norm;
}

double residual_body_neighbour(
    const double strain[], const double stress[], const double dt[], const double weight[],
    int nset, const int index[], const int select[], const double sims[]
) {
    double L2_norm[5] = {0.0};
    for (int m = -2; m < 3; m++) {
        double L2_norm_p;
        for (int k = 0; k < nset; k++) {
            for (int i = index[select[k]]; i < index[select[k] + 1]; i++) {
                // Compute the residual for sim data difference
                L2_norm_p += (stress[i] - sims[i + m]) * (stress[i] - sims[i + m]);
            }
        }
        L2_norm[m + 2] = L2_norm_p;
    }
    return *std::max_element(L2_norm, L2_norm + 5);
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
    double res = residual_body_l2(strain, stress, dt, weight, nset, index, select, sims);
    delete[] sims;
    return res;
}

template <class matlaw, PenaltyFunction pen_func>
double calc_objective_zeroed(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double stress[], double dt[], double weight[], int n, int nset, int index[], int select[]
) {
    double *sims = new double[n]();
    simulate_general<matlaw>(pars, fiber, visco, Tf, Cmax, strain, dt, &sims[0], n);
    double res = residual_body_zeroed(strain, stress, dt, weight, nset, index, select, sims);
    delete[] sims;
    return res * (pen_func(pars, fiber, visco));
}

} // namespace optimization_1d