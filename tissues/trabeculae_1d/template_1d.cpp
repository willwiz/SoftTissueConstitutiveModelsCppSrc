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

void add_drift(
    const double pars[], const double strain[], double stress_inout[], const double dt[],
    const int n
) {
    double drift = pars[3];
    for (int i = 0; i < n; i++) {
        drift = drift + dt[i] * pars[4];
        stress_inout[i] = stress_inout[i] + drift * (1 / sqrt(strain[i]));
    }
}

template <class matlaw>
void simulate_general(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double dt[], double stress_out[], int n
) {
    double vals[1];
    kinematics::deformation1D kin;
    matlaw law(pars, fiber, visco, Tf, Cmax);
    double drift = pars[3];
    for (int i = 1; i < n; i++) {
        kin.precompute(strain[i]);
        drift = drift + dt[i] * pars[4];
        law.stress(kin, dt[i], vals);
        stress_out[i] = vals[0] + drift * (1 / sqrt(strain[i]));
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

template <class matlaw, ResidualNorm norm_func, PenaltyFunction pen_func>
double calc_objective_general(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double stress[], double dt[], double weight[], int n, int nset, int index[], int select[]
) {
    double *sims = new double[n]();
    simulate_general<matlaw>(pars, fiber, visco, Tf, Cmax, strain, dt, &sims[0], n);
    double res =
        residual_body_general<norm_func>(strain, stress, weight, nset, index, select, sims);
    // (void)add_drift(pars, strain, sims, dt, n);
    delete[] sims;
    return res * (pen_func(pars, fiber, visco));
}

} // namespace optimization_1d