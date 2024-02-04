#include "ensemble_template.hpp"
#include "constants.hpp"
#include "thoracic_full_ensemble_model.hpp"
#include "thoracic_smc_ensemble_model.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace ensemble {

template <class matlaw>
void ensemble_simulate(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double dt[],
    double out_stress[], int n
) {
    double vals[ctv::prob_dim];

    kinematics::deformation_ensemble2D kin;
    matlaw law(pars, visco, Tf, Cmax);
    for (int i = 0; i < n; i++) {
        kin.precompute(strain[i]);
        law.stress(kin, dt[i], vals);
        out_stress[i] = vals[0] + vals[3];
    }
}

template <
    class matlaw, ResidualFunction resfunc, HysteresisFunction hystfunc, PenaltyFunction penfunc>
double calc_ensemble_objective(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double stress[],
    double dt[], double deltaCG[], double hysteresis, int n, int skip
) {
    double *sims = new double[n]();
    ensemble_simulate<matlaw>(pars, visco, Tf, Cmax, strain, dt, &sims[0], n);
    double res = resfunc(strain, stress, n, skip, sims, pars[0]);
    double hyst = hystfunc(sims, deltaCG, n, hysteresis);
    delete[] sims;
    return (res + 10.0 * hyst) * (penfunc(pars, visco));
}

inline double top_sided_residual(double sim, double data, double baseline) {
    double difference = sim - data - baseline;
    double weight = (difference > 0) ? difference * difference : 1.0;

    return weight * difference * difference;
}

inline double quadratic_residual(double sim, double data, double baseline) {
    double difference = sim - data - baseline;
    return difference * difference;
}

double residual_body_below(
    double strain[], double stress[], int n, int skip, double sims[], double baseline
) {
    double res = 0.0;
    for (int k = 0; k < n; k++) {
        res = res + top_sided_residual(sims[k], stress[k], baseline) / std::pow(strain[k], 2);
    }
    return res;
}

double residual_body_quadratic(
    double strain[], double stress[], int n, int skip, double sims[], double baseline
) {
    double res = 0.0;
    for (int k = 0; k < n; k++) {
        res = res + quadratic_residual(sims[k], stress[k], baseline) / std::pow(strain[k], 1);
    }
    return res;
}

double penalty_function_null(double pars[], double visco[]) {
    return 1.0;
};

double penalty_body_ensemble2(double pars[], double visco[]) {
    double k_e = pars[1] - pars[11];
    // double b_s = pars[3];
    return 1.0 + 0.1 * k_e * k_e;
}

double penalty_body_ensemble3(double pars[], double visco[]) {
    double k_e = pars[1] - pars[11];
    // double k_sc = pars[2] + pars[4] - pars[12];
    // double k_s = pars[2] - pars[12];
    double b_s = pars[3];
    // double b_c = pars[5] - 2*pars[13];
    double v_c = visco[0] - pars[14];
    double v_s = visco[1] - pars[14];
    return 1.0 + 1.0 * k_e * k_e + 0.1 * b_s * b_s + 10.0 * v_c * v_c + 10.0 * v_s * v_s;
}

double hysteresis_body_null(double sims[], double deltaCG[], int n, double hysteresis) {
    return 0.0;
}

// Hysteresis Calculation
double hysteresis_body(double sims[], double deltaCG[], int n, double hysteresis) {

    double hyst = 0;
    for (int i = 0; i < n; i++) {
        hyst = hyst + sims[i] * deltaCG[i];
    }
    double ds = hyst - hysteresis;

    return ds * ds;
}

void thoracic_smc_ensemble_simulate(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double dt[],
    double stress[], int n
) {
    ensemble_simulate<ThoracicSMCEnsembleVE>(pars, visco, Tf, Cmax, strain, dt, &stress[0], n);
}

void thoracic_full_ensemble_simulate(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double dt[],
    double stress[], int n
) {
    ensemble_simulate<ThoracicFullEnsembleVE>(pars, visco, Tf, Cmax, strain, dt, &stress[0], n);
}

double thoracic_smc_ensemble_residual(
    double pars[], double visco[], double Tf, double Cmax[], double args[], double stress[],
    double dt[], double deltaCG[], double hysteresis, int n, int skip
) {
    return calc_ensemble_objective<
        ThoracicSMCEnsembleVE, residual_body_below, hysteresis_body_null, penalty_body_ensemble2>(
        pars, visco, Tf, Cmax, args, stress, dt, deltaCG, hysteresis, n, skip
    );
}

double thoracic_full_ensemble_residual(
    double pars[], double visco[], double Tf, double Cmax[], double args[], double stress[],
    double dt[], double deltaCG[], double hysteresis, int n, int skip
) {
    return calc_ensemble_objective<
        ThoracicFullEnsembleVE, residual_body_quadratic, hysteresis_body, penalty_body_ensemble3>(
        pars, visco, Tf, Cmax, args, stress, dt, deltaCG, hysteresis, n, skip
    );
}
} // namespace ensemble