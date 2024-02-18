#include "ensemble_template.hpp"
#include "constants.hpp"
#include "thoracic_full_ensemble_model.hpp"
#include <cmath>

namespace thoracic {

template <class matlaw>
void ensemble_simulate(
    const double pars[], const double visco[], double Tf, const double Cmax[],
    const double strain[], const double dt[], double out_stress[], int n
) {
    double vals[4];

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
    const double pars[], const double visco[], double Tf, const double Cmax[],
    const double strain[], const double stress[], const double dt[], const double deltaCG[],
    double hysteresis, int n, int skip
) {
    double *sims = new double[n]();
    ensemble_simulate<matlaw>(pars, visco, Tf, Cmax, strain, dt, &sims[0], n);
    double res = resfunc(strain, stress, n, skip, sims, pars[0]);
    double hyst = hystfunc(sims, deltaCG, n, hysteresis);
    delete[] sims;
    return (res + 10.0 * hyst) * (penfunc(pars, visco));
}

inline double quart_quad_residual(double sim, double data, double strain) {
    double difference = sim - data;
    double r2 = difference * difference;
    double mult = (difference > 0) ? r2 : 1.0;
    return r2 * mult;
}

inline double quadratic_residual(double sim, double data, double strain) {
    double difference = sim - data;
    return difference * difference;
}

double residual_body_below(
    const double strain[], const double stress[], int n, int skip, const double sims[],
    double baseline
) {
    double res = 0.0;
    for (int k = 0; k < n; k++) {
        res = res + quart_quad_residual(sims[k], stress[k], baseline) / std::pow(strain[k], 2);
    }
    return res;
}

double residual_body_quadratic(
    const double strain[], const double stress[], int n, int skip, const double sims[],
    double baseline
) {
    double res = 0.0;
    for (int k = 0; k < n; k++) {
        res = res + quadratic_residual(sims[k], stress[k], baseline) / std::pow(strain[k], 1);
    }
    return res;
}

double penalty_function_null(const double pars[], const double visco[]) {
    return 1.0;
};

double penalty_body_ensemble2(const double pars[], const double visco[]) {
    double k_e = pars[1] - pars[11];
    // double b_s = pars[3];
    return 1.0 + 0.1 * k_e * k_e;
}

double penalty_body_ensemble3(const double pars[], const double visco[]) {
    double k_e = pars[1] - pars[11];
    // double k_sc = pars[2] + pars[4] - pars[12];
    // double k_s = pars[2] - pars[12];
    double b_s = pars[3];
    // double b_c = pars[5] - 2*pars[13];
    double v_c = visco[0] - pars[14];
    double v_s = visco[1] - pars[14];
    return 1.0 + 1.0 * k_e * k_e + 0.1 * b_s * b_s + 10.0 * v_c * v_c + 10.0 * v_s * v_s;
}

double hysteresis_body_null(const double sims[], const double deltaCG[], int n, double hysteresis) {
    return 0.0;
}

// Hysteresis Calculation
double hysteresis_body(const double sims[], const double deltaCG[], int n, double hysteresis) {

    double hyst = 0;
    for (int i = 0; i < n; i++) {
        hyst = hyst + sims[i] * deltaCG[i];
    }
    double ds = hyst - hysteresis;

    return ds * ds;
}

void thoracic_full_ensemble_simulate(
    double pars[], double visco[], double Tf, double Cmax[], double strain[], double dt[],
    double stress[], int n
) {
    ensemble_simulate<ThoracicFullEnsembleVE>(pars, visco, Tf, Cmax, strain, dt, &stress[0], n);
}

double thoracic_full_ensemble_residual(
    const double pars[], const double visco[], double Tf, const double Cmax[],
    const double strain[], const double stress[], const double dt[], const double deltaCG[],
    double hysteresis, int n, int skip
) {
    return calc_ensemble_objective<
        ThoracicFullEnsembleVE, residual_body_quadratic, hysteresis_body, penalty_body_ensemble3>(
        pars, visco, Tf, Cmax, strain, stress, dt, deltaCG, hysteresis, n, skip
    );
}
} // namespace thoracic