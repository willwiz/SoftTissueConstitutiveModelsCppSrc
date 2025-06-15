#define _USE_MATH_DEFINES
#include "api.hpp"
#include "template_1d.hpp"
#include "trabeculae_calcium.hpp"
#include "trabeculae_calcium_diffeq.hpp"
#include "trabeculae_calcium_only.hpp"
#include "trabeculae_passive.hpp"
#include "trabeculae_passive_diffeq.hpp"
#include <cmath>

namespace optimization_1d {

double prior_func4(double pars[], double visco[]) {
    double v = (visco[0] - 0.2) * (visco[0] - 0.2);
    return 1.0 + 4.0 * v;
}

double prior_func7(double pars[], double visco[]) {
    double v = (visco[0] - 0.2) * (visco[0] - 0.2);
    double p = (pars[4] - pars[5]) * (pars[4] - pars[5]) +
               (pars[4] - 2.0 * pars[6]) * (pars[4] - 2.0 * pars[6]) +
               (2.0 * pars[6] - pars[5]) * (2.0 * pars[6] - pars[5]);
    return 1.0 + 10.0 * v + 0.00001 * p;
}

void simulation_function_trabeculae_active(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double dt[], double stress_out[], int n
) {
    simulate_drift<TrabeculaeCalciumOnly>(
        pars, fiber, visco, Tf, drift, Cmax, strain, dt, stress_out, n
    );
}
void simulation_function_trabeculae_passive(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double dt[], double stress_out[], int n
) {
    simulate_drift<TrabeculaePassive>(
        pars, fiber, visco, Tf, drift, Cmax, strain, dt, stress_out, n
    );
}

void simulation_function_trabeculae_calcium(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double dt[], double stress_out[], int n
) {
    simulate_drift<TrabeculaeCalcium>(
        pars, fiber, visco, Tf, drift, Cmax, strain, dt, stress_out, n
    );
}

void simulation_function_trabeculae_passive_diffeq(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double dt[], double stress_out[], int n
) {
    simulate_drift<TrabeculaePassiveDiffEq>(
        pars, fiber, visco, Tf, drift, Cmax, strain, dt, stress_out, n
    );
}

void simulation_function_trabeculae_calcium_diffeq(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double dt[], double stress_out[], int n
) {
    simulate_drift<TrabeculaeCalciumDiffEq>(
        pars, fiber, visco, Tf, drift, Cmax, strain, dt, stress_out, n
    );
}

double objective_function_trabeculae_active(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double stress[], double dt[], double weight[], int n, int nset, int index[],
    int select[]
) {
    return prior_func7(pars, visco) *
           calc_objective_lgres<TrabeculaeCalciumOnly, penalty_function_null>(
               pars, fiber, visco, Tf, Cmax, strain, stress, dt, weight, n, nset, index, select
           );
}

double objective_function_trabeculae_passive(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double stress[], double dt[], double weight[], int n, int nset, int index[],
    int select[]
) {
    return prior_func4(pars, visco) *
           calc_objective_lgres<TrabeculaePassive, penalty_function_null>(
               pars, fiber, visco, Tf, Cmax, strain, stress, dt, weight, n, nset, index, select
           );
}

double objective_function_trabeculae_calcium_diffeq(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double stress[], double dt[], double weight[], int n, int nset, int index[],
    int select[]
) {
    return prior_func7(pars, visco) *
           calc_objective_lgres<TrabeculaeCalcium, penalty_function_null>(
               pars, fiber, visco, Tf, Cmax, strain, stress, dt, weight, n, nset, index, select
           );
}

double objective_function_trabeculae_passive_diffeq(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double stress[], double dt[], double weight[], int n, int nset, int index[],
    int select[]
) {
    return prior_func7(pars, visco) *
           calc_objective_lgres<TrabeculaePassiveDiffEq, penalty_function_null>(
               pars, fiber, visco, Tf, Cmax, strain, stress, dt, weight, n, nset, index, select
           );
}

double objective_function_trabeculae_calcium(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double stress[], double dt[], double weight[], int n, int nset, int index[],
    int select[]
) {
    return prior_func7(pars, visco) *
           calc_objective_lgres<TrabeculaeCalciumDiffEq, penalty_function_null>(
               pars, fiber, visco, Tf, Cmax, strain, stress, dt, weight, n, nset, index, select
           );
}

} // namespace optimization_1d
