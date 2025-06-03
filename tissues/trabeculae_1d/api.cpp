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
    return calc_objective_drift<TrabeculaeCalciumOnly, quadratic_residual, penalty_function_null>(
        pars, fiber, visco, Tf, drift, Cmax, strain, stress, dt, weight, n, nset, index, select
    );
}

double objective_function_trabeculae_passive(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double stress[], double dt[], double weight[], int n, int nset, int index[],
    int select[]
) {
    return calc_objective_drift<TrabeculaePassive, quadratic_residual, penalty_function_null>(
        pars, fiber, visco, Tf, drift, Cmax, strain, stress, dt, weight, n, nset, index, select
    );
}

double objective_function_trabeculae_calcium_diffeq(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double stress[], double dt[], double weight[], int n, int nset, int index[],
    int select[]
) {
    return calc_objective_drift<TrabeculaeCalcium, quadratic_residual, penalty_function_null>(
        pars, fiber, visco, Tf, drift, Cmax, strain, stress, dt, weight, n, nset, index, select
    );
}

double objective_function_trabeculae_passive_diffeq(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double stress[], double dt[], double weight[], int n, int nset, int index[],
    int select[]
) {
    return calc_objective_drift<TrabeculaePassiveDiffEq, quadratic_residual, penalty_function_null>(
        pars, fiber, visco, Tf, drift, Cmax, strain, stress, dt, weight, n, nset, index, select
    );
}

double objective_function_trabeculae_calcium(
    double pars[], double fiber[], double visco[], double Tf, double drift[], double Cmax[],
    double strain[], double stress[], double dt[], double weight[], int n, int nset, int index[],
    int select[]
) {
    return calc_objective_drift<TrabeculaeCalciumDiffEq, quadratic_residual, penalty_function_null>(
        pars, fiber, visco, Tf, drift, Cmax, strain, stress, dt, weight, n, nset, index, select
    );
}

} // namespace optimization_1d
