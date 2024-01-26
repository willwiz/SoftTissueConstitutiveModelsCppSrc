#define _USE_MATH_DEFINES
#include "trabeculae_1d.hpp"
#include "template_1d.hpp"
#include <cmath>

namespace optimization_1d {

void Trabeculae1D::stress(const kinematics::kinematics<1> &kin, const double dt, double stress[]) {
    double p = 0.0;
    double iso[1], smc[1];
    p = m_matrix.stress(kin, iso);
    p = p + muscle.stress(kin, dt, smc);
    stress[0] = iso[0] + smc[0] - p * kin.I_n * kin.Cinv[0];
}

void Trabeculae2Phase1D::stress(
    const kinematics::kinematics<1> &kin, const double dt, double stress[]
) {
    double p = 0.0;
    double iso[1], smc[1], col[1];
    p = m_matrix.stress(kin, iso);
    p = p + muscle.stress(kin, dt, smc);
    p = p + collagen.stress(kin, dt, col);
    stress[0] = iso[0] + smc[0] + col[0] - p * kin.I_n * kin.Cinv[0];
}

void simulation_function_trabeculae_1d(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double dt[], double stress_out[], int n
) {
    simulate_general<Trabeculae1D>(pars, fiber, visco, Tf, Cmax, strain, dt, stress_out, n);
}

double objective_function_trabeculae_stacked_1d(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double stress[], double dt[], double weight[], int n[]
) {
    return calc_objective_general_stacked<Trabeculae1D, quadratic_residual, penalty_function_null>(
        pars, fiber, visco, Tf, Cmax, strain, stress, dt, weight, n
    );
}

double objective_function_trabeculae_1d(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double stress[], double dt[], double weight[], int n, int nset, int index[], int select[]
) {
    return calc_objective_general<Trabeculae1D, quadratic_residual, penalty_function_null>(
        pars, fiber, visco, Tf, Cmax, strain, stress, dt, weight, n, nset, index, select
    );
}

void simulation_function_trabeculae_2phase_1d(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double dt[], double stress_out[], int n
) {
    simulate_general<Trabeculae2Phase1D>(pars, fiber, visco, Tf, Cmax, strain, dt, stress_out, n);
}

double objective_function_trabeculae_2phase_1d(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double stress[], double dt[], double weight[], int n[]
) {
    return calc_objective_general_stacked<
        Trabeculae2Phase1D, quadratic_residual, penalty_function_null>(
        pars, fiber, visco, Tf, Cmax, strain, stress, dt, weight, n
    );
}

} // namespace optimization_1d
