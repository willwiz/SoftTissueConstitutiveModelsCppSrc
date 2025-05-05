#pragma once

namespace optimization_1d {
void simulation_function_trabeculae_passive(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double dt[], double stress_out[], int n
);
void simulation_function_trabeculae_calcium(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double dt[], double stress_out[], int n
);

double objective_function_trabeculae_passive(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double stress[], double dt[], double weight[], int n, int nset, int index[], int select[]
);

double objective_function_trabeculae_calcium(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double stress[], double dt[], double weight[], int n, int nset, int index[], int select[]
);
} // namespace optimization_1d
