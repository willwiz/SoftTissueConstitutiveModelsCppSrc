#pragma once

namespace sim {

void simulate_matrix(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
);

void simulate_elastin(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
);

void simulate_smc(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
);

void simulate_collagen(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
);

void simulate_smc_ve(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
);

void simulate_collagen_ve(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
);

} // namespace sim