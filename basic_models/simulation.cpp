#include "models.hpp"
#include "simulation.hpp"
#include "template.hpp"

/*----------------------------------------------------------------------
 |  This file provides the definitions of the different model forms
 |  which combines the primative constitive model from the other cpp files
 |  in the constitutive_models namespace
 |
 |  These models are mainly used for the collaboration with Alexey Kamenskiy
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

using namespace constitutive_models;

namespace sim {

/*----------------------------------------------------------------------
 |  This provides the main models in the full constitutive model
 -----------------------------------------------------------------------*/
void simulate_matrix(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
) {
    templates::simulate<ModelMatrix>(pars, fiber, caputo, Tf, args, dt, stress, n);
}

void simulate_elastin(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
) {
    templates::simulate<ModelElastin>(pars, fiber, caputo, Tf, args, dt, stress, n);
}

void simulate_smc(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
) {
    templates::simulate<ModelSMC>(pars, fiber, caputo, Tf, args, dt, stress, n);
}

void simulate_collagen(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
) {
    templates::simulate<ModelCollagen>(pars, fiber, caputo, Tf, args, dt, stress, n);
}

void simulate_smc_ve(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
) {
    templates::simulate<ModelSMCVE>(pars, fiber, caputo, Tf, args, dt, stress, n);
}

void simulate_collagen_ve(
    double pars[], double fiber[], double caputo[], double Tf, double args[], double dt[],
    double stress[], int n
) {
    templates::simulate<ModelCollagenVE>(pars, fiber, caputo, Tf, args, dt, stress, n);
}

} // namespace sim