#include "ensemble_template.hpp"
#include "constants.hpp"
#include "thoracic_ensemble_model.hpp"
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

template <class matlaw, NormFunction nfunc, HysteresisFunction hfunc, PenaltyFunction pfunc>
double calc_ensemble_objective(
    const double pars[], const double visco[], double Tf, const double Cmax[],
    const double strain[], const double stress[], const double dt[], const double deltaCG[],
    double hysteresis, int n, int skip
) {
    double *sims = new double[n]();
    ensemble_simulate<matlaw>(pars, visco, Tf, Cmax, strain, dt, &sims[0], n);
    double res = ensemble_residual<nfunc>(strain, stress, n, skip, sims, pars[0]);
    double hyst = hfunc(sims, deltaCG, n, hysteresis);
    delete[] sims;
    return (res + 10.0 * hyst) * (1.0 + pfunc(pars, visco));
}

template <NormFunction nfunc>
double ensemble_residual(
    const double strain[], const double stress[], int n, int skip, const double sims[],
    double baseline
) {
    double res = 0.0;
    for (int k = 0; k < n; k++) {
        res = res + nfunc(sims[k], stress[k], strain[k]);
    }
    return res;
}

} // namespace thoracic