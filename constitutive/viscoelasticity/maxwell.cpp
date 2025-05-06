#include "./maxwell.hpp"
#include <cmath>

namespace constitutive_models {

template <int dim>
double MaxwellVE<dim>::stress(
    const kinematics::kinematics<dim> &kin, const double dt, double stress[]
) {
    // Here, p is the hydrostatic pressure
    double p, visco_p, decay;
    double fvals[dim];
    // get the hyperelastic stress
    p = m_law->stress(kin, fvals);
    // compute the viscoelastic stress
    decay = 1.0 / (1.0 + tau * dt);
    for (int i = 0; i < dim; i++) {
        stress[i] = decay * (store_visco[i] + neta * (fvals[i] - store_hyper[i]));
    }
    visco_p = decay * (store_visco_p + neta * (p - store_hyper_p));
    // update the internal variables
    for (int i = 0; i < dim; i++) {
        store_visco[i] = stress[i];
        store_hyper[i] = fvals[i];
    }
    store_visco_p = visco_p;
    store_hyper_p = p;
    return visco_p;
}

} // namespace constitutive_models