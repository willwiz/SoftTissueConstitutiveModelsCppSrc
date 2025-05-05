#include "./maxwell.hpp"
#include <cmath>

namespace constitutive_models {

template <int dim>
double MaxwellVE<dim>::stress(
    const kinematics::kinematics<dim> &kin, const double dt, double stress[]
) {
    // Here, p is the hydrostatic pressure
    double p, p_out, decay;
    double fvals[dim];
    // get the hyperelastic stress
    p = m_law->stress(kin, fvals);
    // compute the viscoelastic stress
    decay = 1.0 / (1.0 + alpha * dt);
    for (int i = 0; i < dim; i++) {
        stress[i] = decay * (store_visco[i] + (fvals[i] - store_hyper[i]));
    }
    p_out = decay * (store_visco_p + (p - store_hyper_p));
    // update the internal variables
    for (int i = 0; i < dim; i++) {
        store_visco[i] = stress[i];
        store_hyper[i] = fvals[i];
    }
    store_visco_p = p_out;
    store_hyper_p = p;
    return p_out;
}

} // namespace constitutive_models