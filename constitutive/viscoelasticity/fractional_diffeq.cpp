#include "fractional_diffeq.hpp"
#include <cmath>

namespace constitutive_models {

template <int dim>
double FractionalDiffeq<dim>::stress(
    const kinematics::kinematics<dim> &kin, const double dt, double stress[]
) {
    double p;
    double frac[dim];

    p = m_law->stress(kin, dt, frac);
    store.diffeq_iter(frac, dt, stress);
    p = store_p.diffeq_iter(p, dt);
    return p;
}

} // namespace constitutive_models