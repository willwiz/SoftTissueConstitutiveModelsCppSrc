#include "fractional.hpp"
#include <cmath>

namespace constitutive_models {

template <int dim>
double FractionalVE<dim>::stress(
    const kinematics::kinematics<dim> &kin, const double dt, double stress[]
) {
    double p;
    double frac[dim];

    p = m_law->stress(kin, frac);
    store.caputo_iter(frac, dt, stress);
    p = store_p.caputo_iter(p, dt);
    return p;
}

template <int dim>
void FractionalVE3D<dim>::stress(
    const kinematics::kinematics<dim> &kin, double dt, double stress[]
) {
    double frac[dim];

    m_law->stress(kin, frac);
    store.caputo_iter(frac, dt, stress);
}

} // namespace constitutive_models