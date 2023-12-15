#include "../kinematics/kinematics.hpp"
#include "fractional.hpp"
#include <cmath>

namespace constitutive_models
{
  template<int dim>
  FractionalVE<dim>::FractionalVE(MatLaw<dim> &law, const double alpha, const double Tf)
    : store(alpha, Tf, 0.0), store_p(alpha, Tf, 0.0), m_law(&law) {}

  template<int dim>
  FractionalVE<dim>::~FractionalVE() {}

  template<int dim>
  double FractionalVE<dim>::stress(const kinematics::kinematics<dim> &kin, const double dt, double stress[])
  {
    double p;
    double frac[dim];

    p = m_law->stress(kin, frac);
    store.caputo_iter(frac, dt, stress);
    p = store_p.caputo_iter(p, dt);
    return p;
  }



} // namespace constitutive