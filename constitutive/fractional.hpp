#pragma once

#include "../constitutive/caputo.hpp"
#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace constitutive_models {
template <int dim> class FractionalVE {
  private:
    caputo::caputo_init_vec<dim> store;
    caputo::caputo_init_scl store_p;

  public:
    MatLaw<dim> *m_law;
    FractionalVE(MatLaw<dim> &law, const double alpha, const double Tf);
    ~FractionalVE();
    double stress(const kinematics::kinematics<dim> &kin, const double dt, double stress[]);
};

} // namespace constitutive_models
