#pragma once

#include "../../interfaces.hpp"
#include "../../kinematics/kinematics.hpp"
#include "caputo.hpp"

namespace constitutive_models {
template <int dim> class FractionalDiffeq {
  private:
    caputo::caputo_init_vec<dim, 9> store;
    caputo::caputo_init_scl<9> store_p;

  public:
    MatLawTime<dim> *m_law;
    FractionalDiffeq(MatLawTime<dim> &law, const double alpha, const double delta, const double Tf)
        : store(alpha, Tf, delta), store_p(alpha, Tf, delta), m_law(&law) {
    }
    ~FractionalDiffeq() {
    }
    double stress(const kinematics::kinematics<dim> &kin, const double dt, double stress[]);
};

} // namespace constitutive_models
