#pragma once

#include "../../interfaces.hpp"
#include "../../kinematics/kinematics.hpp"
#include "./caputo.hpp"

namespace constitutive_models {
template <int dim> class FractionalVE : public MatLawTime<dim> {
  private:
    caputo::caputo_init_vec<dim> store;
    caputo::caputo_init_scl store_p;

  public:
    MatLaw<dim> *m_law;
    FractionalVE(MatLaw<dim> &law, const double alpha, const double Tf)
        : store(alpha, Tf, 0.0), store_p(alpha, Tf, 0.0), m_law(&law) {
    }
    ~FractionalVE() {
    }
    double stress(const kinematics::kinematics<dim> &kin, const double dt, double stress[]);
};

template <int dim> class FractionalVE3D {
  private:
    caputo::caputo_init_vec<dim> store;

  public:
    MatLaw3D<dim> *m_law;
    FractionalVE3D(MatLaw3D<dim> &law, const double alpha, const double Tf)
        : store(alpha, Tf, 0.0), m_law(&law) {
    }
    ~FractionalVE3D() {
    }
    void stress(const kinematics::kinematics<dim> &kin, const double dt, double stress[]);
};

} // namespace constitutive_models
