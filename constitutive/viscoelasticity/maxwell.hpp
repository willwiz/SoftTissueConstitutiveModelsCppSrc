#pragma once

#include "../../interfaces.hpp"
#include "../../kinematics/kinematics.hpp"

namespace constitutive_models {
template <int dim> class MaxwellVE {
  private:
    double alpha;
    double store_hyper[dim];
    double store_visco[dim];
    double store_hyper_p;
    double store_visco_p;

  public:
    MatLaw<dim> *m_law;
    MaxwellVE(MatLaw<dim> &law, const double alpha)
        : alpha{alpha}, store_hyper{0}, store_visco{0}, store_hyper_p{0}, store_visco_p{0},
          m_law(&law) {};
    ~MaxwellVE() {};
    ;
    double stress(const kinematics::kinematics<dim> &kin, const double dt, double stress[]);
};

} // namespace constitutive_models
