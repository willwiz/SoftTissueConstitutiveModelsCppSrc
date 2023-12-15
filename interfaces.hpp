#pragma once

#include "kinematics/kinematics.hpp"

namespace constitutive_models
{

  template<int dim>
  class MatLaw
  {
  public:
    virtual double stress(const kinematics::kinematics<dim> &kin, double stress[]) = 0;
  };

  template<int dim>
  class MatLawTime
  {
  public:
    virtual void stress(const kinematics::kinematics<dim> &kin, const double dt, double stress[]) = 0;
  };

} // namespace constitutive_models
