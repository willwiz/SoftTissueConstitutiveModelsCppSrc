#pragma once

#define _USE_MATH_DEFINES

#include <cmath>
#include "../interfaces.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/planar_hog.hpp"


namespace sim {

  class PlanarElastinMatrix: public constitutive_models::MatLawTime<4>
  {
  protected:
    constitutive_models::NeoHookean m_matrix;
    constitutive_models::PlanarHog2D m_elastin;

  public:
    PlanarElastinMatrix(double pars[], double fiber[], double visco[], double Tf, double Cmax[]):
      m_matrix(pars[0]),
      m_elastin(pars[1], 0.0, fiber[0], fiber[1])
    {}

    ~PlanarElastinMatrix() {}

    void get_scaled_pars(double pars[]);

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);

  };



}
