#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include "../CTvalues_optimization.hpp"
#include "../interfaces.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/planar_hog.hpp"
#include "../constitutive/struc_hog_2D.hpp"
#include "../constitutive/hog_2D.hpp"
#include "../constitutive/fractional.hpp"


namespace sim {


  class ThoracicCircBase: public constitutive_models::MatLawTime<4>
  {
  protected:
    constitutive_models::NeoHookean m_matrix;
    constitutive_models::PlanarHog2D m_elastin;
    constitutive_models::PlanarHog2D m_muscle;


  public:
    ThoracicCircBase(double pars[], double fiber[]):
      m_matrix(pars[0]),
      m_elastin(pars[1], 0.0, fiber[0], fiber[1]),
      m_muscle(pars[2], pars[3], fiber[0], fiber[2])
    {}

    ThoracicCircBase(double pars[], double fiber[], double Cmax[]):
      m_matrix(pars[0]),
      m_elastin(pars[1], 0.0, fiber[0], fiber[1]),
      m_muscle(pars[2], pars[3], fiber[0], fiber[2], Cmax)
    {}

    ~ThoracicCircBase() {}

    void get_scaled_pars(double pars[]);

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);

  };


  class ThoracicCircVEBase: public ThoracicCircBase
  {
  protected:
    constitutive_models::FractionalVE<4> muscle;

  public:
    ThoracicCircVEBase(double pars[], double fiber[], double visco[], double Tf):
      ThoracicCircBase(pars, fiber),
      muscle(m_muscle, visco[0], Tf)

    {}

    ThoracicCircVEBase(double pars[], double fiber[], double visco[], double Tf, double Cmax[]):
      ThoracicCircBase(pars, fiber, Cmax),
      muscle(m_muscle, visco[0], Tf)

    {}

    ~ThoracicCircVEBase() {}

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);

  };




  class ThoracicCircVEScaled: public ThoracicCircVEBase
  {
  public:
    ThoracicCircVEScaled(double pars[], double fiber[], double visco[], double Tf, double Cmax[]):
      ThoracicCircVEBase(pars, fiber, visco, Tf, Cmax) {}
    ~ThoracicCircVEScaled() {}
  };

}
