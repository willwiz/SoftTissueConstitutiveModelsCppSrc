#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include "../CTvalues_optimization.hpp"
#include "../interfaces.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/planar_hog.hpp"
#include "../constitutive/fractional.hpp"


namespace ensemble {


  class ThoracicSMCEnsembleBase: public constitutive_models::MatLawTime<4>
  {
  protected:
    constitutive_models::NeoHookean m_matrix;
    constitutive_models::PlanarHog2D m_elastin;
    constitutive_models::PlanarHog2D m_muscle;


  public:
    ThoracicSMCEnsembleBase(double pars[]):
      m_matrix(0),
      m_elastin(pars[1], 0.0, 0.0, 0.5),
      m_muscle(pars[2], pars[3], 0.0, 0.5)
    {}

    ThoracicSMCEnsembleBase(double pars[], double Cmax[]):
      m_matrix(0),
      m_elastin(pars[1], 0.0, 0.0, 0.5),
      m_muscle(pars[2], pars[3], 0.0, 0.5, Cmax)
    {}

    ~ThoracicSMCEnsembleBase() {}

    void get_scaled_pars(double pars[]);

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);

  };


  class ThoracicSMCEnsembleVE: public ThoracicSMCEnsembleBase
  {
  protected:
    constitutive_models::FractionalVE<4> muscle;

  public:
    ThoracicSMCEnsembleVE(double pars[], double visco[], double Tf):
      ThoracicSMCEnsembleBase(pars),
      muscle(m_muscle, visco[0], Tf)
    {}

    ThoracicSMCEnsembleVE(double pars[], double visco[], double Tf, double Cmax[]):
      ThoracicSMCEnsembleBase(pars, Cmax),
      muscle(m_muscle, visco[0], Tf)
    {}

    ~ThoracicSMCEnsembleVE() {}

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);

  };


}
