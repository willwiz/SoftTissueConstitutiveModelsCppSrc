#pragma once

#define _USE_MATH_DEFINES

#include <cmath>
#include "../CTvalues_optimization.hpp"
#include "../interfaces.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/isotropic_fung.hpp"
#include "../constitutive/struc_hog_2D.hpp"
#include "../constitutive/fractional.hpp"


namespace sim {

  class ThoracicIsoBase: public constitutive_models::MatLawTime<4>
  {
  protected:
    constitutive_models::NeoHookean m_elastin;
    constitutive_models::fungIso m_muscle;
    constitutive_models::StrucHOG2D m_collagen;

  public:
    ThoracicIsoBase(double pars[], double fiber[]):
      m_elastin(pars[0]),
      m_muscle(pars[1], 0.5),
      m_collagen(pars[2], pars[3], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop)
    {}

    ThoracicIsoBase(double pars[], double fiber[], double Cmax[]):
      m_elastin(pars[0]),
      m_muscle(pars[1], 0.5),
      m_collagen(pars[2], pars[3], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop, Cmax)
    {}

    ~ThoracicIsoBase() {}

    void get_scaled_pars(double pars[]);

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);

  };


  class ThoracicIsoVEBase: public ThoracicIsoBase
  {
  protected:
    constitutive_models::FractionalVE<4> collagen;
    constitutive_models::FractionalVE<4> muscle;

  public:
    ThoracicIsoVEBase(double pars[], double fiber[], double visco[], double Tf):
      ThoracicIsoBase(pars, fiber),
      collagen(m_collagen, visco[0], Tf),
      muscle(m_muscle, visco[1], Tf)
    {}

    ThoracicIsoVEBase(double pars[], double fiber[], double visco[], double Tf, double Cmax[]):
      ThoracicIsoBase(pars, fiber, Cmax),
      collagen(m_collagen, visco[0], Tf),
      muscle(m_muscle, visco[1], Tf)
    {}

    ~ThoracicIsoVEBase() {}

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);

  };




  class ThoracicIsoHE: public ThoracicIsoBase
  {
  public:
    ThoracicIsoHE(double pars[], double fiber[], double visco[], double Tf, double Cmax[]):
      ThoracicIsoBase(pars, fiber) {}
    ~ThoracicIsoHE() {}
  };


  class ThoracicIsoHEScaled: public ThoracicIsoBase
  {
  public:
    ThoracicIsoHEScaled(double pars[], double fiber[], double visco[], double Tf, double Cmax[]):
      ThoracicIsoBase(pars, fiber, Cmax) {}
    ~ThoracicIsoHEScaled() {}
  };

  class ThoracicIsoVE: public ThoracicIsoVEBase
  {
  public:
    ThoracicIsoVE(double pars[], double fiber[], double visco[], double Tf, double Cmax[]):
      ThoracicIsoVEBase(pars, fiber, visco, Tf) {}
    ~ThoracicIsoVE() {}
  };


  class ThoracicIsoVEScaled: public ThoracicIsoVEBase
  {
  public:
    ThoracicIsoVEScaled(double pars[], double fiber[], double visco[], double Tf, double Cmax[]):
      ThoracicIsoVEBase(pars, fiber, visco, Tf, Cmax) {}
    ~ThoracicIsoVEScaled() {}
  };

}
