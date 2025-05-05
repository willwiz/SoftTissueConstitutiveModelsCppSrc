#pragma once

#define _USE_MATH_DEFINES

#include "../../interfaces.hpp"
#include "../CTvalues_optimization.hpp"
#include "../constitutive/fractional.hpp"
#include "../constitutive/hog_2D.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/planar_hog.hpp"
#include "../constitutive/struc_hog_2D.hpp"
#include <cmath>

namespace sim {

class ThoracicDefaultBase : public constitutive_models::MatLawTime<4> {
  protected:
    constitutive_models::NeoHookean m_matrix;
    constitutive_models::PlanarHog2D m_elastin;
    constitutive_models::PlanarHog2D m_muscle;
    constitutive_models::StrucHOG2D m_collagen;

  public:
    ThoracicDefaultBase(double pars[], double fiber[])
        : m_matrix(pars[0]), m_elastin(pars[1], 0.0, fiber[0], 0.5),
          m_muscle(pars[2], pars[3], fiber[0], fiber[1]),
          m_collagen(pars[4], pars[5], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop) {
    }

    ThoracicDefaultBase(double pars[], double fiber[], double Cmax[])
        : m_matrix(pars[0]), m_elastin(pars[1], 0.0, fiber[0], fiber[1]),
          m_muscle(pars[2], pars[3], fiber[0], fiber[2], Cmax),
          m_collagen(
              pars[4], pars[5], fiber[0], fiber[3], -fiber[3], ctv::M_kip, ctv::M_kop, Cmax
          ) {
    }

    ~ThoracicDefaultBase() {
    }

    void get_scaled_pars(double pars[]);

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class ThoracicDefaultVEBase : public ThoracicDefaultBase {
  protected:
    constitutive_models::FractionalVE<4> muscle;
    constitutive_models::FractionalVE<4> collagen;

  public:
    ThoracicDefaultVEBase(double pars[], double fiber[], double visco[], double Tf)
        : ThoracicDefaultBase(pars, fiber), muscle(m_muscle, visco[0], Tf),
          collagen(m_collagen, visco[1], Tf) {
    }

    ThoracicDefaultVEBase(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : ThoracicDefaultBase(pars, fiber, Cmax), muscle(m_muscle, visco[0], Tf),
          collagen(m_collagen, visco[1], Tf) {
    }

    ~ThoracicDefaultVEBase() {
    }

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class ThoracicDefaultVEStandard : public ThoracicDefaultVEBase {
  public:
    ThoracicDefaultVEStandard(
        double pars[], double fiber[], double visco[], double Tf, double Cmax[]
    )
        : ThoracicDefaultVEBase(pars, fiber, visco, Tf) {
    }
    ~ThoracicDefaultVEStandard() {
    }
};

class ThoracicDefaultVEScaled : public ThoracicDefaultVEBase {
  public:
    ThoracicDefaultVEScaled(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : ThoracicDefaultVEBase(pars, fiber, visco, Tf, Cmax) {
    }
    ~ThoracicDefaultVEScaled() {
    }
};

} // namespace sim
