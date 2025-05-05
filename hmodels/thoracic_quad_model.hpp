#pragma once

#include "../../interfaces.hpp"
#include "../CTvalues_optimization.hpp"
#include "../constitutive/fractional.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/quadratic_fiber.hpp"
#include "../constitutive/struc_hog_2D.hpp"

namespace sim {

class ThoracicQuadBase : public constitutive_models::MatLawTime<4> {
  protected:
    constitutive_models::NeoHookean m_elastin;
    constitutive_models::QuadraticFiber m_muscle;
    constitutive_models::StrucHOG2D m_collagen;

  public:
    ThoracicQuadBase(double pars[], double fiber[])
        : m_elastin(pars[0]), m_muscle(pars[1], fiber[0], fiber[1], -fiber[1]),
          m_collagen(pars[2], pars[3], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop) {
    }

    ThoracicQuadBase(double pars[], double fiber[], double Cmax[])
        : m_elastin(pars[0]), m_muscle(pars[1], fiber[0], fiber[1], -fiber[1]),
          m_collagen(
              pars[2], pars[3], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop, Cmax
          ) {
    }

    ~ThoracicQuadBase() {
    }

    void get_scaled_pars(double pars[]);

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class ThoracicQuadVEBase : public ThoracicQuadBase {
  protected:
    constitutive_models::FractionalVE<4> collagen;
    constitutive_models::FractionalVE<4> muscle;

  public:
    ThoracicQuadVEBase(double pars[], double fiber[], double visco[], double Tf)
        : ThoracicQuadBase(pars, fiber), collagen(m_collagen, visco[0], Tf),
          muscle(m_muscle, visco[1], Tf) {
    }

    ThoracicQuadVEBase(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : ThoracicQuadBase(pars, fiber, Cmax), collagen(m_collagen, visco[0], Tf),
          muscle(m_muscle, visco[1], Tf) {
    }

    ~ThoracicQuadVEBase() {
    }

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class ThoracicQuadVEScaled : public ThoracicQuadVEBase {
  public:
    ThoracicQuadVEScaled(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : ThoracicQuadVEBase(pars, fiber, visco, Tf, Cmax) {
    }
    ~ThoracicQuadVEScaled() {
    }
};

} // namespace sim
