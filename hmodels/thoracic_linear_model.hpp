#pragma once

#include "../../interfaces.hpp"
#include "../CTvalues_optimization.hpp"
#include "../constitutive/fractional.hpp"
#include "../constitutive/linear_fiber.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/struc_hog_2D.hpp"

namespace sim {

class ThoracicLinBase : public constitutive_models::MatLawTime<4> {
  protected:
    constitutive_models::NeoHookean m_elastin;
    constitutive_models::LinearFiber m_muscle;
    constitutive_models::StrucHOG2D m_collagen;

  public:
    ThoracicLinBase(double pars[], double fiber[])
        : m_elastin(pars[0]), m_muscle(pars[1], fiber[0], fiber[1], -fiber[1]),
          m_collagen(pars[2], pars[3], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop) {
    }

    ThoracicLinBase(double pars[], double fiber[], double Cmax[])
        : m_elastin(pars[0]), m_muscle(pars[1], fiber[0], fiber[1], -fiber[1]),
          m_collagen(
              pars[2], pars[3], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop, Cmax
          ) {
    }

    ~ThoracicLinBase() {
    }

    void get_scaled_pars(double pars[]);

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class ThoracicLinVEBase : public ThoracicLinBase {
  protected:
    constitutive_models::FractionalVE<4> collagen;
    constitutive_models::FractionalVE<4> muscle;

  public:
    ThoracicLinVEBase(double pars[], double fiber[], double visco[], double Tf)
        : ThoracicLinBase(pars, fiber), collagen(m_collagen, visco[0], Tf),
          muscle(m_muscle, visco[1], Tf) {
    }

    ThoracicLinVEBase(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : ThoracicLinBase(pars, fiber, Cmax), collagen(m_collagen, visco[0], Tf),
          muscle(m_muscle, visco[1], Tf) {
    }

    ~ThoracicLinVEBase() {
    }

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class ThoracicLinVEScaled : public ThoracicLinVEBase {
  public:
    ThoracicLinVEScaled(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : ThoracicLinVEBase(pars, fiber, visco, Tf, Cmax) {
    }
    ~ThoracicLinVEScaled() {
    }
};

} // namespace sim
