#pragma once

#define _USE_MATH_DEFINES

#include "../../CTvalues_optimization.hpp"
#include "../../interfaces.hpp"
#include "../constitutive/hyperelasticity/isotropic_fung.hpp"
#include "../constitutive/hyperelasticity/mooney.hpp"
#include "../constitutive/hyperelasticity/struc_hog_2D.hpp"
#include "../constitutive/viscoelasticity/fractional.hpp"
#include <cmath>

namespace sim {

class ThoracicDefaultFungBase : public constitutive_models::MatLawTime<4> {
  protected:
    constitutive_models::fungIso m_elastin;
    constitutive_models::StrucHOG2D m_collagen;

  public:
    ThoracicDefaultFungBase(double pars[], double fiber[])
        : m_elastin(pars[0], 1.0),
          m_collagen(pars[1], pars[2], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop) {
    }

    ThoracicDefaultFungBase(double pars[], double fiber[], double Cmax[])
        : m_elastin(pars[0], 1.0),
          m_collagen(
              pars[1], pars[2], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop, Cmax
          ) {
    }

    ~ThoracicDefaultFungBase() {
    }

    void get_scaled_pars(double pars[]);

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class ThoracicDefaultFungVEBase : public ThoracicDefaultFungBase {
  protected:
    constitutive_models::FractionalVE<4> collagen;

  public:
    ThoracicDefaultFungVEBase(double pars[], double fiber[], double visco[], double Tf)
        : ThoracicDefaultFungBase(pars, fiber), collagen(m_collagen, visco[0], Tf) {
    }

    ThoracicDefaultFungVEBase(
        double pars[], double fiber[], double visco[], double Tf, double Cmax[]
    )
        : ThoracicDefaultFungBase(pars, fiber, Cmax), collagen(m_collagen, visco[0], Tf) {
    }

    ~ThoracicDefaultFungVEBase() {
    }

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class ThoracicDefaultFungVEScaled : public ThoracicDefaultFungVEBase {
  public:
    ThoracicDefaultFungVEScaled(
        double pars[], double fiber[], double visco[], double Tf, double Cmax[]
    )
        : ThoracicDefaultFungVEBase(pars, fiber, visco, Tf, Cmax) {
    }
    ~ThoracicDefaultFungVEScaled() {
    }
};

} // namespace sim
