#pragma once

#define _USE_MATH_DEFINES

#include "../constitutive/fractional.hpp"
#include "../constitutive/hog_double_2D.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/planar_hog.hpp"
#include "../constitutive/struc_hog_2D.hpp"
#include "../interfaces.hpp"
#include "constants.hpp"
#include <cmath>

namespace thoracic {

class ThoracicBase : public constitutive_models::MatLawTime<4> {
  protected:
    constitutive_models::NeoHookean m_matrix;
    constitutive_models::PlanarHog2D m_elastin;
    constitutive_models::HOGDouble2D m_muscle;
    constitutive_models::StrucHOG2D m_collagen;

  public:
    ThoracicBase(double pars[], double fiber[])
        : m_matrix(pars[0]), m_elastin(pars[1], 0.0, fiber[0], 0.5),
          m_muscle(pars[2], pars[3], fiber[0], fiber[1]),
          m_collagen(pars[4], pars[5], fiber[0], fiber[1], -fiber[1], M_kip, M_kop) {
    }

    ThoracicBase(double pars[], double fiber[], double Cmax[])
        : m_matrix(pars[0]), m_elastin(pars[1], 0.0, fiber[0], 0.5),
          m_muscle(pars[2], pars[3], fiber[0], fiber[1], Cmax),
          m_collagen(pars[4], pars[5], fiber[0], fiber[1], -fiber[1], M_kip, M_kop, Cmax) {
    }

    ~ThoracicBase() {
    }

    void get_scaled_pars(double pars[]);

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class ThoracicVEBase : public ThoracicBase {
  protected:
    constitutive_models::FractionalVE<4> muscle;
    constitutive_models::FractionalVE<4> collagen;

  public:
    ThoracicVEBase(double pars[], double fiber[], double visco[], double Tf)
        : ThoracicBase(pars, fiber), muscle(m_muscle, visco[0], Tf),
          collagen(m_collagen, visco[1], Tf) {
    }

    ThoracicVEBase(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : ThoracicBase(pars, fiber, Cmax), muscle(m_muscle, visco[0], Tf),
          collagen(m_collagen, visco[1], Tf) {
    }

    ~ThoracicVEBase() {
    }

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class ThoracicVE : public ThoracicVEBase {
  public:
    ThoracicVE(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : ThoracicVEBase(pars, fiber, visco, Tf) {
    }
    ~ThoracicVE() {
    }
};

class ThoracicVEScaled : public ThoracicVEBase {
  public:
    ThoracicVEScaled(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : ThoracicVEBase(pars, fiber, visco, Tf, Cmax) {
    }
    ~ThoracicVEScaled() {
    }
};

} // namespace thoracic
