#pragma once

#define _USE_MATH_DEFINES

#include "../constitutive/fractional.hpp"
#include "../constitutive/hog_2D.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/struc_hog_2D.hpp"
#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"
#include "constants.hpp"
#include <cmath>

namespace sim {

class FemoralBase : public constitutive_models::MatLawTime<4> {
  protected:
    constitutive_models::NeoHookean m_matrix;
    constitutive_models::StrucHOG2D m_collagen;
    constitutive_models::Hog2D m_elastin;
    constitutive_models::Hog2D m_muscle;

  public:
    FemoralBase(double pars[], double fiber[]);
    FemoralBase(double pars[], double fiber[], double Cmax[]);
    ~FemoralBase();

    void get_scaled_pars(double pars[]);
    virtual void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class FemoralVEBase : public FemoralBase {
  protected:
    constitutive_models::FractionalVE<4> collagen;
    constitutive_models::FractionalVE<4> muscle;

  public:
    FemoralVEBase(double pars[], double fiber[], double visco[], double Tf);
    FemoralVEBase(double pars[], double fiber[], double visco[], double Tf, double Cmax[]);
    ~FemoralVEBase();
    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class FemoralHE : public FemoralBase {
  public:
    FemoralHE(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : FemoralBase(pars, fiber) {
    }
    ~FemoralHE() {
    }
};

class FemoralHEScaled : public FemoralBase {
  public:
    FemoralHEScaled(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : FemoralBase(pars, fiber, Cmax) {
    }
    ~FemoralHEScaled() {
    }
};

class FemoralVE : public FemoralVEBase {
  public:
    FemoralVE(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : FemoralVEBase(pars, fiber, visco, Tf) {
    }
    ~FemoralVE() {
    }
};

class FemoralVEScaled : public FemoralVEBase {
  public:
    FemoralVEScaled(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : FemoralVEBase(pars, fiber, visco, Tf, Cmax) {
    }
    ~FemoralVEScaled() {
    }
};

} // namespace sim
