#pragma once

#define _USE_MATH_DEFINES

#include "../../constitutive/hyperelasticity/hog_2D.hpp"
#include "../../constitutive/hyperelasticity/neohookean.hpp"
#include "../../constitutive/hyperelasticity/struc_hog_2D.hpp"
#include "../../constitutive/viscoelasticity/fractional.hpp"
#include "../../interfaces.hpp"
#include "../kinematics/kinematics.hpp"
#include "constants.hpp"
#include <cmath>

namespace femoral {

class FemoralBase : public constitutive_models::MatLawTime3D<4> {
  protected:
    constitutive_models::NeoHookean m_matrix;
    constitutive_models::Hog2D m_elastin;
    constitutive_models::Hog2D m_muscle;
    constitutive_models::StrucHOG2D m_collagen;

  public:
    FemoralBase(const double pars[], const double fiber[]);
    FemoralBase(const double pars[], const double fiber[], const double Cmax[]);
    ~FemoralBase();

    void get_scaled_pars(double pars[]);
    virtual void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class FemoralVEBase : public FemoralBase {
  protected:
    constitutive_models::FractionalVE<4> muscle;
    constitutive_models::FractionalVE<4> collagen;

  public:
    FemoralVEBase(const double pars[], const double fiber[], const double visco[], double Tf);
    FemoralVEBase(
        const double pars[], const double fiber[], const double visco[], double Tf,
        const double Cmax[]
    );
    ~FemoralVEBase();
    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class FemoralHE : public FemoralBase {
  public:
    FemoralHE(
        const double pars[], const double fiber[], const double visco[], double Tf,
        const double Cmax[]
    )
        : FemoralBase(pars, fiber) {
    }
    ~FemoralHE() {
    }
};

class FemoralHEScaled : public FemoralBase {
  public:
    FemoralHEScaled(
        const double pars[], const double fiber[], const double visco[], double Tf,
        const double Cmax[]
    )
        : FemoralBase(pars, fiber, Cmax) {
    }
    ~FemoralHEScaled() {
    }
};

class FemoralVE : public FemoralVEBase {
  public:
    FemoralVE(
        const double pars[], const double fiber[], const double visco[], double Tf,
        const double Cmax[]
    )
        : FemoralVEBase(pars, fiber, visco, Tf) {
    }
    ~FemoralVE() {
    }
};

class FemoralVEScaled : public FemoralVEBase {
  public:
    FemoralVEScaled(
        const double pars[], const double fiber[], const double visco[], double Tf,
        const double Cmax[]
    )
        : FemoralVEBase(pars, fiber, visco, Tf, Cmax) {
    }
    ~FemoralVEScaled() {
    }
};

} // namespace femoral
