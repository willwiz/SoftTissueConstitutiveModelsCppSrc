#pragma once

#define _USE_MATH_DEFINES
#include "../CTvalues_optimization.hpp"
#include "../constitutive/fractional.hpp"
#include "../interfaces.hpp"
#include "neohookean.hpp"
#include "planar_hog.hpp"
#include "struc_hog.hpp"
#include <cmath>

namespace sim {
class ModelMatrix : public constitutive_models::MatLawTime<9> {
  protected:
    constitutive_models::NeoHookean m_model;

  public:
    ModelMatrix(double pars[], double fiber[], double visco[], double Tf) : m_model(pars[0]) {
    }

    ~ModelMatrix() {
    }

    void stress(const kinematics::kinematics<9> &kin, const double dt, double stress[]);
};

class ModelElastin : public constitutive_models::MatLawTime<9> {
  protected:
    constitutive_models::PlanarHog3D m_model;

  public:
    ModelElastin(double pars[], double fiber[], double visco[], double Tf)
        : m_model(pars[0], 0, fiber[0], fiber[1]) {
    }

    ~ModelElastin() {
    }

    void stress(const kinematics::kinematics<9> &kin, const double dt, double stress[]);
};
class ModelSMC : public constitutive_models::MatLawTime<9> {
  protected:
    constitutive_models::PlanarHog3D m_model;

  public:
    ModelSMC(double pars[], double fiber[], double visco[], double Tf)
        : m_model(pars[0], pars[1], fiber[0], fiber[1]) {
    }

    ~ModelSMC() {
    }

    void stress(const kinematics::kinematics<9> &kin, const double dt, double stress[]);
};

class ModelCollagen : public constitutive_models::MatLawTime<9> {
  protected:
    constitutive_models::StrucHOG2D m_model;

  public:
    ModelCollagen(double pars[], double fiber[], double visco[], double Tf)
        : m_model(pars[0], pars[1], fiber[0], fiber[1], fiber[1], ctv::M_kip, ctv::M_kop) {
    }

    ~ModelCollagen() {
    }

    void stress(const kinematics::kinematics<9> &kin, const double dt, double stress[]);
};

class ModelSMCVE : public ModelSMC {
  protected:
    constitutive_models::FractionalVE3D<9> ve_model;

  public:
    ModelSMCVE(double pars[], double fiber[], double visco[], double Tf)
        : ModelSMC(pars, fiber, visco, Tf), ve_model(m_model, visco[0], Tf) {
    }
    ~ModelSMCVE() {
    }
    void stress(const kinematics::kinematics<9> &kin, const double dt, double stress[]);
};

class ModelCollagenVE : public ModelCollagen {
  protected:
    constitutive_models::FractionalVE3D<9> ve_model;

  public:
    ModelCollagenVE(double pars[], double fiber[], double visco[], double Tf)
        : ModelCollagen(pars, fiber, visco, Tf), ve_model(m_model, visco[0], Tf) {
    }
    ~ModelCollagenVE() {
    }
    void stress(const kinematics::kinematics<9> &kin, const double dt, double stress[]);
};

} // namespace sim
