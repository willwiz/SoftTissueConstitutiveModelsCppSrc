#pragma once
#include "../../constitutive/hyperelasticity/hog_double_2D.hpp"
#include "../../constitutive/hyperelasticity/neohookean.hpp"
#include "../../constitutive/hyperelasticity/planar_hog.hpp"
#include "../../constitutive/hyperelasticity/struc_hog_2D.hpp"
#include "../../constitutive/viscoelasticity/fractional.hpp"
#include "../../interfaces.hpp"
#include "../kinematics/kinematics.hpp"
#include "constants.hpp"
#include <cmath>

namespace thoracic {

class ThoracicEnsembleBase : public constitutive_models::MatLawTime3D<4> {
  protected:
    constitutive_models::NeoHookean m_matrix;
    constitutive_models::PlanarHog2D m_elastin;
    constitutive_models::HOGDouble2D m_muscle;
    constitutive_models::StrucHOG2D m_collagen;

  public:
    ThoracicEnsembleBase(const double pars[]);
    ThoracicEnsembleBase(const double pars[], const double Cmax[]);

    ~ThoracicEnsembleBase() {
    }

    void get_scaled_pars(double pars[6]);

    void stress(const kinematics::kinematics<4> &kin, double dt, double stress[]);
};

class ThoracicEnsembleVE : public ThoracicEnsembleBase {
  protected:
    constitutive_models::FractionalVE<4> muscle;
    constitutive_models::FractionalVE<4> collagen;

  public:
    ThoracicEnsembleVE(const double pars[], const double visco[], double Tf)
        : ThoracicEnsembleBase(pars), muscle(m_muscle, visco[0], Tf),
          collagen(m_collagen, visco[1], Tf) {
    }

    ThoracicEnsembleVE(const double pars[], const double visco[], double Tf, const double Cmax[])
        : ThoracicEnsembleBase(pars, Cmax), muscle(m_muscle, visco[0], Tf),
          collagen(m_collagen, visco[1], Tf) {
    }

    ~ThoracicEnsembleVE() {
    }

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

} // namespace thoracic
