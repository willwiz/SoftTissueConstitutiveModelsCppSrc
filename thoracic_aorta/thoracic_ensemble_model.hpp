#pragma once
#include "../constitutive/fractional.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/planar_hog.hpp"
#include "../constitutive/struc_hog_2D.hpp"
#include "../interfaces.hpp"
#include <cmath>

namespace thoracic {

class ThoracicFullEnsembleBase : public constitutive_models::MatLawTime<4> {
  protected:
    constitutive_models::NeoHookean m_matrix;
    constitutive_models::PlanarHog2D m_elastin;
    constitutive_models::PlanarHog2D m_muscle;
    constitutive_models::StrucHOG2D m_collagen;

  public:
    ThoracicFullEnsembleBase(double pars[]);
    ThoracicFullEnsembleBase(double pars[], double Cmax[]);

    ~ThoracicFullEnsembleBase() {
    }

    void get_scaled_pars(double pars[]);

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

class ThoracicFullEnsembleVE : public ThoracicFullEnsembleBase {
  protected:
    constitutive_models::FractionalVE<4> muscle;
    constitutive_models::FractionalVE<4> collagen;

  public:
    ThoracicFullEnsembleVE(double pars[], double visco[], double Tf)
        : ThoracicFullEnsembleBase(pars), muscle(m_muscle, visco[0], Tf),
          collagen(m_collagen, visco[1], Tf) {
    }

    ThoracicFullEnsembleVE(double pars[], double visco[], double Tf, double Cmax[])
        : ThoracicFullEnsembleBase(pars, Cmax), muscle(m_muscle, visco[0], Tf),
          collagen(m_collagen, visco[1], Tf) {
    }

    ~ThoracicFullEnsembleVE() {
    }

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress[]);
};

} // namespace thoracic
