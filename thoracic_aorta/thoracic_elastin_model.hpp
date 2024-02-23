#pragma once

#include "../constitutive/neohookean.hpp"
#include "../constitutive/planar_hog.hpp"
#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace thoracic {

class PlanarElastinMatrix : public constitutive_models::MatLawTime<4> {
  protected:
    constitutive_models::NeoHookean m_matrix;
    constitutive_models::PlanarHog2D m_elastin;

  public:
    PlanarElastinMatrix(
        const double pars[], const double fiber[], const double visco[], double Tf,
        const double Cmax[]
    )
        : m_matrix(pars[0]), m_elastin(pars[1], 0.0, 0.0, 0.5) {
    }

    ~PlanarElastinMatrix() {
    }

    void get_scaled_pars(double pars[]);

    void stress(const kinematics::kinematics<4> &kin, double dt, double stress[]);
};

} // namespace thoracic
