#pragma once

#include "../../constitutive/hyperelasticity_1d/neohookean_1d.hpp"
#include "../../constitutive/hyperelasticity_1d/transverse_iso.hpp"
#include "../../constitutive/viscoelasticity/fractional.hpp"
#include "../../constitutive/viscoelasticity/maxwell.hpp"
#include "../../interfaces.hpp"

namespace optimization_1d {
class TrabeculaeCalcium : public constitutive_models::MatLawTime3D<1> {
  protected:
    constitutive_models::NeoHookean1D m_matrix;
    constitutive_models::TransverseIso1D m_muscle;
    constitutive_models::FractionalVE<1> muscle;
    constitutive_models::MaxwellVE<1> calcium;

  public:
    TrabeculaeCalcium(double pars[], double fiber[], double visco[], double Tf)
        : m_matrix(pars[0]), m_muscle(pars[1], pars[2], fiber[0]), muscle(m_muscle, visco[0], Tf),
          calcium(m_muscle, pars[3], visco[1]) {};
    TrabeculaeCalcium(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : m_matrix(pars[0]), m_muscle(pars[1], pars[2], fiber[0], Cmax[0]),
          muscle(m_muscle, visco[0], Tf), calcium(m_muscle, pars[3], visco[1]) {};
    ~TrabeculaeCalcium() {};

    void stress(const kinematics::kinematics<1> &kin, const double dt, double stress[]);
};

} // namespace optimization_1d
