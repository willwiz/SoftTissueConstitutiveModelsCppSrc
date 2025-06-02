#pragma once

#include "../../constitutive/hyperelasticity_1d/neohookean_1d.hpp"
#include "../../constitutive/hyperelasticity_1d/transverse_iso.hpp"
#include "../../constitutive/viscoelasticity/fractional.hpp"
#include "../../constitutive/viscoelasticity/maxwell.hpp"
#include "../../interfaces.hpp"

namespace optimization_1d {

class TrabeculaeCalciumOnly : public constitutive_models::MatLawTime3D<1> {
  protected:
    constitutive_models::NeoHookean1D matrix;
    constitutive_models::TransverseIso1D m_muscle;
    constitutive_models::FractionalVE<1> muscle;
    constitutive_models::TransverseIso1D m_titin;
    constitutive_models::FractionalVE<1> titin;
    constitutive_models::TransverseIso1D m_calcium;
    constitutive_models::MaxwellVE<1> calcium;

  public:
    TrabeculaeCalciumOnly(double pars[7], double fiber[], double visco[], double Tf)
        : matrix(pars[0]), m_muscle(pars[1], pars[4], fiber[0]), muscle(m_muscle, visco[0], Tf),
          m_titin(pars[2], pars[5], fiber[0]), titin(m_muscle, visco[1], Tf),
          m_calcium(pars[3], pars[6], fiber[0]), calcium(m_calcium, 1.0, visco[2]) {};
    TrabeculaeCalciumOnly(double pars[7], double fiber[], double visco[], double Tf, double Cmax[])
        : matrix(pars[0]), m_muscle(pars[1], pars[4], fiber[0], Cmax[0]),
          muscle(m_muscle, visco[0], Tf), m_titin(pars[2], pars[5], fiber[0], Cmax[0]),
          titin(m_muscle, visco[1], Tf), m_calcium(pars[3], pars[6], fiber[0], Cmax[0]),
          calcium(m_calcium, 1.0, visco[2]) {};
    ~TrabeculaeCalciumOnly() {};

    void stress(const kinematics::kinematics<1> &kin, const double dt, double stress[]);
};

} // namespace optimization_1d
