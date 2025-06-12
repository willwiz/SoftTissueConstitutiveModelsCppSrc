#pragma once
#include "../../constitutive/hyperelasticity_1d/neohookean_1d.hpp"
#include "../../constitutive/hyperelasticity_1d/transverse_iso.hpp"
#include "../../constitutive/viscoelasticity/fractional.hpp"
#include "../../interfaces.hpp"
#include <iostream>

namespace optimization_1d {
class TrabeculaePassive : public constitutive_models::MatLawTime3D<1> {
  protected:
    constitutive_models::NeoHookean1D m_matrix;
    constitutive_models::TransverseIso1D m_muscle;
    constitutive_models::FractionalVE<1> muscle;

  public:
    TrabeculaePassive(double pars[], double fiber[], double visco[], double Tf)
        : m_matrix(pars[0]), m_muscle(pars[1], pars[2], fiber[0]),
          muscle(m_muscle, visco[0], Tf) {};
    TrabeculaePassive(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : m_matrix(pars[0]), m_muscle(pars[1], pars[2], fiber[0], Cmax[0]),
          muscle(m_muscle, visco[0], Tf) {
              // std::cout << "TrabeculaePassive initialized with parameters: " << std::endl;
              // std::cout << "Pars are " << pars[0] << " " << pars[1] << " " << pars[2] <<
              // std::endl;
              // std::cout << "Fiber is " << fiber[0] << std::endl;
              // std::cout << "Visco is " << visco[0] << std::endl;
              // std::cout << "Tf is " << Tf << std::endl;
          };
    ~TrabeculaePassive() {};

    void stress(const kinematics::kinematics<1> &kin, const double dt, double stress[]);
};

} // namespace optimization_1d
