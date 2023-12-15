#pragma once

#include "../constitutive/fractional.hpp"
#include "../constitutive_laws_1d/neohookean_1d.hpp"
#include "../constitutive_laws_1d/transverse_hog_1d.hpp"
#include "../constitutive_laws_1d/transverse_iso.hpp"
#include "../interfaces.hpp"

namespace optimization_1d {
class Trabeculae1D : public constitutive_models::MatLawTime<1> {
  protected:
    constitutive_models::NeoHookean1D m_matrix;
    constitutive_models::TransverseIso1D m_muscle;
    constitutive_models::FractionalVE<1> muscle;

  public:
    Trabeculae1D(double pars[], double fiber[], double visco[], double Tf)
        : m_matrix(pars[0]), m_muscle(pars[1], pars[2], fiber[0]), muscle(m_muscle, visco[0], Tf){};
    Trabeculae1D(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : m_matrix(pars[0]), m_muscle(pars[1], pars[2], fiber[0], Cmax[0]),
          muscle(m_muscle, visco[0], Tf){};
    ~Trabeculae1D(){};

    void stress(const kinematics::kinematics<1> &kin, const double dt, double stress[]);
};

class Trabeculae2Phase1D : public constitutive_models::MatLawTime<1> {
  protected:
    constitutive_models::NeoHookean1D m_matrix;
    constitutive_models::TransverseIso1D m_muscle;
    constitutive_models::TransverseHog1D m_collagen;
    constitutive_models::FractionalVE<1> muscle;
    constitutive_models::FractionalVE<1> collagen;

  public:
    Trabeculae2Phase1D(double pars[], double fiber[], double visco[], double Tf)
        : m_matrix(pars[0]), m_muscle(pars[1], pars[2], fiber[0]),
          m_collagen(pars[3], pars[4], fiber[1]), muscle(m_muscle, visco[0], Tf),
          collagen(m_collagen, visco[1], Tf){};
    Trabeculae2Phase1D(double pars[], double fiber[], double visco[], double Tf, double Cmax[])
        : m_matrix(pars[0]), m_muscle(pars[1], pars[2], fiber[0]),
          m_collagen(pars[3], pars[4], fiber[1], Cmax[0]), muscle(m_muscle, visco[0], Tf),
          collagen(m_collagen, visco[1], Tf){};
    ~Trabeculae2Phase1D(){};

    void stress(const kinematics::kinematics<1> &kin, const double dt, double stress[]);
};

void simulation_function_trabeculae_1d(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double dt[], double stress_out[], int n
);

double objective_function_trabeculae_stacked_1d(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double stress[], double dt[], double weight[], int n[]
);

double objective_function_trabeculae_1d(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double stress[], double dt[], double weight[], int n, int nset, int index[], int select[]
);

void simulation_function_trabeculae_2phase_1d(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double dt[], double stress_out[], int n
);

double objective_function_trabeculae_2phase_1d(
    double pars[], double fiber[], double visco[], double Tf, double Cmax[], double strain[],
    double stress[], double dt[], double weight[], int n[]
);
} // namespace optimization_1d
