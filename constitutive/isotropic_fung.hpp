#pragma once

#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

class fungIso : public MatLaw<4> {
  public:
    double m_mu;
    double m_b;
    double m_Emax;
    double exponent;

    fungIso(){};
    fungIso(double mu, double b);
    fungIso(double mu, double b, const double Cmax[]);
    ~fungIso(){};
    void set_pars(double mu, double b);
    void set_pars(double mu, double b, const double Cmax[]);
    double get_scaled_modulus();
    double stress(const kinematics::kinematics<4> &kin, double stress[]);
    void stress(const double args[], double stress[]);
};
} // namespace constitutive_models
