#pragma once

#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

class mooney : public MatLaw<4> {
  public:
    double m_mu;
    double m_b;
    mooney(){};
    mooney(double mu, double b);
    ~mooney(){};
    void set_pars(double mu, double b);
    double get_scaled_modulus();
    double stress(const kinematics::kinematics<4> &kin, double stress[]);
    void stress(const double args[], double stress[]);
};
} // namespace constitutive_models
