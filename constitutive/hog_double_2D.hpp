#pragma once

#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

class HOGDouble2D : public MatLaw<4> {
  public:
    double k1, k2;
    double m4[4], m6[4];
    double E1;
    double E2;

    HOGDouble2D();
    ~HOGDouble2D();
    HOGDouble2D(double k1, double k2, double theta, double alpha);
    HOGDouble2D(double k1, double k2, double theta, double alpha, double Cmax[]);

    void set_pars(double k1, double k2, double theta, double alpha);
    void set_pars(double k1, double k2, double theta, double alpha, double Cmax[]);
    double get_scaled_modulus();
    double stress(const kinematics::kinematics<4> &kin, double stress[4]);
    void stress(double args[4], double stress[4]);
};

} // namespace constitutive_models
