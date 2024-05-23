#pragma once

#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

class Hog2D : public MatLaw<4> {
  public:
    double k, b;
    double m[4];
    double z;
    double slope;

    Hog2D();
    ~Hog2D();
    Hog2D(double k1, double k2, double theta);
    Hog2D(double k1, double k2, double theta, const double Cmax[]);

    void set_pars(double k1, double k2, double theta);
    void set_pars(double k1, double k2, double theta, const double Cmax[]);
    double get_scaled_modulus();
    double stress(const kinematics::kinematics<4> &kin, double stress[4]);
    void stress(const double args[4], double stress[4]);
};

} // namespace constitutive_models
