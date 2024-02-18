#pragma once

#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

class PlanarHog2D : public MatLaw<4> {
  public:
    double k1, k2;
    double m[4];
    double E1;
    double E2;

    PlanarHog2D();
    ~PlanarHog2D();
    PlanarHog2D(double k1, double k2, double theta, double kappa);
    PlanarHog2D(double k1, double k2, double theta, double kappa, const double Cmax[]);

    void set_pars(double theta, double kappa);
    void set_pars(double theta, double kappa, const double Cmax[]);
    double get_scaled_modulus();
    double stress(const kinematics::kinematics<4> &kin, double stress[4]);
    void stress(const double args[4], double stress[4]);
};

} // namespace constitutive_models