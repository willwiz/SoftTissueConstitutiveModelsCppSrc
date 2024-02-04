#pragma once

#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

class PlanarHog3D : public MatLaw3D<9> {
  public:
    double k1;
    double k2;
    double E1;
    double E2;
    double m[9];

    PlanarHog3D();
    ~PlanarHog3D();
    PlanarHog3D(double k1, double k2, double theta, double kappa);

    void set_pars(double theta, double kappa);
    void stress(const kinematics::kinematics<9> &kin, double stress[9]);
};

} // namespace constitutive_models