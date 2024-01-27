#pragma once

#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

class StrucHOG2D : public MatLaw3D<9> {
  public:
    double k1, k2;
    double A, B, C;
    double H4[9], H6[9];
    double E1;
    double E2;

    StrucHOG2D();
    ~StrucHOG2D();
    StrucHOG2D(
        double k1, double k2, double theta, double alpha, double beta, double kip, double kop
    );

    void set_pars(
        double k1, double k2, double theta, double alpha, double beta, double kip, double kop
    );

    void stress(const kinematics::kinematics<9> &kin, double stress[9]);
};

} // namespace constitutive_models
