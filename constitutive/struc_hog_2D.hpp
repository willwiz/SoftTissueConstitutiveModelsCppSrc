#pragma once

#include "../interfaces.hpp"
#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

class StrucHOG2D : public MatLaw<4> {
  public:
    double k, b;
    double H4[4], H6[4];
    double H_33;
    double z;
    double low_slope;
    double high_slope;

    StrucHOG2D();
    ~StrucHOG2D();
    StrucHOG2D(
        double k1, double k2, double theta, double alpha, double beta, double kip, double kop
    );
    StrucHOG2D(
        double k1, double k2, double theta, double alpha, double beta, double kip, double kop,
        const double Cmax[]
    );

    void set_pars(
        double k1, double k2, double theta, double alpha, double beta, double kip, double kop
    );
    void set_pars(
        double k1, double k2, double theta, double alpha, double beta, double kip, double kop,
        const double Cmax[]
    );
    double get_scaled_modulus();
    double stress(const kinematics::kinematics<4> &kin, double stress[4]);
    void stress(const double args[4], double stress[4]);
};

} // namespace constitutive_models
