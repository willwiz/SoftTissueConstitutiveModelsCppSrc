#pragma once

#include "../constitutive/DoubleE.hpp"
#include "../constitutive/fractional.hpp"
#include "../kinematics/kinematics.hpp"

namespace caputo_test {
class DoubleEVE : public constitutive_models::MatLawTime<4> {
  protected:
    constitutive_models::DoubleE2D HE;
    constitutive_models::FractionalVE<4> VE;

  public:
    DoubleEVE() : HE(0.0, 0.0, 1.0, 1.0, 1.0, 0.0), VE(HE, 0.1, 1.0) {
    }

    DoubleEVE(double pars[], double visco[], double Tf)
        : HE(pars[0], pars[1], pars[2], pars[3], pars[4], pars[5]), VE(HE, visco[0], Tf) {
    }

    ~DoubleEVE() {
    }

    void stress(const kinematics::kinematics<4> &kin, const double dt, double stress_out[]);
};

void double_model_sim(
    double pars[], double visco[], double Tf, double strain[], double dt[], double stress_out[],
    int n
);

} // namespace caputo_test
