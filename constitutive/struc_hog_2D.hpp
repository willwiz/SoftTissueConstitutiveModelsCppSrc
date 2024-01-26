#pragma once

#include "../kinematics/kinematics.hpp"
#include "../interfaces.hpp"

namespace constitutive_models {

  class StrucHOG2D: public MatLaw<4>
  {
  public:
    double k1, k2;
    double A, B, C;
    double m4[4], m6[4], H4[4], H6[4];
    double E1;
    double E2;

    StrucHOG2D();
    ~StrucHOG2D();
    StrucHOG2D(double k1, double k2, double theta, double alpha, double beta, double kip,
      double kop);
    StrucHOG2D(double k1, double k2, double theta, double alpha, double beta, double kip,
      double kop, double Cmax[]);

    void set_pars(double k1, double k2, double theta, double alpha, double beta, double kip,
      double kop);
    void set_pars(double k1, double k2, double theta, double alpha, double beta, double kip,
      double kop, double Cmax[]);
    double get_scaled_modulus();
    double stress(const kinematics::kinematics<4> &kin, double stress[4]);
    void stress(double args[4], double stress[4]);
  };

}
