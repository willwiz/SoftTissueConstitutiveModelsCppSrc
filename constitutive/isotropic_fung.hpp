#pragma once

#include "../kinematics/kinematics.hpp"
#include "../interfaces.hpp"

namespace constitutive_models {

  class fungIso: public MatLaw<4>
  {
    public:
      double m_mu;
      double m_b;
      double m_Emax;
      double exponent;

      fungIso () {};
      fungIso(double mu, double b);
      fungIso(double mu, double b, double Cmax[]);
      ~fungIso () {};
      void set_pars(double mu, double b);
      void set_pars(double mu, double b, double Cmax[]);
      double get_scaled_modulus();
      double stress(const kinematics::kinematics<4> &kin, double stress[]);
      void stress(double args[], double stress[]);
  };
}
