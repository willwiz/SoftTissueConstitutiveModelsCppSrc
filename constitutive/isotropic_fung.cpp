#include <cmath>
#include "../kinematics/tensor_algebra.hpp"
#include "../kinematics/kinematics.hpp"
#include "isotropic_fung.hpp"

namespace constitutive_models {

/*----------------------------------------------------------------------
 |  Standard Neohookean, no much to say
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

  // Fung form
  fungIso::fungIso(double mu, double b): m_Emax(3.0)
  {
    this -> set_pars(mu, b);
  };

  fungIso::fungIso(double mu, double b, double Cmax[])
  {
    this -> set_pars(mu, b, Cmax);
  }


  void fungIso::set_pars(double mu, double b) {
    m_mu = mu;
    m_b  = b;
  }


  void fungIso::set_pars(double mu, double b, double Cmax[]) {
    this->set_pars(mu, b);
    double det = Cmax[0]*Cmax[3] - Cmax[1]*Cmax[1];
    double I_n = 1 / det;
    m_Emax = Cmax[0] + Cmax[3] + I_n;
  }


  double fungIso::get_scaled_modulus()
  {
    return m_mu * exp(-m_b*m_Emax);
  }


  double fungIso::stress(const kinematics::kinematics<4> &kin, double stress[]){

    exponent = m_mu*exp(m_b*(kin.I_1 - m_Emax));

    for (int i = 0; i < 4; i++)
    {
      stress[i] = exponent*id2d[i];
    }
    return exponent;
  }

  void fungIso::stress(double args[], double stress[]){

    kinematics::deformation2D kin(args);
    double p = this->stress(kin, stress);
    for (int i = 0; i < 4; i++)
    {
      stress[i] = stress[i] - p*kin.I_n*kin.Cinv[i];
    }
  }
}
