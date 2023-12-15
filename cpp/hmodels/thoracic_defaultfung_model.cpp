#define _USE_MATH_DEFINES

#include <cmath>
#include "thoracic_defaultfung_model.hpp"

/*----------------------------------------------------------------------
 |  This file provides the definitions of the different model forms
 |  which combines the primative constitive model from the other cpp files
 |  in the constitutive_models namespace
 |
 |  These models are mainly used for the collaboration with Alexey Kamenskiy
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/


using namespace constitutive_models;

namespace sim {


/*----------------------------------------------------------------------
 |  This provides the main models in the full constitutive model
 -----------------------------------------------------------------------*/

    void ThoracicDefaultFungBase::get_scaled_pars(double pars[])
    {
        pars[0] = m_elastin.m_mu;
        pars[1] = m_collagen.get_scaled_modulus();
        pars[2] = m_collagen.k2;
    }


    void ThoracicDefaultFungBase::stress(const kinematics::kinematics<4> &kin, const double dt, double stress[])
    {
        double p = 0.0;
        double el[4], col[4];
        p = p + m_elastin.stress(kin, el);
        p = p + m_collagen.stress(kin, col);
        for (int j = 0; j < ctv::prob_dim; j++)
        {
            stress[j] = el[j] + col[j] - p * kin.I_n*kin.Cinv[j];
        }
    }


    void ThoracicDefaultFungVEBase::stress(const kinematics::kinematics<4> &kin, const double dt, double stress[])
    {
        double p = 0.0;
        double el[4], col[4];
        p = p + m_elastin.stress(kin, el);
        p = p + collagen.stress(kin, dt, col);
        for (int j = 0; j < ctv::prob_dim; j++)
        {
            stress[j] = el[j] + col[j] - p * kin.I_n*kin.Cinv[j];
        }
    }

}