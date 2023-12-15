#define _USE_MATH_DEFINES

#include <cmath>
#include "femoral_model.hpp"

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
    FemoralBase::FemoralBase(double pars[], double fiber[])
        :   m_matrix(pars[0]),
            m_collagen(pars[1], pars[2], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop),
            m_elastin(pars[3], pars[4], fiber[0]),
            m_muscle(pars[5], pars[6], fiber[0] + M_PI_2)
    {}
    FemoralBase::FemoralBase(double pars[], double fiber[], double Cmax[])
        :   m_matrix(pars[0]),
            m_collagen(pars[1], pars[2], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop, Cmax),
            m_elastin(pars[3], pars[4], fiber[0], Cmax),
            m_muscle(pars[5], pars[6], fiber[0] + M_PI_2, Cmax)
    {}

    FemoralBase::~FemoralBase() {}


    void FemoralBase::get_scaled_pars(double pars[])
    {
        pars[0] = m_matrix.mu;
        pars[1] = m_collagen.get_scaled_modulus();
        pars[2] = m_collagen.k2;
        pars[3] = m_elastin.get_scaled_modulus();
        pars[4] = m_elastin.k2;
        pars[5] = m_muscle.get_scaled_modulus();
        pars[6] = m_muscle.k2;
    }


    void FemoralBase::stress(const kinematics::kinematics<4> &kin, const double dt, double stress[])
    {
        double p = 0.0;
        double iso[4], el[4], smc[4], col[4];
        p = m_matrix.stress(kin, iso);
        (void) m_elastin.stress(kin, el);
        (void) m_muscle.stress(kin, smc);
        p = p + m_collagen.stress(kin, col);
        for (int j = 0; j < ctv::prob_dim; j++)
        {
            stress[j] = iso[j] + col[j] + el[j] + smc[j] - p * kin.I_n*kin.Cinv[j];
        }
    }





    FemoralVEBase::FemoralVEBase(double pars[], double fiber[], double visco[], double Tf):
        FemoralBase(pars, fiber),
        collagen(m_collagen, visco[0], Tf),
        muscle(m_muscle, visco[1], Tf)
    {}

    FemoralVEBase::FemoralVEBase(double pars[], double fiber[], double visco[], double Tf, double Cmax[]):
        FemoralBase(pars, fiber, Cmax),
        collagen(m_collagen, visco[0], Tf),
        muscle(m_muscle, visco[1], Tf)
    {}
    FemoralVEBase::~FemoralVEBase() {}


    void FemoralVEBase::stress(const kinematics::kinematics<4> &kin, const double dt, double stress[])
    {
        double p = 0.0;
        double iso[4], el[4], smc[4], col[4];
        p = m_matrix.stress(kin, iso);
        (void) m_elastin.stress(kin, el);
        p = p + muscle.stress(kin, dt, smc);
        p = p + collagen.stress(kin, dt, col);
        for (int j = 0; j < ctv::prob_dim; j++)
        {
            stress[j] = iso[j] + col[j] + el[j] + smc[j] - p * kin.I_n*kin.Cinv[j];
        }
    }
}