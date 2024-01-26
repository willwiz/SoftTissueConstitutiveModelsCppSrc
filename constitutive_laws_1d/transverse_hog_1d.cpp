#include <algorithm>
#include <cmath>
#include <iostream>
#include "../kinematics/tensor_algebra.hpp"
#include "../kinematics/kinematics.hpp"
#include "transverse_hog_1d.hpp"

/*----------------------------------------------------------------------
 |  This file provides the holzapfel class of models in 2D. Single and
 |  dual family models are available
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

namespace constitutive_models
{

    /*******************************************************************************
     * Standard holzapfel ogden with one fiber family
     *
     * COMMENTS:
     *
     *******************************************************************************/
    TransverseHog1D::TransverseHog1D() : E1{1}, E2{1} {};
    TransverseHog1D::~TransverseHog1D(){};

    TransverseHog1D::TransverseHog1D(double k1, double k2, double kappa)
        : k1(k1), k2(k2), fiber{1 - kappa + kappa / 3.0}, iso{kappa / 3.0}, E1{1}, E2{1}
    {
        this->set_pars(kappa);
    };

    TransverseHog1D::TransverseHog1D(double k1, double k2, double kappa, double Cmax)
        : k1(k1), k2(k2), fiber{1 - kappa + kappa / 3.0}, iso{kappa / 3.0}
    {
        this->set_pars(kappa, Cmax);
    };

    inline double calc_I4(double lambda, double iso, double fiber)
    {
        return lambda * fiber + 2 * iso / std::sqrt(lambda) - 1.0;
    }

    void TransverseHog1D::set_pars(double kappa)
    {
    }

    // Scaled version
    void TransverseHog1D::set_pars(double kappa, double Cmax)
    {
        this->set_pars(kappa);
        E1 = calc_I4(Cmax, iso, fiber);
        k1 = k1 / E1;
        E2 = E1 * E1;
    }

    // Stress functions
    double TransverseHog1D::stress(const kinematics::kinematics<1> &kin, double stress[])
    {

        double I_4 = kin.C[0] * fiber + 2 * iso * kin.I_n - 1.0;
        double dWd4 = k1 * I_4 * exp(k2 * (I_4 * I_4 - E2));

        stress[0] = dWd4 * fiber;
        // std::cout << kin.C[0] << " " << kin.I_n << "\n";
        // std::cout << dWd4 << " " << I_4 << " " << E2 << " " << E1 << "\n";

        return dWd4 * iso;
    }

    double TransverseHog1D::stress(double args)
    {
        double stress[1];
        kinematics::deformation1D kin(args);
        double p = this->stress(kin, stress);
        return stress[0] - p * kin.I_n * kin.Cinv[0];
    }

    /*------------------------------------------------------------------------------
     |  THE END
     -----------------------------------------------------------------------------*/

}