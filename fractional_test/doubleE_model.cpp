#include "doubleE_model.hpp"

namespace caputo_test
{
    void DoubleEVE::stress(const kinematics::kinematics<4> &kin, const double dt, double stress[])
    {
        double p;
        double mat[4];
        p = VE.stress(kin, dt, mat);
        for (int j = 0; j < 4; j++)
        {
            stress[j] = mat[j] - p * kin.I_n*kin.Cinv[j];
        }
    }


    void double_model_sim(
        double pars[], double visco[], double Tf, double strain[], double dt[],
        double stress_out[], int n
    ) {
        DoubleEVE model(pars, visco, Tf);
        kinematics::deformation2D kin;
        for (int i = 0; i < n; i++)
        {
            kin.precompute(&strain[4*i]);
            model.stress(kin, dt[i], &stress_out[4*i]);
        }

    }

} // namespace caputo_test
