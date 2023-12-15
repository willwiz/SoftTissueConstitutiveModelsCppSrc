// #include <stdio.h>
// #include <iostream>
// #include <stdlib.h>
// #include <algorithm>
#include "objective_template.hpp"
#include "objective.hpp"
#include "models.hpp"

namespace residuals {

  using namespace sim;
/* ------------------------------------------------------------------------------
 |  This file provides the definitions for calculating the residuals for
 |  the different model forms
 |
 |  These models are mainly used for the collaboration with Alexey Kamenskiy
 |
 |  Author: Will Zhang
 |  Dependencies: None
 ----------------------------------------------------------------------------- */



/*******************************************************************************
 * Calculating the residual for the hyperelastic functions
 *
 * COMMENTS:
 * The main form is predetermined, the forms differs by whether the model is
 * scaled.
*******************************************************************************/

  double planar_elastin_matrix_residual(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual_general<PlanarElastinMatrix, residual_body_below, hysteresis_body_null, penalty_body_kappa, calculate_viscopart_body_null>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


  double thoracic_circ_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual_general<ThoracicCircVEScaled, residual_body_below2, hysteresis_body_null, penalty_body_kappa_2, calculate_viscopart_body_null>(
      pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


  double thoracic_default_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual_general<ThoracicDefaultVEScaled, residual_body, hysteresis_body, penalty_body_kappa_3, calculate_viscopart_body_null>(
      pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }
















  double femoral_residual_HE(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual<FemoralHE, quadratic_residual>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


  double femoral_residual_HE_scaled(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual<FemoralHEScaled, quadratic_residual>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


  double thoracic_residual_HE(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual<ThoracicIsoHE, quadratic_residual>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


  double thoracic_residual_HE_scaled(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual<ThoracicIsoHEScaled, quadratic_residual>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


/*******************************************************************************
 * Calculating the residual for the viscoelastic models
 *
 * COMMENTS:
 * The main form is predetermined, the forms differs by whether the model is
 * scaled.
*******************************************************************************/

  double femoral_residual_VE(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual<FemoralVE, quadratic_residual>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


  double femoral_residual_VE_hyst(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual_hyst<FemoralVE, quadratic_residual, penalty_body_4>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


  double femoral_residual_VE_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual_hyst_relax<FemoralVE, penalty_body_4, calculate_viscopart_body>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


  double femoral_residual_VE_scaled_hyst(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual_hyst<FemoralVEScaled, quadratic_residual, penalty_body_4>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


  double femoral_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual_hyst_relax<FemoralVEScaled, penalty_body_4, calculate_viscopart_body>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }



// The thoracic arteries


  double thoracic_defaultfung_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual_hyst_relax<ThoracicDefaultFungVEScaled, penalty_body_default, calculate_viscopart_body1>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }

  double thoracic_iso_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual_hyst_relax<ThoracicIsoVEScaled, penalty_body_1_2,calculate_viscopart_body2>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


  double thoracic_lin_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual_hyst_relax<ThoracicLinVEScaled, penalty_body_1_2,calculate_viscopart_body3>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }


  double thoracic_quad_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    return calc_residual_hyst_relax<ThoracicQuadVEScaled, penalty_body_1_2,calculate_viscopart_body3>(pars, fiber, visco, Tf, Cmax, \
      args, stress, dt, weights, deltaCG, hysteresis, alphas, index, select, n, dim, \
      nprot, skip);
  }



/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/
}

