#pragma once

namespace residuals {

  // Exported functions
  double planar_elastin_matrix_residual(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double femoral_residual_HE(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double femoral_residual_HE_scaled(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double femoral_residual_VE(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double femoral_residual_VE_hyst(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double femoral_residual_VE_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double femoral_residual_VE_scaled_hyst(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double femoral_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double thoracic_default_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double thoracic_defaultfung_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double thoracic_iso_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double thoracic_lin_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double thoracic_quad_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double thoracic_circ_residual_VE_scaled_hyst_relax(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);
}
