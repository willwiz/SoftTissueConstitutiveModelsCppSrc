#pragma once

namespace caputo {
class caputo_init {
  protected:
    double bek[9];
    double e2[9];
    double dt, C0, K0, K1;
    int N = 9;

  public:
    double alpha, Tf, delta;
    double betas[9];
    double taus[9];
    double beta0;
    caputo_init();
    caputo_init(double alpha, double Tf, double delta);
    ~caputo_init();
    void set_pars(double alpha, double Tf, double delta);
    void update_dt(double dt);
    void update_dt_lin(double dt);
};

class caputo_init_scl : public caputo_init {
  public:
    double Q[9];
    double df;
    double f_prev;
    caputo_init_scl();
    caputo_init_scl(double alpha, double Tf, double delta);
    ~caputo_init_scl();
    double caputo_iter(double fn, double dt);
    double diffeq_iter(double fn, double dt);
};

template <int dim> class caputo_init_vec : public caputo_init {
  public:
    double Q[9 * dim];
    double df[dim];
    double f_prev[dim];

    caputo_init_vec();

    caputo_init_vec(double alpha, double Tf, double delta);

    ~caputo_init_vec();

    void caputo_iter(const double fn[], double dt, double v[]);
    void diffeq_iter(const double fn[], double dt, double v[]);
};

typedef caputo_init_vec<4> caputo_init_4;

double interpolate1D_newton_linear(double p1, double p2, double t);
double extrapolate1D_newton_linear(double p1, double p2, double t);

double interpolate_caputo_parameter_arr(double alpha, const double arr[100]);
double interpolate_caputo_parameter_beta(double alpha, const double arr[100]);
double interpolate_caputo_parameter_taus(double alpha, const double arr[100]);

extern const double betam[10][100];
extern const double taum[9][100];

} // namespace caputo
