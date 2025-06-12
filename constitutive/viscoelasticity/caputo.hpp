#pragma once

namespace caputo {

template <int n_prony> class caputo_init {
  protected:
    double bek[n_prony];
    double e2[n_prony];
    double dt, C0, K0, K1;

  public:
    double alpha, Tf, delta;
    double betas[n_prony];
    double taus[n_prony];
    double beta0;
    caputo_init() : dt{}, betas{}, taus{} {};
    caputo_init(double alpha, double Tf, double delta);
    ~caputo_init() {};
    void set_pars(double alpha, double Tf, double delta);
    void update_dt(double dt);
    void update_dt_lin(double dt);
};

template <int n_prony> class caputo_init_scl : public caputo_init<n_prony> {
  public:
    double Q[n_prony];
    double df;
    double f_prev;
    caputo_init_scl() : caputo_init<n_prony>(), Q{}, f_prev{} {};
    caputo_init_scl(double alpha, double Tf, double delta)
        : caputo_init<n_prony>(alpha, Tf, delta), Q{}, f_prev() {};
    ~caputo_init_scl() {};
    double caputo_iter(double fn, double dt);
    double diffeq_iter(double fn, double dt);
};

template <int dim, int n_prony> class caputo_init_vec : public caputo_init<n_prony> {
  public:
    double Q[n_prony * dim];
    double df[dim];
    double f_prev[dim];

    caputo_init_vec() : caputo_init<n_prony>(), Q{}, f_prev{} {};

    caputo_init_vec(double alpha, double Tf, double delta)
        : caputo_init<n_prony>(alpha, Tf, delta), Q{}, f_prev{} {};

    ~caputo_init_vec() {};

    void caputo_iter(const double fn[], const double dt, double v[]);
    void diffeq_iter(const double fn[], const double dt, double v[]);
};

typedef caputo_init_vec<4, 9> caputo_init_4;

double interpolate1D_newton_linear(double p1, double p2, double t);
double extrapolate1D_newton_linear(double p1, double p2, double t);

double interpolate_caputo_parameter_arr(double alpha, const double arr[100]);
double interpolate_caputo_parameter_beta(double alpha, const double arr[100]);
double interpolate_caputo_parameter_taus(double alpha, const double arr[100]);

} // namespace caputo
