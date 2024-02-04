#pragma once

#define _USE_MATH_DEFINES

namespace ctv {

extern const double p_fiber = 0.001;
extern const double p_alpha = 1.0;
extern const double p_elastin = 1e-3;
extern const double p_collagen = 1e-1;
extern const double b_visco = 0.01;
extern const double b_modulus = 0.01;

extern const double w_hyst = 0.01;
extern const double w_visco = 10.0;

extern const int prob_dim = 4;

extern const double M_kip = 0.155;
extern const double M_kop = 0.424285714285714;

extern const double ideal_alpha = M_PI_4;

} // namespace ctv
