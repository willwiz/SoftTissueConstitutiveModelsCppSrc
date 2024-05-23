#pragma once

#define _USE_MATH_DEFINES

namespace femoral {

extern const double M_p_fiber = 0.001;
extern const double M_p_alpha = 0.1;
extern const double M_p_elastin = 1e-3;
extern const double M_p_collagen = 1e-1;
extern const double M_b_visco = 0.01;
extern const double M_b_modulus = 0.01;

extern const double M_w_hyst = 0.001;
extern const double M_w_visco = 10.0;

extern const double M_kip = 0.155;
extern const double M_kop = 0.424285714285714;

extern const double M_ideal_alpha = M_PI_4;

} // namespace femoral
