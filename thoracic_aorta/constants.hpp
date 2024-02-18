#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

namespace thoracic {

const double M_p_fiber = 0.001;
const double M_p_alpha = 1.0;
const double M_p_elastin = 1e-3;
const double M_p_collagen = 1e-1;
const double M_b_visco = 0.01;
const double M_b_modulus = 0.01;

const double M_w_hyst = 0.1;
const double M_w_visco = 10.0;

const double M_kip = 0.155;
const double M_kop = 0.424285714285714;

const double M_ideal_alpha = M_PI_4;

} // namespace thoracic
