#pragma once

#define _USE_MATH_DEFINES

namespace ensemble {

const double p_fiber = 0.001;
const double p_alpha = 1.0;
const double p_elastin = 1e-3;
const double p_collagen = 1e-1;
const double b_visco = 0.01;
const double b_modulus = 0.01;

const double w_hyst = 0.01;
const double w_visco = 10.0;

const double M_kip = 0.155;
const double M_kop = 0.424285714285714;

const double M_ideal_alpha = M_PI_4;

} // namespace ensemble
