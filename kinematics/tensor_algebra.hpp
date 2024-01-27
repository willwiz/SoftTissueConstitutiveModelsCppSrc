#pragma once

namespace constitutive_models {

extern const double id2d[4];
extern const double id3d[9];

inline double ddot2D(const double a[4], const double b[4]);
inline double ddot(const double a[], const double b[], const int dim);

void addto(const double a[], double b[], int dim);

} // namespace constitutive_models
