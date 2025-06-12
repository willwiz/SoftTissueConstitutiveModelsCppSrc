#pragma once
namespace res {
void regress_linear(
    const double a_data[], const double b_data[], const double dt[], int start, int end,
    double res[]
);
void frechet_residual_onesided(
    const double a_data[], const double a_dt[], int a_start, int a_end, const double b_data[],
    const double b_dt[], int b_start, int b_end, double res[]
);

void window_residual_onesided(
    const double a_data[], const double b_data[], const int start, const int end, double res[]
);
} // namespace res
