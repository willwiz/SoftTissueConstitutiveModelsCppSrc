#include "onesided.hpp"
#include <cmath>
namespace res {

void regress_linear(
    const double a_data[], const double b_data[], const double dt[], int start, int end,
    double res[]
) {
    double time = 0.0;
    double t2 = 0.0;
    double ty = 0.0;
    double slope;
    for (int i = start; i < end; i++) {
        time = time + dt[i];
        t2 = t2 + time * time;
        ty = ty + time * (a_data[i] - b_data[i]);
    }
    slope = ty / t2;
    time = 0.0;
    for (int i = start; i < end; i++) {
        time = time + dt[i];
        res[i] = a_data[i] - slope * time;
    }
}

void frechet_residual_onesided(
    const double a_data[], const double a_dt[], int a_start, int a_end, const double b_data[],
    const double b_dt[], int b_start, int b_end, double res[]
) {
    double t_a = -a_dt[a_start];
    double t_b = 0.0;
    double old_res;
    double new_res;
    int j = b_start;
    for (int i = a_start; i < a_end; i++) {
        t_a = t_a + a_dt[i];
        old_res = (t_a - t_b) * (t_a - t_b) + (a_data[i] - b_data[j]) * (a_data[i] - b_data[j]);
        for (int k = j; k < b_end; k++) {
            new_res = (t_a - t_b - b_dt[k]) * (t_a - t_b - b_dt[k]) +
                      (a_data[i] - b_data[k]) * (a_data[i] - b_data[k]);
            if (old_res < new_res) break;
            j = k;
            t_b += b_dt[k];
            old_res = new_res;
        }
        res[i] = a_data[i] - b_data[j];
    }
}

void window_residual_onesided(
    const double a_data[], const double b_data[], const int start, const int end, double res[]
) {
    // double res_i;
    for (int i = start; i < end; i++) {
        // res_i[i] = a_data[i] - b_data[i];
        // if (std::abs(a_data[i] - b_data[i + 1]) < std::abs(res_i)) {
        //     res_i = a_data[i] - b_data[i + 1];
        // }
        // res[i] = res_i;
        res[i] = a_data[i] - b_data[i];
    }
}
} // namespace res
