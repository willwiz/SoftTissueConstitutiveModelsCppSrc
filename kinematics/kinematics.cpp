#define _USE_MATH_DEFINES

#include "kinematics.hpp"
#include <cmath>
#include <iostream>
#include <stdio.h>

namespace kinematics {

/*----------------------------------------------------------------------
 |  This file provides the structure for storing the info for deformation
 |  gradient. It is helpful for reducing recomputation.
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

template <int dim> kinematics<dim>::kinematics() : C{}, Cinv{} {
}

template <int dim> kinematics<dim>::~kinematics() {
}

deformation1D::deformation1D() : kinematics<1>() {
}
deformation1D::deformation1D(const double args) : kinematics<1>() {
    this->precompute(args);
}

void deformation1D::precompute(const double args) {
    this->det = args;
    this->I_n = 1 / std::sqrt(args);
    this->Cinv[0] = 1 / args;
    this->C[0] = args;
    this->I_1 = args + I_n + I_n;
    this->I_1m3 = I_1 - 3.0;
}

deformation2D::deformation2D() {
}
deformation2D::deformation2D(const double args[4]) : kinematics<4>() {
    this->precompute(args);
}

void deformation2D::precompute(const double args[4]) {

    this->det = args[0] * args[3] - args[1] * args[2];
    this->I_n = 1 / det;
    this->Cinv[0] = args[3] * I_n;
    this->Cinv[1] = -args[1] * I_n;
    this->Cinv[2] = -args[2] * I_n;
    this->Cinv[3] = args[0] * I_n;

    for (int i = 0; i < 4; i++) {
        this->C[i] = args[i];
    }

    this->I_1 = args[0] + args[3] + I_n;
    this->I_1m3 = I_1 - 3.0;
}

deformation_ensemble2D::deformation_ensemble2D() {
}
deformation_ensemble2D::deformation_ensemble2D(double eb_strain) : kinematics<4>() {
    this->precompute(eb_strain);
}

void deformation_ensemble2D::precompute(double eb_strain) {

    this->det = eb_strain * eb_strain;
    this->I_n = 1 / det;
    this->C[0] = eb_strain;
    this->C[1] = 0;
    this->C[2] = 0;
    this->C[3] = eb_strain;
    this->Cinv[0] = 1.0 / eb_strain;
    this->Cinv[1] = 0;
    this->Cinv[2] = 0;
    this->Cinv[3] = 1.0 / eb_strain;

    // for (int i = 0; i < 4; i++)
    // {
    //   this->C33Cinv[i] = I_n*Cinv[i];
    // }

    this->I_1 = eb_strain + eb_strain + I_n;
    this->I_1m3 = I_1 - 3.0;
}

deformation3D::deformation3D() {
}
deformation3D::deformation3D(const double args[9]) : kinematics<9>() {
    this->precompute(args);
}

double determinant_3D_tensor(const double tensor[9]) {
    return -tensor[2] * tensor[2] * tensor[4] + 2 * tensor[1] * tensor[2] * tensor[5] -
           tensor[0] * tensor[5] * tensor[5] - tensor[1] * tensor[1] * tensor[8] +
           tensor[0] * tensor[4] * tensor[8];
}

void deformation3D::precompute(const double args[9]) {

    det = determinant_3D_tensor(args);
    I_n = 1 / det;
    Cinv[0] = I_n * (-args[5] * args[5] + args[4] * args[8]);
    Cinv[1] = I_n * (args[2] * args[5] - args[1] * args[8]);
    Cinv[2] = I_n * (-args[2] * args[4] + args[1] * args[5]);
    Cinv[3] = Cinv[1];
    Cinv[4] = I_n * (-args[2] * args[2] + args[0] * args[8]);
    Cinv[5] = I_n * (args[1] * args[2] - args[0] * args[5]);
    Cinv[6] = Cinv[2];
    Cinv[7] = Cinv[5];
    Cinv[8] = I_n * (-args[1] * args[1] + args[0] * args[4]);

    for (int i = 0; i < 9; i++) {
        C[i] = args[i];
    }

    I_1 = args[0] + args[4] + args[8];
    I_1m3 = I_1 - 3.0;
}
} // namespace kinematics