/*
 * Copyright (c) 2013, Marc Lebrun <marc.debrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file LibMatrix.cpp.
 * @brief Tools for matrix manipulation, based on ccmath functions by Daniel A. Atkinson (public code).
 *
 * @author Marc Lebrun <marc.debrun.ik@gmail.com>
 **/

#include "LibMatrix.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>

 using namespace std;

//! Invert (in place) a symmetric real matrix, V -> Inv(V).
//! The input matrix V is symmetric (V[i,j] = V[j,i]).
int inverseMatrix(
  std::vector<double> &io_mat,
  const unsigned int p_N){

  //! Initializations
  double z;
  unsigned int p, q, r, s, t, j, k;

  for (j = 0, p = 0; j < p_N ; j++, p += p_N + 1) {
    for (q = j * p_N; q < p ; q++) {
      io_mat[p] -= io_mat[q] * io_mat[q];
    }

    if (io_mat[p] <= 0.05) {
      return EXIT_FAILURE;
    }

    io_mat[p] = sqrtl(io_mat[p]);

    for (k = j + 1, q = p + p_N; k < p_N ; k++, q += p_N) {
      for (r = j * p_N, s = k * p_N, z = 0; r < p; r++, s++) {
        z += io_mat[r] * io_mat[s];
      }

      io_mat[q] -= z;
      io_mat[q] /= io_mat[p];
    }
  }

  transposeMatrix(io_mat, p_N);

  for (j = 0, p = 0; j < p_N; j++, p += p_N + 1) {
    io_mat[p] = 1 / io_mat[p];

    for (q = j, t = 0; q < p; t += p_N + 1, q += p_N) {
      for (s = q, r = t, z = 0; s < p; s += p_N, r++) {
        z -= io_mat[s] * io_mat[r];
      }

      io_mat[q] = z * io_mat[p];
    }
  }

  for (j = 0, p = 0; j < p_N; j++, p += p_N + 1) {
    for (q = j, t = p - j; q <= p; q += p_N, t++) {
      for (k = j, r = p, s = q, z = 0; k < p_N; k++, r++, s++) {
        z += io_mat[r] * io_mat[s];
      }

      io_mat[t] = io_mat[q] = z;
    }
  }

  return EXIT_SUCCESS;
}


//! Transpose a real square matrix in place io_mat -> io_mat^t.
void transposeMatrix(
  std::vector<double> &io_mat,
  const unsigned int p_N){

  for (unsigned int i = 0; i < p_N - 1; i++) {
    unsigned int p = i * (p_N + 1) + 1;
    unsigned int q = i * (p_N + 1) + p_N;

    for (unsigned int j = 0; j < p_N - 1 - i; j++, p++, q += p_N) {
      const double s = io_mat[p];
      io_mat[p] = io_mat[q];
      io_mat[q] = s;
    }
  }
}


//! Compute the covariance matrix.
void covarianceMatrix(
  std::vector<float> const& i_patches,
  std::vector<double> &o_covMat,
  const unsigned int p_nb,
  const unsigned int p_N){

  const double coefNorm = 1.f / (double) (p_nb-1);

  for (unsigned int i = 0; i < p_N; i++) {
    const float* Pi = &i_patches[i * p_nb];
    double* cMi = &o_covMat[i * p_N];
    double* cMj = &o_covMat[i];

    for (unsigned int j = 0; j < i + 1; j++) {
      const float* Pj = &i_patches[j * p_nb];
      double val = 0.f;

      for (unsigned int k = 0; k < p_nb; k++) {
        val += (double) Pi[k] * (double) Pj[k];
      }

      cMj[j * p_N] = cMi[j] = val * coefNorm;
    }
  }
}

//! Compute the covariance matrix.
void covarianceMatrix(
  std::vector<double> const& i_patches,
  std::vector<double> &o_covMat,
  const unsigned int p_nb,
  const unsigned int p_N){

  const double coefNorm = 1.f / (double) (p_nb-1);

  for (unsigned int i = 0; i < p_N; i++) {
    const double* Pi = &i_patches[i * p_nb];
    double* cMi = &o_covMat[i * p_N];
    double* cMj = &o_covMat[i];

    for (unsigned int j = 0; j < i + 1; j++) {
      const double* Pj = &i_patches[j * p_nb];
      double val = 0.f;

      for (unsigned int k = 0; k < p_nb; k++) {
        val += (double) Pi[k] * (double) Pj[k];
      }

      cMj[j * p_N] = cMi[j] = val * coefNorm;
    }
  }
}


//! Multiply two matrix A * B. (all matrices stored in row order).
void productMatrix(
  std::vector<double> &o_mat,
  std::vector<double> const& i_A,
  std::vector<double> const& i_B,
  const unsigned int p_n,
  const unsigned int p_m,
  const unsigned int p_l){

  if (o_mat.size() != p_n * p_l) {
    o_mat.resize(p_n * p_l);
  }
  std::vector<double> q0(p_m);
  double* const pq0 = &q0[0];

  for (unsigned int i = 0; i < p_l; i++) {
    const double* pB = &i_B[i];

    for (unsigned int k = 0; k < p_m; k++, pB += p_l) {
      pq0[k] = *pB;
    }

    double* pO = &o_mat[i];
    const double* pA = &i_A[0];
    for (unsigned int j = 0; j < p_n; j++, pO += p_l, pA += p_m) {
      double z = 0;

      for (unsigned int k = 0; k < p_m; k++) {
        z += pA[k] * pq0[k];
      }

      *pO = z;
    }
  }
}



















