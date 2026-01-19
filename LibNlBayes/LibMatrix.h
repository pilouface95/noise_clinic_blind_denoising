/*
 * Copyright (c) 2014, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LIB_MATRIX_H_INCLUDED
#define LIB_MATRIX_H_INCLUDED

#include <vector>

/**
 * @brief Invert (in place) a symmetric real matrix, V -> Inv(V).
 *       The input matrix V is symmetric (V[i,j] = V[j,i]).
 *
 * @param io_mat = array containing a symmetric input matrix. This is converted
 *       to the inverse matrix;
 * @param p_N = dimension of the system (dim(v)=n*n).
 *
 * @return EXIT_FAILURE if input matrix not positive definite, otherwise EXIT_SUCCESS.
 **/
int inverseMatrix(
  std::vector<double> &io_mat,
  const unsigned int p_N);


/**
 * @brief Transpose a real square matrix in place io_mat -> io_mat~.
 *
 * @param io_mat : pointer to an array of p_N by p_N input matrix io_mat. This is
 *       overloaded by the transpose of io_mat;
 * @param p_N : dimension (dim(io_mat) = p_N * p_N).
 *
 * @return none.
 **/
void transposeMatrix(
  std::vector<double> &io_mat,
  const unsigned int p_N);


/**
 * @brief Compute the covariance matrix.
 *
 * @param i_patches: set of patches of size (nb x N);
 * @param o_covMat: will contain patches' * patches;
 * @param p_N : size of a patch;
 * @param p_nb: number of similar patches in the set of patches.
 *
 * @return none.
 **/
void covarianceMatrix(
  std::vector<float> const& i_patches,
  std::vector<double> &o_covMat,
  const unsigned int p_nb,
  const unsigned int p_N);


/**
 * @brief Compute the covariance matrix.
 *
 * @param i_patches: set of patches of size (nb x N);
 * @param o_covMat: will contain patches' * patches;
 * @param p_N : size of a patch;
 * @param p_nb: number of similar patches in the set of patches.
 *
 * @return none.
 **/
void covarianceMatrix(
  std::vector<double> const& i_patches,
  std::vector<double> &o_covMat,
  const unsigned int p_nb,
  const unsigned int p_N);


/**
 * @brief Multiply two matrix A * B. (all matrices stored in row order).
 *
 * @param o_mat : array containing n by l product matrix at exit;
 * @param i_A : input array containing n by m matrix;
 * @param i_B : input array containing m by l matrix;
 * @param p_n : nb lines of A
 * @param p_l : nb columns of B
 * @param p_m : nb columns of A (equals nb lines of B).
 *
 * @return  none.
 *
 * @note The right and output term may point to the same vector.
 */
void productMatrix(
  std::vector<double> &o_mat,
  std::vector<double> const& i_A,
  std::vector<double> const& i_B,
  const unsigned int p_n,
  const unsigned int p_m,
  const unsigned int p_l);


#endif // LIB_MATRIX_H_INCLUDED
