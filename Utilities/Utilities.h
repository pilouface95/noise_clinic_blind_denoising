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

#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

#include <vector>

/**
 * @brief Convenient function to use the sort function provided by the vector library.
 **/
bool comparaisonFirst(
  std::pair<float, unsigned int> const& i_pair1,
  std::pair<float, unsigned int> const& i_pair2);

/**
 * @brief Get (a,b) such that n = a * b, with a and b maximum, (a < b).
 *
 * @param i_n : number to divide;
 * @param o_a : first divisor;
 * @param o_b : second divisor.
 *
 * @return none.
 **/
void getDivisors(
  const unsigned int i_n,
  unsigned int &o_a,
  unsigned int &o_b);


/**
 * @brief Return the PPCM of (a, b).
 *
 * @param i_a : first number,
 * @param i_b : second number.
 *
 * @return PPCM(a, b).
 **/
unsigned int getPpcm(
  const unsigned int i_a,
  const unsigned int i_b);


/**
 * @brief Get the average intensity range of a 3D group.
 *
 * @param i_values : 3D group on which the mean is computed;
 * @param i_offSet : offset to get the first values;
 * @param p_nb : total number of pixel to average;
 * @param p_off : offset value of minimum intensity;
 * @param p_max : maximum value of intensity.
 *
 * @return the average intensity range.
 **/
unsigned int getMean(
  std::vector<float> const& i_values,
  const unsigned int p_offSet,
  const unsigned int p_nb,
  const unsigned int p_off,
  const unsigned int p_max);


/**
 * @brief Obtain and substract the baricenter of io_group3d.
 *
 * @param io_group3d(p_rows x p_cols) : data to center;
 * @param o_baricenter(p_cols): will contain the baricenter of io_group3d;
 * @param p_rows, p_cols: size of io_group3d.
 *
 * @return none.
 **/
void centerData(
  std::vector<float> const& i_group3d,
  std::vector<double> &o_group3d,
  std::vector<float> &o_baricenter,
  const unsigned p_rows,
  const unsigned p_cols);


/**
 * @brief Affine interpolation of sparse values in a range [0, offset + max].
 *
 * @param i_xValues : abscissa;
 * @param i_yValues : ordinate;
 * @param o_values : will get the interpolation of y for all x \in [0, offset + max];
 * @param p_offset : offset to apply to abscissa;
 * @param p_max : maximum that abscissa must reach.
 *
 * @return none.
 **/
void interpValues(
  std::vector<float> const& i_xValues,
  std::vector<float> const& i_yValues,
  std::vector<float> &o_values,
  const unsigned int p_offset,
  const unsigned int p_max);


/**
 * @brief return (a, b) = ((y2 - y1) / (x2 - x1), (y2 * x1 - y1 * x2) / (x1 - x2)).
 *
 * @param i_x1 : abscissa;
 * @param i_x2 : abscissa;
 * @param i_y1 : ordinate;
 * @param i_y2 : ordinate;
 * @param o_a : will be (y2 - y1) / (x2 - x1);
 * @param o_b : will be (y2 * x1 - y1 * x2) / (x1 - x2).
 *
 * @return none.
 **/
void getAffineCoef(
  const float i_x1,
  const float i_x2,
  const float i_y1,
  const float i_y2,
  float &o_a,
  float &o_b);


/**
 * @brief Write low and high-frequency average noise curves into a specified txt file,
 *       depending on the current scale.
 *
 * @param i_means : array of estimated means for each channel and each bin;
 * @param i_stds : array of standard deviation for every frequencies, channels and bins;
 * @param p_sP : size of patches;
 * @param p_nBins : number of bins;
 * @param p_path : path to the output folder;
 * @param p_currentScale : current scale of the multi-scale algorithm;
 * @param p_nChnls : number of channels.
 *
 * @return EXIT_FAILURE in case of problems.
 **/
int writeNoiseCurves(
  std::vector<float> const& i_means,
  std::vector<float> const& i_stds,
  const unsigned int p_sP,
  const unsigned int p_nBins,
  const char* p_path,
  const unsigned int p_currentScale,
  const unsigned int p_nChnls);


/**
 * @brief Return the median of an array of any size.
 *
 * @param i_values : array of values.
 *
 * @return median(i_values).
 **/
float getMedian(
  std::vector<float> const& i_values);


/**
 * @brief Return the median of an array of any size.
 *
 * @param i_values : array of values.
 *
 * @return median(i_values).
 **/
float getMedian(
  std::vector<double> const& i_values);


/**
 * @brief Compute the weight for the distance computation according to
 *        the number of channels.
 *
 * @param o_weight : will contain the weight for each channel;
 * @param p_nbChnls : number of channels;
 * @param p_useAllChnls : if true, all channels will be used, with different weights.
 *
 * @return none.
 **/
void getWeight(
  std::vector<float> &o_weight,
  const unsigned int p_nbChnls,
  const bool p_useAllChnls);


/**
 * @brief Obtain random permutations for a tabular.
 *
 * @param o_table : will contain a random sequence of [1, N];
 * @param p_N : size of the sequence.
 *
 * @return none.
 **/
void getRandPerm(
  std::vector<unsigned int> &o_table,
  const unsigned int p_N);


#endif // UTILITIES_H_INCLUDED
