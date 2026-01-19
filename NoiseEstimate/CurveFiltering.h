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
#ifndef CURVEFILTERING_H_INCLUDED
#define CURVEFILTERING_H_INCLUDED

#include <vector>

/**
 * @brief Filter/smooth noise curves.
 *
 * @param i_meanBins : contains for each channel and each bin, the mean;
 * @param io_stdsBins : contains for each channel and each bin the frequencies of the noise.
 *       Will be modified.
 * @param p_sP : size of patches;
 * @param p_nBins : number of bins;
 * @param p_nChannels : number of channels;
 * @param p_rangeOffset : offset for intensity range;
 * @param p_rangeMax : maximum of intensity range. i \isin [-offset, max] -> [0, max + offset].
 *
 * @return none.
 **/
void filterNoiseCurves(
  std::vector<float> const& i_meanBins,
  std::vector<float> &io_stdsBins,
  const unsigned int p_sP,
  const unsigned int p_nBins,
  const unsigned int p_nChannels,
  const unsigned int p_rangeOffset,
  const unsigned int p_rangeMax);


/**
 * @brief Smooth control point values according to the whole interpolated curve.
 *
 * @param i_X : contains control point abscissa;
 * @param o_Y : contains control point ordinate which will be smoothed;
 * @param i_Yall : contains all ordinate;
 * @param p_size : size of the smooth filter;
 * @param p_offset : offset to apply to abscissa in array.
 *
 * @return none.
 **/
void smoothValues(
  std::vector<float> const& i_X,
  std::vector<float> &o_Y,
  std::vector<float> const& i_Yall,
  const unsigned int p_size,
  const unsigned int p_offset);


/**
 * @brief Filter coefficient of a matrix by taking the median of the closest 4-neighbours.
 *
 * @param io_mat : matrix to filter;
 * @param p_sP : size of the matrix (p_sP \times p_sP).
 *
 * @return none.
 **/
void applyMatrixMedianFilter(
  std::vector<float> &io_mat,
  const unsigned int p_sP);


/**
 * @brief Return the median of 3 values.
 *
 * @param i_a : first value,
 * @param i_b : second value,
 * @param i_c : third value.
 *
 * @return median(i_a, i_b, i_c).
 **/
float getMedian3(
  const float i_a,
  const float i_b,
  const float i_c);


/**
 * @brief Return the median of 4 values.
 *
 * @param i_a : first value,
 * @param i_b : second value,
 * @param i_c : third value,
 * @param i_d : fourth value.
 *
 * @return median(i_a, i_b, i_c, i_d).
 **/
float getMedian4(
  const float i_a,
  const float i_b,
  const float i_c,
  const float i_d);


/**
 * @brief Return the median of 4 values.
 *
 * @param i_a : first value,
 * @param i_b : second value,
 * @param i_c : third value,
 * @param i_d : fourth value,
 * @param i_e : fifth value.
 *
 * @return median(i_a, i_b, i_c, i_d, i_e).
 **/
float getMedian5(
  const float i_a,
  const float i_b,
  const float i_c,
  const float i_d,
  const float i_e);


#endif // CURVEFILTERING_H_INCLUDED
