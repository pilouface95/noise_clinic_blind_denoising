/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
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
 * @file Utilities.cpp
 * @brief Set of useful functions.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "Utilities.h"

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>

using namespace std;


//! Convenient function to use the sort function provided by the vector library.
bool comparaisonFirst(
  std::pair<float, unsigned int> const& i_pair1,
  std::pair<float, unsigned int> const& i_pair2){

  return i_pair1.first < i_pair2.first;
}


//! Get (a,b) such that n = a * b, with a and b maximum, a < b
void getDivisors(
  const unsigned int i_n,
  unsigned int &o_a,
  unsigned int &o_b){

  o_a = (unsigned int) sqrtf(float(i_n));
  while (i_n % o_a != 0) {
    o_a--;
  }
  o_b = i_n / o_a;

  if (o_b < o_a) {
    int c = o_b;
    o_b = o_a;
    o_a = c;
  }
}


//! Return the PPCM of (a, b).
unsigned int getPpcm(
  const unsigned int i_a,
  const unsigned int i_b){

  if (i_a == 0 || i_b == 0) {
    return 0;
  }

  int a = i_a;
  int b = i_b;
  while (a != b) {
    if (a < b) {
      a += i_a;
    }
    else {
      b += i_b;
    }
  }

  return a;
}


//! Get the average intensity range of a 3D group.
unsigned int getMean(
  std::vector<float> const& i_values,
  const unsigned int p_offSet,
  const unsigned int p_nb,
  const unsigned int p_off,
  const unsigned int p_max){

  //! Compute the mean
  float mean = 0.f;
  const float* iV = &i_values[p_offSet];
  for (unsigned int k = 0; k < p_nb; k++) {
    mean += iV[k];
  }
  mean /= float(p_nb);
  mean += p_off;

  //! Return the intensity
  return (mean < 0 ? 0 : (mean > (float) (p_max + p_off) ? (p_max + p_off) : (unsigned int) mean));
}


//! Obtain and substract the baricenter of io_group3d.
void centerData(
  std::vector<float> const& i_group3d,
  std::vector<double> &o_group3d,
  std::vector<float> &o_baricenter,
  const unsigned p_rows,
  const unsigned p_cols){

  const float inv = 1.f / (float) p_rows;
  for (unsigned j = 0; j < p_cols; j++) {
    const float* iG = &i_group3d[j * p_rows];
    float sum = 0.f;
    float* iB = &o_baricenter[j];
    double* oG = &o_group3d[j * p_rows];

    for (unsigned i = 0; i < p_rows; i++) {
      sum += iG[i];
    }
    *iB = sum * inv;

    for (unsigned i = 0; i < p_rows; i++) {
      oG[i] = iG[i] - *iB;
    }
  }
}


//! Affine interpolation of sparse values in a range [0, offset + max].
void interpValues(
  std::vector<float> const& i_xValues,
  std::vector<float> const& i_yValues,
  std::vector<float> &o_values,
  const unsigned int p_offset,
  const unsigned int p_max){

  //! Initializations
  const unsigned int N     = i_xValues.size();
  const unsigned int range = p_offset + p_max + 1;
  const float off          = (float) p_offset;
  o_values.resize(range);

  //! Weird case
  if (N == 1) {
    const float value = i_yValues[0];
    for (unsigned int r = 0; r < range; r++) {
      o_values[r] = value;
    }
    return;
  }

  //! Initialization
  unsigned int i = 0;
  float a, b;
  float* oV = &o_values[0];

  //! First bin
  getAffineCoef(i_xValues[0] + off, i_xValues[1] + off, i_yValues[0], i_yValues[1], a, b);
  const float xFirst = i_xValues[0];
  while (i <= xFirst) {
    oV[i] = a * (float) i + b;
    i++;
  }

  //! medium bins
  for (unsigned int n = 1; n < N - 1; n++) {
    const float xMax = i_xValues[n];

    if (xMax < range - 1) {
      getAffineCoef(i_xValues[n - 1] + off, i_xValues[n] + off, i_yValues[n - 1], i_yValues[n], a, b);

      while (i <= xMax) {
        oV[i] = a * float(i) + b;
        i++;
      }
    }
  }

  //! Last bin
  getAffineCoef(i_xValues[N - 2] + off, i_xValues[N - 1] + off, i_yValues[N - 2], i_yValues[N - 1], a, b);
  while (i < range) {
    oV[i] = a * float(i) + b;
    i++;
  }
}


//! Return (a, b) = ((y2 - y1) / (x2 - x1), (y2 * x1 - y1 * x2) / (x1 - x2)).
void getAffineCoef(
  const float i_x1,
  const float i_x2,
  const float i_y1,
  const float i_y2,
  float &o_a,
  float &o_b){

  o_a = (i_y2 - i_y1) / (i_x2 - i_x1);
  o_b = (i_y2 * i_x1 - i_y1 * i_x2) / (i_x1 - i_x2);
}


//! Write low and high-frequency average noise curves into a specified txt file,
//! depending on the current scale.
int writeNoiseCurves(
  std::vector<float> const& i_means,
  std::vector<float> const& i_stds,
  const unsigned int p_sP,
  const unsigned int p_nBins,
  const char* p_path,
  const unsigned int p_currentScale,
  const unsigned int p_nChnls){

  //! Initializations
  ostringstream oss_cs;
  oss_cs << p_currentScale;

  string fileName_H = p_path;
  fileName_H  += "_noiseCurves_H_"  + oss_cs.str() + ".txt";
  //
  string fileName_L = p_path;
  fileName_L  += "_noiseCurves_L_"  + oss_cs.str() + ".txt";

  const unsigned int sH  = (unsigned int) floor((float) p_sP * 0.5f);
  const unsigned int sP2 = p_sP * p_sP;

  //! Create files
  ofstream file_H((char*) fileName_H.c_str(), ios::out | ios::trunc);
  ofstream file_L((char*) fileName_L.c_str(), ios::out | ios::trunc);

  //! Check created OK
  if (!file_L) {
    fprintf(stderr, "Couldn't create the file: %s\n", fileName_L.c_str());
    return EXIT_FAILURE;
  }

  if (!file_H) {
    fprintf(stderr, "Couldn't create the file: %s\n", fileName_H.c_str());
    return EXIT_FAILURE;
  }

  //! Write data
  for (unsigned int b = 0; b < p_nBins; b++) {

    //! Means
    for (unsigned int c = 0; c < p_nChnls; c++) {
      float mean = i_means[c * p_nBins + b];
      file_L << mean << " ";
      file_H << mean << " ";
    }

    //! Low Frequency
    for (unsigned int c = 0; c < p_nChnls; c++) {
      float lf = 0.f;
      const float* iS = &i_stds[c * p_nBins * sP2 + b * sP2];
      for (unsigned int p = 0; p < sH; p++) {
        for (unsigned int q = 0; q < sH; q++) {
          lf += iS[p * p_sP + q];
        }
      }
      file_L << lf / (float) (sH * sH) << " ";
    }

    //! High Frequency
    for (unsigned int c = 0; c < p_nChnls; c++) {
      float hf = 0.f;
      const float* iS = &i_stds[c * p_nBins * sP2 + b * sP2];
      for (unsigned int p = 0; p < p_sP; p++) {
        for (unsigned int q = 0; q < p_sP; q++) {
          if (p >= sH || q >= sH) {
            hf += iS[p * p_sP + q];
          }
        }
      }
      file_H << hf / (float) (sP2 - sH * sH) << " ";
    }
    file_L << "\n";
    file_H << "\n";
  }
  file_L.close();
  file_H.close();

  return EXIT_SUCCESS;
}


//! Return the median of an array of any size
float getMedian(
  std::vector<float> const& i_values){

  //! Get size
  const unsigned int N = i_values.size();
  vector<float> valTmp = i_values;

  //! Weird cases
  if (N <= 0) {
      cout << "Error: not enough samples to compute median" << endl;
      return 0;
  }
  if (N == 1) {
    return i_values[0];
  }

  //! Sort the values
  const unsigned int ind = (N % 2 ? (N + 1) / 2 : N / 2 + 1);
  partial_sort(valTmp.begin(), valTmp.begin() + ind, valTmp.end());

  //! Return the result
  if (N % 2 == 0) {
    return (valTmp[N / 2 - 1] + valTmp[N / 2]) * 0.5f;
  }
  else {
    return valTmp[(N - 1) / 2];
  }
}


//! Return the median of an array of any size (double)
float getMedian(
  std::vector<double> const& i_values){

  //! Get size
  const unsigned int N = i_values.size();
  vector<double> valTmp = i_values;

  //! Weird cases
  if (N <= 0) {
      cout << "Error: not enough samples to compute median: " << N << endl;
      return 0;
  }
  if (N == 1) {
    return i_values[0];
  }

  //! Sort the values
  const unsigned int ind = (N % 2 ? (N + 1) / 2 : N / 2 + 1);
  partial_sort(valTmp.begin(), valTmp.begin() + ind, valTmp.end());

  //! Return the result
  if (N % 2 == 0) {
    return (valTmp[N / 2 - 1] + valTmp[N / 2]) * 0.5f;
  }
  else {
    return valTmp[(N - 1) / 2];
  }
}


//! Compute the weight for the distance computation according to
//! the number of channels.
void getWeight(
  std::vector<float> &o_weight,
  const unsigned int p_nbChnls,
  const bool p_useAllChnls){

  //! Initialization
  o_weight.resize(p_nbChnls);

  //! Compute weight depending on the number of channels
  switch (p_nbChnls) {
    case 0:
      o_weight[0] = 1.f;
      break;
    case 3:
      if (p_useAllChnls) {
        o_weight[0] = sqrtf(0.5f);
        o_weight[1] = sqrtf(0.25f);
        o_weight[2] = sqrtf(0.25f);
      }
      else {
        o_weight[0] = 1.f;
        o_weight[1] = 0.f;
        o_weight[2] = 0.f;
      }
      break;
    case 4:
      if (p_useAllChnls) {
        o_weight[0] = sqrtf(1.f / 3.f);
        o_weight[1] = sqrtf(1.f / 6.f);
        o_weight[2] = sqrtf(1.f / 6.f);
        o_weight[3] = sqrtf(1.f / 3.f);
      }
      else {
        o_weight[0] = 1.f;
        o_weight[1] = 0.f;
        o_weight[2] = 0.f;
        o_weight[3] = 0.f;
      }
      break;
    default:
      break;
  }
}


//! Obtain random permutations for a tabular.
void getRandPerm(
  std::vector<unsigned int> &o_table,
  const unsigned int p_N){

  //! Initializations
  o_table.resize(p_N);
  vector<unsigned int> tmp(p_N + 1, 0);
  tmp[1] = 1;
  srand((time(NULL)));

  for (unsigned int i = 2; i < p_N + 1; i++) {
    unsigned int j = rand() % i + 1;
    tmp[i] = tmp[j];
    tmp[j] = i;
  }

  for (unsigned int n = 0; n < p_N; n++) {
    o_table[n] = tmp[n + 1] - 1;
  }
}




















