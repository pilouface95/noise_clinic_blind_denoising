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
 * @file LibImages.cpp
 * @brief Usefull functions on images
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "CurveFiltering.h"
#include "../Utilities/Utilities.h"

 using namespace std;


//! Filter/smooth noise curves.
void filterNoiseCurves(
  std::vector<float> const& i_meanBins,
  std::vector<float> &io_stdsBins,
  const unsigned int p_sP,
  const unsigned int p_nBins,
  const unsigned int p_nChannels,
  const unsigned int p_rangeOffset,
  const unsigned int p_rangeMax){

  //! Initializations
  const unsigned int sP2        = p_sP * p_sP;
  const unsigned int nbFilter   = 5;
  const unsigned int sizeFilter = 10;

  //! Channel independant
  for (unsigned int c = 0; c < p_nChannels; c++) {

    //! Get the means
    vector<float> meanBins(p_nBins);
    const float* iM = &i_meanBins[c * p_nBins];
    for (unsigned int n = 0; n < p_nBins; n++) {
      meanBins[n] = iM[n];
    }

    //! Loop of smooth filtering
    for (unsigned int nbF = 0; nbF < nbFilter; nbF++) {

      //! Get coefficients along the bins
      for (unsigned int k = 0; k < sP2; k++) {

        //! Get frequency values
        vector<float> values(p_nBins);
        float* iS = &io_stdsBins[c * p_nBins * sP2 + k];
        for (unsigned int n = 0; n < p_nBins; n++) {
          values[n] = iS[n * sP2] * iS[n * sP2];
        }

        //! Interpolation of the values
        vector<float> allValues;
        interpValues(meanBins, values, allValues, p_rangeOffset, p_rangeMax);

        //! Values can't be negative
        for (unsigned int k = 0; k < allValues.size(); k++) {
          allValues[k] = std::max(0.f, allValues[k]);
        }

        //! Filter the values
        smoothValues(meanBins, values, allValues, sizeFilter, p_rangeOffset);

        //! Save values
        for (unsigned int n = 0; n < p_nBins; n++) {
          iS[n * sP2] = sqrtf(values[n]);
        }
      }

      //! Filter matrix by matrix
      for (unsigned int n = 0; n < p_nBins; n++) {

        //! Get values
        vector<float> mat(sP2);
        float* iS = &io_stdsBins[c * p_nBins * sP2 + n * sP2];
        for (unsigned int k = 0; k < sP2; k++) {
          mat[k] = iS[k] * iS[k];
        }

        //! Filter by median the matrix
        applyMatrixMedianFilter(mat, p_sP);

        //! Save values
        for (unsigned int k = 0; k < sP2; k++) {
          iS[k] = sqrtf(mat[k]);
        }
      }
    }
  }
}


//! Smooth control point values according to the whole interpolated curve.
void smoothValues(
  std::vector<float> const& i_X,
  std::vector<float> &o_Y,
  std::vector<float> const& i_Yall,
  const unsigned int p_size,
  const unsigned int p_offset){

  //! Initializations
  const unsigned int nb = i_X.size();

  //! Loop over control points
  for (unsigned int n = 0; n < nb; n++) {

    //! Get the min/max index
    const unsigned int minInd = std::max((int) (i_X[n] - (int) p_size + (int) p_offset), 0);
    const unsigned int maxInd = std::min((int) (i_X[n] + p_size + p_offset), (int)i_Yall.size() - 1);

    //! Get the value
    float mean = 0.f;
    const unsigned int nbInd = maxInd - minInd + 1;
    const float* iY = &i_Yall[minInd];
    for (unsigned int k = 0; k < nbInd; k++) {
      mean += iY[k];
    }
    o_Y[n] = mean / float(nbInd);
  }
}


//! Filter coefficient of a matrix by taking the median of the closest 4-neighbours.
void applyMatrixMedianFilter(
  std::vector<float> &io_mat,
  const unsigned int p_sP){

  //! Top left corner (i == 0, j == 0)
  io_mat[0] = getMedian3(io_mat[0],
                         io_mat[1],
                         io_mat[p_sP]);

  //! Top right corner (i == 0, j == p_sP - 1)
  io_mat[p_sP - 1] = getMedian3(io_mat[p_sP - 1],
                                io_mat[p_sP - 2],
                                io_mat[p_sP + p_sP - 1]);

  //! Bottom left corner (i == p_sP - 1, j == 0)
  io_mat[(p_sP - 1) * p_sP] = getMedian3(io_mat[(p_sP - 1) * p_sP],
                                         io_mat[(p_sP - 2) * p_sP],
                                         io_mat[(p_sP - 1) * p_sP + 1]);

  //! Bottom right corner (i == p_sP - 1, j == p_sP - 1)
  io_mat[(p_sP - 1) * p_sP + p_sP - 1] = getMedian3(io_mat[(p_sP - 1) * p_sP + p_sP - 1],
                                                    io_mat[(p_sP - 2) * p_sP + p_sP - 1],
                                                    io_mat[(p_sP - 1) * p_sP + p_sP - 2]);

  for (unsigned int j = 1; j < p_sP - 1; j++) {

    //! Top line except corners
    io_mat[j] = getMedian4(io_mat[j],
                           io_mat[j - 1],
                           io_mat[j + 1],
                           io_mat[j + p_sP]);

    //! Bottom line except corners
    io_mat[(p_sP - 1) * p_sP + j] = getMedian4(io_mat[(p_sP - 1) * p_sP + j],
                                               io_mat[(p_sP - 1) * p_sP + j - 1],
                                               io_mat[(p_sP - 1) * p_sP + j + 1],
                                               io_mat[(p_sP - 2) * p_sP + j]);

    //! Left column except corners
    io_mat[j * p_sP] = getMedian4(io_mat[j * p_sP],
                                  io_mat[(j - 1) * p_sP],
                                  io_mat[(j + 1) * p_sP],
                                  io_mat[j * p_sP + 1]);

    //! Right column except corners
    io_mat[j * p_sP + p_sP - 1] = getMedian4(io_mat[j * p_sP + p_sP - 1],
                                             io_mat[(j - 1) * p_sP + p_sP - 1],
                                             io_mat[(j + 1) * p_sP + p_sP - 1],
                                             io_mat[j * p_sP + p_sP - 2]);
  }

  //! Center of the matrix
  for (unsigned int i = 1; i < p_sP - 1; i++) {
    for (unsigned int j = 1; j < p_sP - 1; j++) {
      io_mat[i * p_sP + j] = getMedian5(io_mat[i * p_sP + j],
                                        io_mat[i * p_sP + j + 1],
                                        io_mat[i * p_sP + j - 1],
                                        io_mat[(i - 1) * p_sP + j],
                                        io_mat[(i + 1) * p_sP + j]);
    }
  }
}


//! Return the median of 3 values.
float getMedian3(
  const float i_a,
  const float i_b,
  const float i_c){

  if (i_a < i_b) {
    if (i_a < i_c) {
      return std::min(i_b, i_c);
    }
    else {
      return i_a;
    }
  }
  else {
    if (i_b < i_c) {
      return std::min(i_a, i_c);
    }
    else {
      return i_b;
    }
  }
}


//! Return the median of 4 values.
float getMedian4(
  const float i_a,
  const float i_b,
  const float i_c,
  const float i_d){

  if (i_a < i_b) {
    if (i_a < i_c) {
      if (i_a < i_d) {
        if (i_b < i_c) {
          if (i_b < i_d) {
            return (i_b + std::min(i_c, i_d)) / 2.f;
          }
          else {
            return (i_b + i_d) / 2.f;
          }
        }
        else {
          if (i_c < i_d) {
            return (i_c + std::min(i_b, i_d)) / 2.f;
          }
          else {
            return (i_c + i_d) / 2.f;
          }
        }
      }
      else {
        return (i_a + std::min(i_b, i_c)) / 2.f;
      }
    }
    else {
      if (i_a < i_d) {
        return (i_a + std::min(i_b, i_d)) / 2.f;
      }
      else {
        return (i_a + std::max(i_c, i_d)) / 2.f;
      }
    }
  }
  else {
    if (i_b < i_c) {
      if (i_b < i_d) {
        if (i_a < i_c) {
          if (i_a < i_d) {
            return (i_a + std::min(i_c, i_d)) / 2.f;
          }
          else {
            return (i_a + i_d) / 2.f;
          }
        }
        else {
          if (i_c < i_d) {
            return (i_c + std::min(i_a, i_d)) / 2.f;
          }
          else {
            return (i_c + i_d) / 2.f;
          }
        }
      }
      else {
        return (i_b + std::min(i_a, i_c)) / 2.f;
      }
    }
    else {
      if (i_b < i_d) {
        return (i_b + std::min(i_a, i_d)) / 2.f;
      }
      else {
        return (i_b + std::max(i_c, i_d)) / 2.f;
      }
    }
  }
}


//! Return the median of 5 values.
float getMedian5(
  const float i_a,
  const float i_b,
  const float i_c,
  const float i_d,
  const float i_e){

  if (i_a < i_b) {
    if (i_a < i_c) {
      if (i_a < i_d) {
        if (i_a < i_e) {
          if (i_b < i_c) {
            if (i_b < i_d) {
              if (i_b < i_e) {
                return std::min(std::min(i_c, i_d), i_e);
              }
              else {
                return i_b;
              }
            }
            else {
              if (i_e < i_b) {
                return std::max(i_d, i_e);
              }
              else {
                return i_b;
              }
            }
          }
          else {
            if (i_b < i_d) {
              if (i_b < i_e) {
                return i_b;
              }
              else {
                return std::max(i_c, i_e);
              }
            }
            else {
              if (i_c < i_d) {
                if (i_c < i_e) {
                  return std::min(std::min(i_b, i_d), i_e);
                }
                else {
                  return i_c;
                }
              }
              else {
                if (i_e < i_c) {
                  return std::max(i_d, i_e);
                }
                else {
                  return i_c;
                }
              }
            }
          }
        }
        else {
          return std::min(std::min(i_b, i_c), i_d);
        }
      }
      else {
        if (i_a < i_e) {
          return std::min(std::min(i_b, i_c), i_e);
        }
        else {
          return i_a;
        }
      }
    }
    else {
      if (i_a < i_d) {
        if (i_a < i_e) {
          return std::min(std::min(i_b, i_d), i_e);
        }
        else {
          return i_a;
        }
      }
      else {
        if (i_a < i_e) {
          return i_a;
        }
        else {
          return std::max(std::max(i_c, i_d), i_e);
        }
      }
    }
  }
  else {
    if (i_b < i_c) {
      if (i_b < i_d) {
        if (i_b < i_e) {
          if (i_a < i_c) {
            if (i_a < i_d) {
              if (i_a < i_e) {
                return std::min(std::min(i_c, i_d), i_e);
              }
              else {
                return i_a;
              }
            }
            else {
              if (i_a < i_e) {
                return i_a;
              }
              else {
                return std::max(i_d, i_e);
              }
            }
          }
          else {
            if (i_c < i_d) {
              if (i_c < i_e) {
                return std::min(std::min(i_a, i_d), i_e);
              }
              else {
                return i_c;
              }
            }
            else {
              if (i_c < i_e) {
                return i_c;
              }
              else {
                return std::max(i_d, i_e);
              }
            }
          }
        }
        else {
          return std::min(std::min(i_a, i_c), i_d);
        }
      }
      else {
        if (i_b < i_e) {
          return std::min(std::min(i_a, i_c), i_e);
        }
        else {
          return i_b;
        }
      }
    }
    else {
      if (i_b < i_d) {
        if (i_b < i_e) {
          return std::min(std::min(i_a, i_d), i_e);
        }
        else {
          return i_b;
        }
      }
      else {
        if (i_b < i_e) {
          return i_b;
        }
        else {
          return std::max(std::max(i_c, i_d), i_e);
        }
      }
    }
  }



  return 0.f;
}



























