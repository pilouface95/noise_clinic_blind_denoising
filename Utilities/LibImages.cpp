/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Copyright (c) 2013, Miguel Colom <miguel.colom@cmla.ens-cachan.fr>
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
 * @brief Useful functions on images
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 * @author Miguel Colom <miguel.colom@cmla.ens-cachan.fr>
 **/

#include "LibImages.h"
#include "Utilities.h"
#include "CImage.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

 using namespace std;

//! Load image, check the number of channels.
int loadImage(
  char* p_name,
  std::vector<float> &o_im,
  ImageSize &o_imSize,
  const bool p_verbose){

  //! read input image
  if (p_verbose) {
    cout << endl << "Read input image...";
  }

  //! Initialization
  CImage imTmp;
  imTmp.load(p_name);
  o_imSize.width     = imTmp.get_width();
  o_imSize.height    = imTmp.get_height();
  o_imSize.nChannels = imTmp.get_num_channels();
  o_imSize.wh        = o_imSize.width * o_imSize.height;
  o_imSize.whc       = o_imSize.wh * o_imSize.nChannels;
  o_imSize.nBits     = imTmp.get_bits_per_channel();

  //! Test if image is really a color image and exclude the alpha channel
  if (o_imSize.nChannels > 2 && imTmp.channels_are_equal()) {
    o_imSize.nChannels = 1;
    o_imSize.whc       = o_imSize.wh;
  }

  //! Some image informations
  if (p_verbose) {
    cout << "image size :" << endl;
    cout << " - width          = " << o_imSize.width << endl;
    cout << " - height         = " << o_imSize.height << endl;
    cout << " - nb of channels = " << o_imSize.nChannels << endl;
  }

  //! Get values
  o_im.resize(o_imSize.whc);
  const float factor = (o_imSize.nBits == 16 ? 1.f / 256.f : 1.f);
  for (unsigned int c = 0; c < o_imSize.nChannels; c++) {
    const float* iI = imTmp.get_channel(c);
    float* oI = &o_im[c * o_imSize.wh];

    for (unsigned int k = 0; k < o_imSize.wh; k++) {
      oI[k] = iI[k] * factor;
    }
  }

  return EXIT_SUCCESS;
}


//! write image.
int saveImage(
  char* p_name,
  std::vector<float> const& i_im,
  const ImageSize &p_imSize){

  //! Initialization
  const float minVal = 0.f;
  const float maxVal = (p_imSize.nBits == 8 ? 255.f : 65535.f);

  //! Allocate Memory
  CImage imTmp(p_imSize.width, p_imSize.height,
               p_imSize.nBits, p_imSize.nChannels);

  //! Save values
  const float factor = (p_imSize.nBits == 16 ? 256.f : 1.f);
  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
    const float* iI = &i_im[c * p_imSize.wh];
    float* oI = imTmp.get_channel(c);

    for (unsigned int k = 0; k < p_imSize.wh; k++) {
      oI[k] = std::min(maxVal, std::max(minVal, factor * iI[k]));
    }
  }

  //! Save image
  imTmp.save(p_name, p_imSize.nBits);

  return EXIT_SUCCESS;
}


//! Save a set of difference (noise removed) images.
int saveNoiseRemovalImages(
  std::vector<std::vector<float> > &io_imDiffs,
  const ImageSize &p_imSize,
  const char* p_path){

  //! Initializations
  const unsigned int Nb = io_imDiffs.size();

  //! Sizes check
  for (unsigned int n = 0; n < Nb; n++) {
    if (io_imDiffs[n].size() != p_imSize.whc) {
      cout << "Size of " << n << " difference image is not consistant !!!" << endl;
      return EXIT_FAILURE;
    }
  }

  //! Get standard deviation and absolute maximum
  float sd  = 0.f;
  float max = 0.f;
  const float norm = 1.f / ((float) p_imSize.whc);
  for (unsigned int n = 0; n < Nb; n++) {
    const float* iI = &io_imDiffs[n][0];
    float mean = 0.f;
    float sq = 0.f;

    for (unsigned int k = 0; k < p_imSize.whc; k++) {
      max = std::max(max, fabsf(iI[k]));
      mean += iI[k];
      sq += iI[k] * iI[k];
    }
    mean *= norm;

    sd = std::max(sqrtf(sq * norm - mean * mean) * (1 << n), sd);
  }

  //! For each image
  for (unsigned int n = 0; n < Nb; n++) {

    //! Get the final name
    string imName = p_path;
    imName.erase(imName.end() - 4, imName.end());
    ostringstream oss_n;
    oss_n << n;
    imName += "_diff_" + oss_n.str() + ".png";

    //! Enhance the difference and add an offset
    float* iI = &io_imDiffs[n][0];
    const float coeff = ((max + sd) * max / (2.f * sd) > 120.f ? 120.f / (max + sd) : max / (2.f * sd));
    for (unsigned int k = 0; k < p_imSize.whc; k++) {
      iI[k] = (iI[k] + sd) * coeff + 128.f;
    }

    //! Save the image
    if (saveImage((char*) imName.c_str(), io_imDiffs[n], p_imSize) != EXIT_SUCCESS) {
      return EXIT_SUCCESS;
    }
  }

  return EXIT_SUCCESS;
}


//! Transform the color space of an image, from RGB to YUV, or vice-versa.
void transformColorSpace(
  std::vector<float> &io_im,
  const ImageSize &p_imSize,
  const bool p_isForward){

  //! If the image as only one channel, do nothing
  if (p_imSize.nChannels == 1) {
    return;
  }

  //! Initialization
  vector<float> imTmp(p_imSize.whc);

  //! Only 3 channels
  if (p_imSize.nChannels == 3) {
    //! Initializations
    float* iTr       = &imTmp[0];
    const float* iIr = &io_im[0];
    float* iTg       = &imTmp[p_imSize.wh];
    const float* iIg = &io_im[p_imSize.wh];
    float* iTb       = &imTmp[p_imSize.wh * 2];
    const float* iIb = &io_im[p_imSize.wh * 2];

    //! RGB to YUV
    if (p_isForward) {
      const float a = 1.f / sqrtf(3.f);
      const float b = 1.f / sqrtf(2.f);
      const float c = 2.f * a * sqrtf(2.f);

      for (unsigned int k = 0; k < p_imSize.wh; k++) {
        iTr[k] = a * (iIr[k] + iIg[k] + iIb[k]);
        iTg[k] = b * (iIr[k] - iIb[k]);
        iTb[k] = c * (0.25f * iIr[k] - 0.5f * iIg[k] + 0.25f * iIb[k]);
      }
    }

    //! YUV to RGB
    else {
      const float a = 1.f / sqrtf(3.f);
      const float b = 1.f / sqrtf(2.f);
      const float c = a / b;

      for (unsigned int k = 0; k < p_imSize.wh; k++) {
        iTr[k] = a * iIr[k] + b * iIg[k] + c * 0.5f * iIb[k];
        iTg[k] = a * iIr[k] - c * iIb[k];
        iTb[k] = a * iIr[k] - b * iIg[k] + c * 0.5f * iIb[k];
      }
    }
  }

  //! 4 channels
  else {
    //! Initializations
    float* iTgr       = &imTmp[0];
    const float* iIgr = &io_im[0];
    float* iTr        = &imTmp[p_imSize.wh];
    const float* iIr  = &io_im[p_imSize.wh];
    float* iTb        = &imTmp[p_imSize.wh * 2];
    const float* iIb  = &io_im[p_imSize.wh * 2];
    float* iTgb       = &imTmp[p_imSize.wh * 3];
    const float* iIgb = &io_im[p_imSize.wh * 3];

    //! RGB to YUV
    if (p_isForward) {
      const float a = 0.5f;
      const float b = 1.f / sqrtf(2.f);

      for (unsigned int k = 0; k < p_imSize.wh; k++) {
        iTgr[k] = a * (iIgr[k] + iIr[k] + iIb[k] + iIgb[k]);
        iTr [k] = b * (iIr[k] - iIb[k]);
        iTb [k] = a * (-iIgr[k] + iIr[k] + iIb[k] - iIgb[k]);
        iTgb[k] = b * (-iIgr[k] + iIgb[k]);
      }
    }

    //! YUV to RGB
    else {
      const float a = 0.5f;
      const float b = 1.f / sqrtf(2.f);

      for (unsigned int k = 0; k < p_imSize.wh; k++) {
        iTgr[k] = a * iIgr[k] - a * iIb[k] - b * iIgb[k];
        iTr [k] = a * iIgr[k] + b * iIr[k] + a * iIb [k];
        iTb [k] = a * iIgr[k] - b * iIr[k] + a * iIb [k];
        iTgb[k] = a * iIgr[k] - a * iIb[k] + b * iIgb[k];
      }
    }
  }

  io_im = imTmp;
}


//! Add boundaries by symetry around an image to increase its size until a multiple of p_N.
int enhanceBoundaries(
  std::vector<float> const& i_im,
  std::vector<float> &o_imSym,
  const ImageSize &p_imSize,
  ImageSize &o_imSizeSym,
  const unsigned p_N){

  //! Check the size of the input image
  if (i_im.size() != p_imSize.whc) {
    cout << "Sizes are not consistent !!!" << endl;
    return EXIT_FAILURE;
  }

  //! Initialization
  unsigned int dh = 0;
  unsigned int dw = 0;

  //! Obtain size of boundaries
  if (p_imSize.height % p_N != 0) {
    dh = p_N - p_imSize.height % p_N;
  }
  if (p_imSize.width % p_N != 0) {
    dw = p_N - p_imSize.width % p_N;
  }

  //! Update size of o_imSym
  o_imSizeSym.height    = p_imSize.height + dh;
  o_imSizeSym.width     = p_imSize.width  + dw;
  o_imSizeSym.nChannels = p_imSize.nChannels;
  o_imSizeSym.wh        = o_imSizeSym.width * o_imSizeSym.height;
  o_imSizeSym.whc       = o_imSizeSym.wh * o_imSizeSym.nChannels;
  o_imSizeSym.nBits     = p_imSize.nBits;
  o_imSym.resize(o_imSizeSym.whc);

  //! Particular case
  if (dh == 0 && dw == 0) {
    o_imSym = i_im;
    return EXIT_SUCCESS;
  }

  //! Determine Left, Right, Up and Down boundaries
  const unsigned int dhU = dh / 2;
  const unsigned int dhD = dh - dhU;
  const unsigned int dwL = dw / 2;
  const unsigned int dwR = dw - dwL;

  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {

    //! For convenience
    const unsigned int wi  = p_imSize.width;
    const unsigned int ws  = o_imSizeSym.width;
    const unsigned int dci = c * p_imSize.wh;
    const unsigned int dcs = c * o_imSizeSym.wh;

    //! Center of the image
    for (unsigned int i = 0; i < p_imSize.height; i++) {
      float* oS       = &o_imSym[dcs + (i + dhU) * ws + dwL];
      const float* iI = &i_im   [dci + i * wi];

      for (unsigned int j = 0; j < p_imSize.width; j++) {
        oS[j] = iI[j];
      }
    }

    //! Up
    for (unsigned int i = 0; i < dhU; i++) {
      const float* iI = &i_im   [dci + (dhU - i) * wi];
      float* oS       = &o_imSym[dcs + i * ws + dwL];

      for (unsigned int j = 0; j < p_imSize.width; j++) {
        oS[j] = iI[j];
      }
    }

    //! Down
    for (unsigned int i = 0; i < dhD; i++) {
      const float *iI = &i_im   [dci + (p_imSize.height - i - 2) * wi];
      float* oS       = &o_imSym[dcs + (p_imSize.height + i + dhU) * ws + dwL];

      for (unsigned int j = 0; j < p_imSize.width; j++) {
        oS[j] = iI[j];
      }
    }

    //! Left
    for (unsigned int i = 0; i < o_imSizeSym.height; i++) {
      float* oS = &o_imSym[dcs + i * ws];

      for (int j = 0; j < (int) dwL; j++) {
        oS[j] = oS[2 * dwL - j];
      }
    }

    //! Right
    for (unsigned int i = 0; i < o_imSizeSym.height; i++) {
      float* oS = &o_imSym[dcs + i * ws + dwL];

      for (int j = 0; j < (int) dwR; j++) {
        oS[j + wi] = oS[wi - j - 2];
      }
    }
  }

  return EXIT_SUCCESS;
}


//! Remove boundaries added with enhanceBoundaries
int removeBoundaries(
  std::vector<float> &o_im,
  std::vector<float> const& i_imSym,
  const ImageSize &p_imSize,
  const ImageSize &p_imSizeSym){

  //! Check if sizes are consistent
  if (i_imSym.size() != p_imSizeSym.whc) {
    cout << "Sizes are not consistent !!! " << endl;
    return EXIT_FAILURE;
  }
  if (p_imSize.width > p_imSizeSym.width || p_imSize.height > p_imSizeSym.height
    || p_imSize.nChannels != p_imSizeSym.nChannels) {
    cout << "Inner image is too big" << endl;
    return EXIT_FAILURE;
  }
  if (o_im.size() != p_imSize.whc) {
    o_im.resize(p_imSize.whc);
  }

  //! Obtain size of boundaries
  const unsigned int dhU = (p_imSizeSym.height - p_imSize.height) / 2;
  const unsigned int dwL = (p_imSizeSym.width  - p_imSize.width ) / 2;

  //! Get the inner image
  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
    const float* iS = &i_imSym[c * p_imSizeSym.wh + dhU * p_imSizeSym.width + dwL];
    float* oI       = &o_im   [c * p_imSize.wh];

    for (unsigned int i = 0; i < p_imSize.height; i++) {
      for (unsigned int j = 0; j < p_imSize.width; j++) {
        oI[i * p_imSize.width + j] = iS[i * p_imSizeSym.width + j];
      }
    }
  }

  return EXIT_SUCCESS;
}


//! Split an image into 4 sub-scale images, 4 times smaller by averaging a 2x2 square of pixels
void meanSplit(
  std::vector<float> const& i_im,
  std::vector<std::vector<float> >& o_subIm,
  const ImageSize &i_imSize,
  ImageSize &o_imSize){

  //! Built size of sub-scale images
  o_imSize.width      = i_imSize.width  / 2;
  o_imSize.height     = i_imSize.height / 2;
  o_imSize.nChannels  = i_imSize.nChannels;
  o_imSize.wh         = o_imSize.width * o_imSize.height;
  o_imSize.whc        = o_imSize.wh * o_imSize.nChannels;

  //! Initializations
  o_subIm.resize(4);
  for (unsigned int n = 0; n < 4; n++) {
    o_subIm[n].resize(o_imSize.whc);
  }
  const unsigned int wi = i_imSize.width;
  const unsigned int wo = o_imSize.width;

  //! Average pixels in 4 different ways
  for (unsigned int c = 0; c < i_imSize.nChannels; c++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < o_imSize.height; i++) {

      //! Initializations
      const float* iI0 = &i_im[c * i_imSize.wh + 2 * i * wi];
      const float* iI1 = iI0 + 1;
      const float* iI2 = iI0 + wi;
      const float* iI3 = iI0 + wi + 1;
      float* oS0       = &o_subIm[0][c * o_imSize.wh + i * wo];
      float* oS1       = &o_subIm[1][c * o_imSize.wh + i * wo];
      float* oS2       = &o_subIm[2][c * o_imSize.wh + i * wo];
      float* oS3       = &o_subIm[3][c * o_imSize.wh + i * wo];

      for (unsigned int j = 0; j < o_imSize.width; j++) {

        //! Image 0: top-left
        float sum = 0.f;
        for (unsigned int p = 0; p < 2; p++) {
          for (unsigned int q = 0; q < 2; q++) {
            sum += iI0[2 * j + p * wi + q];
          }
        }
        oS0[j] = sum * 0.25f;

        //! Image 1: top-right
        if (j == o_imSize.width - 1) {
          oS1[j] = (iI1[2 * j] + iI1[2 * j + wi]) * 0.5f;
        }
        else {
          float sum = 0.f;
          for (unsigned int p = 0; p < 2; p++) {
            for (unsigned int q = 0; q < 2; q++) {
              sum += iI1[2 * j + p * wi + q];
            }
          }
          oS1[j] = sum * 0.25f;
        }

        //! Image 2: bottom-left
        if (i == o_imSize.height- 1) {
          oS2[j] = (iI2[2 * j] + iI2[2 * j + 1]) * 0.5f;
        }
        else {
          float sum = 0.f;
          for (unsigned int p = 0; p < 2; p++) {
            for (unsigned int q = 0; q < 2; q++) {
              sum += iI2[2 * j + p * wi + q];
            }
          }
          oS2[j] = sum * 0.25f;
        }

        //! Image 3: bottom-right
        if (i == o_imSize.height - 1 && j == o_imSize.width - 1) {
          oS3[j] = iI3[2 * j];
        }
        else if (j == o_imSize.width - 1) {
          oS3[j] = (iI3[2 * j] + iI3[2 * j + wi]) * 0.5f;
        }
        else if (i == o_imSize.height- 1) {
          oS3[j] = (iI3[2 * j] + iI3[2 * j + 1]) * 0.5f;
        }
        else {
          float sum = 0.f;
          for (unsigned int p = 0; p < 2; p++) {
            for (unsigned int q = 0; q < 2; q++) {
              sum += iI3[2 * j + p * wi + q];
            }
          }
          oS3[j] = sum * 0.25f;
        }
      }
    }
  }
}


//! Reconstruct an image splitted by te meanSplit function.
void meanBuilt(
  std::vector<std::vector<float> > const& i_subIm,
  std::vector<float> &o_im,
  const ImageSize &i_imSize,
  ImageSize &o_imSize){

  //! Obtain final size
  o_imSize.width     = i_imSize.width * 2;
  o_imSize.height    = i_imSize.height * 2;
  o_imSize.nChannels = i_imSize.nChannels;
  o_imSize.wh        = o_imSize.width * o_imSize.height;
  o_imSize.whc       = o_imSize.wh * o_imSize.nChannels;
  o_im.resize(o_imSize.whc);

  //! Reconstruct the image from sub-images
  for (unsigned int c = 0; c < i_imSize.nChannels; c++) {
    const unsigned int dci = c * i_imSize.wh;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < o_imSize.height; i++) {
      const unsigned int di1 = dci + (unsigned) floor(float(i    ) * 0.5f) * i_imSize.width;
      const unsigned int di2 = dci + (unsigned) floor(float(i - 1) * 0.5f) * i_imSize.width;
      const float *iS0       = &i_subIm[0][di1];
      const float *iS1       = &i_subIm[1][di1];
      const float *iS2       = &i_subIm[2][di2];
      const float *iS3       = &i_subIm[3][di2];
      float* oI              = &o_im[c * o_imSize.wh + i * o_imSize.width];

      for (unsigned int j = 0; j < o_imSize.width; j++) {
        const unsigned int dj1 = (unsigned) floor(float(j    ) * 0.5f);
        const unsigned int dj2 = (unsigned) floor(float(j - 1) * 0.5f);

        //! Top left pixel
        if (i == 0 && j == 0) {
          oI[j] = iS0[dj1];
        }

        //! First row
        else if (i == 0) {
          oI[j] = (iS0[dj1] + iS1[dj2]) * 0.5f;
        }

        //! First column
        else if (j == 0) {
          oI[j] = (iS0[dj1] + iS2[dj1]) * 0.5f;
        }

        //! Reconstruct the image
        else {
          oI[j] = (iS0[dj1] + iS1[dj2] + iS2[dj1] + iS3[dj2]) * 0.25f;
        }
      }
    }
  }
}


//! Compose a mosaic of sub-images, sticked by symetry.
void constructMosaic(
  std::vector<std::vector<float> > const& i_subIm,
  std::vector<float> &o_im,
  const ImageSize &i_imSize,
  ImageSize &o_imSize){

  //! Particular case
  if (i_subIm.size() == 1) {
    o_im = i_subIm[0];
    o_imSize = i_imSize;
    return;
  }

  //! Initialization
  const unsigned int N = (unsigned int) sqrtf((float) i_subIm.size());
  o_imSize.width       = N * i_imSize.width;
  o_imSize.height      = N * i_imSize.height;
  o_imSize.nChannels   = i_imSize.nChannels;
  o_imSize.wh          = o_imSize.width * o_imSize.height;
  o_imSize.whc         = o_imSize.wh * o_imSize.nChannels;
  o_im.resize(o_imSize.whc);

  //! Construct the mosaic
  for (unsigned int c = 0; c < o_imSize.nChannels; c++) {
    for (unsigned int p = 0; p < N; p++) {
      for (unsigned int q = 0; q < N; q++) {
        const float* iSn = &i_subIm[p * N + q][c * i_imSize.wh];
        float* oI        = &o_im[c * o_imSize.wh + p * i_imSize.height * o_imSize.width + q * i_imSize.width];

        for (unsigned int i = 0; i < i_imSize.height; i++) {
          const unsigned int di = (p % 2 == 0 ? i : i_imSize.height - 1 - i) * o_imSize.width;

          for (unsigned int j = 0; j < i_imSize.width; j++) {
            const unsigned int dj = (q % 2 == 0 ? j : i_imSize.width - 1 - j);

            oI[di + dj] = iSn[i * i_imSize.width + j];
          }
        }
      }
    }
  }
}


//! Upgrade in-place a mosaic: from a mosaic of 4^n sub-images, built a mosaic with 4^(n-1) sub-images.
void upgradeMosaic(
  std::vector<float> &io_mosaic,
  std::vector<ImageSize> &p_imSizes,
  const unsigned int p_currentScale){

  //! Declarations
  vector<vector<float> > imSplit;
  const unsigned int Nb = 1 << (2 * p_currentScale);
  const unsigned int N = Nb / 4;

  //! Split the mosaic in Nb sub-images
  splitMosaic(io_mosaic, imSplit, p_imSizes[0], p_imSizes[p_currentScale], Nb);

  //! Reconstruct sub-images of the previous scale
  for (unsigned int n = 0; n < N; n++) {
    vector<vector<float> > imTmp(4);
    for (unsigned int k = 0; k < 4; k++) {
      imTmp[k] = imSplit[n * 4 + k];
    }
    meanBuilt(imTmp, io_mosaic, p_imSizes[p_currentScale], p_imSizes[p_currentScale - 1]);
    imSplit.push_back(io_mosaic);
  }

  //! Keep only the new sub-images
  imSplit.erase(imSplit.begin(), imSplit.begin() + Nb);

  //! Build the new mosaic with the new sub-images
  constructMosaic(imSplit, io_mosaic, p_imSizes[p_currentScale - 1], p_imSizes[0]);
}


//! Split a mosaic obtained with the constructMosaic function into sub-images.
void splitMosaic(
  std::vector<float> const& i_im,
  std::vector<std::vector<float> >& o_subIm,
  const ImageSize &i_imSize,
  ImageSize &o_imSize,
  const unsigned int p_N){

  //! Particular case
  if (p_N == 1) {
    o_imSize = i_imSize;
    o_subIm.clear();
    o_subIm.push_back(i_im);
    return;
  }

  //! Built size of sub-scale images
  const unsigned int N = (unsigned int) sqrtf((float) p_N);
  o_imSize.width       = i_imSize.width  / N;
  o_imSize.height      = i_imSize.height / N;
  o_imSize.nChannels   = i_imSize.nChannels;
  o_imSize.wh          = o_imSize.width * o_imSize.height;
  o_imSize.whc         = o_imSize.wh * o_imSize.nChannels;
  o_subIm.resize(p_N);
  for (unsigned int n = 0; n < p_N; n++) {
    o_subIm[n].resize(o_imSize.whc);
  }

  //! Split the mosaic
  for (unsigned int c = 0; c < o_imSize.nChannels; c++) {
    for (unsigned int p = 0; p < N; p++) {
      for (unsigned int q = 0; q < N; q++) {
      const float* iI = &i_im[c * i_imSize.wh + p * o_imSize.height * i_imSize.width + q * o_imSize.width];
      float* oSn      = &o_subIm[p * N + q][c * o_imSize.wh];

      for (unsigned int i = 0; i < o_imSize.height; i++) {
        const unsigned int di = (p % 2 == 0 ? i : o_imSize.height - 1 - i) * i_imSize.width;

        for (unsigned int j = 0; j < o_imSize.width; j++) {
            const unsigned int dj = (q % 2 == 0 ? j : o_imSize.width - 1 - j);

            oSn[i * o_imSize.width + j] = iI[di + dj];
          }
        }
      }
    }
  }
}


//! Add boundaries around image to avoid boundary problems for the big loop over pixels.
void setBoundaries(
  std::vector<float> const& i_im,
  std::vector<float> &o_im,
  const ImageSize &i_imSize,
  ImageSize &o_imSize,
  const unsigned int p_border){

  //! Get the new size
  o_imSize.width     = i_imSize.width + 2 * p_border;
  o_imSize.height    = i_imSize.height + 2 * p_border;
  o_imSize.nChannels = i_imSize.nChannels;
  o_imSize.wh        = o_imSize.width * o_imSize.height;
  o_imSize.whc       = o_imSize.wh * o_imSize.nChannels;
  o_im.assign(o_imSize.whc, 0.f);

  //! For convenience
  const unsigned wI = i_imSize.width;
  const unsigned wO = o_imSize.width;
  const unsigned hO = o_imSize.height;

  //! Get the new (bigger) image
  for (unsigned c = 0; c < i_imSize.nChannels; c++) {
    const float* iI = &i_im[c * i_imSize.wh];

    //! Inner image
    float* oI = &o_im[c * o_imSize.wh + p_border * wO + p_border];
    for (unsigned i = 0; i < i_imSize.height; i++) {
      for (unsigned j = 0; j < i_imSize.width; j++) {
        oI[i * wO + j] = iI[i * wI + j];
      }
    }

    //! Left and right
    oI = &o_im[c * o_imSize.wh + p_border * wO];
    for (unsigned i = 0; i < i_imSize.height; i++) {
      for (unsigned j = 0; j < p_border; j++) {
        oI[i * wO + j] = iI[i * wI + p_border - j];
        oI[i * wO + wO - p_border + j] = iI[i * wI + wI - j - 2];
      }
    }

    //! Top and bottom
    oI = &o_im[c * o_imSize.wh];
    for (unsigned i = 0; i < p_border; i++) {
      for (unsigned j = 0; j < o_imSize.width; j++) {
        oI[i * wO + j] = oI[(2 * p_border - i) * wO + j];
        oI[(i + hO - p_border) * wO + j] = oI[(hO - p_border - i - 1) * wO + j];
      }
    }
  }
}


//! Divide an image into sub-parts of equal sizes with a boundaries around each.
int subDivide(
  std::vector<float> const& i_im,
  std::vector<std::vector<float> > &o_imSub,
  const ImageSize &i_imSize,
  ImageSize &o_imSizeSub,
  const unsigned int p_boundaries,
  const unsigned int p_nbParts){

  //! Add boundaries around the image
  unsigned int a, b;
  getDivisors(p_nbParts, a, b);
  const unsigned int d = getPpcm(a, b);
  ImageSize imSizeBig;
  vector<float> imBig;
  if (enhanceBoundaries(i_im, imBig, i_imSize, imSizeBig, d) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Get the size of sub-images
  vector<float> imNew;
  ImageSize imSizeNew;
  setBoundaries(imBig, imNew, imSizeBig, imSizeNew, p_boundaries);
  o_imSizeSub.height    = imSizeBig.height / a + 2 * p_boundaries;
  o_imSizeSub.width     = imSizeBig.width  / b + 2 * p_boundaries;
  o_imSizeSub.nChannels = imSizeBig.nChannels;
  o_imSizeSub.wh        = o_imSizeSub.width * o_imSizeSub.height;
  o_imSizeSub.whc       = o_imSizeSub.wh * o_imSizeSub.nChannels;
  o_imSub.resize(p_nbParts);
  for (unsigned int n = 0; n < p_nbParts; n++) {
    o_imSub[n].resize(o_imSizeSub.whc);
  }

  //! For convenience
  const unsigned int ws = o_imSizeSub.width;
  const unsigned int wn = imSizeNew.width;

  //! Get all sub-images with the boundaries
  for (unsigned int c = 0; c < i_imSize.nChannels; c++) {
    for (unsigned int ni = 0; ni < a; ni++) {
      for (unsigned int nj = 0; nj < b; nj++) {
        const unsigned int di = ni * (o_imSizeSub.height - 2 * p_boundaries);
        const unsigned int dj = nj * (o_imSizeSub.width  - 2 * p_boundaries);
        float* oS = &o_imSub[ni * b + nj][c * o_imSizeSub.wh];
        const float* iN = &imNew[c * imSizeNew.wh + di * wn + dj];

        for (unsigned int i = 0; i < o_imSizeSub.height; i++) {
          for (unsigned int j = 0; j < o_imSizeSub.width; j++) {
            oS[i * ws + j] = iN[i * wn + j];
          }
        }
      }
    }
  }

  return EXIT_SUCCESS;
}


//! Build an image from a set of sub-parts obtained with subDivide.
int subBuild(
  std::vector<std::vector<float> > const& i_imSub,
  std::vector<float> &o_im,
  const ImageSize &p_imSizeSub,
  const ImageSize &p_imSize,
  const unsigned int p_boundaries){

  //! Initializations
  const unsigned int nbParts = i_imSub.size();
  unsigned int a, b;
  getDivisors(nbParts, a, b);
  const unsigned int hs = p_imSizeSub.height;
  const unsigned int ws = p_imSizeSub.width;
  const unsigned int hn = (hs - 2 * p_boundaries) * a;
  const unsigned int wn = (ws - 2 * p_boundaries) * b;
  //ImageSize imSizeNew   = {wn, hn, p_imSize.nChannels, -1, wn * hn, wn * hn * p_imSize.nChannels};
  ImageSize imSizeNew(wn, hn, p_imSize.nChannels, -1,
                      wn * hn, wn * hn * p_imSize.nChannels);


  o_im.resize(p_imSize.whc);
  vector<float> imNew(imSizeNew.whc);

  //! Get the image
  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
    for (unsigned int ni = 0; ni < a; ni++) {
      for (unsigned int nj = 0; nj < b; nj++) {
        const unsigned int di = ni * (hs - 2 * p_boundaries);
        const unsigned int dj = nj * (ws - 2 * p_boundaries);
        const float* iS = &i_imSub[ni * b + nj][c * p_imSizeSub.wh + p_boundaries * ws + p_boundaries];
        float* iN = &imNew[c * imSizeNew.wh + di * wn + dj];
        for (unsigned int i = 0; i < hs - 2 * p_boundaries; i++) {
          for (unsigned int j = 0; j < ws - 2 * p_boundaries; j++) {
            iN[i * wn + j] = iS[i * ws + j];
          }
        }
      }
    }
  }

  //! Extract image
  if (removeBoundaries(o_im, imNew, p_imSize, imSizeNew) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


//! Save the difference image and compute the RMSE of it.
int saveDiff(
  std::vector<float> const& i_imRef,
  std::vector<float> const& i_imFinal,
  const ImageSize &p_imSize,
  const char* p_imName){

  //! Check the size
  if (i_imRef.size() != i_imFinal.size()) {
    cout << "Images aren't of the same size !!!" << endl;
    return EXIT_FAILURE;
  }
  const unsigned int whc = i_imRef.size();

  //! Compute the difference
  vector<float> imDiff(whc);
  double mean = 0.f;
  double std  = 0.f;
  float max   = 0.f;
  for (unsigned int k = 0; k < whc; k++) {
    imDiff[k] = i_imRef[k] - i_imFinal[k];
    mean += (double) imDiff[k];
    std  += (double) imDiff[k] * imDiff[k];
    max = std::max(max, fabsf(imDiff[k]));
  }

  //! Compute the rmse
  mean /= (double) whc;
  const float rmse = sqrtf((float) (std / (double) whc - mean * mean));
  cout << "RMSE: " << rmse << endl;

  //! Enhance the difference image
  const float coeff = ((max + rmse) * max / (2.f * rmse) > 120.f ? 120.f / (max + rmse) : max / (2.f * rmse));
  for (unsigned int k = 0; k < whc; k++) {
    imDiff[k] = (imDiff[k] + rmse) * coeff + 128.f;
  }

  //! Reescale to 16 bits if the input image is 16 bits
  if (p_imSize.nBits == 16) {
    for (unsigned int i = 0; i < p_imSize.whc; i++) {
      imDiff[i] *= 256.0;
    }
  }

  //! Save the difference image
  if (saveImage((char*) p_imName, imDiff, p_imSize) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Save the rmse
  string fileName = p_imName;
  fileName.erase(fileName.end() - 4, fileName.end());
  fileName += "_rmse.txt";
  ofstream file((char*) fileName.c_str(), ios::out | ios::trunc);

  //! Check the opening
  if (!file) {
    fprintf(stderr, "Couldn't create the file: %s\n", fileName.c_str());
    return EXIT_FAILURE;
  }

  //! Write measure
  file << "RMSE = " << rmse << endl;
  file.close();

  return EXIT_SUCCESS;
}


// Image size information class
ImageSize::ImageSize() {
}

ImageSize::ImageSize(const unsigned width, const unsigned height,
                     const unsigned nChannels, const unsigned nBits,
                     const unsigned wh, const unsigned whc) {
  this->width = width;
  this->height = height;
  this->nChannels = nChannels;
  this->nBits = nBits;
  this->wh = wh;
  this->whc = whc;
}
