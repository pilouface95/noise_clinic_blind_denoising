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
 * @file MultiScalesDenoising.cpp
 * @brief Functions used for the multi scales Denoising method.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "MultiScalesDenoising.h"
#include "../LibNlBayes/NlBayes.h"

using namespace std;

//! Run a multi-scales denoising.
int runMultiScalesDenoising(
  std::vector<float> const& i_imNoisy,
  std::vector<float> &o_imFinal,
  std::vector<std::vector<float> > &o_imDiffs,
  const ImageSize &p_imSize,
  const MsdParams &p_params){

  //! Parameters initializations
  const unsigned int pow = 1 << p_params.nbScales;

  //! Check the size
  if (p_imSize.whc != i_imNoisy.size()) {
    cout << "runMultiScalesDenoising - error : i_imNoisy doesn't have the good size !!!" << endl;
    cout << p_imSize.whc << ", " << i_imNoisy.size() << endl;
    return EXIT_FAILURE;
  }

  //! Adapt size of the image
  ImageSize imSizeNew;
  vector<float> imNoisyNew;
  if (enhanceBoundaries(i_imNoisy, imNoisyNew, p_imSize, imSizeNew, pow) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Multi-scale denoising
  vector<float> imFinalNew;
  vector<vector<float> > imDiffsNew;
  MsdParams params = p_params;

  if (applyMultiScaleDenoising(imNoisyNew, imFinalNew, imDiffsNew, imSizeNew, params) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Remove boundaries
  if (removeBoundaries(o_imFinal, imFinalNew, p_imSize, imSizeNew) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }
  o_imDiffs.resize(imDiffsNew.size());
  for (unsigned int n = 0; n < imDiffsNew.size(); n++) {
    if (removeBoundaries(o_imDiffs[n], imDiffsNew[n], p_imSize, imSizeNew) != EXIT_SUCCESS) {
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}


//! Apply a mosaic multi-scales denoising.
int applyMultiScaleDenoising(
  const std::vector<float> &i_imNoisy,
  std::vector<float> &o_imFinal,
  std::vector<std::vector<float> > &o_imDiffs,
  const ImageSize &p_imSize,
  MsdParams &io_params){

  //! Declarations
  vector<vector<float> > imDownSampled;
  vector<vector<float> > imDiff(io_params.nbScales, vector<float> (p_imSize.whc, 0.f));
  vector<float> imMosaic = i_imNoisy;
  vector<ImageSize> imSizeDown(io_params.nbScales);

  //! Check the size
  if (p_imSize.whc != i_imNoisy.size()) {
    cout << "applyMultiScalesDenoising - error : i_imNoisy doesn't have the good size !!!" << endl;
    cout << p_imSize.whc << ", " << i_imNoisy.size() << endl;
    return EXIT_FAILURE;
  }

  //! Initialization
  o_imFinal.resize(p_imSize.whc);
  imDownSampled.push_back(i_imNoisy);
  imSizeDown[0] = p_imSize;
  o_imDiffs.resize(io_params.nbScales);

  //! Obtain down sampled and differences images for all scales
  for (unsigned int s = 1; s < io_params.nbScales; s++) {

    const unsigned int Nb = 1 << (2 * (s - 1));
    vector<vector<float> > imDiffTmp(Nb);
    vector<vector<float> > imDownSampledTmp(4 * Nb);

    //! For each sub-images
    for (unsigned int n = 0; n < Nb; n++) {

      //! Down Sampling
      vector<vector<float> > imDownTmp;
      meanSplit(imDownSampled[n], imDownTmp, imSizeDown[s - 1], imSizeDown[s]);

      //! Save downSampled image
      for (unsigned int k = 0; k < imDownTmp.size(); k++) {
        imDownSampledTmp[n * 4 + k] = imDownTmp[k];
      }

      //! Up sample the down sampled images
      vector<float> imTmp;
      meanBuilt(imDownTmp, imTmp, imSizeDown[s], imSizeDown[s - 1]);

      //! Compute the difference
      float *iDSn = &imDownSampled[n][0];
      for (unsigned int k = 0; k < imSizeDown[s - 1].whc; k++) {
        imTmp[k] = iDSn[k] - imTmp[k];
      }

      //! Save the difference
      imDiffTmp[n] = imTmp;
    }

    //! Compute the difference mosaic
    constructMosaic(imDiffTmp, imDiff[s], imSizeDown[s - 1], imSizeDown[0]);

    //! Keep only the latest sub-images
    imDownSampled = imDownSampledTmp;

    //! Compute the noisy mosaic only for the latest scale
    if (s == io_params.nbScales - 1) {
      constructMosaic(imDownSampled, imMosaic, imSizeDown[s], imSizeDown[0]);
    }
  }

  //! Estimate the noise, denoise and reconstruct images for all scales
  for (int s = io_params.nbScales - 1; s >= 0; s--) {

    //! Parameters
    io_params.currentScale = s;
    if (io_params.verbose) {
      cout << "Current scale : " << s << endl;
    }

    //! Denoise the Image at This Scale
    vector<float> imDenoised;
    if (runNlBayes(imMosaic, imDenoised, p_imSize, io_params) != EXIT_SUCCESS) {
      return EXIT_FAILURE;
    }

    //! Get the difference (noise removed)
    o_imDiffs[s].resize(p_imSize.whc);
    float* iNr = &o_imDiffs[s][0];
    for (unsigned int k = 0; k < p_imSize.whc; k++) {
      iNr[k] = imDenoised[k] - imMosaic[k];
    }
    unsigned int ns = s;
    while (ns > 0) {
      upgradeMosaic(o_imDiffs[s], imSizeDown, ns--);
    }

    //! Up sample the denoised result
    if (s > 0) {
      upgradeMosaic(imDenoised, imSizeDown, s);
    }

    //! Add the saved details of this scale
    const float* iDs = &imDiff[s][0];
    for (unsigned int k = 0; k < p_imSize.whc; k++) {
      imMosaic[k] = iDs[k] + imDenoised[k];
    }
  }

  //! Get the final denoised image
  o_imFinal = imMosaic;

  return EXIT_SUCCESS;
}










