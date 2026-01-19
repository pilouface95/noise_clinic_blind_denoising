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
 * @file NlBayes.cpp
 * @brief NlBayes denoising algorithm.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 * @author Miguel Colom <miguel.colom@cmla.ens-cachan.fr>
 **/

#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <cmath>

#include "NlBayes.h"
#include "LibMatrix.h"
#include "../Utilities/Utilities.h"
#include "../NoiseEstimate/RunNoiseEstimate.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

//! Initialize Parameters of the NL-Bayes algorithm.
int initializeNlbParameters(
  nlbParams &o_paramStep1,
  nlbParams &o_paramStep2,
  MsdParams &io_msdParams,
  std::vector<float> const& i_imNoisy,
  const ImageSize &p_imSize){

  //! Verbose
  o_paramStep1.verbose = io_msdParams.verbose;
  o_paramStep2.verbose = io_msdParams.verbose;

  //! Is finest scale ?
  o_paramStep1.isFinestScale = (io_msdParams.currentScale == 0);
  o_paramStep2.isFinestScale = (io_msdParams.currentScale == 0);

  //! Size of patches
  o_paramStep1.sizePatch = 4;
  o_paramStep2.sizePatch = 4;

  //! Number of similar patches
  o_paramStep1.nSimilarPatches = 2 * o_paramStep1.sizePatch * o_paramStep1.sizePatch;
  o_paramStep2.nSimilarPatches = 2 * o_paramStep2.sizePatch * o_paramStep2.sizePatch;

  //! Maximum number of similar patches taken during the second step
  o_paramStep1.maxSimilarPatches = 16 * o_paramStep1.nSimilarPatches;
  o_paramStep2.maxSimilarPatches = 16 * o_paramStep2.nSimilarPatches;

  //! Size of the search window around the reference patch (must be odd)
  o_paramStep1.sizeSearchWindow = o_paramStep1.nSimilarPatches / 2;
  o_paramStep1.sizeSearchWindow |= 1;
  o_paramStep2.sizeSearchWindow = o_paramStep2.nSimilarPatches / 2;
  o_paramStep2.sizeSearchWindow |= 1;

  //! Parameter used to determine similar patches
  o_paramStep1.tau = 3.f * (float) (o_paramStep1.sizePatch * o_paramStep1.sizePatch);
  o_paramStep2.tau = 3.f * (float) (o_paramStep2.sizePatch * o_paramStep2.sizePatch * p_imSize.nChannels);

  //! Is first step?
  o_paramStep1.isFirstStep = true;
  o_paramStep2.isFirstStep = false;

  //! Compute the covariance matrix of the noise
  if (computeNoiseCovarianceMatrix(i_imNoisy, p_imSize, o_paramStep1, io_msdParams) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }
  if (computeNoiseCovarianceMatrix(i_imNoisy, p_imSize, o_paramStep2, io_msdParams) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Threshold to decide if a 3D group belongs to an homogeneous area
  o_paramStep1.kLow  = (o_paramStep1.isFinestScale ? 4.f : 2.f);
  o_paramStep1.kHigh = (o_paramStep1.isFinestScale ? 8.f : 4.f);
  o_paramStep1.gamma = 1.0f;

  //! Use Area
  o_paramStep1.useArea = true;

  //! Use random patches
  o_paramStep1.useRandomPatches = true;
  o_paramStep2.useRandomPatches = true;

  return EXIT_SUCCESS;
}


//! Main function to process the whole NL-Bayes algorithm.
int runNlBayes(
  std::vector<float> const& i_imNoisy,
  std::vector<float> &o_imFinal,
  const ImageSize &p_imSize,
  MsdParams &io_msdParams){

  //! Check
  const unsigned int chnls = p_imSize.nChannels;
  if (chnls != 1 && chnls != 3 && chnls != 4) {
    cout << "runNlBayes - error : Number of channels not valid. Must be 1, 3 or 4." << endl;
    return EXIT_FAILURE;
  }
  if (i_imNoisy.size() != p_imSize.whc) {
    cout << "runNlBayes - error : Size of the image not conform." << endl;
    return EXIT_FAILURE;
  }

  //! Initialization
  std::vector<float> imBasic(i_imNoisy.size());
  const float alpha = 1.5f;

  //! Number of available cores
#ifdef _OPENMP
  const unsigned int nbParts = omp_get_max_threads();
#else
  const unsigned int nbParts = 1;
#endif

  //! Parameters Initialization
  nlbParams paramStep1, paramStep2;
  if (initializeNlbParameters(paramStep1, paramStep2, io_msdParams, i_imNoisy, p_imSize) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! First step
  if (paramStep1.verbose) {
    cout << "First step" << endl;
  }

  //! RGB to YUV
  vector<float> imNoisy = i_imNoisy;
  transformColorSpace(imNoisy, p_imSize, true);

  //! Divide the noisy image into sub-images in order to easier parallelize the process
  vector<vector<float> > imNoisySub(nbParts), imBasicSub(nbParts), imFinalSub(nbParts);
  ImageSize imSizeSub;
  if (subDivide(imNoisy, imNoisySub, p_imSize, imSizeSub, alpha * paramStep1.sizeSearchWindow, nbParts) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Call denoising
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
            shared(imNoisySub, imBasicSub, imFinalSub, imSizeSub) \
            firstprivate (paramStep1)
#endif
  for (int n = 0; n < (int) nbParts; n++) {
    processNlBayes(imNoisySub[n], imBasicSub[n], imFinalSub[n], imSizeSub, paramStep1);
  }

  //! Get the basic estimate
  if (subBuild(imBasicSub, imBasic, imSizeSub, p_imSize, alpha * paramStep1.sizeSearchWindow) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! YUV to RGB
  transformColorSpace(imBasic, p_imSize, false);

  //! 2nd Step
  if (paramStep2.verbose) {
    cout << "Second step" << endl;
  }

  //! Divide the noisy and basic images into sub-images in order to easier parallelize the process
  if (subDivide(i_imNoisy, imNoisySub, p_imSize, imSizeSub, alpha * paramStep2.sizeSearchWindow, nbParts) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }
  if (subDivide(imBasic, imBasicSub, p_imSize, imSizeSub, alpha * paramStep2.sizeSearchWindow, nbParts) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Call denoising
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
            shared(imNoisySub, imBasicSub, imFinalSub, imSizeSub) \
            firstprivate (paramStep2)
#endif
  for (int n = 0; n < (int) nbParts; n++) {
    processNlBayes(imNoisySub[n], imBasicSub[n], imFinalSub[n], imSizeSub, paramStep2);
  }

  //! Get the final result
  if (subBuild(imFinalSub, o_imFinal, imSizeSub, p_imSize, alpha * paramStep2.sizeSearchWindow) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


//! Generic step of the NL-Bayes denoising (could be the first or the second).
void processNlBayes(
  std::vector<float> const& i_imNoisy,
  std::vector<float> &io_imBasic,
  std::vector<float> &o_imFinal,
  const ImageSize &p_imSize,
  nlbParams &p_params){

  //! Parameters initialization
  const unsigned int sW  = p_params.sizeSearchWindow;
  const unsigned int sP  = p_params.sizePatch;
  const unsigned int sP2 = sP * sP;
  const unsigned int sPC = sP2 * p_imSize.nChannels;
  const unsigned int nSP = p_params.maxSimilarPatches;

  //! Allocate Sizes
  p_params.sigma.resize(p_imSize.nChannels);
  if (p_params.isFirstStep) {
    io_imBasic.resize(p_imSize.whc);
    p_params.barycenter .resize(sPC);
    if (p_params.useArea) {
      p_params.meanPatches.resize(p_imSize.nChannels);
      p_params.noiseTrace .resize(p_imSize.nChannels);
      p_params.patchTrace .resize(p_imSize.nChannels);
      p_params.noiseTrace .resize(p_imSize.nChannels);
    }
  }
  else {
    o_imFinal.resize(p_imSize.whc);
  }

  //! weight sum per pixel
  std::vector<float> weight(i_imNoisy.size(), 0.f);

  //! Used matrices during Bayes' estimate
  std::vector< std::vector<float> > group3d(p_imSize.nChannels, std::vector<float> (nSP * sP2));
  std::vector<float> group3dNoisy(nSP * sPC), group3dBasic(nSP * sPC);
  std::vector<unsigned int> index(nSP);

  matParams mat;
  mat.barycenter.resize(p_params.isFirstStep ? sP2 : sPC);
  mat.covMat    .resize(p_params.isFirstStep ? sP2 * sP2 : sPC * sPC);
  mat.covMatTmp .resize(p_params.isFirstStep ? sP2 * sP2 : sPC * sPC);
  mat.group3d   .resize(p_params.isFirstStep ? nSP * sP2 : sW * sW * sPC);

  //! Only pixels of the center of the image must be processed (not the boundaries)
  std::vector<bool> mask(p_imSize.wh, false);
  for (unsigned int i = sW; i < p_imSize.height - sW; i++) {
    for (unsigned int j = sW; j < p_imSize.width - sW; j++) {
      mask[i * p_imSize.width + j] = true;
    }
  }

  //! Loop over pixels
  for (unsigned int ij = 0; ij < p_imSize.wh; ij++) {

    //! Only non-seen patches are processed
    if (mask[ij]) {

      //! Get the average sigma corresponding to the reference patch
      p_params.indRef = ij;
      getSigma(ij, i_imNoisy, p_imSize, p_params);

      //! Search for similar patches around the reference one
      if (p_params.isFirstStep) {
        estimateSimilarPatchesStep1(i_imNoisy, group3d, index, p_imSize, p_params);
      }
      else {
        estimateSimilarPatchesStep2(i_imNoisy, io_imBasic, group3dNoisy, group3dBasic, index, p_imSize, p_params);
      }

      //! Compute the Bayes estimation
      if (p_params.isFirstStep) {
        computeBayesEstimateStep1(group3d, mat, p_params);

        //! If desired, compute the Area trick
        if (p_params.useArea) {
          computeSoftHomogeneousAreaStep1(group3d, index, p_params, p_imSize);
        }
      }
      else {
        computeBayesEstimateStep2(group3dNoisy, group3dBasic, mat, p_imSize, p_params);
      }

      //! Aggregation
      if (p_params.isFirstStep) {
        computeAggregationStep1(io_imBasic, weight, mask, group3d, index, p_imSize, p_params);
      }
      else {
        computeAggregationStep2(o_imFinal, weight, mask, group3dBasic, index, p_imSize, p_params);
      }
    }
  }

  //! Final Weighted aggregation
  computeWeightedAggregation(i_imNoisy, p_params.isFirstStep ? io_imBasic : o_imFinal, weight, p_params, p_imSize);
}


//! Estimate the best similar patches to a reference one.
void estimateSimilarPatchesStep1(
  std::vector<float> const& i_im,
  std::vector<std::vector<float> > &o_group3d,
  std::vector<unsigned int> &o_index,
  const ImageSize &p_imSize,
  nlbParams &io_params){

  //! Initialization
  const unsigned int sW    = io_params.sizeSearchWindow;
  const unsigned int sP    = io_params.sizePatch;
  const unsigned int w     = p_imSize.width;
  const unsigned int wh    = p_imSize.wh;

  //! Will contain distances and position of patches in the search windows.
  vector<pair<float, unsigned int> > distance(sW * sW);

  //! Compute weight according to channels
  vector<float> weightTable;
  getWeight(weightTable, p_imSize.nChannels, true);

  //! Get the threshold value according to the weighting of the distance
  float sigmaMean = 0.f;
  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
    sigmaMean += io_params.sigma[c] * io_params.sigma[c] * weightTable[c] * weightTable[c];
  }
  const float tau = io_params.tau * sigmaMean;

  //! Compute distance between patches
  computeDistance(distance, io_params, weightTable, i_im, p_imSize);

  //! Get the similar patches
  const unsigned int nSimP = getSimilarPatches(o_index, distance, io_params, tau);

  //! Register similar patches into the 3D group
  //! and update the barycenter for the large set of patches (homogeneous areas)
  const float normBary  = 1.f / (float) nSimP;
  const float normPatch = 1.f / (float) (nSimP * sP * sP);
  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
    float* iG = &o_group3d[c][0];
    float mean = 0.f;
    float* iB = &io_params.barycenter[c * sP * sP];
    const float* iI = &i_im[c * wh];

    for (unsigned p = 0, k = 0; p < sP; p++) {
      for (unsigned q = 0; q < sP; q++) {
        float bary = 0.f;

        for (unsigned int n = 0; n < nSimP; n++, k++) {
          float imValue = iI[o_index[n] + p * w + q];
          iG[k] = imValue;
          mean += imValue;
          bary += imValue;
        }
        iB[p * sP + q] = bary * normBary;
      }
    }

    //! Kept the mean of all patches
    if (io_params.useArea) {
      io_params.meanPatches[c] = mean * normPatch;
    }
  }

  //! Save the current number of patches
  io_params.nSimP = nSimP;
}


//! Keep from all near patches the similar ones to the reference patch for the second step.
void estimateSimilarPatchesStep2(
  std::vector<float> const& i_imNoisy,
  std::vector<float> const& i_imBasic,
  std::vector<float> &o_group3dNoisy,
  std::vector<float> &o_group3dBasic,
  std::vector<unsigned int> &o_index,
  const ImageSize &p_imSize,
  nlbParams &io_params){

  //! Initialization
  const unsigned int w     = p_imSize.width;
  const unsigned int chnls = p_imSize.nChannels;
  const unsigned int wh    = p_imSize.wh;
  const unsigned int sP    = io_params.sizePatch;
  const unsigned int sW    = io_params.sizeSearchWindow;

  //! Will contain distances and position of patches in the search windows.
  vector<pair<float, unsigned int> > distance(sW * sW);

  //! Compute weight according to channels
  vector<float> weightTable(p_imSize.nChannels, 1.f);

  //! Compute distance between patches
  computeDistance(distance, io_params, weightTable, i_imBasic, p_imSize);

  //! Get similar patches
  const unsigned int nSimP = getSimilarPatches(o_index, distance, io_params, io_params.tau);

  //! Save similar patches into 3D groups
  for (unsigned int c = 0; c < chnls; c++) {
    for (unsigned int p = 0; p < sP; p++) {
      for (unsigned int q = 0; q < sP; q++) {
        float* oGn       = &o_group3dNoisy[c * sP * sP * nSimP + (p * sP + q) * nSimP];
        float* oGb       = &o_group3dBasic[c * sP * sP * nSimP + (p * sP + q) * nSimP];
        const float* iIn = &i_imNoisy[c * wh + p * w + q];
        const float* iIb = &i_imBasic[c * wh + p * w + q];

        for (unsigned int n = 0; n < nSimP; n++) {
          oGn[n] = iIn[o_index[n]];
          oGb[n] = iIb[o_index[n]];
        }
      }
    }
  }

  //! Save the current number of similar patches
  io_params.nSimP = nSimP;
}


//! Compute the Bayes estimation for the first step.
void computeBayesEstimateStep1(
  std::vector<std::vector<float> > &io_group3d,
  matParams &i_mat,
  nlbParams &p_params){

  //! Initializations
  const unsigned int chnls = io_group3d.size();
  const unsigned int sP2   = p_params.sizePatch * p_params.sizePatch;
  const unsigned int nSimP = p_params.nSimP;

  //! For each 3D group
  for (unsigned int c = 0; c < chnls; c++) {

    float* iGc = &io_group3d[c][0];
    float max  = iGc[0];
    float min  = iGc[0];
    for (unsigned int k = 0; k < sP2 * nSimP; k++) {
      max = std::max(max, iGc[k]);
      min = std::min(min, iGc[k]);
    }
    min -= p_params.sigma[c];
    max += p_params.sigma[c];

    //! Compute the average intensity of the reference block
    const unsigned int intRef = getMean(io_group3d[c], 0, sP2 * nSimP, p_params.rangeOffset, p_params.rangeMax);

    //! Get the right covariance matrix of the noise
    const float* iCn = &p_params.noiseCovMatrix[intRef][c * sP2 * sP2];

    //! Center data around the barycenter
    centerData(io_group3d[c], i_mat.group3d, i_mat.barycenter, nSimP, sP2);

    //! Compute the covariance matrix of the set of similar patches
    covarianceMatrix(io_group3d[c], i_mat.covMat, nSimP, sP2);
    double* iCp = &i_mat.covMat[0];
    double* iCt = &i_mat.covMatTmp[0];
    for (unsigned int k = 0; k < sP2 * sP2; k++) {
      iCt[k] = iCp[k] - (double) iCn[k];
    }

    //! Compute the difference between the trace of covariance matrices
    float noiseTrace = 0.f;
    float patchTrace = 0.f;

    for (unsigned int k = 0; k < sP2; k++) {
      noiseTrace += iCn[k * sP2 + k];
      patchTrace += iCp[k * sP2 + k];
      iCp[k * sP2 + k] = std::max((double) iCn[k * sP2 + k], iCp[k * sP2 + k]);
    }

    if (p_params.useArea) {
      p_params.patchTrace[c] = patchTrace / ((float) sP2);
      p_params.noiseTrace[c] = noiseTrace / ((float) sP2);
    }

    //! Bayes' Filtering
    if (inverseMatrix(i_mat.covMat, sP2) == EXIT_SUCCESS) {
      productMatrix(i_mat.covMat, i_mat.covMatTmp, i_mat.covMat, sP2, sP2, sP2);
      productMatrix(i_mat.group3d, i_mat.covMat   , i_mat.group3d, sP2, sP2, nSimP);
    }

    //! Add barycenter
    for (unsigned int j = 0; j < sP2; j++) {
      float*  oG = &io_group3d[c][j * nSimP];
      const double* iG = &i_mat.group3d[j * nSimP];
      const float bary = i_mat.barycenter[j];

      for (unsigned int i = 0; i < nSimP; i++) {
        oG[i] = std::min(max, std::max(bary + (float) iG[i], min));
      }
    }
  }
}


//! Compute the Bayes estimation for the second step.
void computeBayesEstimateStep2(
  std::vector<float> &i_group3dNoisy,
  std::vector<float> &io_group3dBasic,
  matParams &i_mat,
  const ImageSize &p_imSize,
  nlbParams &p_params){

  //! Initializations
  const unsigned int sP2   = p_params.sizePatch * p_params.sizePatch;
  const unsigned int sPC   = sP2 * p_imSize.nChannels;
  const unsigned int nSimP = p_params.nSimP;

  vector<float> min(p_imSize.nChannels);
  vector<float> max(p_imSize.nChannels);
  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
    const float* iGb = &io_group3dBasic[c * sP2 * nSimP];

    for (unsigned int k = 0; k < sP2 * nSimP; k++) {
      min[c] = std::min(min[c], iGb[k]);
      max[c] = std::max(max[c], iGb[k]);
    }
    min[c] -= p_params.sigma[c];
    max[c] += p_params.sigma[c];
  }

  //! covariance matrix of the "oracle"
  for (unsigned int k = 0; k < sPC; k++) {
    float* iGn = &i_group3dNoisy [k * nSimP];
    float* iGb = &io_group3dBasic[k * nSimP];
    double* iG = &i_mat.group3d  [k * nSimP];
    double val = 0.f;
    for (unsigned int n = 0; n < nSimP; n++) {
      val += iGn[n];
    }
    val /= (double) (nSimP);
    i_mat.barycenter[k] = val;
    for (unsigned int n = 0; n < nSimP; n++) {
      iG[n] = iGb[n] - val;
    }
  }

  covarianceMatrix(i_mat.group3d, i_mat.covMat, nSimP, sPC);

  //! Construct the current covMatNoise + covMatPatch
  i_mat.covMatTmp = i_mat.covMat;
  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {

    //! Compute the average intensity of the reference block
    const unsigned int intRef = getMean(i_group3dNoisy, c * sP2 * nSimP,
                                        sP2 * nSimP, p_params.rangeOffset, p_params.rangeMax);

    //! Get the right covariance matrix of the noise
    const float* iCn = &p_params.noiseCovMatrix[intRef][c * sP2 * sP2];

    //! Constructing Cp + Cn
    double* iCp = &i_mat.covMat[c * sP2 * (sPC + 1)];
    for (unsigned int i = 0; i < sP2; i++) {
      for (unsigned int j = 0; j < sP2; j++) {
        iCp[i * sPC + j] += (double) iCn[i * sP2 + j];
      }
    }
  }

  //! centering the current noisy data
  for (unsigned int k = 0; k < sPC; k++) {
    const float bary = i_mat.barycenter[k];
    double* iG       = &i_mat.group3d [k * nSimP];
    const float* iGn = &i_group3dNoisy[k * nSimP];

    for (unsigned int n = 0; n < nSimP; n++) {
      iG[n] = iGn[n] - bary;
    }
  }

  //! Bayes' estimate
  if (inverseMatrix(i_mat.covMat, sPC) == EXIT_SUCCESS) {
    productMatrix(i_mat.covMat, i_mat.covMatTmp, i_mat.covMat, sPC, sPC, sPC);
    productMatrix(i_mat.group3d, i_mat.covMat   , i_mat.group3d, sPC, sPC, nSimP);
  }

  //! Add barycenter
  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
    const float maxi = max[c];
    const float mini = min[c];

    for (unsigned int k = 0; k < sP2; k++) {
      const float bary = i_mat.barycenter[c * sP2 + k];
      float* oGb  = &io_group3dBasic[(c * sP2 + k) * nSimP];
      double* iGb = &i_mat.group3d  [(c * sP2 + k) * nSimP];

      for (unsigned int i = 0; i < nSimP; i++) {
        oGb[i] = std::max(std::min(bary + (float) iGb[i], maxi), mini);
      }
    }
  }
}


//! Aggregate estimates of all similar patches contained in the 3D group durint the first step.
void computeAggregationStep1(
  std::vector<float> &io_im,
  std::vector<float> &io_weight,
  std::vector<bool> &io_mask,
  std::vector<std::vector<float> > const& i_group3d,
  std::vector<unsigned int> const& i_index,
  const ImageSize &p_imSize,
  const nlbParams &p_params){

  //! For convenience
  const unsigned int w     = p_imSize.width;
  const unsigned int sP    = p_params.sizePatch;
  const unsigned int nSimP = p_params.nSimP;

  //! Aggregate estimates
  for (unsigned int n = 0; n < nSimP; n++) {
    const unsigned int ind = i_index[n];

    for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
      const float* iG = &i_group3d[c][n];
      float* iI       = &io_im    [ind + c * p_imSize.wh];
      float* iW       = &io_weight[ind + c * p_imSize.wh];

      //! accumulates each estimator of the pixels, with the same weight (=1 currently)
      for (unsigned int p = 0; p < sP; p++) {
        for (unsigned int q = 0; q < sP; q++) {
          iI[p * w + q] += iG[(p * sP + q) * nSimP];
          iW[p * w + q]++;
        }
      }
    }

    //! Use Paste Trick
    io_mask[ind    ] = false;
    io_mask[ind - 1] = false;
    io_mask[ind + 1] = false;
    io_mask[ind - w] = false;
    io_mask[ind + w] = false;
  }
}


//! Aggregate estimates of all similar patches contained in the 3D group durint the second step.
void computeAggregationStep2(
  std::vector<float> &io_im,
  std::vector<float> &io_weight,
  std::vector<bool> &io_mask,
  std::vector<float> const& i_group3d,
  std::vector<unsigned int> const& i_index,
  const ImageSize &p_imSize,
  const nlbParams &p_params){

  //! For convenience
  const unsigned int w     = p_imSize.width;
  const unsigned int sP    = p_params.sizePatch;
  const unsigned int nSimP = p_params.nSimP;

  //! Aggregate estimates
  for (unsigned int n = 0; n < nSimP; n++) {
    const unsigned int ind = i_index[n];

    for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
      const float* iG = &i_group3d[c * sP * sP * nSimP + n];
      float* iI       = &io_im    [ind + c * p_imSize.wh];
      float* iW       = &io_weight[ind + c * p_imSize.wh];

      //! Save Values
      for (unsigned int p = 0; p < sP; p++) {
        for (unsigned int q = 0; q < sP; q++) {
          iI[p * w + q] += iG[(p * sP + q) * nSimP];
          iW[p * w + q]++;
        }
      }
    }

    //! Apply Paste Trick (test because of the \tau parameter)
    if (n < p_params.nSimilarPatches || !p_params.isFinestScale) {
      io_mask[ind    ] = false;
      io_mask[ind - 1] = false;
      io_mask[ind + 1] = false;
      io_mask[ind - w] = false;
      io_mask[ind + w] = false;
    }
  }
}


//! Compute the final weighted aggregation.
void computeWeightedAggregation(
  std::vector<float> const& i_imNoisy,
  std::vector<float> &io_imOutput,
  std::vector<float> const& i_weight,
  const nlbParams &p_params,
  const ImageSize &p_imSize){

  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
    const float* iW = &i_weight   [c * p_imSize.wh];
    const float* iN = &i_imNoisy  [c * p_imSize.wh];
    float* iO       = &io_imOutput[c * p_imSize.wh];

    for (unsigned int k = 0; k < p_imSize.wh; k++) {

      //! To avoid weighting problem (particularly near boundaries of the image)
      if (iW[k] > 0.f) {
        iO[k] /= iW[k];
      }
      else {
        iO[k] = iN[k];
      }
    }
  }
}


//! Compute a soft thresholding on homogeneous area for the first step.
void computeSoftHomogeneousAreaStep1(
  std::vector<std::vector<float> > &io_group3d,
  std::vector<unsigned int> const& i_index,
  const nlbParams &p_params,
  ImageSize const& p_imSize){
  //! Parameters
  const unsigned int sP2   = p_params.sizePatch * p_params.sizePatch;
  const unsigned int nSimP = p_params.nSimP;

  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
    const float mean       = p_params.meanPatches[c];
    const float patchTrace = p_params.patchTrace [c];
    const float noiseTrace = p_params.noiseTrace [c];
    const float kL         = p_params.kLow  * noiseTrace;
    const float kH         = p_params.kHigh * noiseTrace;

    //! Check the difference to the barycenter
    float diff = 0.f;
    const float* iB = &p_params.barycenter[c * sP2];
    for (unsigned int k = 0; k < sP2; k++) {
      diff += (iB[k] - mean) * (iB[k] - mean);
    }
    diff /= (float) sP2;
    //const bool doIt = true;
    const bool doIt = (diff < p_params.gamma * noiseTrace);

    if (patchTrace <= kL && doIt) {
      float* iG = &io_group3d[c][0];

      for (unsigned int k = 0; k < nSimP * sP2; k++) {
        iG[k] = mean;
      }
    }

    else if (patchTrace > kL && patchTrace <= kH && doIt) {
      const float alpha = (patchTrace - kL) / (kH - kL);
      float* iG         = &io_group3d[c][0];

      for (unsigned int k = 0; k < nSimP * sP2; k++) {
        iG[k] = alpha * iG[k] + (1.f - alpha) * mean;
      }
    }

    else {
      //! Do nothing
    }
  }
}


//! Compute the noise covariance matrix.
int computeNoiseCovarianceMatrix(
  std::vector<float> const& i_imNoisy,
  const ImageSize &p_imSize,
  nlbParams &io_params,
  const MsdParams &p_msdParams){

  //! Color space transform (if needed)
  vector<float> imNoisy = i_imNoisy;
  if (io_params.isFirstStep) {
    transformColorSpace(imNoisy, p_imSize, true);
  }

  //! Get the intensity range
  if (io_params.isFirstStep) {
    float min = imNoisy[0];
    float max = imNoisy[0];
    for (unsigned int k = 0; k < p_imSize.whc; k++) {
      const float value = imNoisy[k];
      min = std::min(min, value);
      max = std::max(max, value);
    }
    io_params.rangeOffset = (unsigned int) fabs(std::min(0.f, floor(min)));
    io_params.rangeMax    = (unsigned int) ceil(max);
  }
  else {
    io_params.rangeOffset = 0;
    io_params.rangeMax    = 255;
  }

  //! To compute the execution time
  const clock_t uptimeStart = clock() / (CLOCKS_PER_SEC / 1000);
  struct timespec start, finish;
  clock_gettime(CLOCK_MONOTONIC, &start);

  //! Estimate covariance matrices for each channel and intensity
  int res = runNoiseEstimate(imNoisy, p_imSize, io_params, p_msdParams);

  //! Detect the amount of time
  const clock_t uptimeEnd = clock() / (CLOCKS_PER_SEC / 1000);
  clock_gettime(CLOCK_MONOTONIC, &finish);
  cout << "time elapsed (ms) = " << (uptimeEnd - uptimeStart) << endl;
  double elapsed = (finish.tv_sec - start.tv_sec) * 1000;
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000.0;
  cout << "real elapsed (ms) = " << elapsed / 1.0 << endl;

  return res;
}


//! Get the average sigma value for the current reference patch.
void getSigma(
  const unsigned int p_ij,
  std::vector<float> const& i_im,
  const ImageSize &p_imSize,
  nlbParams &io_params){

  //! For convenience
  const unsigned int sP  = io_params.sizePatch;
  const unsigned int sP2 = sP * sP;
  const unsigned int max = io_params.rangeMax;
  const unsigned int off = io_params.rangeOffset;

  //! Get the average intensity
  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {

    float mean = 0.f;
    for (unsigned int p = 0; p < sP; p++) {
      const float* iI = &i_im[c * p_imSize.wh + p * p_imSize.width];
      for (unsigned int q = 0; q < sP; q++) {
        mean += iI[q];
      }
    }
    mean /= (float) sP2;
    mean += off;
    const unsigned int intRef = (mean < 0.f? 0 : (mean > (float) (max + off) ? (max + off) : (unsigned int) mean));

    //! Update sigma value
    const float* iCn = &io_params.noiseCovMatrix[intRef][c * sP2 * sP2];
    float tr = 0.f;
    for (unsigned int k = 0; k < sP2; k++) {
      tr += iCn[k * sP2 + k];
    }
    io_params.sigma[c] = sqrtf(std::max(tr, 0.f) / (float) sP2);
  }
}


//! Compute the distance between patches of the search window and the
//! reference patch.
void computeDistance(
  std::vector<std::pair<float, unsigned int> > &o_distance,
  const nlbParams &p_params,
  std::vector<float> const& p_weightTable,
  std::vector<float> const& i_im,
  const ImageSize &p_imSize){

  //! For convenience
  const unsigned int sW  = p_params.sizeSearchWindow;
  const unsigned int sP  = p_params.sizePatch;
  const unsigned int w   = p_imSize.width;
  const unsigned int wh  = p_imSize.wh;
  const unsigned int ij  = p_params.indRef;
  const unsigned int ind = ij - (sW - 1) * (w + 1) / 2;

  for (unsigned int i = 0; i < sW; i++) {
    for (unsigned int j = 0; j < sW; j++) {
      const unsigned int k = i * w + j + ind;

      //! center is at distance 0
      if(k == ij) {
        o_distance[i * sW + j] = std::make_pair(0.f, k);
        continue;
      }

      float diff = 0.f;

      for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
        const float weight = p_weightTable[c];
        const float* iIr = &i_im[c * wh + ij];
        const float* iIg = &i_im[c * wh + k];

        for (unsigned int p = 0; p < sP; p++) {
          for (unsigned int q = 0; q < sP; q++) {
            const float tmpValue = iIr[p * w + q] - iIg[p * w + q];
            diff += weight * tmpValue * tmpValue;
          }
        }
      }

      o_distance[i * sW + j] = std::make_pair(diff, k);
    }
  }
}


//! Get similar patches according to their distance.
unsigned int getSimilarPatches(
  std::vector<unsigned int> &o_index,
  std::vector<std::pair<float, unsigned int> > &io_distance,
  const nlbParams &p_params,
  const float p_threshold){

  //! Keep only the N2 best similar patches
  if (p_params.useRandomPatches) {
    sort(io_distance.begin(), io_distance.end(), comparaisonFirst);
  }
  else {
    partial_sort(io_distance.begin(), io_distance.begin() + p_params.nSimilarPatches, io_distance.end(), comparaisonFirst);
  }
  const float threshold = std::max(p_threshold, io_distance[p_params.nSimilarPatches - 1].first);

  //! Check the total number of eligible patches
  unsigned int nSimP = 0;
  if (p_params.useRandomPatches && threshold == p_threshold) {
    unsigned int nbTotal = 0;
    while (nbTotal < io_distance.size() && io_distance[nbTotal].first <= threshold) {
      nbTotal++;
    }

    //! Get the wanted number of patches, randomly
    vector<unsigned int> index;
    getRandPerm(index, nbTotal - 1);
    nSimP = std::min(p_params.maxSimilarPatches, nbTotal);
    o_index[0] = p_params.indRef;
    for (unsigned int n = 0; n < nSimP - 1; n++) {
      o_index[n + 1] = io_distance[index[n] + 1].second;
    }
  }
  else {
    //! Register position of similar patches
    //! automatically register the first p_params.nSimilarPatches, then perform the loop
    nSimP = p_params.nSimilarPatches;
    for (unsigned int n = 0; n < nSimP; n++) {
      o_index[n] = io_distance[n].second;
    }
  }

  return nSimP;
}






