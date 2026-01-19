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

#ifndef NL_BAYES_H_INCLUDED
#define NL_BAYES_H_INCLUDED

#include <vector>

#include "../Utilities/LibImages.h"
#include "../MultiScalesDenoising/MultiScalesDenoising.h"


/**
 * @brief Structures of parameters dedicated to NL-Bayes process.
 *
 * @param sizePatch : Size of side of the patches (the patch will have sizePatch^2 pixels);
 * @param nSimilarPatches : Number of similar patches wanted;
 * @param sizeSearchWindows : Size of the search window around the reference patch;
 * @param tau: parameter used to determine similar patches;
 * @param isFirstStep: true if the first step of the algorithm is wanted;
 * @param kLow : low threshold to decide for homogeneous area;
 * @param kHigh : high threshold to decide for homogeneous area;
 * @param useArea : if true, use the homogeneous area trick;
 * @param noiseCovMatrix : estimated covariance matrix of the noise, for each intensity;
 * @param meanPatches : will contain mean of all patches of the current 3d group;
 * @param isFinestScale : true if the current scale is the finest scale;
 * @param diffTrace : contains Tr(C_P) - Tr(C_N) for each channel;
 * @param noiseTrace : contains Tr(C_N) for each channel;
 * @param rangeOffset : offset to apply to the mean intensity of a patch to get the index
 *       of intensity for the noiseCovMatrix;
 * @param rangeMax : maximum value of the intensity range for the noiseCovMatrix;
 * @param sigma : will contain for all channels and the current intensity an estimation of the
 *                standard deviation of the noise, as the normalized trace of the covariance
 *                matrix of the noise;
 * @param indRef : current index of the reference patch;
 * @param barycenter : will contain the barycenter of the 3D group;
 * @param gamma : parameter to determine if a 3D group belongs to an homogeneous area;
 * @param useRandomPatches : if true, select a random number of patches from the set of valid patches;
 * @param maxSimilarPatches : maximal number of similar patches allowed;
 * @param nSimP : current number of patches in the current 3D group;
 * @param verbose : if true, print some informations.
 **/
struct nlbParams
{
  unsigned int sizePatch;
  unsigned int nSimilarPatches;
  unsigned int sizeSearchWindow;
  float tau;
  bool isFirstStep;
  float kLow;
  float kHigh;
  bool useArea;
  std::vector<std::vector<float> > noiseCovMatrix;
  std::vector<float> meanPatches;
  bool isFinestScale;
  std::vector<float> patchTrace;
  std::vector<float> noiseTrace;
  unsigned int rangeOffset;
  unsigned int rangeMax;
  std::vector<float> sigma;
  unsigned int indRef;
  std::vector<float> barycenter;
  float gamma;
  bool useRandomPatches;
  unsigned int maxSimilarPatches;
  unsigned int nSimP;
  bool verbose;
};


/**
 * @brief Structure containing usefull matrices for the Bayes estimations.
 *
 * @param barycenter: allocated memory. Used to contain the barycenter of io_group3dBasic;
 * @param covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 * @param covMatTmp: allocated memory. Used to process the Bayes estimate.
 **/
struct matParams
{
  std::vector<float> barycenter;
  std::vector<double> covMat;
  std::vector<double> covMatTmp;
  std::vector<double> group3d;
};


/**
 * @brief Initialize Parameters of the NL-Bayes algorithm.
 *
 * @param o_paramStep1 : will contain the nlbParams for the first step of the algorithm;
 * @param o_paramStep2 : will contain the nlbParams for the second step of the algorithm;
 * @param io_msdParams : see MsdParams. The noiseCurve will be updated;
 * @param i_imNoisy : noisy image;
 * @param p_imSize: size of the image.
 *
 * @return EXIT_FAILURE in case of problems, EXIT_SUCCESS otherwise.
 **/
int initializeNlbParameters(
  nlbParams &o_paramStep1,
  nlbParams &o_paramStep2,
  MsdParams &io_msdParams,
  std::vector<float> const& i_imNoisy,
  const ImageSize &p_imSize);


/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * @param i_imNoisy: contains the noisy image;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_imSize: size of the image;
 * @param io_msdParams : see MsdParams. NoiseCurve will be updated in the coarsest scale.
 *
 * @return EXIT_FAILURE if something wrong happens during the whole process.
 **/
int runNlBayes(
  std::vector<float> const& i_imNoisy,
  std::vector<float> &o_imFinal,
  const ImageSize &p_imSize,
  MsdParams &io_msdParams);


/**
 * @brief Generic step of the NL-Bayes denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy image;
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_imSize: size of i_imNoisy;
 * @param p_params: see nlbParams.
 *
 * @return none.
 **/
void processNlBayes(
  std::vector<float> const& i_imNoisy,
  std::vector<float> &io_imBasic,
  std::vector<float> &o_imFinal,
  const ImageSize &p_imSize,
  nlbParams &p_params);


/**
 * @brief Estimate the best similar patches to a reference one.
 *
 * @param i_im: contains the noisy image on which distances are processed;
 * @param o_group3d: will contain values of similar patches;
 * @param o_index: will contain index of similar patches;
 * @param p_imSize: size of the image;
 * @param p_params: see nlbParams for more explanation. Will contain the mean,
 *                  barycenter and the number of similar patches of the 3D group.
 *
 * @return none.
 **/
void estimateSimilarPatchesStep1(
  std::vector<float> const& i_im,
  std::vector<std::vector<float> > &o_group3d,
  std::vector<unsigned int> &o_index,
  const ImageSize &p_imSize,
  nlbParams &io_params);


/**
 * @brief Keep from all near patches the similar ones to the reference patch for the second step.
 *
 * @param i_imNoisy: contains the original noisy image;
 * @param i_imBasic: contains the basic estimation;
 * @param o_group3dNoisy: will contain similar patches for all channels of i_imNoisy;
 * @param o_group3dBasic: will contain similar patches for all channels of i_imBasic;
 * @param o_index: will contain index of similar patches;
 * @param p_imSize: size of images;
 * @param p_params: see nlbParams for more explanations. Will contain the number of
 *                  similar patches.
 *
 * @return none.
 **/
void estimateSimilarPatchesStep2(
  std::vector<float> const& i_imNoisy,
  std::vector<float> const& i_imBasic,
  std::vector<float> &o_group3dNoisy,
  std::vector<float> &o_group3dBasic,
  std::vector<unsigned int> &o_index,
  const ImageSize &p_imSize,
  nlbParams &io_params);


/**
 * @brief Compute the Bayes estimation for the first step.
 *
 * @param io_group3d: contains all similar patches. Will contain estimates for all similar patches;
 * @param i_mat: see matParams;
 * @param p_params: see nlbParams.
 *
 * @return none.
 **/
void computeBayesEstimateStep1(
  std::vector<std::vector<float> > &io_group3d,
  matParams &i_mat,
  nlbParams &p_params);


/**
 * @brief Compute the Bayes estimation for the second step.
 *
 * @param i_group3dNoisy: contains all similar patches in the noisy image;
 * @param io_group3dBasic: contains all similar patches in the basic image. Will contain estimates
 *       for all similar patches;
 * @param i_mat: see matParams
 * @param p_imSize: size of the image;
 * @param p_params: see nlbParams.
 *
 * @return none.
 **/
void computeBayesEstimateStep2(
  std::vector<float> &i_group3dNoisy,
  std::vector<float> &io_group3dBasic,
  matParams &i_mat,
  const ImageSize &p_imSize,
  nlbParams &p_params);


/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group durint the first step.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_imSize: size of io_im;
 * @param p_params: see nlbParams for more explanation.
 *
 * @return none.
 **/
void computeAggregationStep1(
  std::vector<float> &io_im,
  std::vector<float> &io_weight,
  std::vector<bool> &io_mask,
  std::vector<std::vector<float> > const& i_group3d,
  std::vector<unsigned int> const& i_index,
  const ImageSize &p_imSize,
  const nlbParams &p_params);


/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group during the second step.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_imSize: size of io_im;
 * @param p_params: see nlbParams.
 *
 * @return none.
 **/
void computeAggregationStep2(
  std::vector<float> &io_im,
  std::vector<float> &io_weight,
  std::vector<bool> &io_mask,
  std::vector<float> const& i_group3d,
  std::vector<unsigned int> const& i_index,
  const ImageSize &p_imSize,
  const nlbParams &p_params);


/**
 * @brief Compute the final weighted aggregation.
 *
 * @param i_imNoisy : noisy image;
 * @param io_imOutput: will contain the aggregated estimation;
 * @param i_weight associated weight for each estimate of pixels;
 * @param p_params : see nlbParams;
 * @param p_imSize : size of images.
 *
 * @return none.
 **/
void computeWeightedAggregation(
  std::vector<float> const& i_imNoisy,
  std::vector<float> &io_imOutput,
  std::vector<float> const& i_weight,
  const nlbParams &p_params,
  const ImageSize &p_imSize);


/**
 * @brief Compute a soft thresholding on homogeneous area for the first step.
 *
 * @param io_group3d : 3D group of similar patches;
 * @param i_index : index of position of the patches contained in i_group3d;
 * @param p_params : see nlbParams;
 * @param p_imSize : size of the image.
 *
 * @return none.
 **/
void computeSoftHomogeneousAreaStep1(
  std::vector<std::vector<float> > &io_group3d,
  std::vector<unsigned int> const& i_index,
  const nlbParams &p_params,
  ImageSize const& p_imSize);


/**
 * @brief Compute the noise covariance matrix.
 *
 * @param i_imNoisy : input noisy image;
 * @param p_imSize : size of the input image;
 * @param io_params : see nlbParams;
 * @param p_msdParams : see MsdParams.
 *
 * @return EXIT_FAILURE in case of problems, EXIT_SUCCESS otherwise.
 **/
int computeNoiseCovarianceMatrix(
  std::vector<float> const& i_imNoisy,
  const ImageSize &p_imSize,
  nlbParams &io_params,
  const MsdParams &p_msdParams);


/**
 * @brief Get the average sigma value for the current reference patch.
 *
 * @param p_ij : position of the reference patch;
 * @param p_currChnl : current channel;
 * @param i_im : input image;
 * @param p_imSize : size of i_im;
 * @param io_params : see nlbParams. Will update the sigma value.
 *
 * @return none.
 **/
void getSigma(
  const unsigned int p_ij,
  std::vector<float> const& i_im,
  const ImageSize &p_imSize,
  nlbParams &io_params);


/**
 * @brief Compute the distance between patches of the search window and the
 *        reference patch. If needed, substract a planar regression before.
 *
 * @param o_distance : will contain distance and position of patches;
 * @param p_params : see nlbParams;
 * @param p_weightTable : weight of channels;
 * @param i_im : image on which distance is computed;
 * @param i_plan : plane to substract before computing the distance;
 * @param p_imSize : size of the image,
 * @param p_usePlan : if true, substract the plane.
 *
 * @return none.
 **/
void computeDistance(
  std::vector<std::pair<float, unsigned int> > &o_distance,
  const nlbParams &p_params,
  std::vector<float> const& p_weightTable,
  std::vector<float> const& i_im,
  const ImageSize &p_imSize);


/**
 * @brief Get similar patches according to their distance.
 *
 * @param o_index : will contain position of similar patches;
 * @param io_distance : contains position and distance of all patches;
 * @param p_params : see nlbParams,
 * @param p_threshold : threshold of maximal similarity.
 *
 * @return the number of similar patches.
 **/
unsigned int getSimilarPatches(
  std::vector<unsigned int> &o_index,
  std::vector<std::pair<float, unsigned int> > &io_distance,
  const nlbParams &p_params,
  const float p_threshold);


#endif // NL_BAYES_H_INCLUDED
