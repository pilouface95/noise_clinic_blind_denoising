/*
 * Copyright (c) 2013, Miguel Colom <miguel.colom@cmla.ens-cachan.fr>
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

#ifndef RUNNOISEESTIMATE_H_INCLUDED
#define RUNNOISEESTIMATE_H_INCLUDED

#include <vector>
#include <map>

#include "../Utilities/LibImages.h"
#include "../LibNlBayes/NlBayes.h"
#include "../MultiScalesDenoising/MultiScalesDenoising.h"
#include "CHistogram.cpp"


//! Maps blocks_dct with their sparse distances
typedef std::map<float *, float> TMapSD;
typedef CHistogram<float*> histo_t;


/**
 * @brief Structure of parameters dedicated to the noise estimate algorithm.
 *
 * @param sizePatch : size of patches;
 * @param r1 : minimum distance between two patches to compute their MSE;
 * @param r2 : maximum distance between two patches to compute their MSE;
 * @param percentile : percentage of blocks which will be used for the estimation;
 * @param numBins : total number of bins;
 * @param numBlocks : total number of validated blocks in the image;
 * @param imSize : size of the image;
 * @param currentChannel : current channel to work on;
 * @param correctionFactor : empiric factor to correct the noise estimation.
 **/
struct neParams{
  float sizePatch;
  unsigned int r1;
  unsigned int r2;
  float percentile;
  unsigned int numBins;
  unsigned int numBlocks;
  unsigned int currentChannel;
  float correctionFactor;
};


/**
 * @brief Ponomarenko-diff noise estimation algorithm.
 *
 * @param i_im : input noisy image;
 * @param p_imSize : size of the image;
 * @param io_nlbParams : see nlbParams. The covariance matrices will be updated.
 * @param p_msdParams : see MsdParams.
 *
 * @return EXIT_FAILURE in case of problems, EXIT_SUCCESS otherwise.
 **/
int runNoiseEstimate(
  std::vector<float> const &i_im,
  const ImageSize& p_imSize,
  nlbParams &io_nlbParams,
  const MsdParams &p_msdParams);


/**
 * @brief Transpose a matrix.
 *
 * @param o_mat : will contain the transpose of i_mat;
 * @param i_mat : input matrix (must be squared of size p_sP x p_sP);
 * @param p_sP : side of i_mat.
 *
 * @return none.
 **/
void transposeMatrix(
  std::vector<float> &o_mat,
  std::vector<float> const& i_mat,
  const unsigned int p_sP);


/**
 * @brief Get the DCT normalized matrix.
 *
 * @param o_mat : will contain the matrix of the 2D normalized DCT;
 * @param p_rows : number of lines wanted for the 2D normalized DCT;
 * @param p_cols : number of columns wanted for the 2D normalized DCT;
 *
 * @return none.
 **/
void getNormalized2dDctMatrix(
  std::vector<float> &o_mat,
  const unsigned int p_rows,
  const unsigned int p_cols);


/**
 * @brief Estimate by linear interpolation one covariance matrix of the noise for
 *       each intensity between [0, max + offset].
 *
 * @param i_meanBins : for each bin and chnl, contains the mean;
 * @param i_matrixBins : for each bin and chnl, contains the estimated covariance
 *       matrix of the noise;
 * @param o_covMat : will contain for each channel and each intensity value the
 *       covariance matrix of the noise;
 * @param p_sizePatch : size of blocks;
 * @param p_numBins : total number of bins;
 * @param p_rangeOffset : offset (minimum value of the intensity);
 * @param p_rangeMax : maximum value of the intensity.
 *
 * @return none.
 **/
int estimatedStdMatrices(
  std::vector<std::vector<float> > const& i_meanBins,
  std::vector<std::vector<float> > const& i_matrixBins,
  std::vector<std::vector<float> > &o_covMat,
  const unsigned int p_sizePatch,
  const unsigned int p_numBins,
  const unsigned int p_rangeOffset,
  const unsigned int p_rangeMax);


/**
 * @brief From matrices M in frequency space, obtain the Dt M D matrix in the spatial space;
 *       i.e. the covariance matrix of the noise.
 *
 * @param i_vMeans : contains for each channel and each bin the mean;
 * @param i_vStds : contains for each channel, bin and frequency the estimated
 *       standard deviation of the noise (i.e. the M matrices);
 * @param o_meanBins : will contain for each channel and each bin the mean;
 * @param o_matrixBins : for each channel and bin will contain the DtMD matrix;
 * @param p_params : see neParams.
 *
 * @return none.
 **/
void getMatrixBins(
  std::vector<float> const& i_vMeans,
  std::vector<float> const& i_vStds,
  std::vector<std::vector<float> > &o_meanBins,
  std::vector<std::vector<float> > &o_matrixBins,
  neParams const& p_params,
  const ImageSize &p_imsize);


#endif // RUNNOISEESTIMATE_H_INCLUDED
