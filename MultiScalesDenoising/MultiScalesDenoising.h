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
#ifndef MULTISCALESDENOISING_H_INCLUDED
#define MULTISCALESDENOISING_H_INCLUDED


#include <vector>
#include "../Utilities/LibImages.h"


/**
 * @brief Parameters used in the multi-scales denoising algorithm.
 *
 * @param currentScale : current scale of the multi-scale algorithm;
 * @param nbScales: number of scales which will be applied to this algorithm;
 * @param verbose: if true, print some informations;
 * @param path : path to the output folder;
 * @param stdFactor : will modify the estimated standard deviation of the noise.
 **/
struct MsdParams
{
  unsigned int currentScale;
  unsigned int nbScales;
  bool verbose;
  char* path;
  float stdFactor;
};


/**
 * @brief Run a multi-scales denoising.
 *
 * @param i_imNoisy: contains the original noisy image;
 * @param o_imFinal: will contain the final denoised image for each scale;
 * @param o_imDiffs : will contain for each scale the noise that has been removed;
 * @param p_imSize: size of i_imNoisy;
 * @param p_params : see msdParams.
 *
 * @return EXIT_FAILURE if anything wrong happens, EXIT_SUCCESS otherwise.
 **/
int runMultiScalesDenoising(
  std::vector<float> const& i_imNoisy,
  std::vector<float> &o_imFinal,
  std::vector<std::vector<float> > &o_imDiffs,
  const ImageSize &p_imSize,
  const MsdParams &p_params);


/**
 * @brief Apply a mosaic multi-scales denoising.
 *
 * @param io_imNoisy: original noisy image to denoise;
 * @param o_imFinal: will contain denoised image;
 * @param o_imDiffs : will contain for each scale the noise that has been removed;
 * @param p_imSize: size of i_imNoisy;
 * @param p_params : see msdParams.
 *
 * @return EXIT_FAILURE in case of problems, EXIT_SUCCESS otherwise.
 **/
int applyMultiScaleDenoising(
  const std::vector<float> &i_imNoisy,
  std::vector<float> &o_imFinal,
  std::vector<std::vector<float> > &o_imDiffs,
  const ImageSize &p_imSize,
  MsdParams &io_params);


#endif // MULTISCALESDENOISING_H_INCLUDED
