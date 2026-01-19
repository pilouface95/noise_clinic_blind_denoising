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


/**
 * @file main.cpp
 * @brief main function.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 * @author Miguel Colom <miguel.colom@cmla.ens-cachan.fr>
 **/

#include <iostream>
#include <stdlib.h>
#include <vector>

#include "Utilities/LibImages.h"
#include "MultiScalesDenoising/MultiScalesDenoising.h"

using namespace std;

int main(
  int argc,
  char **argv){

  //! Check if there is the right number of input parameters
  if (argc < 6) {
    cout << "Usage: NoiseClinic imNoisy imDenoised imDiff nbScales factor verbose" << endl;
    return EXIT_FAILURE;
  }
  // /home/ibal/Images/input/Noisy/Lena_BW.png /home/ibal/Images/output/NoiseClinic/Denoised.png /home/ibal/Images/output/NoiseClinic/Diff.png 1 1 1
  // /home/ibal/Bureau/single.png /home/ibal/Images/output/NoiseClinic/Denoised.png /home/ibal/Images/output/NoiseClinic/Diff.png 3 1 1

  //! Declarations
  vector<float> imNoisy, imFinal;
  vector<vector<float> > imDiffs;
  ImageSize imSize;
  MsdParams params;
  string output             = argv[2];
  output.erase(output.end() - 4, output.end());
  params.path               = (char*) output.c_str();
  params.nbScales           = atoi(argv[4]);
  params.stdFactor          = atof(argv[5]);
  params.verbose            = (bool) atoi(argv[6]);

  //! Read Image
  if (loadImage(argv[1], imNoisy, imSize, params.verbose) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Denoise the image
  if (runMultiScalesDenoising(imNoisy, imFinal, imDiffs, imSize, params) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Save differences for each scale
  if (saveNoiseRemovalImages(imDiffs, imSize, argv[2]) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Save the final difference and compute the RMSE over this image
  if (saveDiff(imFinal, imNoisy, imSize, argv[3]) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Save the final denoised image
  return saveImage(argv[2], imFinal, imSize);
}
