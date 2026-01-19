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
#ifndef LIB_IMAGES_H_INCLUDED
#define LIB_IMAGES_H_INCLUDED

#include <vector>


/**
 * @brief Class containing size informations of an image.
 *
 * @param width     : width of the image;
 * @param height    : height of the image;
 * @param nChannels : number of channels in the image;
 * @param nBits     : number of bits/channel
 * @param wh        : equal to width * height. Provided for convenience;
 * @param whc       : equal to width * height * nChannels. Provided for convenience.
 **/
class ImageSize {
  public:
    unsigned width;
    unsigned height;
    unsigned nChannels;
    unsigned nBits;
    unsigned wh;
    unsigned whc;

	ImageSize();
	ImageSize(const unsigned width, const unsigned height,
		  const unsigned nChannels, const unsigned nBits,
		  const unsigned wh, const unsigned whc);
};


/**
 * @brief Load image, check the number of channels.
 *
 * @param p_name : name of the image to read;
 * @param o_im : vector which will contain the image : R, G and B concatenated;
 * @param o_imSize : will contain the size of the image;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_SUCCESS if the image has been loaded, EXIT_FAILURE otherwise.
 **/
int loadImage(
  char* p_name,
  std::vector<float> &o_im,
  ImageSize &o_imSize,
  const bool p_verbose);


/**
 * @brief write image.
 *
 * @param p_name : path+name+extension of the image;
 * @param i_im : vector which contains the image;
 * @param p_imSize : size of the image;
 *
 * @return EXIT_SUCCESS if the image has been saved, EXIT_FAILURE otherwise
 **/
int saveImage(
  char* p_name,
  std::vector<float> const& i_im,
  const ImageSize &p_imSize);


/**
 * @brief Save a set of difference (noise removed) images.
 *
 * @param i_imDiffs : difference images to save. Need an offset to be visualized.
 * @param p_imSize : size of images;
 * @param p_path : name of the generic image (the extension must be removed).
 *
 * @return EXIT_FAILURE in case of problems.
 **/
int saveNoiseRemovalImages(
  std::vector<std::vector<float> > &io_imDiffs,
  const ImageSize &p_imSize,
  const char* p_path);


/**
 * @brief Transform the color space of an image, from RGB to YUV, or vice-versa.
 *
 * @param io_im: image on which the transform will be applied;
 * @param p_imSize: size of io_im;
 * @param p_isForward: if true, go from RGB to YUV, otherwise go from YUV to RGB.
 *
 * @return none.
 **/
void transformColorSpace(
  std::vector<float> &io_im,
  const ImageSize &p_imSize,
  const bool p_isForward);


/**
* @brief Add boundaries by symetry around an image to increase its size until a multiple of p_N.
*
* @param i_im : image to symetrize;
* @param o_imSym : will contain i_im with symetrized boundaries;
* @param p_imSize : size of i_im;
* @param p_imSizeSym : will contain the size of o_imSym;
* @param p_N : parameter of the boundaries.
*
* @return EXIT_FAILURE in case of problems.
**/
int enhanceBoundaries(
  std::vector<float> const& i_im,
  std::vector<float> &o_imSym,
  const ImageSize &p_imSize,
  ImageSize &o_imSizeSym,
  const unsigned p_N);


/**
* @brief Remove boundaries added with enhanceBoundaries
*
* @param o_im : will contain the inner image;
* @param i_imSym : contains i_im with symetrized boundaries;
* @param p_imSize: size of o_im;
* @param p_imSizeSym : size of i_imSym.
*
* @return none.
**/
int removeBoundaries(
  std::vector<float> &o_im,
  std::vector<float> const& i_imSym,
  const ImageSize &p_imSize,
  const ImageSize &p_imSizeSym);


/**
* @brief Split an image into 4 sub-scale images, 4 times smaller by averaging a 2x2 square of pixels
*
* @param i_im : original image;
* @param o_subIm : will contain the 4 sub-scale images;
* @param i_imSize : size of i_im;
* @param o_imSize : will contain the size of sub-scales images.
*
* @return none.
**/
void meanSplit(
  std::vector<float> const& i_im,
  std::vector<std::vector<float> >& o_subIm,
  const ImageSize &i_imSize,
  ImageSize &o_imSize);


/**
* @brief Reconstruct an image splitted by te meanSplit function.
*
* @param i_subIm : contains the 4 sub-images;
* @param o_im : will contain the reconstructed image of i_subIm;
* @param i_imSize : size of sub-image;
* @param o_imSize : will contain the image of o_im.
*
* @return none.
**/
void meanBuilt(
  std::vector<std::vector<float> > const& i_subIm,
  std::vector<float> &o_im,
  const ImageSize &i_imSize,
  ImageSize &o_imSize);


/**
* @brief Compose a mosaic of sub-images, sticked by symetry.
*
* @param i_subIm : set of sub-images;
* @param o_im : will contain a mosaic of all sub-images;
* @param i_imSize : size of sub-images;
* @param o_imSize : will contain the size of the mosaic image.
*
* @return none.
**/
void constructMosaic(
  std::vector<std::vector<float> > const& i_subIm,
  std::vector<float> &o_im,
  const ImageSize &i_imSize,
  ImageSize &o_imSize);


/**
* @brief Upgrade in-place a mosaic: from a mosaic of 4^n sub-images, built a mosaic
*      with 4^(n-1) sub-images.
*
* @param io_mosaic : contain the 4^n sub-images mosaic. Will contain the 4^(n-1) mosaic;
* @param p_imSizes : contain for each scale the size of the sub-images contained in the mosaics;
* @param p_currentScale: current scale of the multi-scale algorithm;
* @param p_n : number of sub-images in the input mosaic (must be a power of 4).
*
* @return none.
**/
void upgradeMosaic(
  std::vector<float> &io_mosaic,
  std::vector<ImageSize> &p_imSizes,
  const unsigned int p_currentScale);


/**
* @brief Split a mosaic obtained with the constructMosaic function into sub-images.
*
* @param i_im : contains the mosaic;
* @param o_subIm : will contain the set of sub-images;
* @param i_imSize : size of i_im;
* @param o_imSize : will contain the size of sub-images;
* @param p_N : number of sub-images contained in the mosaic.
*
* @return none.
**/
void splitMosaic(
  std::vector<float> const& i_im,
  std::vector<std::vector<float> >& o_subIm,
  const ImageSize &i_imSize,
  ImageSize &o_imSize,
  const unsigned p_N);


/**
 * @brief Add boundaries around image to avoid boundary problems for the big loop over pixels.
 *
 * @param i_im : image to enwrap;
 * @param o_im : image which contains i_im in its center, and symetry in its borders;
 * @param i_imSize : size of i_im;
 * @param o_imSize : will contain the size of o_im;
 * @param p_border : size of the border.
 *
 * @return none.
 **/
void setBoundaries(
  std::vector<float> const& i_im,
  std::vector<float> &o_im,
  const ImageSize &i_imSize,
  ImageSize &o_imSize,
  const unsigned int p_border);


/**
 * @brief Divide an image into sub-parts of equal sizes with a boundaries around each.
 *
 * @param i_im : image to divide;
 * @param o_imSub : will contain sub-images;
 * @param i_imSize : size of i_im;
 * @param o_imSizeSub : will contain the size of o_imSub;
 * @param p_boundaries : boundaries to add around each sub-images;
 * @param p_nbParts : number of sub-images.
 *
 * @return EXIT_FAILURE in case of problems, EXIT_SUCCESS otherwise.
 **/
int subDivide(
  std::vector<float> const& i_im,
  std::vector<std::vector<float> > &o_imSub,
  const ImageSize &i_imSize,
  ImageSize &o_imSizeSub,
  const unsigned int p_boundaries,
  const unsigned int p_nbParts);


/**
 * @brief Build an image from a set of sub-parts obtained with subDivide.
 *
 * @param i_imSub : set of sub-images;
 * @param o_im : will contain the reconstructed image;
 * @param p_imSizeSub : size of sub-images;
 * @param p_imSize : size that o_im will have;
 * @param p_boundaries : boundaries that were added around each sub-images.
 *
 * @return EXIT_FAILURE in case of problems, EXIT_SUCCESS otherwise.
 **/
int subBuild(
  std::vector<std::vector<float> > const& i_imSub,
  std::vector<float> &o_im,
  const ImageSize &p_imSizeSub,
  const ImageSize &p_imSize,
  const unsigned int p_boundaries);


/**
 * @brief Save the difference image and compute the RMSE of it.
 *
 * @param i_imRef : reference image;
 * @param i_imFinal : denoised image;
 * @param p_imSize : size of images;
 * @param p_imName : name of the difference image to save.
 *
 * @return EXIT_FAILURE in case of problems, EXIT_SUCCESS otherwise.
 **/
int saveDiff(
  std::vector<float> const& i_imRef,
  std::vector<float> const& i_imFinal,
  const ImageSize &p_imSize,
  const char* p_imName);



#endif // LIB_IMAGES_H_INCLUDED
