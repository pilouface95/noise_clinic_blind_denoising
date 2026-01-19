/*
 *
 * This file is part of the Multiscale SFD Noise Estimation algorithm.
 *
 * Copyright(c) 2014 Miguel Colom.
 * miguel@mcolom.info
 *
 * This file may be licensed under the terms of of the
 * GNU General Public License Version 2 (the ``GPL'').
 *
 * Software distributed under the License is distributed
 * on an ``AS IS'' basis, WITHOUT WARRANTY OF ANY KIND, either
 * express or implied. See the GPL for the specific language
 * governing rights and limitations.
 *
 * You should have received a copy of the GPL along with this
 * program. If not, go to http://www.gnu.org/licenses/gpl.html
 * or write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 */

#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fftw3.h>
#include <iostream>
#include <malloc.h>
#include <xmmintrin.h>
#include <x86intrin.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "operations.cpp"
#include "RunNoiseEstimate.h"
#include "../Utilities/Utilities.h"
#include "CurveFiltering.h"
#include "../Utilities/CImage.h"


#define GET_BLOCK_PTR(x, y, MAT, w, Nx, Ny) (&MAT[(y*(Nx-w+1) + x)*w*w])

using namespace std;

// NDEBUG = NO Debug
#define NDEBUG

#ifndef NDEBUG
float *GET_BLOCK_PTR(int x, int y, float *MAT, int w, int Nx, int Ny) {
  assert(x < Nx-w+1);
  assert(y < Ny-w+1);
  return &MAT[(y*(Nx-w+1) + x)*w*w];
}
#else
#define GET_BLOCK_PTR(x, y, MAT, w, Nx, Ny) (&MAT[(y*(Nx-w+1) + x)*w*w])
#endif


typedef map<float *, float> TMapSD; // Maps blocks_dct with their sparse distances
typedef CHistogram<float*> histo_t;



/**
 * @brief Debug function to print the contents of a 2D block
 *
 * @param ptr : block pointer
 * @param w : block side
 *
 **/
void debug_print_block(float *ptr, int w) {
  float mean = 0;
  //
  for (int y = 0; y < w; y++) {
    for (int x = 0; x < w; x++) {
      printf("%f\t", ptr[y*w+x]);
      mean += ptr[y*w+x];
    }
    printf("\n");
  }
  mean /= w*w;
  printf("Mean: %f\n", mean);
}


/**
 * @brief Corrects hat_sigma to obtain tilde_sigma
 *
 * @param hat_sigma : biased STD of the noise
 * @param w : block side
 *
 * @return Corrected tilde_sigma
 *
 **/
float get_tilde_sigma(float hat_sigma, int w) {
  float K = (w == 8 ? 1.651094 : 1.791106);
  return K * hat_sigma;
}


/**
 * @brief Build a mask for valide pixel. If mask(i, j) = true, the pixels will not be used.
 *
 * @param i_im : noisy image;
 * @param o_mask : will contain the mask for all pixels in the image size;
 * @param p_imSize : size of the image;
 * @param p_sizePatch : size of a patch.
 *
 * @return number of valid blocks.
 *
 **/
unsigned buildMask(CImage &i_im, int *o_mask,
                    unsigned Nx, unsigned Ny, unsigned w,
                    unsigned num_channels) {
  unsigned count  = 0;

  for (unsigned ij = 0; ij < Nx*Ny; ij++) {
    const unsigned j = ij / Nx;
    const unsigned i = ij - j * Nx;

    //! Look if the pixel is not to close to boundaries of the image
    if (i < Nx - w + 1 && j < Ny - w + 1) {
      for (unsigned c = 0; c < num_channels; c++) {
        float *u = i_im.get_channel(c);

        //! Look if the square 2x2 of pixels is constant
        int invalid_pixel = (c == 0 ? 1 : o_mask[ij]);

        // Try to validate pixel
        if (fabs(u[ij] - u[ij + 1]) > 0.001f) {
                invalid_pixel = 0;
        }
        else
        if (fabs(u[ij + 1] - u[ij + Nx]) > 0.001f) {
                invalid_pixel = 0;
        }
        else
        if (fabs(u[ij + Nx] - u[ij + Nx + 1]) > 0.001f) {
                invalid_pixel = 0;
        }
        o_mask[ij] = invalid_pixel;
      }
    }
    else {
            o_mask[ij] = 1; // Not valid
    }

    if (o_mask[ij] == 0)
        count++;
  }
  return count;
}


/**
 * @brief Computes the sparse distance between two blocks
 *
 * @param A : pointer to the first block
 * @param B : pointer to the second block
 * @param indices : indices to the signal coefficients
 * @param num_coeffs: number of signal coefficients
 * @param w : block side
 *
 * @return Sparse distance between blocks A and B
 *
 **/
double compute_sparse_distance(const float *A, const float *B, const int *indices, int num_coeffs, int w) {
  double SD = 0;
  for (int i = w*w-num_coeffs; i < w*w; i++) {
    const double a = A[indices[i]];
    const double b = B[indices[i]];
    SD += fabs(a - b) * max(fabs(a), fabs(b));
  }

  return SD;
}


/**
 * @brief Update the histo position according to the num_coeffs bigger difference values.
 *
 * @param : i_block_A : first DCT block;
 * @param : i_block_B : second DCT block;
 * @param : o_histo_coeffs : histogram of coefficients to update;
 * @param : p_w2 : wxw size of the blocks;
 * @param : p_num_coeffs : number of coefficients to update.
 *
 * @return none.
 */
void updateHisto(
  const float* i_block_A,
  const float* i_block_B,
  int* o_histo_coeffs,
  const int p_w2,
  const int p_num_coeffs) {

  int indices_dists[p_num_coeffs];
  float distances[p_num_coeffs];

  for (int i = 0; i < p_w2; i++) {
    float diff = fabs(i_block_A[i] - i_block_B[i]);
    int pos = i;

    if (i < p_num_coeffs) {
      for (int j = 0; j < i; j++) {
        if (distances[j] > diff) {
          const float t = distances[j];
          distances[j] = diff;
          diff = t;
          const int n = pos;
          pos = indices_dists[j];
          indices_dists[j] = n;
        }
      }
      distances[i] = diff;
      indices_dists[i] = pos;
    }
    else {
      if (distances[0] < diff) {
        distances[0] = diff;
        indices_dists[0] = pos;
        for (int j = 0; j < p_num_coeffs - 1; j++) {
          if (distances[j] < distances[j + 1]) {
            break;
          }
          const float t = distances[j];
          distances[j] = distances[j + 1];
          distances[j + 1] = t;
          const int n = indices_dists[j];
          indices_dists[j] = indices_dists[j + 1];
          indices_dists[j + 1] = n;
        }
      }
    }
  }

  for (int i = 0; i < p_num_coeffs; i++) {
    o_histo_coeffs[indices_dists[i]]++;
  }
}


/**
 * @brief Optimized version of updateHisto when p_num_coeffs = 4.
 *
 * @param : i_block_A : first DCT block;
 * @param : i_block_B : second DCT block;
 * @param : o_histo_coeffs : histogram of coefficients to update;
 * @param : p_w2 : wxw size of the blocks.
 *
 * @return none.
 */
void updateHisto4(
  const float* i_block_A,
  const float* i_block_B,
  int* o_histo_coeffs,
  const int p_w2) {

  //! Sort the 4 first values in ascending orders.
  float diff0 = fabs(i_block_A[0] - i_block_B[0]);
  int pos0 = 0;
  float diff1 = fabs(i_block_A[1] - i_block_B[1]);
  int pos1 = 1;
  if (diff1 < diff0) {
    const float d = diff1;
    diff1 = diff0;
    diff0 = d;
    const int p = pos1;
    pos1 = pos0;
    pos0 = p;
  }
  float diff2 = fabs(i_block_A[2] - i_block_B[2]);
  int pos2 = 2;
  if (diff2 < diff0) {
    const float d = diff2;
    diff2 = diff1;
    diff1 = diff0;
    diff0 = d;
    const int p = pos2;
    pos2 = pos1;
    pos1 = pos0;
    pos0 = p;
  }
  else if (diff2 < diff1) {
    const float d = diff2;
    diff2 = diff1;
    diff1 = d;
    const int p = pos2;
    pos2 = pos1;
    pos1 = p;
  }
  float diff3 = fabs(i_block_A[3] - i_block_B[3]);
  int pos3 = 3;
  if (diff3 < diff0) {
    const float d = diff3;
    diff3 = diff2;
    diff2 = diff1;
    diff1 = diff0;
    diff0 = d;
    const int p = pos3;
    pos3 = pos2;
    pos2 = pos1;
    pos1 = pos0;
    pos0 = p;
  }
  else if (diff3 < diff1) {
    const float d = diff3;
    diff3 = diff2;
    diff2 = diff1;
    diff1 = d;
    const int p = pos3;
    pos3 = pos2;
    pos2 = pos1;
    pos1 = p;
  }
  else if (diff3 < diff2) {
    const float d = diff3;
    diff3 = diff2;
    diff2 = d;
    const int p = pos3;
    pos3 = pos2;
    pos2 = p;
  }

  //! Update the 4 bigger values
  for (int i = 4; i < p_w2; i++) {
    const float diff = fabs(i_block_A[i] - i_block_B[i]);
    const int pos = i;
    if (diff < diff0) {
      continue;
    }
    else if (diff < diff1) {
      diff0 = diff;
      pos0 = pos;
      continue;
    }
    else if (diff < diff2) {
      diff0 = diff1;
      diff1 = diff;
      pos0 = pos1;
      pos1 = pos;
      continue;
    }
    else if (diff < diff3) {
      diff0 = diff1;
      diff1 = diff2;
      diff2 = diff;
      pos0  = pos1;
      pos1  = pos2;
      pos2  = pos;
      continue;
    }
    else {
      diff0 = diff1;
      diff1 = diff2;
      diff2 = diff3;
      diff3 = diff;
      pos0  = pos1;
      pos1  = pos2;
      pos2  = pos3;
      pos3  = pos;
    }
  }

  //! Update the histogram
  o_histo_coeffs[pos0]++;
  o_histo_coeffs[pos1]++;
  o_histo_coeffs[pos2]++;
  o_histo_coeffs[pos3]++;
}


void updateHisto4(
  const float* i_diff,
  int* o_histo_coeffs,
  const int p_w2) {

  //! Sort the 4 first values in ascending orders.
  float diff0 = i_diff[0];
  int pos0 = 0;
  float diff1 = i_diff[1];
  int pos1 = 1;
  if (diff1 < diff0) {
    const float d = diff1;
    diff1 = diff0;
    diff0 = d;
    const int p = pos1;
    pos1 = pos0;
    pos0 = p;
  }
  float diff2 = i_diff[2];
  int pos2 = 2;
  if (diff2 < diff0) {
    const float d = diff2;
    diff2 = diff1;
    diff1 = diff0;
    diff0 = d;
    const int p = pos2;
    pos2 = pos1;
    pos1 = pos0;
    pos0 = p;
  }
  else if (diff2 < diff1) {
    const float d = diff2;
    diff2 = diff1;
    diff1 = d;
    const int p = pos2;
    pos2 = pos1;
    pos1 = p;
  }
  float diff3 = i_diff[3];
  int pos3 = 3;
  if (diff3 < diff0) {
    const float d = diff3;
    diff3 = diff2;
    diff2 = diff1;
    diff1 = diff0;
    diff0 = d;
    const int p = pos3;
    pos3 = pos2;
    pos2 = pos1;
    pos1 = pos0;
    pos0 = p;
  }
  else if (diff3 < diff1) {
    const float d = diff3;
    diff3 = diff2;
    diff2 = diff1;
    diff1 = d;
    const int p = pos3;
    pos3 = pos2;
    pos2 = pos1;
    pos1 = p;
  }
  else if (diff3 < diff2) {
    const float d = diff3;
    diff3 = diff2;
    diff2 = d;
    const int p = pos3;
    pos3 = pos2;
    pos2 = p;
  }

  //! Update the 4 bigger values
  for (int i = 4; i < p_w2; i++) {
    const float diff = i_diff[i];
    const int pos = i;
    if (diff < diff0) {
      continue;
    }
    else if (diff < diff1) {
      diff0 = diff;
      pos0 = pos;
      continue;
    }
    else if (diff < diff2) {
      diff0 = diff1;
      diff1 = diff;
      pos0 = pos1;
      pos1 = pos;
      continue;
    }
    else if (diff < diff3) {
      diff0 = diff1;
      diff1 = diff2;
      diff2 = diff;
      pos0  = pos1;
      pos1  = pos2;
      pos2  = pos;
      continue;
    }
    else {
      diff0 = diff1;
      diff1 = diff2;
      diff2 = diff3;
      diff3 = diff;
      pos0  = pos1;
      pos1  = pos2;
      pos2  = pos3;
      pos3  = pos;
    }
  }

  //! Update the histogram
  o_histo_coeffs[pos0]++;
  o_histo_coeffs[pos1]++;
  o_histo_coeffs[pos2]++;
  o_histo_coeffs[pos3]++;
}


/**
 * @brief Gets the coordinates of the most similar block to the given block
 *
 * @param xB: X coordinate of the given block
 * @param yB: Y coordinate of the given block
 * @param Nx : image width
 * @param Ny : image height
 * @param blocks_dct : pointer to the DCT transformed blocks
 * @param r1 : minimal distance for the candidate blocks
 * @param r2 : maximal distance for the candidate blocks
 * @param w : block side
 *
 * @return The coordinates of the most similar block, compute with the SD
 *
 **/
double get_most_similar_block(int xB, int yB,
                              int Nx, int Ny,
                              float *blocks_dct,
                              int r1, int r2,
                              int w) {

  const int num_coeffs = w*w/4;
  int indices[w*w];

  bool block_found = false;
  float* absDiff = (float*) memalign(16, w*w*sizeof(float));

  // Reference block. Compare with all block_A candidates
  const float *block_B = GET_BLOCK_PTR(xB, yB, blocks_dct, w, Nx, Ny);

  // Preselection
  int histo_coeffs[w*w];
  for (int i = 0; i < w*w; i++)
    histo_coeffs[i] = 0;

  for (int y = yB+r1; y<yB+r2 && y<Ny-w+1; y++)
    for (int x = xB+r1; x<xB+r2 && x<Nx-w+1; x++) {

      // Check spatial distance
      const int r = max(x-xB, y-yB);
      if (r > r1 && r < r2) {
        const float *block_A = GET_BLOCK_PTR(x, y, blocks_dct, w, Nx, Ny);

        // Precompute the absolute difference
        for (int i = 0; i < w*w; i+=4){
          const __m128 a = _mm_loadu_ps(block_A + i);
          const __m128 b = _mm_loadu_ps(block_B + i);
          _mm_storeu_ps(absDiff + i, _mm_max_ps(_mm_sub_ps(_mm_setzero_ps(), a-b), a-b));
        }

        // Sorting indices for the distances
        if (num_coeffs == 4) {
          //updateHisto4(block_A, block_B, histo_coeffs, w*w);
          updateHisto4(absDiff, histo_coeffs, w*w);
        }
        else {
          updateHisto(block_A, block_B, histo_coeffs, w*w, num_coeffs);
        }
      }
    } // for x

  free(absDiff);
  argsort(histo_coeffs, indices, w*w);

  double min_dist = FLT_MAX;
  for (int y = yB+r1; y<yB+r2 && y<Ny-w+1; y++)
    for (int x = xB+r1; x<xB+r2 && x<Nx-w+1; x++) {

      // Check spatial distance
      const int r = max(x-xB, y-yB);

      if (r > r1 && r < r2) {
        const float *block_A = GET_BLOCK_PTR(x, y, blocks_dct, w, Nx, Ny);

        const float distance = compute_sparse_distance(
                                   block_A, block_B,
                                   indices,
                                   num_coeffs, w);

        if (distance < min_dist) {
          min_dist = distance;

          block_found = true;
        }
      } // if r_ok
    } // for x
  return block_found ? min_dist : -1;
}


/**
 * @brief Computes the final unbiased STD for the list of given similar blocks
 *
 * @param blocks_ptr : list of blocks with minimal sparse distance
 * @param indices_SD : the sparse distance associated to each block
 * @param w : block side
 * @param K : quantile for the list of similar block
 * @param ci : horizontal frequency where to measure the noise
 * @param cj : vertical frequency where to measure the noise
 *
 * @return Unbiased tilda_sigmal of the noise
 *
 **/
float compute_std(float **blocks_ptr,
                  int *indices_SD, int w, int K,
                  int ci, int cj) {
  double *coeffs = new double[K];

  // Obtain the coefficients according to indices_SD
  for (int k = 0; k < K; k++) {
    float *block = blocks_ptr[indices_SD[k]];
    float coef_ci_cj = block[cj*w + ci];
    coeffs[k] = coef_ci_cj;
  }

  // MAD estimator for the STD
  double hat_std = compute_MAD<double>(coeffs, K);
  if (hat_std < 0) hat_std = 0;

  //delete[] squared_coeffs;
  delete[] coeffs;

  //float tilde_sigma = get_tilde_sigma(hat_std, w);
  float tilde_sigma = get_tilde_sigma(hat_std, w);
  //
  if (tilde_sigma < 0) tilde_sigma = 0;
  return tilde_sigma;
}


/**
 * @brief Reads all images blocks taking care of the image boundaries
 *
 * @param D : list of read blocks (output)
 * @param u : input image
 * @param Nx: image width
 * @param Ny: image height
 * @param w : block side
 *
 * @return Number of read blocks
 *
 **/
int read_all_blocks(float *D,
                    float *u,
                    int Nx, int Ny,
                    int w) {
  int q = 0;
  // Locate a valid pixel in the image
  for (int y = 0; y < Ny-w+1; y++) {
    for (int x = 0; x < Nx-w+1; x++) {
      // Read block
      float *block = GET_BLOCK_PTR(x, y, D, w, Nx, Ny);
      for (int j = 0; j < w; j++) {
        for (int i = 0; i < w; i++) {
          block[j*w+i] = u[(y+j)*Nx+(x+i)];
        }
      }
      q++;
    }
  }
  return q;
}


/**
 * @brief Returns the mean associated to the given bin
 *
 * @param mean_method : method to compute the mean (1: mean of means, 2: median of means)
 * @param K : number of elements to compute the mean intensity
 * @param indices_SD : SD associated to each block in the bin
 * @param histo : histogram object
 *
 * @return The mean of the bin
 *
 **/
float get_bin_mean(int mean_method, int K, int *indices_SD,
                   int bin, histo_t &histo) {
  float *means_bin = histo.get_datal_bin(bin);
  float bin_mean;

  float *values = new float[K+1];
  for (int i = 0; i <= K; i++)
    values[i] = means_bin[indices_SD[i]];

  switch (mean_method) {
    case 1: { // mean of means
      bin_mean = 0.0;
      for (int i = 0; i <= K; i++) {
        bin_mean += values[i];
      }
      bin_mean /= (K+1);
      break;
    }
    case 2: { // median of means
      bin_mean = median<float>(values, K+1);
      break;
    }
    default: {
      PRINT_ERROR("Unknown mean method: %d\n", mean_method);
      exit(-1);
    }
  }
  delete[] values;
  return bin_mean;
}


/**
 * @brief In-place Normalization of the FFTW to get a orthonormal 2D DCT-II
 *
 * @param blocks : input/output list of transformed blocks
 * @param w : block side
 * @param num_blocks : number of blocks in the list
 *
 **/
void normalize_FFTW(float *blocks, int w, int num_blocks) {
  const float ONE_DIV_2w = 1.0 / (2.0*w);
  const float ONE_DIV_SQRT_2 = 1.0 / sqrtf(2);

  // Divide all coefficients by 2*w
  //#pragma omp parallel for shared(blocks)
  for (int i = 0; i < num_blocks*w*w; i++)
    blocks[i] *= ONE_DIV_2w;

  #ifdef _OPENMP
  #pragma omp parallel for shared(blocks) schedule(static)
  #endif
  for (int b = 0; b < num_blocks; b++) {
    // {(i, j)} with i == 0 or j == 0
    for (int j = 0; j < w; j++) {
      int i = 0;
      blocks[b*w*w+j*w+i] *= ONE_DIV_SQRT_2;
    }
    for (int i = 0; i < w; i++) {
      int j = 0;
      blocks[b*w*w+j*w+i] *= ONE_DIV_SQRT_2;
    }
  }
}


/**
 * @brief Maps the pointers of the blocks to their sparse distance
 *
 * @param SDs : (output) map whichs links blocks with their sparse distance
 * @param blocks_dct : list of blocks
 * @param w : block side
 * @param r1 : minimal distance for similar block candidates
 * @param r2 : maximal distance for similar block candidates
 * @param Nx : image width
 * @param Ny : image height
 *
 **/
void map_blocks(TMapSD &SDs,
                float *blocks_dct,
                int w,
                int r1, int r2,
                int Nx, int Ny) {
  // For each block B, get the most similar block A
  int num_blocks = (Nx-w+1)*(Ny-w+1);
  vector<float*> v_pointers(num_blocks);
  vector<float> v_SDs(num_blocks);

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (int yB = 0; yB < Ny-w+1; yB++) {
    for (int xB = 0; xB < Nx-w+1; xB++) {

      const double SD = get_most_similar_block(xB, yB,
                                         Nx, Ny,
                                         blocks_dct,
                                         r1, r2,
                                         w);
      const bool block_found = (SD != -1);

      // Map the SD with the block pointer
      float *ptr = GET_BLOCK_PTR(xB, yB, blocks_dct, w, Nx, Ny);

      // Store and do mapping
      int pos = yB*(Nx-w+1) + xB;
      v_pointers[pos] = ptr;
      v_SDs[pos] = block_found ? SD : -1;
    } // for xB
  } // for yB

  // Do the mapping
  for (int xByB = 0; xByB < num_blocks; xByB++) {
    float *ptr = v_pointers[xByB];
    const float SD = v_SDs[xByB];
    SDs[ptr] = SD;
  }
}


/**
 * @brief Computes all the means given the whole list of blocks in the image
 *
 * @param means : (output) list of means of the blocks
 * @param blocks_dct : list of DTC transformed blocks
 * @param w : blocks side
 * @param Nx : image width
 * @param Ny : image height
 *
 * @return Number of blocks
 *
 **/
int compute_means(float *means, float *blocks_dct,
                  int w, int Nx, int Ny) {
  int count = 0;
  for (int yB = 0; yB < Ny-w+1; yB++)
    for (int xB = 0; xB < Nx-w+1; xB++)
     means[count++] = *GET_BLOCK_PTR(xB, yB, blocks_dct, w, Nx, Ny) / w;
  return count;

  // D[q*w*w+j*w+i] = u[(y+j)*Nx+(x+i)];
}


/**
 * @brief Saves a noise curve into a text file
 *
 * @param vstds : list of STDs
 * @param vmeans : list means corresponding to each STD
 * @param num_bins : number of estimations
 * @param ci : horizontal frequency
 * @param cj : vertical frequency
 * @param num_channels : number of channels
 * @param w : block side
 * @param output_format : format for the output file name
 *
 **/
void save_noise_curve(float *vstds, float *vmeans,
                        int num_bins,
                        int ci, int cj,
                        int num_channels,
                        int w,
                        const char *output_format) {
  FILE *fp;
  char filename[1024];

  sprintf(filename, output_format, ci, cj);
  fp = fopen(filename, "wt");

  int ij = w*cj + ci;

  for (int bin = 0; bin < num_bins; bin++) {
    // Means
    for (int ch = 0; ch < num_channels; ch++)
        fprintf(fp, "%f  ", vmeans[num_bins*w*w*ch + num_bins*ij + bin]);

    // Standard deviations
    for (int ch = 0; ch < num_channels; ch++)
        fprintf(fp, "%f  ", vstds[num_bins*w*w*ch + num_bins*ij + bin]);
    //
    fprintf(fp, "\n");
  }
  fclose(fp);
}


/**
 * @brief Splits the blocks into bins according to their mean intensity
 *
 * @param Nx : image width
 * @param Ny : image height
 * @param w : block side
 * @param num_bins : number of bins to split the blocks into
 * @param num_blocks : number of blocks to split
 * @param blocks : blocks that will be split
 * @param mask_all : mask to avoid blocks with completely saturated pixels
 * @param SDs: map to link block pointers with the sparse distances
 *
 **/
histo_t create_histogram(int Nx, int Ny, int w, int num_bins, int num_blocks,
                         float *blocks,
                         float *means_histo, float *means,
                         int *mask_all, TMapSD &SDs) {
  // Block list for the histogram
  float **blocks_ptr = new float*[num_blocks];

  int num_blocks_histo = 0;
  for (int y = 0; y < Ny-w+1; y++) {
    for (int x = 0; x < Nx-w+1; x++) {
      bool mask_ok = (mask_all == NULL || mask_all[y*Nx+x] == 0);

      float *block_ptr = GET_BLOCK_PTR(x, y, blocks, w, Nx, Ny);
      bool valid_block = (SDs[block_ptr] > 0); // avoid blocks
      // with SD = -1 (without a candidate)

      if (mask_ok && valid_block) {
        blocks_ptr[num_blocks_histo] = block_ptr;
        if (means_histo != NULL)
          means_histo[num_blocks_histo] = means[y*(Nx-w+1)+x];
        //
        num_blocks_histo++;
      }
    }
  }

  // Create histogram according to the means
  histo_t histo = histo_t(num_bins,
                          blocks_ptr,
                          means_histo,
                          num_blocks_histo,
                          true);

  return histo;
}


/**
 * @brief Gets the sorting indices of the blocks depending on their SD
 *
 * @param indices_SD : (output) sorting indices
 * @param histo : histogram object
 * @param bin : number of bin
 * @param SDs: map whichs links blocks with their sparse distance
 *
 **/
void get_indices_SD(int *indices_SD, histo_t &histo, int bin, TMapSD &SDs) {
  int elems_bin = histo.get_num_elements_bin(bin);

  float **block_ptr_bin = histo.get_data_bin(bin);

  //
  float *array_SD = new float[elems_bin];
  for (int i = 0; i < elems_bin; i++) {
    float *ptr = block_ptr_bin[i];
    array_SD[i] = SDs[ptr];
  }

  argsort(array_SD, indices_SD, elems_bin);

  delete[] array_SD;
}


/**
 * @brief Estimates and write the noise estimation for all channels
 *
 * @param blocks_DCT : the DCT transformed blocks pointer assigned memory
 * @param blocks : the blocks pointer assigned memory
 * @param means : the means of the blocks
 * @param Nx : image width
 * @param Ny : image height
 * @param fft_plan : FFTW plan used to transform the blocks
 * @param num_blocks: number of blocks
 * @param num_bins: number of bins
 * @param p: quantile of the blocks on their sparse distance
 * @param mean_method: method to compute the mean; mean (1)/median(2) of means
 * @param vmeans : list means corresponding to each STD
 * @param vstds : list of STDs
 * @param w : block side
 * @param r1 : minimal distance for the candidate blocks
 * @param r2 : maximal distance for the candidate blocks
 * @param num_channels : number of channels
 * @param input : input image
 * @param mask_all : mask to avoid blocks with completely saturated pixels
 *
 **/
void process_all(float *blocks_dct, float *blocks,
                 float *means,
                 int Nx, int Ny,
                 fftwf_plan fft_plan,
                 int num_blocks, int num_bins, float p,
                 int mean_method, float *vmeans, float *vstds,
                 int w, int r1, int r2,
                 int num_channels,
                 CImage &input,
                 int *mask_all) {
  float *means_histo = new float[Nx*Ny];

  for (int ch = 0; ch < num_channels; ch++) {
    float *u = input.get_channel(ch);

    // Read all overlapping the blocks
    int total_blocks = read_all_blocks(blocks, u, Nx, Ny, w);

    // Transform blocks with FFTW
    fftwf_execute_r2r(fft_plan, blocks, blocks_dct);
    normalize_FFTW(blocks_dct, w, total_blocks);

    // Map the SDs to the blocks_dct pointers
    TMapSD SDs;
    map_blocks(SDs, blocks_dct,
               w,
               r1, r2,
               Nx, Ny);

    // Compute means
    compute_means(means, blocks_dct,
                                  w, Nx, Ny);

    // Split blocks according to their mean intensity
    histo_t histo = create_histogram(Nx, Ny, w, num_bins, num_blocks,
                                     blocks_dct,
                                     means_histo, means,
                                     mask_all, SDs);

    // For each bin...
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int bin = 0; bin < num_bins; bin++) {
      int elems_bin = histo.get_num_elements_bin(bin);
      int *indices_SD = new int[elems_bin];
      get_indices_SD(indices_SD, histo, bin, SDs);

      // Process all coefficients in the bin
      for (int ij = 0; ij < w*w; ij++) {
        const unsigned cj = ij / w;
        const unsigned ci = ij - cj * w;

        float **block_ptr_bin = histo.get_data_bin(bin);
        int K = elems_bin * p;

        float tilde_sigma = compute_std(block_ptr_bin,
                                        indices_SD, w, K,
                                        ci, cj);
        float bin_mean = get_bin_mean(mean_method, K, indices_SD, bin, histo);

        // Store results
        vmeans[num_bins*ch + bin] = bin_mean;
        vstds[ch*w*w*num_bins + bin*w*w +ij] = tilde_sigma;
      } // for ij

      delete[] indices_SD;
    } // for bin
  } // for ch

  delete[] means_histo;
}


/**
 * @brief From here, it's a wrapper to join the Noise Estimation algorithm
 * developped by Miguel Colom to the NL-Bayes multi-scale denoising.
 *
 **/

//! Ponomarenko-diff noise estimation algorithm.
int runNoiseEstimate(
  std::vector<float> const &i_im,
  const ImageSize& p_imSize,
  nlbParams &io_nlbParams,
  const MsdParams &p_msdParams){

  //! Parameters initialization
  neParams params;
  params.sizePatch        = io_nlbParams.sizePatch;
  params.r1               = 4;
  params.r2               = 14;
  params.percentile       = 0.005f;
  params.correctionFactor = 0.5f * 7.460923461f; // Work only for sP = 4. Needs to be tabulated for other values.
  params.numBins          = std::max(p_imSize.wh / 42000.f, 1.f);
  const int meanMethod    = 2;

  //! For convenience
  const unsigned int sP  = params.sizePatch;
  const unsigned int sP2 = sP * sP;

//! Parallelization config
#ifdef _OPENMP
  omp_set_num_threads(omp_get_num_procs());
#endif

  //! Number of overlapping blocks
  const int numBlocks = (p_imSize.width - sP + 1) * (p_imSize.height - sP + 1);

  //! Build the mask of valid pixels
  vector<int> maskAll(p_imSize.wh);
  CImage im(p_imSize.width, p_imSize.height, p_imSize.nBits, p_imSize.nChannels);
  for (unsigned int c = 0; c < p_imSize.nChannels; c++) {
    const float* iI = &i_im[c * p_imSize.wh];
    float* oI = im.get_channel(c);

    for (unsigned int k = 0; k < p_imSize.wh; k++) {
      oI[k] = iI[k];
    }
  }
  buildMask(im, &maskAll[0], p_imSize.width, p_imSize.height, sP, p_imSize.nChannels);

  //! Init FFTW threads
  fftwf_init_threads();

  int nbTable[2] = {sP, sP};
  int nembed [2] = {sP, sP};

#ifdef _OPENMP
  fftwf_plan_with_nthreads(omp_get_num_procs());
#endif

  fftwf_r2r_kind kindTable[2] = {FFTW_REDFT10, FFTW_REDFT10};

  //! Blocks from the image
  vector<float> blocks(numBlocks * sP2);

  //! Means of the blocks from the image
  vector<float> means(numBlocks * sP2);

  //! DCT blocks
  vector<float> blocksDct (numBlocks * sP2);

  //! Arrays to store the final means and noise estimations
  vector<float> vMeans(p_imSize.nChannels * params.numBins);
  vector<float> vStds (p_imSize.nChannels * sP2 * params.numBins);

  fftwf_plan fft_plan = fftwf_plan_many_r2r(2, nbTable, numBlocks, &blocks[0],
                                            nembed, 1, sP2, &blocksDct[0], nembed,
                                            1, sP2, kindTable, FFTW_ESTIMATE);

  //! Process all coefficients
  process_all(&blocksDct[0], &blocks[0], &means[0], p_imSize.width, p_imSize.height,
              fft_plan, numBlocks, params.numBins, params.percentile, meanMethod, &vMeans[0],
              &vStds[0], sP, params.r1, params.r2, p_imSize.nChannels, im, &maskAll[0]);

  //! Save low and high-frequencies noise curve into txt file.
  if (!io_nlbParams.isFirstStep) {
    if (writeNoiseCurves(vMeans, vStds, sP, params.numBins, p_msdParams.path,
                         p_msdParams.currentScale, p_imSize.nChannels) != EXIT_SUCCESS) {
      return EXIT_FAILURE;
    }
  }

  //! Filter curves
  filterNoiseCurves(vMeans, vStds, sP, params.numBins, p_imSize.nChannels,
                    io_nlbParams.rangeOffset, io_nlbParams.rangeMax);

  //! Apply correction factor
  for (unsigned int k = 0; k < p_imSize.nChannels * sP2 * params.numBins; k++) {
    vStds[k] *= p_msdParams.stdFactor;
  }

  //! Look for NaN
  for (unsigned int k = 0; k < p_imSize.nChannels * sP2 * params.numBins; k++) {
    if (vStds[k] != vStds[k]) {
      cout << "NaN detected in SD of the estimated noise curve" << endl;
      return EXIT_FAILURE;
    }
    if (vStds[k] < 0.f) {
      cout << "Negative values detected in SD of the estimated noise curve : " << vStds[k] << endl;
      vStds[k] = 0.f;
    }
  }

  //! Compute for each matrix the Dt M D product to get back to the spatial space.
  vector<vector<float> > matrixBins, meanBins;
  getMatrixBins(vMeans, vStds, meanBins, matrixBins, params, p_imSize);

  //! Interpolation and filtering of the estimated matrices
  if (estimatedStdMatrices(meanBins, matrixBins, io_nlbParams.noiseCovMatrix, sP, params.numBins,
                              io_nlbParams.rangeOffset, io_nlbParams.rangeMax) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! FFTW Cleanup
  fftwf_destroy_plan(fft_plan);
  fftwf_cleanup_threads();
  fftwf_cleanup();

  return EXIT_SUCCESS;
  //return EXIT_FAILURE;
 }


//! From matrices M in frequency space, obtain the Dt M D matrix in the spatial space,
//! i.e. the covariance matrix of the noise.
void getMatrixBins(
  std::vector<float> const& i_vMeans,
  std::vector<float> const& i_vStds,
  std::vector<std::vector<float> > &o_meanBins,
  std::vector<std::vector<float> > &o_matrixBins,
  neParams const& p_params,
  const ImageSize &p_imSize){

  //! For convenience
  const unsigned int nChnls = p_imSize.nChannels;
  const unsigned int sP     = p_params.sizePatch;
  const unsigned int sP2    = sP * sP;
  const unsigned int nBins  = p_params.numBins;

  //! Sizes initialization
  o_meanBins.resize(nChnls);
  for (unsigned int c = 0; c < nChnls; c++) {
    o_meanBins[c].resize(nBins);
  }
  o_matrixBins.resize(nChnls * nBins);
  for (unsigned int n = 0; n < nChnls * nBins; n++) {
    o_matrixBins[n].resize(sP2 * sP2);
  }

  //! Get the normalized 2D DCT matrix
  vector<float> DCT;
  getNormalized2dDctMatrix(DCT, sP, sP);

  //! Compute DCT transposed
  vector<float> DCTt;
  transposeMatrix(DCTt, DCT, sP2);

  //! Get for each bin the matrices
  for (unsigned int c = 0; c < nChnls; c++) {
    float* oM = &o_meanBins[c][0];
    const float* iM = &i_vMeans[c * nBins];

    for (unsigned int bin = 0; bin < nBins; bin++) {

      //! Get the mean
      oM[bin] = iM[bin];

      //! Get the diagonal estimated matrix mat = Dt * S * D
      const float* iS = &i_vStds[c * nBins * sP2 + bin * sP2];
      for (unsigned int i = 0; i < sP2; i++) {
        float* oS = &o_matrixBins[c * nBins + bin][i * sP2];
        const float* iDt = &DCTt[i * sP2];

        for (unsigned int j = 0; j < sP2; j++) {
          const float* iD = &DCTt[j * sP2];
          float val = 0.f;

          for (unsigned int k = 0; k < sP2; k++) {
            val += iDt[k] * iD[k] * iS[k] * iS[k];
          }
          oS[j] = val;
        }
      }
    }
  }
}


//! Estimate by linear interpolation one covariance matrix of the noise for
//! each intensity between [0, max + offset].
int estimatedStdMatrices(
  std::vector<std::vector<float> > const& i_meanBins,
  std::vector<std::vector<float> > const& i_matrixBins,
  std::vector<std::vector<float> > &o_covMat,
  const unsigned int p_sizePatch,
  const unsigned int p_numBins,
  const unsigned int p_rangeOffset,
  const unsigned int p_rangeMax){

  //! Initializations
  const unsigned int chnls = i_meanBins.size();
  const unsigned int sP2   = p_sizePatch * p_sizePatch;
  const unsigned int nBins = p_numBins;
  const unsigned int range = p_rangeOffset + p_rangeMax + 1;
  o_covMat.assign(range, vector<float> (chnls * sP2 * sP2, 0.f));

  //! Check sizes
  if (nBins == 0) {
    cout << "estimatedStdMatrices - error : at least one bin is needed" << endl;
    return EXIT_FAILURE;
  }
  for (unsigned int c = 0; c < chnls; c++) {
    if (i_meanBins[c].size() != nBins) {
      cout << "estimatedStdMatrices - error : meanBins doesn't have the right number of bins" << endl;
      return EXIT_FAILURE;
    }
    //! Check bins
    const float* iM = &i_meanBins[c][0];
    for (unsigned int n = 0; n < nBins - 1; n++) {
      if (iM[n + 1] < iM[n]) {
        cout << "estimatedStdMatrices - error : means must be sorted" << endl;
        return EXIT_FAILURE;
      }
    }
  }
  if (i_matrixBins.size() != chnls * nBins) {
    cout << "estimatedStdMatrices - error : number of bins in matrix doesn't fit" << endl;
    for (unsigned int n = 0; n < chnls * nBins; n++) {
      if (i_matrixBins[n].size() != sP2 * sP2) {
        cout << "estimatedStdMatrices - error : size of patches doesn't fit with matrix sizes" << endl;
        return EXIT_FAILURE;
      }
    }
  }

  //! For each channels...
  for (unsigned int c = 0; c < chnls; c++) {

    //! Get coefficients along the bins
    for (unsigned int k = 0; k < sP2 * sP2; k++) {
      vector<float> values(nBins);

      for (unsigned int n = 0; n < nBins; n++) {
        values[n] = i_matrixBins[c * nBins + n][k];
      }

      //! Interpolation of the values
      vector<float> allValues;
      interpValues(i_meanBins[c], values, allValues, p_rangeOffset, p_rangeMax);

      //! Save values
      for (unsigned int r = 0; r < range; r++) {
        o_covMat[r][c * sP2 * sP2 + k] = allValues[r];
      }
    }
  }

  //! Normal exit
  return EXIT_SUCCESS;
}


//! Get the DCT normalized matrix
void getNormalized2dDctMatrix(
  std::vector<float> &o_mat,
  const unsigned int p_rows,
  const unsigned int p_cols){

  //! Initialization
  const float norm     = 2.f / (sqrtf(float(p_rows)) * sqrtf(float(p_cols)));
  const float sqrt2    = sqrtf(2.f);
  const float invSqrt2 = 1.f / sqrt2;
  const float coefRows = M_PI / float(p_rows);
  const float coefCols = M_PI / float(p_cols);
  o_mat.resize(p_rows * p_rows * p_cols * p_cols);

  //! Get the matrix
  for (unsigned int a = 0; a < p_rows; a++) {
    for (unsigned int b = 0; b < p_cols; b++) {
      for (unsigned int m = 0; m < p_rows; m++) {
        float* iM = &o_mat[(a * p_cols + b) * p_cols * p_rows + m * p_cols];
        for (unsigned int n = 0; n < p_cols; n++) {
          iM[n]  = cos(coefRows * a * (m + 0.5f)) * cos(coefCols * b * (n + 0.5f));
          iM[n] *= (a == 0 ? invSqrt2 : 1.f) * (b == 0 ? invSqrt2 : 1.f) * norm;
        }
      }
    }
  }
}


//! Transpose a matrix
void transposeMatrix(
  std::vector<float> &o_mat,
  std::vector<float> const& i_mat,
  const unsigned int p_sP){

  o_mat.resize(p_sP * p_sP);

  for (unsigned int i = 0; i < p_sP; i++) {
    const float* iM = &i_mat[i * p_sP];
    float*       oM = &o_mat[i];

    for (unsigned int j = 0; j < p_sP; j++) {
      oM[j * p_sP] = iM[j];
    }
  }
}



