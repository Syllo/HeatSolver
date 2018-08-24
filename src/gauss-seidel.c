/*
 * Copyright (c) 2018 Maxime Schmitt <max.schmitt@unistra.fr>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <tgmath.h>

#include "gauss-seidel.h"

#define ceild(n, d) (((n) + (d)-1) / (d))
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

void solve_gauss_seidel(size_t sizeX, size_t sizeY,
                        double t[restrict sizeX][sizeY], double stop_criteria,
                        double dx, double dy, double *restrict error_end,
                        size_t *num_iterations) {
  assert(sizeX >= 3 && sizeY >= 3);
  double factorX = 0.5 * dy * dy / ((dx * dx) + (dy * dy));
  double factorY = 0.5 * dx * dx / ((dx * dx) + (dy * dy));

  stop_criteria *= stop_criteria;
  do {
    double max_error = 0.;
    for (size_t i = 1; i < sizeX - 1; ++i) {
      for (size_t j = 1; j < sizeY - 1; ++j) {
        double old_val = t[i][j];
        t[i][j] = (t[i - 1][j] + t[i + 1][j]) * factorX +
                  (t[i][j - 1] + t[i][j + 1]) * factorY;
        double error = old_val - t[i][j];
        error *= error;
        max_error = max(max_error, error);
      }
    }

    *error_end = max_error;

    (*num_iterations)++;
    if (!finite(*error_end))
      break;
  } while (*error_end > stop_criteria);
  *error_end = sqrt(*error_end);
}

#define S1(i, j)                                                               \
  t[i][j] = (t[i - 1][j] + t[i + 1][j]) * factorX +                            \
            (t[i][j - 1] + t[i][j + 1]) * factorY;

void solve_gauss_seidel_parallel(size_t sizeX, size_t sizeY,
                                 double t[restrict sizeX][sizeY],
                                 double stop_criteria, double dx, double dy,
                                 double *restrict error_end,
                                 size_t *num_iterations) {
  assert(sizeX >= 3 && sizeY >= 3);
  double factorX = 0.5 * dy * dy / ((dx * dx) + (dy * dy));
  double factorY = 0.5 * dx * dx / ((dx * dx) + (dy * dy));

  stop_criteria *= stop_criteria;
  double max_error;
#pragma omp parallel shared(num_iterations, error_end, stop_criteria, sizeX,   \
                            sizeY, t, max_error, factorX, factorY)
  {
    do {
#pragma omp barrier
#pragma omp single
      max_error = 0.;
      for (size_t t1 = 2; t1 <= sizeX + sizeY - 4; t1++) {
        size_t lbp = sizeX > t1 + 1 ? 1 : t1 - sizeX + 2;
        size_t ubp = min(sizeY - 2, t1 - 1);
#pragma omp for reduction(max : max_error)
        for (size_t t2 = lbp; t2 <= ubp; t2++) {
          double old_val = t[t1 - t2][t2];
          S1((t1 - t2), t2);
          double error = old_val - t[t1 - t2][t2];
          error *= error;
          max_error = max(max_error, error);
        }
      }
#pragma omp master
      (*num_iterations)++;
      if (!finite(max_error))
        break;
    } while (max_error > stop_criteria);
  }
  *error_end = sqrt(max_error);
}

void solve_gauss_seidel_parallel_tiled(size_t sizeX, size_t sizeY,
                                       double t[restrict sizeX][sizeY],
                                       double stop_criteria, double dx,
                                       double dy, double *restrict error_end,
                                       size_t *num_iterations) {
  assert(sizeX >= 3 && sizeY >= 3);
  double factorX = 0.5 * dy * dy / ((dx * dx) + (dy * dy));
  double factorY = 0.5 * dx * dx / ((dx * dx) + (dy * dy));

  stop_criteria *= stop_criteria;
  double max_error;
#pragma omp parallel shared(num_iterations, error_end, stop_criteria, sizeX,   \
                            sizeY, t, max_error, factorX, factorY)
  {
    do {
#pragma omp barrier
#pragma omp single
      max_error = 0.;

      for (size_t t1 = 0; t1 <= (sizeX + sizeY - 4) / 32; t1++) {
        size_t lbp = 32 * t1 + 2 >= sizeX ? ceild(32 * t1 - sizeX + 2, 32) : 0;
        size_t ubp = min((sizeY - 2) / 32, t1);
#pragma omp for reduction(max : max_error)
        for (size_t t2 = lbp; t2 <= ubp; t2++) {
          for (size_t t3 = max(1, 32 * t1 - 32 * t2);
               t3 <= min(sizeX - 2, 32 * t1 - 32 * t2 + 31); t3++) {
            for (size_t t4 = max(1, 32 * t2);
                 t4 <= min(sizeY - 2, 32 * t2 + 31); t4++) {
              double old_val = t[t3][t4];
              S1(t3, t4);
              double error = old_val - t[t3][t4];
              error *= error;
              max_error = max(max_error, error);
            }
          }
        }
      }
#pragma omp master
      (*num_iterations)++;
      if (!finite(max_error))
        break;
    } while (max_error > stop_criteria);
  }
  *error_end = sqrt(max_error);
}
