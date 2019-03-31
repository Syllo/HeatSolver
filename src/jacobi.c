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
#include <stdint.h>
#include <tgmath.h>

#include "jacobi.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

void solve_jacobi(size_t sizeX, size_t sizeY, double t[restrict sizeX][sizeY],
                  double tNext[restrict sizeX][sizeY], double stop_criteria,
                  double dx, double dy, double *restrict error_end,
                  size_t *num_iterations) {

  double rev_dx2 = 0.5 * dy * dy / ((dx * dx) + (dy * dy));
  double rev_dy2 = 0.5 * dx * dx / ((dx * dx) + (dy * dy));

  assert(sizeX >= 3 && sizeY >= 3);

  const size_t iteration_stop =
      *num_iterations != 0 ? *num_iterations : SIZE_MAX;
  stop_criteria =
      *num_iterations != 0 ? -HUGE_VAL : stop_criteria * stop_criteria;
  *num_iterations = 0;

  double max_error;
  do {
    max_error = 0.;
    for (size_t i = 1; i < sizeX - 1; ++i) {
      for (size_t j = 1; j < sizeY - 1; ++j) {
        tNext[i][j] = (t[i - 1][j] + t[i + 1][j]) * rev_dx2 +
                      (t[i][j - 1] + t[i][j + 1]) * rev_dy2;
        double error = t[i][j] - tNext[i][j];
        error *= error;
        max_error = max(error, max_error);
      }
    }
    (*num_iterations)++;
    if (!finite(max_error))
      break;
    void *tmp = tNext;
    tNext = t;
    t = tmp;
  } while (max_error > stop_criteria && *num_iterations != iteration_stop);
  *error_end = sqrt(max_error);
}

void solve_jacobi_vectorized(size_t sizeX, size_t sizeY,
                             double t[restrict sizeX][sizeY],
                             double tNext[restrict sizeX][sizeY],
                             double stop_criteria, double dx, double dy,
                             double *restrict error_end,
                             size_t *num_iterations) {
  assert(sizeX >= 3 && sizeY >= 3);

  double rev_dx2 = 0.5 * (dy * dy / ((dx * dx) + (dy * dy)));
  double rev_dy2 = 0.5 * (dx * dx / ((dx * dx) + (dy * dy)));

  const size_t iteration_stop =
      *num_iterations != 0 ? *num_iterations : SIZE_MAX;
  stop_criteria =
      *num_iterations != 0 ? -HUGE_VAL : stop_criteria * stop_criteria;
  *num_iterations = 0;

  double max_error;
  do {
    max_error = 0.;
    for (size_t i = 1; i < sizeX - 1; ++i) {
#pragma omp simd reduction(max : max_error)
      for (size_t j = 1; j < sizeY - 1; ++j) {
        tNext[i][j] = (t[i - 1][j] + t[i + 1][j]) * rev_dx2 +
                      (t[i][j - 1] + t[i][j + 1]) * rev_dy2;
        double error = t[i][j] - tNext[i][j];
        error *= error;
        max_error = max(error, max_error);
      }
    }
    (*num_iterations)++;
    if (!finite(max_error))
      break;
    void *tmp = tNext;
    tNext = t;
    t = tmp;
  } while (max_error > stop_criteria && *num_iterations != iteration_stop);
  *error_end = sqrt(max_error);
}

void solve_jacobi_parallel(size_t sizeX, size_t sizeY,
                           double t[restrict sizeX][sizeY],
                           double tNext[restrict sizeX][sizeY],
                           double stop_criteria, double dx, double dy,
                           double *restrict error_end, size_t *num_iterations) {
  assert(sizeX >= 3 && sizeY >= 3);

  double rev_dx2 = 0.5 * (dy * dy / ((dx * dx) + (dy * dy)));
  double rev_dy2 = 0.5 * (dx * dx / ((dx * dx) + (dy * dy)));

  const size_t iteration_stop =
      *num_iterations != 0 ? *num_iterations : SIZE_MAX;
  stop_criteria =
      *num_iterations != 0 ? -HUGE_VAL : stop_criteria * stop_criteria;
  *num_iterations = 0;

  double max_error;
#pragma omp parallel default(none)                                             \
    shared(sizeX, sizeY, max_error, stop_criteria, num_iterations, rev_dx2,    \
           rev_dy2) firstprivate(t, tNext)
  {
    size_t local_num_iteration = 0;
    do {
#pragma omp barrier
#pragma omp single
      max_error = 0.;
#pragma omp for reduction(max : max_error)
      for (size_t i = 1; i < sizeX - 1; ++i) {
#pragma omp simd reduction(max : max_error)
        for (size_t j = 1; j < sizeY - 1; ++j) {
          tNext[i][j] = (t[i - 1][j] + t[i + 1][j]) * rev_dx2 +
                        (t[i][j - 1] + t[i][j + 1]) * rev_dy2;
          double error = t[i][j] - tNext[i][j];
          error *= error;
          max_error = max(error, max_error);
        }
      }
#pragma omp master
      (*num_iterations)++;
      local_num_iteration++;
      if (!finite(max_error))
        break;
      void *tmp = tNext;
      tNext = t;
      t = tmp;
    } while (max_error > stop_criteria &&
             local_num_iteration != iteration_stop);
  }
  *error_end = sqrt(max_error);
}

#define S1(i, j)                                                               \
  tNext[i][j] = (t[i - 1][j] + t[i + 1][j]) * rev_dx2 +                        \
                (t[i][j - 1] + t[i][j + 1]) * rev_dy2;

void solve_jacobi_parallel_tiled(size_t sizeX, size_t sizeY,
                                 double t[restrict sizeX][sizeY],
                                 double tNext[restrict sizeX][sizeY],
                                 double stop_criteria, double dx, double dy,
                                 double *restrict error_end,
                                 size_t *num_iterations) {
  assert(sizeX >= 3 && sizeY >= 3);
  double rev_dx2 = 0.5 * (dy * dy / ((dx * dx) + (dy * dy)));
  double rev_dy2 = 0.5 * (dx * dx / ((dx * dx) + (dy * dy)));

  const size_t iteration_stop =
      *num_iterations != 0 ? *num_iterations : SIZE_MAX;
  stop_criteria =
      *num_iterations != 0 ? -HUGE_VAL : stop_criteria * stop_criteria;
  *num_iterations = 0;

  double max_error;
#pragma omp parallel shared(sizeX, sizeY, max_error, stop_criteria,            \
                            num_iterations) firstprivate(t, tNext)
  {
    size_t local_num_iteration = 0;
    do {
#pragma omp barrier
#pragma omp single
      max_error = 0.;
#pragma omp for reduction(max : max_error)
      for (size_t t1 = 0; t1 <= (sizeX - 2) / 32; t1++) {
        for (size_t t2 = 0; t2 <= (sizeY - 2) / 32; t2++) {
          for (size_t t3 = max(1, 32 * t1); t3 <= min(sizeX - 2, 32 * t1 + 31);
               t3++) {
            size_t lbv = max(1, 32 * t2);
            size_t ubv = min(sizeY - 2, 32 * t2 + 31);
#pragma omp simd reduction(max : max_error)
            for (size_t t4 = lbv; t4 <= ubv; t4++) {
              S1(t3, t4);
              double error = t[t3][t4] - tNext[t3][t4];
              error *= error;
              max_error = max(error, max_error);
            }
          }
        }
      }
#pragma omp master
      (*num_iterations)++;
      local_num_iteration++;
      if (!finite(max_error))
        break;
      void *tmp = tNext;
      tNext = t;
      t = tmp;
    } while (max_error > stop_criteria &&
             local_num_iteration != iteration_stop);
  }
  *error_end = sqrt(max_error);
}
