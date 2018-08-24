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

#include <stdbool.h>
#include <stdio.h>
#include <sys/stat.h>

#include "heatSolver.h"

void init_heat_values(size_t sizeX, size_t sizeY, double t[sizeX][sizeY],
                      double (*f)(size_t, size_t)) {
  for (size_t i = 0; i < sizeX; ++i) {
    for (size_t j = 0; j < sizeY; ++j) {
      t[i][j] = f(i, j);
    }
  }
}

void print_grid_to_file(const char *fileName, size_t sizeX, size_t sizeY,
                        double dx, double dy, double tab[sizeX][sizeY]) {
  FILE *writeFile = fopen(fileName, "w");
  if (writeFile == NULL) {
    perror("Fopen:");
    return;
  }
  for (size_t i = 0; i < sizeX; ++i) {
    for (size_t j = 0; j < sizeY; ++j) {
      fprintf(writeFile, "%.10e %.10e %.10e\n", (double)i * dx, (double)j * dy,
              tab[i][j]);
    }
  }
  fclose(writeFile);
}

const uint8_t versionAvailable[numSolvers] = {
    [gaussSeidel] = bareVersion | parallelVersion | tiledVersion,
    [jacobi] = bareVersion | vectorVersion | parallelVersion | tiledVersion,
    [overRelaxation] =
        bareVersion | vectorVersion | parallelVersion | tiledVersion,
};

void (*const gauss_seidel_functions[])(size_t sizeX, size_t sizeY,
                                       double t[sizeX][sizeY],
                                       double stop_criteria, double dx,
                                       double dy, double *error,
                                       size_t *num_iterations) = {
    [bareVersion] = solve_gauss_seidel,
    [vectorVersion] = NULL,
    [parallelVersion] = solve_gauss_seidel_parallel,
    [tiledVersion] = solve_gauss_seidel_parallel_tiled,
};

void (*const jacobi_functions[])(size_t sizeX, size_t sizeY,
                                 double t[sizeX][sizeY],
                                 double tNext[sizeX][sizeY],
                                 double stop_criteria, double dx, double dy,
                                 double *error_end, size_t *num_iterations) = {
    [bareVersion] = solve_jacobi,
    [vectorVersion] = solve_jacobi_vectorized,
    [parallelVersion] = solve_jacobi_parallel,
    [tiledVersion] = solve_jacobi_parallel_tiled,
};

void (*const over_relaxation_functions[])(
    size_t sizeX, size_t sizeY, double t[sizeX][sizeY],
    double tNext[sizeX][sizeY], double relaxation_factor, double stop_criteria,
    double dx, double dy, double *error_end, size_t *num_iterations) = {
    [bareVersion] = solve_over_relaxation,
    [vectorVersion] = solve_over_relaxation_vectorized,
    [parallelVersion] = solve_over_relaxation_parallel,
    [tiledVersion] = solve_over_relaxation_parallel_tiled,
};

const char *solverName[] = {
    [gaussSeidel] = "Gauss-Seidel",
    [jacobi] = "Jacobi",
    [overRelaxation] = "Over-Relaxation",
};

const char *solverVersionName[] = {
    [bareVersion] = "Bare",
    [vectorVersion] = "Vectorized",
    [parallelVersion] = "Parallel",
    [tiledVersion] = "Parallel and Tiled",
};
