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

#ifndef _HEAT_SOLVER_H_
#define _HEAT_SOLVER_H_

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "gauss-seidel.h"
#include "jacobi.h"
#include "over-relaxation.h"

void init_heat_values(size_t sizeX, size_t sizeY, double t[sizeX][sizeY],
                      double (*f)(size_t, size_t));

void print_grid_to_file(const char *fileName, size_t sizeX, size_t sizeY,
                        double dx, double dy, double tab[sizeX][sizeY]);

enum solverType {
  gaussSeidel = 0,
  jacobi,
  overRelaxation,
  numSolvers,
};

enum solverVersion {
  bareVersion = 1,
  vectorVersion = 2,
  parallelVersion = 4,
  tiledVersion = 8,
};

extern const char *solverName[];
extern const char *solverVersionName[];

extern const uint8_t versionAvailable[numSolvers];

extern void (*const gauss_seidel_functions[])(size_t sizeX, size_t sizeY,
                                              double t[sizeX][sizeY],
                                              double stop_criteria, double dx,
                                              double dy, double *error,
                                              size_t *num_iterations);

extern void (*const jacobi_functions[])(size_t sizeX, size_t sizeY,
                                        double t[sizeX][sizeY],
                                        double tNext[sizeX][sizeY],
                                        double stop_criteria, double dx,
                                        double dy, double *error_end,
                                        size_t *num_iterations);

extern void (*const over_relaxation_functions[])(
    size_t sizeX, size_t sizeY, double t[sizeX][sizeY],
    double tNext[sizeX][sizeY], double relaxation_factor, double stop_criteria,
    double dx, double dy, double *error_end, size_t *num_iterations);

static inline bool has_version(enum solverType type,
                               enum solverVersion version) {
  return !!(versionAvailable[type] & version);
}

#endif // _HEAT_SOLVER_H_
