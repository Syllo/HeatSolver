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

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include "heatSolver.h"
#include "initialize.h"
#include "time_measurement.h"

static struct option opt_options[] = {
    {"gauss-seidel", no_argument, 0, 'g'},
    {"jacobi", no_argument, 0, 'j'},
    {"over-relaxation", no_argument, 0, 'r'},
    {"bare", no_argument, 0, 'b'},
    {"vector", no_argument, 0, 'v'},
    {"parallel", no_argument, 0, 'p'},
    {"parallel-tiled", no_argument, 0, 't'},
    {"domain-x", required_argument, 0, 'n'},
    {"domain-y", required_argument, 0, 'm'},
    {"discrete-x", required_argument, 0, 'x'},
    {"discrete-y", required_argument, 0, 'y'},
    {"output", optional_argument, 0, 'o'},
    {"help", no_argument, 0, 'h'},
    {"max-error", required_argument, 0, 'e'},
    {"stop-iteration", required_argument, 0, 'i'},
    {0, 0, 0, 0}};

static const char options[] = ":gjrbvptn:m:x:y:o:he:i:";

static const char help_string[] =
    "Options:"
    "\n  -g --gauss-seidel    : Use Gauss-Seidel iterative solver"
    "\n  -j --jacobi          : Use Jacobi iterative solver"
    "\n  -r --over-relaxation : Use over-relaxation iterative solver"
    "\n  -b --bare            : Use non-optimized version of the solvers"
    "\n  -v --vector          : Use vectorized version of the solvers"
    "\n  -p --parallel        : Use parallel version of the solvers"
    "\n  -t --parallel-tiled  : Use parallel and tiled version of the solvers"
    "\n  -m --domain-x        : Set the size of the domain (e.g. 10.)"
    "\n  -n --domain-y        : Set the size of the domain"
    "\n  -x --discrete-x      : Set the discretization of the domain (e.g. "
    "1000)"
    "\n  -y --discrete-y      : Set the discretization of the domain"
    "\n  -o --output          : Select the output file name"
    "\n  -h --help            : Print this help"
    "\n  -e --max-error       : Set the halt condition of the solver (e.g. "
    "1e-3)'"
    "\n  -i --stop-iteration  : Set the halt condition of the solver to a "
    "specific iteration of the solver";

#define originalNx 500
#define originalNy 500
#define originalLx 1.
#define originalLy 1.
#define originalError 1e-2

int main(int argc, char **argv) {
  size_t Nx = originalNx;
  size_t Ny = originalNy;
  double Lx = originalLx;
  double Ly = originalLy;
  double error_criteria = originalError;
  size_t num_iterations = 0;

  enum solverType solver_type = gaussSeidel;
  enum solverVersion solver_version = bareVersion;
  char *output_filename = NULL;

  while (true) {
    int sscanf_return;
    int optchar = getopt_long(argc, argv, options, opt_options, NULL);
    if (optchar == -1)
      break;
    switch (optchar) {
    case 'g':
      solver_type = gaussSeidel;
      break;
    case 'j':
      solver_type = jacobi;
      break;
    case 'r':
      solver_type = overRelaxation;
      break;
    case 'b':
      solver_version = bareVersion;
      break;
    case 'v':
      solver_version = vectorVersion;
      break;
    case 'p':
      solver_version = parallelVersion;
      break;
    case 't':
      solver_version = tiledVersion;
      break;
    case 'n':
      sscanf_return = sscanf(optarg, "%lf", &Lx);
      if (sscanf_return == EOF || sscanf_return == 0 || Lx < 0.) {
        fprintf(stderr,
                "Please enter a positive floating point number for the domain "
                "size instead of \"-%c %s\"\n",
                optchar, optarg);
        Lx = originalLx;
      }
      break;
    case 'e':
      sscanf_return = sscanf(optarg, "%lf", &error_criteria);
      if (sscanf_return == EOF || sscanf_return == 0 || error_criteria < 0.) {
        fprintf(
            stderr,
            "Please enter a positive floating point number for the max error "
            "instead of \"-%c %s\"\n",
            optchar, optarg);
        error_criteria = originalError;
      }
      break;
    case 'm':
      sscanf_return = sscanf(optarg, "%lf", &Ly);
      if (sscanf_return == EOF || sscanf_return == 0 || Ly < 0.) {
        fprintf(stderr,
                "Please enter a positive floating point number for the domain "
                "size instead of \"-%c %s\"\n",
                optchar, optarg);
        Ly = originalLy;
      }
      break;
    case 'i':
      sscanf_return = sscanf(optarg, "%zu", &num_iterations);
      if (sscanf_return == EOF || sscanf_return == 0) {
        fprintf(stderr,
                "Please enter a positive integer for the stop iteration number"
                "size instead of \"-%c %s\"\n",
                optchar, optarg);
        num_iterations = 0;
      }
      break;
    case 'x':
      sscanf_return = sscanf(optarg, "%zu", &Nx);
      if (sscanf_return == EOF || sscanf_return == 0) {
        fprintf(stderr,
                "Please enter a positive floating point number for the domain "
                "size instead of \"-%c %s\"\n",
                optchar, optarg);
        Nx = originalNx;
      }
      break;
    case 'y':
      sscanf_return = sscanf(optarg, "%zu", &Ny);
      if (sscanf_return == EOF || sscanf_return == 0) {
        fprintf(stderr,
                "Please enter a positive floating point number for the domain "
                "size instead of \"-%c %s\"\n",
                optchar, optarg);
        Ny = originalNy;
      }
      break;
    case 'o':
      if (optarg[0] != '-')
        output_filename = optarg;
      else {
        output_filename = "gridData.dat";
        optind--;
      }
      break;
    case 'h':
      printf("Usage: %s <options>\n%s\n", argv[0], help_string);
      return EXIT_SUCCESS;
      break;
    case ':':
      if (optopt == 'o') {
        output_filename = "gridData.dat";
      } else {
        fprintf(stderr, "Option %c requires an argument\n", optopt);
        exit(EXIT_FAILURE);
      }
      break;
    default:
      fprintf(stderr, "Unrecognized option %c\n", optopt);
      exit(EXIT_FAILURE);
      break;
    }
  }
  if (!has_version(solver_type, solver_version)) {
    fprintf(stderr, "The solver \"%s\" does not support the \"%s\" version\n",
            solverName[solver_type], solverVersionName[solver_version]);
    exit(EXIT_FAILURE);
  }

  double Dx = Lx / (double)Nx;
  double Dy = Ly / (double)Ny;
  fprintf(stdout, "Domain (%e , %e) points (%zu,%zu) (Dx: %e , Dy: %e)\n", Lx,
          Ly, Nx, Ny, Dx, Dy);
  double (*init_fun)(size_t, size_t) =
      initialize_borders_fun(Nx, Ny, 400., 800., 600., 900.);
  double(*heatVals)[Ny] = malloc(sizeof(double[Nx][Ny]));
  if (heatVals == NULL) {
    perror("Malloc:");
    return EXIT_FAILURE;
  }
  init_heat_values(Nx, Ny, heatVals, init_fun);
  double(*heatValsTmp)[Ny] = NULL;

  double error = 0.;
  time_measure start, end;

  fprintf(stdout, "Running the solver \"%s\" with \"%s\" version\n",
          solverName[solver_type], solverVersionName[solver_version]);

  switch (solver_type) {
  case gaussSeidel:
    get_current_time(&start);
    gauss_seidel_functions[solver_version](Nx, Ny, heatVals, error_criteria, Dx,
                                           Dy, &error, &num_iterations);
    get_current_time(&end);
    break;
  case jacobi:
    heatValsTmp = malloc(sizeof(double[Nx][Ny]));
    memcpy(heatValsTmp, heatVals, sizeof(double[Nx][Ny]));
    get_current_time(&start);
    jacobi_functions[solver_version](Nx, Ny, heatVals, heatValsTmp,
                                     error_criteria, Dx, Dy, &error,
                                     &num_iterations);
    get_current_time(&end);
    break;
  case overRelaxation: {
    double wOpt = 2. / (1. + sin(M_PI * Dx));
    fprintf(stdout, "Relaxation factor %f\n", wOpt);
    heatValsTmp = malloc(sizeof(double[Nx][Ny]));
    memcpy(heatValsTmp, heatVals, sizeof(double[Nx][Ny]));
    get_current_time(&start);
    over_relaxation_functions[solver_version](Nx, Ny, heatVals, heatValsTmp,
                                              wOpt, error_criteria, Dx, Dy,
                                              &error, &num_iterations);
    get_current_time(&end);
  } break;
  default:
    get_current_time(&start);
    get_current_time(&end);
    break;
  }

  fprintf(stdout, "Kernel time %.4fs\n", measuring_difftime(start, end));
  fprintf(stdout, "Error at the end: %.4e\n", error);
  fprintf(stdout, "Num iterations: %zu\n", num_iterations);
  if (output_filename) {
    double(*heatValsEnd)[Ny] = num_iterations % 2 == 0 ? heatVals : heatValsTmp;
    if (heatValsEnd == NULL)
      heatValsEnd = heatVals;
    print_grid_to_file(output_filename, Nx, Ny, Dx, Dy, heatValsEnd);
  }

  free(heatVals);
  free(heatValsTmp);
  return EXIT_SUCCESS;
}
