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

#include "initialize.h"

static size_t initialize_sizeX;
static size_t initialize_sizeY;
static double initialize_left_val;
static double initialize_right_val;
static double initialize_top_val;
static double initialize_bottom_val;

double (*initialize_borders_fun(size_t sizeX, size_t sizeY, double left_val,
                                double right_val, double top_val,
                                double bottom_val))(size_t, size_t) {
  initialize_left_val = left_val;
  initialize_right_val = right_val;
  initialize_top_val = top_val;
  initialize_bottom_val = bottom_val;
  initialize_sizeX = sizeX;
  initialize_sizeY = sizeY;
  return initialize_border_index;
}

double initialize_border_index(size_t indeX, size_t indexY) {
  if (indeX == 0)
    return initialize_bottom_val;
  if (indeX == initialize_sizeX - 1)
    return initialize_top_val;
  if (indexY == 0)
    return initialize_left_val;
  if (indexY == initialize_sizeY - 1)
    return initialize_right_val;
  return 1.;
}
