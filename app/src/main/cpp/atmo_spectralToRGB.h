/*
 * Copyright (c) 2025 Jeff Boody
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

#ifndef atmo_spectralToRGB_H
#define atmo_spectralToRGB_H

#include "libcc/math/cc_mat3d.h"
#include "libcc/math/cc_vec3d.h"

#define ATMO_SPECTRAL_TO_RGB_MIN 360
#define ATMO_SPECTRAL_TO_RGB_MAX 830

#define ATMO_SPECTRAL_TO_RGB_NORMALIZE_NONE 0
#define ATMO_SPECTRAL_TO_RGB_NORMALIZE_SUM  1
#define ATMO_SPECTRAL_TO_RGB_NORMALIZE_PEAK 2

void atmo_spectralToRGB_getXYZ(int lambda,
                               cc_vec3d_t* xyz);
void atmo_spectralToRGB_getRGB(int normalize, int lambda,
                               cc_vec3d_t* rgb);
void atmo_spectrlToRGB_getMinv(cc_mat3d_t* Minv);

#endif
