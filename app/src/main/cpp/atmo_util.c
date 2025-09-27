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

#include <math.h>
#include <stdlib.h>

#define LOG_TAG "atmo"
#include "libcc/cc_log.h"
#include "atmo_util.h"

uint32_t atmo_clampu(uint32_t v,
                     uint32_t min, uint32_t max)
{
	ASSERT(min < max);

	if(v < min)
	{
		v = min;
	}
	else if(v > max)
	{
		v = max;
	}
	return v;
}

double atmo_clampd(double v, double min, double max)
{
	ASSERT(min < max);

	if(v < min)
	{
		v = min;
	}
	else if(v > max)
	{
		v = max;
	}
	return v;
}

double atmo_signd(double x)
{
	if(x < 0.0)
	{
		return -1.0;
	}

	return 1.0;
}

double atmo_maxd(double a, double b)
{
	return (a > b) ? a : b;
}

double atmo_deg2rad(double x)
{
	return x*M_PI/180.0;
}
