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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LOG_TAG "atmo"
#include "libcc/jsmn/cc_jsmnStream.h"
#include "libcc/math/cc_pow2n.h"
#include "libcc/math/cc_ray3d.h"
#include "libcc/math/cc_sphere3d.h"
#include "libcc/math/cc_vec3d.h"
#include "libcc/math/cc_vec4d.h"
#include "libcc/math/cc_vec4f.h"
#include "libcc/cc_log.h"
#include "libcc/cc_memory.h"
#include "texgz/texgz_tex.h"
#include "texgz/texgz_png.h"
#include "atmo_solver.h"
#include "atmo_spectralIrradiance.h"
#include "atmo_spectralToRGB.h"

// radius of the planet and atmospheric boundary
// use doubles for radius for numerical stability
#define ATMO_RP 6360000.0
#define ATMO_RA 6460000.0

// Rayleigh and Mie scale heights
//
// Clear day: Mie scale height: 1200 m
//            For clear days, we use the standard values.
// Hazy day: Mie scale height: 1500 m to 2000 m
//           For hazy days, we increase the Mie scale height
//           to account for a greater distribution of
//           aerosols in the atmosphere.
// Foggy day: Mie scale height: 500 m to 800 m
//            For foggy days, we decrease the Mie scale
//            height because fog particles are concentrated
//            closer to the ground.
#define ATMO_DENSITY_SCALE_HEIGHT_RAYLEIGH 8000.0
#define ATMO_DENSITY_SCALE_HEIGHT_MIE      1200.0

// also requires changes to shaders
#define ATMO_USE_MODIFIED_RAYLEIGH_PHASE_FUNCTION

// Mie asymmetry factor
//
// Clear day: 0.7 to 0.8
//            Clear conditions have less scattering, so a
//            higher positive value indicates more forward
//            scattering.
// Hazy day: 0.3 to 0.5
//           Haze involves smaller particles, leading to
//           more isotropic scattering.
// Foggy day: -0.1 to 0.1
//            Fog consists of larger water droplets,
//            resulting in more uniform scattering in all
//            directions.
#define ATMO_PHASE_G_MIE 0.8

// Rayleigh and Mie scattering coefficients represent the
// probability of light being scattered as it travels
// through a medium.
//
// Clear day: Mie scattering coefficient: 2e-6
//            For clear days, we use the standard values.
// Hazy day:  Mie scattering coefficient: 4e-6 to 8e-6
//            For hazy days, we increase the Mie scattering
//            coefficient to account for a greater
//            concentration of aerosols in the atmosphere.
// Foggy day: Mie scattering coefficient: 2e-5 to 5e-5
//            For foggy days, we significantly increase the
//            Mie scattering coefficient because fog
//            consists of larger water droplets that scatter
//            light more strongly.
#define ATMO_BETA_S_R_RAYLEIGH 5.802e-6
#define ATMO_BETA_S_G_RAYLEIGH 13.558e-6
#define ATMO_BETA_S_B_RAYLEIGH 33.1e-6
#define ATMO_BETA_S_MIE        3.996e-6
#define ATMO_BETA_A_MIE        4.40e-6
#define ATMO_BETA_A_R_OZONE    0.65e-6
#define ATMO_BETA_A_G_OZONE    1.881e-6
#define ATMO_BETA_A_B_OZONE    0.085e-6

// transmittance numerical integration steps
#define ATMO_TRANSMITTANCE_STEPS 30

#define ATMO_STEP_THRESH 10.0

// gathering direction numerical integration steps of the
// spherical coordinate system for M (theta) and N (phi)
#define ATMO_GATHER_M_STEPS (180/30)
#define ATMO_GATHER_N_STEPS (360/30)

// multiple scattering events
#define ATMO_K 5

// scattering texture size
#define ATMO_TEXTURE_WIDTH  32
#define ATMO_TEXTURE_HEIGHT 256
#define ATMO_TEXTURE_DEPTH  32

// height parameterization
#define ATMO_PARAM_HEIGHT_LINEAR 0
#define ATMO_PARAM_HEIGHT_POWER  1

// view-zenith angle parameterization
// WARNING: ATMO_PARAM_PHI_BODARE is buggy
#define ATMO_PARAM_PHI_LINEAR         0
#define ATMO_PARAM_PHI_POWER          1
#define ATMO_PARAM_PHI_BODARE         2
#define ATMO_PARAM_PHI_WEIGHTED_POWER 3

// weighted power parameters
#define ATMO_PARAM_PHI_WEIGHTED_POWER_PU  2.0
#define ATMO_PARAM_PHI_WEIGHTED_POWER_PL  2.0
#define ATMO_PARAM_PHI_WEIGHTED_POWER_PS  2.0
#define ATMO_PARAM_PHI_WEIGHTED_POWER_WL1 (20.0/32.0)
#define ATMO_PARAM_PHI_WEIGHTED_POWER_WS0 (4.0/32.0)
#define ATMO_PARAM_PHI_WEIGHTED_POWER_WS1 (4.0/32.0)

// sun-zenith angle parameterization
#define ATMO_PARAM_DELTA_LINEAR 0
#define ATMO_PARAM_DELTA_POWER  1
#define ATMO_PARAM_DELTA_BODARE 2

// select parameterization
// requires corresponding change in atmo_solver.c
#define ATMO_PARAM_HEIGHT ATMO_PARAM_HEIGHT_POWER
#define ATMO_PARAM_PHI    ATMO_PARAM_PHI_WEIGHTED_POWER
#define ATMO_PARAM_DELTA  ATMO_PARAM_DELTA_POWER

// sampling modes
#define ATMO_SAMPLE_MODE_NEAREST 0
#define ATMO_SAMPLE_MODE_LINEAR  1

// select sampling mode
#define ATMO_SAMPLE_MODE ATMO_SAMPLE_MODE_LINEAR

// experimental change which attempts to adjust the
// attenuation of multiple scattering passes which have
// much lower intensity than expected
// #define ATMO_ADJUST_ATTENUATION

/***********************************************************
* private - compute                                        *
***********************************************************/

static void
atmo_solver_computeII(atmo_solverParam_t* param)
{
	ASSERT(param);

	cc_mat3d_t Minv;
	atmo_spectrlToRGB_getMinv(&Minv);

	int min = ATMO_SPECTRAL_TO_RGB_MIN;
	if(ATMO_SPECTRAL_IRRADIANCE_MIN > min)
	{
		min = ATMO_SPECTRAL_IRRADIANCE_MIN;
	}

	int max = ATMO_SPECTRAL_TO_RGB_MAX;
	if(ATMO_SPECTRAL_IRRADIANCE_MAX < max)
	{
		max = ATMO_SPECTRAL_IRRADIANCE_MAX;
	}

	int i;
	cc_vec3d_t xyz0;
	cc_vec3d_t xyz1;

	// normalization of xyz essentially disables HDR color
	#ifdef ATMO_NORMALIZE_XYZ
	atmo_spectralToRGB_getXYZ(min, &xyz0);
	double y0 = xyz0.y;
	double y1;

	// compute normalization factor for XYZ
	double si0 = atmo_spectralIrradiance_get((double) min);
	double si1;
	double K = 0.0;
	for(i = min + 1; i <= max; ++i)
	{
		// Solar irradiance values (W / (m^2 * nm)).
		si1 = atmo_spectralIrradiance_get((double) i);
		atmo_spectralToRGB_getXYZ(i, &xyz1);
		y1 = xyz1.y;

		// apply trapezoidal rule
		K += 0.5*(si0*y0 + si1*y1);

		si0 = si1;
		y0  = y1;
	}
	LOGI("K=%lf", K);
	#endif

	// integrate spectral irradiance and xyz
	double si;
	cc_vec3d_t xyz  = { 0 };
	cc_vec3d_t tmp;
	for(i = min; i <= max; ++i)
	{
		// Solar irradiance values (W / (m^2 * nm)).
		si = atmo_spectralIrradiance_get((double) i);

		atmo_spectralToRGB_getXYZ(i, &xyz1);
		cc_vec3d_muls(&xyz1, si);

		if(i > min)
		{
			// apply trapezoidal rule
			cc_vec3d_addv_copy(&xyz0, &xyz1, &tmp);
			cc_vec3d_muls(&tmp, 0.5);
			cc_vec3d_addv(&xyz, &tmp);
		}

		cc_vec3d_copy(&xyz1, &xyz0);
	}

	#ifdef ATMO_NORMALIZE_XYZ
	cc_vec3d_muls(&xyz, 1.0/K);
	#endif

	cc_vec3d_t II;
	cc_mat3d_mulv_copy(&Minv, &xyz, &II);

	param->II_r = II.x;
	param->II_g = II.y;
	param->II_b = II.z;

	LOGI("II: r=%f, g=%f, b=%f",
	     param->II_r, param->II_g, param->II_b);
}

static uint32_t
atmo_clampu(uint32_t v, uint32_t min, uint32_t max)
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

static double
atmo_clampd(double v, double min, double max)
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

static double getUHeight(atmo_solverParam_t* param, double h)
{
	ASSERT(param);

	double u;
	#if ATMO_PARAM_HEIGHT == ATMO_PARAM_HEIGHT_POWER
	u = pow(h/(param->Ra - param->Rp), 1.0/2.0);
	#else
	u = h/(param->Ra - param->Rp);
	#endif

	return atmo_clampd(u, 0.0, 1.0);
}

static double getHeightU(atmo_solverParam_t* param, double u)
{
	ASSERT(param);

	double h;
	#if ATMO_PARAM_HEIGHT == ATMO_PARAM_HEIGHT_POWER
	h = pow(u, 2.0)*(param->Ra - param->Rp);
	#else
	h = u*(param->Ra - param->Rp);
	#endif

	double Ha = param->Ra - param->Rp;

	return atmo_clampd(h, 0.0, Ha);
}

static double
getHeightP(atmo_solverParam_t* param, cc_vec3d_t* P)
{
	ASSERT(param);
	ASSERT(P);

	double h  = cc_vec3d_mag(P) - param->Rp;
	double Ha = param->Ra - param->Rp;

	return atmo_clampd(h, 0.0, Ha);
}

static void
getZenithP(const cc_vec3d_t* P, cc_vec3d_t* Zenith)
{
	ASSERT(P);
	ASSERT(Zenith);

	cc_vec3d_normalize_copy(P, Zenith);
}

static double atmo_signd(double x)
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

static double
getCosPhiV(atmo_solverParam_t* param, double h, double u,
           double v)
{
	ASSERT(param);

	double Rp = param->Rp;

	double cos_phi;
	#if ATMO_PARAM_PHI == ATMO_PARAM_PHI_POWER
	double vv = 2.0*v - 1.0;
	cos_phi = atmo_signd(vv)*pow(fabs(vv), 3.0);
	#elif ATMO_PARAM_PHI == ATMO_PARAM_PHI_BODARE
	double ch = -sqrt(h*(2.0*Rp + h))/(Rp + h);
	if(v > 0.5)
	{
		cos_phi = ch + pow(v - 0.5, 5.0)*(1.0 - ch);
	}
	else
	{
		cos_phi = ch - pow(v, 5.0)*(1.0 + ch);
	}
	#elif ATMO_PARAM_PHI == ATMO_PARAM_PHI_WEIGHTED_POWER
	double hypH      = Rp + h;
	double oppH      = Rp;
	double adjH      = sqrt(hypH*hypH - oppH*oppH);
	double cos_phi_H = -adjH/hypH;

	double PU  = ATMO_PARAM_PHI_WEIGHTED_POWER_PU;
	double PL  = ATMO_PARAM_PHI_WEIGHTED_POWER_PL;
	double PS  = ATMO_PARAM_PHI_WEIGHTED_POWER_PS;
	double WL1 = ATMO_PARAM_PHI_WEIGHTED_POWER_WL1;
	double WS0 = ATMO_PARAM_PHI_WEIGHTED_POWER_WS0;
	double WS1 = ATMO_PARAM_PHI_WEIGHTED_POWER_WS1;
	double WS  = WS1*u + WS0;
	double WL  = WL1*u;
	double WU  = 1.0 - WL - WS;
	double epsilon = 0.00001;
	if(v >= WL + WS)
	{
		cos_phi = pow((v - (WL + WS))/WU, PU);
	}
	else if(v >= WS)
	{
		cos_phi = -cos_phi_H*pow((v - WS)/
		                         atmo_maxd(WL, epsilon), PL) +
		           cos_phi_H;
	}
	else
	{
		cos_phi = (-1.0 - cos_phi_H)*pow(1.0 - v/WS, PS) +
		          cos_phi_H;
	}
	#else
	cos_phi = 2.0*v - 1.0;
	#endif

	return atmo_clampd(cos_phi, -1.0, 1.0);
}

static double
getVCosPhi(atmo_solverParam_t* param, double h,
           double cos_phi, double u)
{
	ASSERT(param);

	double Rp = param->Rp;

	double v;
	#if ATMO_PARAM_PHI == ATMO_PARAM_PHI_POWER
	v = 0.5*(1.0 + atmo_signd(cos_phi)*
	               pow(fabs(cos_phi), 1.0/3.0));
	#elif ATMO_PARAM_PHI == ATMO_PARAM_PHI_BODARE
	double ch = -sqrt(h*(2.0*Rp + h))/(Rp + h);
	if(cos_phi > ch)
	{
		v = 0.5*pow((cos_phi - ch)/(1.0 - ch), 0.2) + 0.5;
	}
	else
	{
		v = 0.5*pow((ch - cos_phi)/(1.0 + ch), 0.2);
	}
	#elif ATMO_PARAM_PHI == ATMO_PARAM_PHI_WEIGHTED_POWER
	double hypH      = Rp + h;
	double oppH      = Rp;
	double adjH      = sqrt(hypH*hypH - oppH*oppH);
	double cos_phi_H = -adjH/hypH;

	double PU      = ATMO_PARAM_PHI_WEIGHTED_POWER_PU;
	double PL      = ATMO_PARAM_PHI_WEIGHTED_POWER_PL;
	double PS      = ATMO_PARAM_PHI_WEIGHTED_POWER_PS;
	double WL1     = ATMO_PARAM_PHI_WEIGHTED_POWER_WL1;
	double WS0     = ATMO_PARAM_PHI_WEIGHTED_POWER_WS0;
	double WS1     = ATMO_PARAM_PHI_WEIGHTED_POWER_WS1;
	double WS      = WS1*u + WS0;
	double WL      = WL1*u;
	double WU      = 1.0 - WL - WS;
	double epsilon = 0.00001;
	if(cos_phi >= 0.0)
	{
		v = WU*pow(cos_phi, 1.0/PU) + (WL + WS);
	}
	else if(cos_phi >= cos_phi_H)
	{
		v = WL*pow((cos_phi - cos_phi_H)/
		           atmo_maxd(-cos_phi_H, epsilon),
		           1.0/PL) + WS;
	}
	else
	{
		v = WS*(1.0 - pow((cos_phi - cos_phi_H)/
		                  (-1.0 - cos_phi_H), 1.0/PS));
	}
	#else
	v = (cos_phi + 1.0)/2.0;
	#endif

	return atmo_clampd(v, 0.0, 1.0);
}

static double getCosDeltaW(double w)
{
	double cos_delta;
	#if ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_POWER
	double ww = 2.0*w - 1.0;
	cos_delta = atmo_signd(ww)*pow(fabs(ww), 3.0);
	#elif ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_BODARE
	cos_delta = tan((2.0*w - 1.0 + 0.26)*0.75)/
	            tan(1.26*0.75);
	#else
	cos_delta = 2.0*w - 1.0;
	#endif

	return atmo_clampd(cos_delta, -1.0, 1.0);
}

static double getWCosDelta(double cos_delta)
{
	double w;
	#if ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_POWER
	w = 0.5*(1.0 + atmo_signd(cos_delta)*
	               pow(fabs(cos_delta), 1.0/3.0));
	#elif ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_BODARE
	w =  0.5*(atan(atmo_maxd(cos_delta, -0.1975)*
	               tan(1.26*1.1))/1.1 +
	           (1.0 - 0.26));
	#else
	w = (cos_delta + 1.0)/2.0;
	#endif

	return atmo_clampd(w, 0.0, 1.0);
}

static cc_vec4d_t*
getDataK(atmo_solverParam_t* param,
         uint32_t k, cc_vec4d_t* data)
{
	ASSERT(param);
	ASSERT(data);

	// k is base-1
	uint32_t i   = k - 1;
	uint32_t w   = param->texture_width;
	uint32_t h   = param->texture_height;
	uint32_t d   = param->texture_depth;
	uint32_t idx = i*w*h*d;

	return &data[idx];
}

static cc_vec4f_t*
getDataKf(atmo_solverParam_t* param,
          uint32_t k, cc_vec4f_t* data)
{
	ASSERT(param);
	ASSERT(data);

	// k is base-1
	uint32_t i   = k - 1;
	uint32_t w   = param->texture_width;
	uint32_t h   = param->texture_height;
	uint32_t d   = param->texture_depth;
	uint32_t idx = i*w*h*d;

	return &data[idx];
}

static void
getData(atmo_solverParam_t* param,
        uint32_t k, uint32_t x,
        uint32_t y, uint32_t z,
        cc_vec4d_t* data, cc_vec4d_t* val)
{
	ASSERT(param);
	ASSERT(data);
	ASSERT(val);

	// k is base-1
	uint32_t i   = k - 1;
	uint32_t w   = param->texture_width;
	uint32_t h   = param->texture_height;
	uint32_t d   = param->texture_depth;
	uint32_t idx = x + y*w + z*w*h + i*w*h*d;

	val->r = data[idx].r;
	val->g = data[idx].g;
	val->b = data[idx].b;
	val->a = data[idx].a;
}

static void
setData(atmo_solverParam_t* param,
        uint32_t k, uint32_t x,
        uint32_t y, uint32_t z,
        cc_vec4d_t* data, cc_vec4d_t* val)
{
	ASSERT(param);
	ASSERT(data);
	ASSERT(val);

	// k is base-1
	uint32_t i   = k - 1;
	uint32_t w   = param->texture_width;
	uint32_t h   = param->texture_height;
	uint32_t d   = param->texture_depth;
	uint32_t idx = x + y*w + z*w*h + i*w*h*d;

	data[idx].r = val->r;
	data[idx].g = val->g;
	data[idx].b = val->b;
	data[idx].a = val->a;
}

// Rayleigh density function
static double densityR(atmo_solverParam_t* param, double h)
{
	ASSERT(param);

	return exp(-h/param->density_scale_height_rayleigh);
}

// Mie density function
static double densityM(atmo_solverParam_t* param, double h)
{
	ASSERT(param);

	return exp(-h/param->density_scale_height_mie);
}

// Ozone density function
static double densityO(atmo_solverParam_t* param, double h)
{
	ASSERT(param);

	double p = 1.0 - fabs(h - 25000.0)/15000.0;
	if(p < 0.0)
	{
		return 0.0;
	}

	return p;
}

// Rayleigh/Mie optical depth
static void
opticalDepth(atmo_solverParam_t* param, cc_vec3d_t* P1,
             cc_vec3d_t* P2, cc_vec3d_t* out)
{
	ASSERT(param);
	ASSERT(P1);
	ASSERT(P2);
	ASSERT(out);

	// initialize integration
	cc_vec3d_t P;
	cc_vec3d_t step;
	cc_vec3d_copy(P1, &P);
	cc_vec3d_subv_copy(P2, P1, &step);
	cc_vec3d_muls(&step, 1.0/param->transmittance_steps);
	double ds  = cc_vec3d_mag(&step);
	double h   = getHeightP(param, &P);
	double pR0 = densityR(param, h);
	double pM0 = densityM(param, h);
	double pO0 = densityO(param, h);

	if(ds < ATMO_STEP_THRESH)
	{
		out->x = 0.0;
		out->y = 0.0;
		out->z = 0.0;
		return;
	}

	// integrate optical depth
	int    i;
	double pR1;
	double pM1;
	double pO1;
	double tR = 0.0;
	double tM = 0.0;
	double tO = 0.0;
	for(i = 0; i < param->transmittance_steps; ++i)
	{
		cc_vec3d_addv(&P, &step);

		h   = getHeightP(param, &P);
		pR1 = densityR(param, h);
		pM1 = densityM(param, h);
		pO1 = densityO(param, h);

		// apply trapesoidal rule
		tR += 0.5*(pR0 + pR1)*ds;
		tM += 0.5*(pM0 + pM1)*ds;
		tO += 0.5*(pO0 + pO1)*ds;

		pR0 = pR1;
		pM0 = pM1;
		pO0 = pO1;
	}

	// extinction coefficients (scattering + absorption)
	cc_vec3d_t bR =
	{
		.x = param->beta_s_r_rayleigh,
		.y = param->beta_s_g_rayleigh,
		.z = param->beta_s_b_rayleigh,
	};

	double bM = param->beta_s_mie + param->beta_a_mie;

	cc_vec3d_t bO =
	{
		.x = param->beta_a_r_ozone,
		.y = param->beta_a_g_ozone,
		.z = param->beta_a_b_ozone,
	};

	// apply Rayleigh/Mie scattering coefficient
	out->x = bR.x*tR + bM*tM + bO.x*tO;
	out->y = bR.y*tR + bM*tM + bO.y*tO;
	out->z = bR.z*tR + bM*tM + bO.z*tO;
}

// compute the points Pa and Pb on the viewing vector
static int
computePaPb(atmo_solverParam_t* param, cc_vec3d_t* P0,
            cc_vec3d_t* V, double Ro, cc_vec3d_t* Pa,
            cc_vec3d_t* Pb)
{
	ASSERT(param);
	ASSERT(P0);
	ASSERT(V);
	ASSERT(Pa);
	ASSERT(Pb);

	// initialize spheres
	// optionally include a radius offset to ensure that the
	// viewing vector does not intersect at P0
	cc_sphere3d_t sphereP;
	cc_sphere3d_t sphereA;
	cc_sphere3d_load(&sphereP, 0.0, 0.0, 0.0, param->Rp - Ro);
	cc_sphere3d_load(&sphereA, 0.0, 0.0, 0.0, param->Ra + Ro);

	// initialize ray
	cc_ray3d_t ray;
	cc_ray3d_load(&ray, P0->x, P0->y, P0->z, V->x, V->y, V->z);

	// intersect ray
	int   countA;
	int   countP;
	double nearA = 0.0;
	double farA  = 0.0;
	double nearP = 0.0;
	double farP  = 0.0;
	countA = cc_ray3d_intersect(&ray, &sphereA,
	                            &nearA, &farA);
	countP = cc_ray3d_intersect(&ray, &sphereP,
	                            &nearP, &farP);

	// check if ray intersects atmosphere
	if(countA)
	{
		cc_ray3d_getpoint(&ray, nearA, Pa);

		// check if ray intersects planet
		if(countP)
		{
			cc_ray3d_getpoint(&ray, nearP, Pb);
			return 0;
		}

		cc_ray3d_getpoint(&ray, farA, Pb);
		return 1;
	}

	cc_vec3d_copy(P0, Pa);
	cc_vec3d_copy(P0, Pb);
	return 0;
}

// compute the point Pc on the light vector
static int
computePc(atmo_solverParam_t* param, cc_vec3d_t* P,
          cc_vec3d_t* L, cc_vec3d_t* Pc)
{
	ASSERT(param);
	ASSERT(P);
	ASSERT(L);
	ASSERT(Pc);

	cc_vec3d_t Pa;
	cc_vec3d_t Sun;
	cc_vec3d_muls_copy(L, -1.0, &Sun);
	double Ro = 10.0;
	return computePaPb(param, P, &Sun, Ro, &Pa, Pc);
}

#ifdef ATMO_USE_MODIFIED_RAYLEIGH_PHASE_FUNCTION
// modified Rayleigh phase function
static double
atmo_phaseR(atmo_solverParam_t* param, double cos_theta)
{
	ASSERT(param);

	return 0.8*(1.4 + 0.5*cos_theta);
}
#else
// standard Rayleigh phase function
static double
atmo_phaseR(atmo_solverParam_t* param, double cos_theta)
{
	ASSERT(param);

	return 0.75*(1.0 + cos_theta*cos_theta);
}
#endif

// Mie phase function
static double
atmo_phaseM(atmo_solverParam_t* param, double cos_theta)
{
	ASSERT(param);

	double g  = param->phase_g_mie;
	double g2 = g*g;
	double n1 = 3.0*(1.0 - g2);
	double n2 = 1.0 + cos_theta*cos_theta;
	double d1 = 2.0*(2.0 + g2);
	double d2 = pow(1.0 + g2 - 2.0*g*cos_theta, 1.5);
	return (n1/d1)*(n2/d2);
}

// factored single-scattered intensity
static void
fIS1(atmo_solverParam_t* param, double h, double phi,
     double delta, cc_vec4d_t* fis1)
{
	ASSERT(param);
	ASSERT(fis1);

	// initialize fis1
	fis1->r = 0.0;
	fis1->g = 0.0;
	fis1->b = 0.0;
	fis1->a = 0.0;

	// canonical form of the scattering intensity
	// parameterization for P0, V and L
	cc_vec3d_t P0 =
	{
		.z = h + param->Rp,
	};
	cc_vec3d_t V =
	{
		.x = sin(phi),
		.z = cos(phi),
	};
	cc_vec3d_t L =
	{
		.x = -sin(delta),
		.z = -cos(delta),
	};

	// compute ray-sphere intersection
	// include a ray offset for the viewing vector
	// to ensure that the ray does not intersect at P0
	double Ro = 10.0;
	cc_vec3d_t Pa;
	cc_vec3d_t Pb;
	cc_vec3d_t Pc;
	computePaPb(param, &P0, &V, Ro, &Pa, &Pb);
	if(computePc(param, &Pa, &L, &Pc) == 0)
	{
		// ray intersects planet
		return;
	}

	// initialize integration
	cc_vec3d_t P;
	cc_vec3d_t step;
	cc_vec4d_t fx0;
	cc_vec4d_t fx1;
	cc_vec3d_t tPPc;
	cc_vec3d_t tPaP;
	cc_vec3d_copy(&Pa, &P);
	cc_vec3d_subv_copy(&Pb, &Pa, &step);
	cc_vec3d_muls(&step, 1.0/param->transmittance_steps);
	double ds = cc_vec3d_mag(&step);
	double pR = densityR(param, h);
	double pM = densityM(param, h);
	opticalDepth(param, &P, &Pc, &tPPc);
	opticalDepth(param, &Pa, &P, &tPaP);
	fx0.r = pR*exp(-tPPc.x -tPaP.x);
	fx0.g = pR*exp(-tPPc.y -tPaP.y);
	fx0.b = pR*exp(-tPPc.z -tPaP.z);
	fx0.a = pM*exp(-(tPPc.x + tPPc.y + tPPc.z)/3.0
	               -(tPaP.y + tPaP.y + tPaP.z)/3.0);

	if(ds < ATMO_STEP_THRESH)
	{
		return;
	}

	// integrate factored single-scattered intensity
	int i;
	for(i = 0; i < param->transmittance_steps; ++i)
	{
		cc_vec3d_addv(&P, &step);

		if(computePc(param, &P, &L, &Pc) == 0)
		{
			// PPc intersects planet
			fx0.r = 0.0;
			fx0.g = 0.0;
			fx0.b = 0.0;
			fx0.a = 0.0;
			continue;
		}

		h  = getHeightP(param, &P);
		pR = densityR(param, h);
		pM = densityM(param, h);
		opticalDepth(param, &P, &Pc, &tPPc);
		opticalDepth(param, &Pa, &P, &tPaP);
		fx1.r = pR*exp(-tPPc.x -tPaP.x);
		fx1.g = pR*exp(-tPPc.y -tPaP.y);
		fx1.b = pR*exp(-tPPc.z -tPaP.z);
		fx1.a = pM*exp(-(tPPc.x + tPPc.y + tPPc.z)/3.0
		               -(tPaP.x + tPaP.y + tPaP.z)/3.0);

		// apply trapesoidal rule
		fis1->r += 0.5*(fx1.r + fx0.r)*ds;
		fis1->g += 0.5*(fx1.g + fx0.g)*ds;
		fis1->b += 0.5*(fx1.b + fx0.b)*ds;
		fis1->a += 0.5*(fx1.a + fx0.a)*ds;

		cc_vec4d_copy(&fx1, &fx0);
	}

	// apply Rayleigh/Mie scattering coefficient
	fis1->r *= param->beta_s_r_rayleigh/(4.0*M_PI);
	fis1->g *= param->beta_s_g_rayleigh/(4.0*M_PI);
	fis1->b *= param->beta_s_b_rayleigh/(4.0*M_PI);
	fis1->a *= param->beta_s_mie/(4.0*M_PI);
}

static void
fISk_sample(atmo_solverParam_t* param, uint32_t k,
            cc_vec3d_t* P, cc_vec3d_t* V, cc_vec3d_t* L,
            cc_vec4d_t* data, cc_vec4d_t* fisk)
{
	ASSERT(param);
	ASSERT(P);
	ASSERT(V);
	ASSERT(L);
	ASSERT(data);
	ASSERT(fisk);

	// initialize fisk
	fisk->r = 0.0;
	fisk->g = 0.0;
	fisk->b = 0.0;
	fisk->a = 0.0;

	cc_vec3d_t Zenith;
	getZenithP(P, &Zenith);

	cc_vec3d_t Sun;
	cc_vec3d_muls_copy(L, -1.0, &Sun);

	// compute u,v,w
	double h         = getHeightP(param, P);
	double cos_phi   = cc_vec3d_dot(&Zenith, V);
	double cos_delta = cc_vec3d_dot(&Zenith, &Sun);
	double u         = getUHeight(param, h);
	double v         = getVCosPhi(param, h, cos_phi, u);
	double w         = getWCosDelta(cos_delta);
	double width1    = (double) (param->texture_width  - 1);
	double height1   = (double) (param->texture_height - 1);
	double depth1    = (double) (param->texture_depth  - 1);

	// compute x0,y0,z0
	uint32_t x0;
	uint32_t y0;
	uint32_t z0;
	x0 = atmo_clampu((uint32_t) (u*width1),
	                 0, param->texture_width  - 1);
	y0 = atmo_clampu((uint32_t) (v*height1),
	                 0, param->texture_height - 1);
	z0 = atmo_clampu((uint32_t) (w*depth1),
	                 0, param->texture_depth  - 1);

	#if ATMO_SAMPLE_MODE == ATMO_SAMPLE_MODE_LINEAR
	// compute x1,y1,z1
	uint32_t x1;
	uint32_t y1;
	uint32_t z1;
	x1 = atmo_clampu(x0 + 1, 0, param->texture_width  - 1);
	y1 = atmo_clampu(y0 + 1, 0, param->texture_height - 1);
	z1 = atmo_clampu(z0 + 1, 0, param->texture_depth  - 1);

	// sampling coordinates
	double uu = u*width1  - ((double) x0);
	double vv = v*height1 - ((double) y0);
	double ww = w*depth1  - ((double) z0);

	// sample corners
	cc_vec4d_t fisk000;
	cc_vec4d_t fisk001;
	cc_vec4d_t fisk010;
	cc_vec4d_t fisk011;
	cc_vec4d_t fisk100;
	cc_vec4d_t fisk101;
	cc_vec4d_t fisk110;
	cc_vec4d_t fisk111;
	getData(param, k, x0, y0, z0, data, &fisk000);
	getData(param, k, x0, y0, z1, data, &fisk001);
	getData(param, k, x0, y1, z0, data, &fisk010);
	getData(param, k, x0, y1, z1, data, &fisk011);
	getData(param, k, x1, y0, z0, data, &fisk100);
	getData(param, k, x1, y0, z1, data, &fisk101);
	getData(param, k, x1, y1, z0, data, &fisk110);
	getData(param, k, x1, y1, z1, data, &fisk111);

	// interpolate x
	cc_vec4d_t fiskx00;
	cc_vec4d_t fiskx01;
	cc_vec4d_t fiskx10;
	cc_vec4d_t fiskx11;
	cc_vec4d_lerp(&fisk000, &fisk100, uu, &fiskx00);
	cc_vec4d_lerp(&fisk001, &fisk101, uu, &fiskx01);
	cc_vec4d_lerp(&fisk010, &fisk110, uu, &fiskx10);
	cc_vec4d_lerp(&fisk011, &fisk111, uu, &fiskx11);

	// interpolate y
	cc_vec4d_t fiskxy0;
	cc_vec4d_t fiskxy1;
	cc_vec4d_lerp(&fiskx00, &fiskx10, vv, &fiskxy0);
	cc_vec4d_lerp(&fiskx01, &fiskx11, vv, &fiskxy1);

	// interpolate z
	cc_vec4d_lerp(&fiskxy0, &fiskxy1, ww, fisk);

	#else
	getData(param, k, x0, y0, z0, data, fisk);
	#endif
}

// factored multiple-scattered gathered intensity step
static void
fGk_step(atmo_solverParam_t* param, uint32_t k,
         cc_vec3d_t* P, cc_vec3d_t* V, cc_vec3d_t* L,
         cc_vec4d_t* data, double s, double xj, double yi,
         cc_vec4d_t* fgk)
{
	ASSERT(param);
	ASSERT(P);
	ASSERT(V);
	ASSERT(L);
	ASSERT(data);
	ASSERT(fgk);

	// compute omega
	// xj => omega spherical angle theta in (0, pi)
	// yi => omega spherical angle phi   in (0, 2*pi)
	cc_vec3d_t omega =
	{
		.x = sin(xj)*cos(yi),
		.y = sin(xj)*sin(yi),
		.z = cos(xj),
	};

	// compute phase
	double cos_theta;
	double FR;
	double FM;
	cc_vec3d_t minus_omega;
	cc_vec3d_t minus_V;
	cc_vec3d_muls_copy(&omega, -1.0, &minus_omega);
	cc_vec3d_muls_copy(V,      -1.0, &minus_V);
	cos_theta = cc_vec3d_dot(&minus_omega, &minus_V);
	FR = atmo_phaseR(param, cos_theta);
	FM = atmo_phaseM(param, cos_theta);

	// sample fisk
	cc_vec4d_t fisk;
	fISk_sample(param, k, P, &omega, L, data, &fisk);

	// compute sin_theta for domega
	// xj => omega spherical angle theta in (0, pi)
	double domega_sin_xj = sin(xj);

	// add factored multiple-scattered gathered intensity
	fgk->r += FR*fisk.r*domega_sin_xj;
	fgk->g += FR*fisk.g*domega_sin_xj;
	fgk->b += FR*fisk.b*domega_sin_xj;
	fgk->a += FM*fisk.a*domega_sin_xj;
}

// factored multiple-scattered gathered intensity
static void
fGk(atmo_solverParam_t* param, uint32_t k, cc_vec3d_t* P,
    cc_vec3d_t* V, cc_vec3d_t* L, cc_vec4d_t* data,
    cc_vec4d_t* fgk)
{
	ASSERT(param);
	ASSERT(P);
	ASSERT(V);
	ASSERT(L);
	ASSERT(data);
	ASSERT(fgk);

	// initalize 2D trapezoidal rule edges
	// x => omega spherical angle theta in (0, pi)
	// y => omega spherical angle phi   in (0, 2*pi)
	double x0 = 0.0;
	double xn = M_PI;
	double y0 = 0.0;
	double ym = 2.0*M_PI;

	// compute 2D trapezoidal rule step size
	int    m  = param->gather_m_steps;
	int    n  = param->gather_n_steps;
	double dx = (xn - x0)/((double) n);
	double dy = (ym - y0)/((double) m);

	// initialize fgk
	fgk->r = 0.0;
	fgk->g = 0.0;
	fgk->b = 0.0;
	fgk->a = 0.0;

	// apply 2D trapezoidal rule for corners
	fGk_step(param, k, P, V, L, data, 1.0, x0, y0, fgk);
	fGk_step(param, k, P, V, L, data, 1.0, xn, y0, fgk);
	fGk_step(param, k, P, V, L, data, 1.0, x0, ym, fgk);
	fGk_step(param, k, P, V, L, data, 1.0, xn, ym, fgk);

	// apply 2D trapezoidal rule for edges
	int    i;
	int    j;
	double yi;
	double xj;
	for(j = 1; j < n; ++j)
	{
		xj = ((double) j)*dx;
		fGk_step(param, k, P, V, L, data, 2.0, xj, y0, fgk);
		fGk_step(param, k, P, V, L, data, 2.0, xj, ym, fgk);
	}
	for(i = 1; i < m; ++i)
	{
		yi = ((double) i)*dy;
		fGk_step(param, k, P, V, L, data, 2.0, x0, yi, fgk);
		fGk_step(param, k, P, V, L, data, 2.0, xn, yi, fgk);
	}

	// apply 2D trapezoidal rule for center
	for(i = 1; i < m; ++i)
	{
		yi = ((double) i)*dy;
		for(j = 1; j < n; ++j)
		{
			xj = ((double) j)*dx;
			fGk_step(param, k, P, V, L, data, 4.0, xj, yi, fgk);
		}
	}

	// apply 2D trapezoidal rule scale
	cc_vec4d_muls(fgk, 0.25*dx*dy);
}

// factored multiple-scattered intensity
static void
fISk(atmo_solverParam_t* param, uint32_t k,
     double h, double phi, double delta, cc_vec4d_t* data,
     cc_vec4d_t* fisk)
{
	ASSERT(param);
	ASSERT(data);
	ASSERT(fisk);

	// initialize fisk
	fisk->r = 0.0;
	fisk->g = 0.0;
	fisk->b = 0.0;
	fisk->a = 0.0;

	// canonical form of the scattering intensity
	// parameterization for P0, V and L
	cc_vec3d_t P0 =
	{
		.z = h + param->Rp,
	};
	cc_vec3d_t V =
	{
		.x = sin(phi),
		.z = cos(phi),
	};
	cc_vec3d_t L =
	{
		.x = -sin(delta),
		.z = -cos(delta),
	};

	// compute ray-sphere intersection
	// include a ray offset for the viewing vector
	// to ensure that the ray does not intersect at P0
	double Ro = 10.0;
	cc_vec3d_t Pa;
	cc_vec3d_t Pb;
	computePaPb(param, &P0, &V, Ro, &Pa, &Pb);

	// initialize integration
	cc_vec3d_t P;
	cc_vec3d_t step;
	cc_vec4d_t fx0;
	cc_vec4d_t fx1;
	cc_vec3d_t tPaP;
	cc_vec4d_t fgk;
	cc_vec3d_copy(&Pa, &P);
	cc_vec3d_subv_copy(&Pb, &Pa, &step);
	cc_vec3d_muls(&step, 1.0/param->transmittance_steps);
	double ds = cc_vec3d_mag(&step);
	double pR = densityR(param, h);
	double pM = densityM(param, h);
	fGk(param, k - 1, &P0, &V, &L, data, &fgk);
	opticalDepth(param, &Pa, &P, &tPaP);
	fx0.r = fgk.r*pR*exp(-tPaP.x);
	fx0.g = fgk.g*pR*exp(-tPaP.y);
	fx0.b = fgk.b*pR*exp(-tPaP.z);
	fx0.a = fgk.a*pM*exp(-(tPaP.x + tPaP.y + tPaP.z)/3.0);

	if(ds < ATMO_STEP_THRESH)
	{
		return;
	}

	// integrate factored multiple-scattered intensity
	int i;
	for(i = 0; i < param->transmittance_steps; ++i)
	{
		cc_vec3d_addv(&P, &step);

		h  = getHeightP(param, &P);
		pR = densityR(param, h);
		pM = densityM(param, h);
		fGk(param, k - 1, &P, &V, &L, data, &fgk);
		opticalDepth(param, &Pa, &P, &tPaP);
		fx1.r = fgk.r*pR*exp(-tPaP.x);
		fx1.g = fgk.g*pR*exp(-tPaP.y);
		fx1.b = fgk.b*pR*exp(-tPaP.z);
		fx1.a = fgk.a*pM*exp(-(tPaP.x + tPaP.y + tPaP.z)/3.0);

		// apply trapesoidal rule
		fisk->r += 0.5*(fx1.r + fx0.r)*ds;
		fisk->g += 0.5*(fx1.g + fx0.g)*ds;
		fisk->b += 0.5*(fx1.b + fx0.b)*ds;
		fisk->a += 0.5*(fx1.a + fx0.a)*ds;

		cc_vec4d_copy(&fx1, &fx0);
	}

	// apply Rayleigh/Mie scattering coefficient
	fisk->r *= param->beta_s_r_rayleigh/(4.0*M_PI);
	fisk->g *= param->beta_s_g_rayleigh/(4.0*M_PI);
	fisk->b *= param->beta_s_b_rayleigh/(4.0*M_PI);
	fisk->a *= param->beta_s_mie/(4.0*M_PI);
}

/***********************************************************
* private - utility                                        *
***********************************************************/

static void atmo_solver_deleteImages(atmo_solver_t* self)
{
	ASSERT(self);

	if(self->image_array)
	{
		int i;
		for(i = 0; i < self->param.k; ++i)
		{
			vkk_image_delete(&self->image_array[i]);
		}
		FREE(self->image_array);
		self->image_array = NULL;
	}
}

static int
atmo_solver_newImages(atmo_solver_t* self,
                      cc_vec4d_t* data)
{
	ASSERT(self);
	ASSERT(data);

	vkk_image_t* img;
	cc_vec4f_t*  datak;

	atmo_solverParam_t* param = &self->param;

	size_t count = param->k*
	               param->texture_width*
	               param->texture_height*
	               param->texture_depth;

	cc_vec4f_t*  dataf;
	dataf = (cc_vec4f_t*)
	        CALLOC(count, sizeof(cc_vec4f_t));
	if(dataf == NULL)
	{
		LOGE("CALLOC failed");
		return 0;
	}

	// copy data
	int i;
	for(i = 0; i < count; ++i)
	{
		dataf[i].x = (float) data[i].x;
		dataf[i].y = (float) data[i].y;
		dataf[i].z = (float) data[i].z;
		dataf[i].w = (float) data[i].w;
	}

	self->image_array = (vkk_image_t**)
	                    CALLOC(param->k,
	                           sizeof(vkk_image_t*));
	if(self->image_array == NULL)
	{
		LOGE("CALLOC failed");
		goto fail_image_array;
	}

	for(i = 0; i < param->k; ++i)
	{
		// k is base-1
		datak = getDataKf(&self->param, i + 1, dataf);

		img = vkk_image_new(self->engine,
		                    param->texture_width,
		                    param->texture_height,
		                    param->texture_depth,
		                    VKK_IMAGE_FORMAT_RGBAF16,
		                    0, VKK_STAGE_FS,
		                    (const void*) datak);
		if(img == NULL)
		{
			goto fail_image;
		}

		self->image_array[i] = img;
	}

	FREE(dataf);

	// success
	return 1;

	// failure
	fail_image:
	{
		for(i = 0; i < param->k; ++i)
		{
			vkk_image_delete(&self->image_array[i]);
		}
		FREE(self->image_array);
		self->image_array = NULL;
	}
	fail_image_array:
		FREE(dataf);
	return 0;
}

static int
atmo_solver_exportData(atmo_solver_t* self,
                       cc_vec4d_t* data)
{
	ASSERT(self);
	ASSERT(data);

	atmo_solverParam_t* param = &self->param;

	cc_jsmnStream_t* jsmn = cc_jsmnStream_new();
	if(jsmn == NULL)
	{
		return 0;
	}
	cc_jsmnStream_beginObject(jsmn);
	cc_jsmnStream_key(jsmn, "%s", "Rp");
	cc_jsmnStream_double(jsmn, param->Rp);
	cc_jsmnStream_key(jsmn, "%s", "Ra");
	cc_jsmnStream_double(jsmn, param->Ra);
	cc_jsmnStream_key(jsmn, "%s", "density_scale_height_rayleigh");
	cc_jsmnStream_double(jsmn, param->density_scale_height_rayleigh);
	cc_jsmnStream_key(jsmn, "%s", "density_scale_height_mie");
	cc_jsmnStream_double(jsmn, param->density_scale_height_mie);
	cc_jsmnStream_key(jsmn, "%s", "phase_g_mie");
	cc_jsmnStream_double(jsmn, param->phase_g_mie);
	cc_jsmnStream_key(jsmn, "%s", "beta_s_r_rayleigh");
	cc_jsmnStream_double(jsmn, param->beta_s_r_rayleigh);
	cc_jsmnStream_key(jsmn, "%s", "beta_s_g_rayleigh");
	cc_jsmnStream_double(jsmn, param->beta_s_g_rayleigh);
	cc_jsmnStream_key(jsmn, "%s", "beta_s_b_rayleigh");
	cc_jsmnStream_double(jsmn, param->beta_s_b_rayleigh);
	cc_jsmnStream_key(jsmn, "%s", "beta_s_mie");
	cc_jsmnStream_double(jsmn, param->beta_s_mie);
	cc_jsmnStream_key(jsmn, "%s", "beta_a_mie");
	cc_jsmnStream_double(jsmn, param->beta_a_mie);
	cc_jsmnStream_key(jsmn, "%s", "beta_a_r_ozone");
	cc_jsmnStream_double(jsmn, param->beta_a_r_ozone);
	cc_jsmnStream_key(jsmn, "%s", "beta_a_g_ozone");
	cc_jsmnStream_double(jsmn, param->beta_a_g_ozone);
	cc_jsmnStream_key(jsmn, "%s", "beta_a_b_ozone");
	cc_jsmnStream_double(jsmn, param->beta_a_b_ozone);
	cc_jsmnStream_key(jsmn, "%s", "II_r");
	cc_jsmnStream_double(jsmn, param->II_r);
	cc_jsmnStream_key(jsmn, "%s", "II_g");
	cc_jsmnStream_double(jsmn, param->II_g);
	cc_jsmnStream_key(jsmn, "%s", "II_b");
	cc_jsmnStream_double(jsmn, param->II_b);
	cc_jsmnStream_key(jsmn, "%s", "transmittance_steps");
	cc_jsmnStream_int(jsmn, param->transmittance_steps);
	cc_jsmnStream_key(jsmn, "%s", "gather_m_steps");
	cc_jsmnStream_int(jsmn, param->gather_m_steps);
	cc_jsmnStream_key(jsmn, "%s", "gather_n_steps");
	cc_jsmnStream_int(jsmn, param->gather_n_steps);
	cc_jsmnStream_key(jsmn, "%s", "k");
	cc_jsmnStream_int(jsmn, param->k);
	cc_jsmnStream_key(jsmn, "%s", "texture_width");
	cc_jsmnStream_int(jsmn, param->texture_width);
	cc_jsmnStream_key(jsmn, "%s", "texture_height");
	cc_jsmnStream_int(jsmn, param->texture_height);
	cc_jsmnStream_key(jsmn, "%s", "texture_depth");
	cc_jsmnStream_int(jsmn, param->texture_depth);
	cc_jsmnStream_end(jsmn);
	cc_jsmnStream_export(jsmn, "atmo-param.json");
	cc_jsmnStream_delete(&jsmn);

	uint32_t k;
	cc_vec4d_t* datak;
	char fname[256];
	FILE* f;
	size_t size = sizeof(cc_vec4d_t)*param->texture_width*
	              param->texture_height*param->texture_depth;
	for(k = 1; k <= param->k; ++k)
	{
		snprintf(fname, 256, "atmo-data-k%u.dat", k);
		datak = getDataK(param, k, data);
		f = fopen(fname, "w");
		if(f)
		{
			if(fwrite((const void*) datak, size, 1, f) != 1)
			{
				LOGE("invalid %s", fname);
				fclose(f);
				return 0;
			}
			fclose(f);
		}
	}

	return 1;
}

#ifdef ATMO_SOLVER_DEBUG_DATA

// https://64.github.io/tonemapping/
static void
atmo_uncharted2TonemapPartial(cc_vec3d_t* x, cc_vec3d_t* y)
{
	ASSERT(x);
	ASSERT(y);

	double A = 0.15;
	double B = 0.50;
	double C = 0.10;
	double D = 0.20;
	double E = 0.02;
	double F = 0.30;
	y->r = ((x->r*(A*x->r+C*B)+D*E)/(x->r*(A*x->r+B)+D*F))-E/F;
	y->g = ((x->g*(A*x->g+C*B)+D*E)/(x->g*(A*x->g+B)+D*F))-E/F;
	y->b = ((x->b*(A*x->b+C*B)+D*E)/(x->b*(A*x->b+B)+D*F))-E/F;
}

// https://64.github.io/tonemapping/
static void atmo_uncharted2Filmic(cc_vec3d_t* v)
{
	ASSERT(v);

	double exposure_bias = 2.0;

	cc_vec3d_t x;
	cc_vec3d_muls_copy(v, exposure_bias, &x);

	cc_vec3d_t xp = { 0 };
	atmo_uncharted2TonemapPartial(&x, &xp);

	cc_vec3d_t W =
	{
		.r = 11.2,
		.g = 11.2,
		.b = 11.2,
	};
	cc_vec3d_t Wp = { 0 };
	atmo_uncharted2TonemapPartial(&W, &Wp);

	cc_vec3d_t white_scale =
	{
		.r = 1.0/Wp.r,
		.g = 1.0/Wp.g,
		.b = 1.0/Wp.b,
	};

	cc_vec3d_mulv_copy(&xp, &white_scale, v);
}

static void atmo_gamma(cc_vec3d_t* color)
{
	ASSERT(color);

	color->r = pow(color->r, 1.0/2.2);
	color->g = pow(color->g, 1.0/2.2);
	color->b = pow(color->b, 1.0/2.2);
}

static int
atmo_solver_plotDensity(atmo_solverParam_t* param)
{
	ASSERT(param);

	FILE* f = fopen("plot_density.dat", "w");
	if(f == NULL)
	{
		LOGE("fopen failed");
		return 0;
	}

	// output density
	int   i;
	int   nh = 100;
	double Ha = param->Ra - param->Rp;
	for(i = 0; i < nh; ++i)
	{
		double h;
		h = Ha*((double) i)/((double) (nh - 1));

		fprintf(f, "%lf, %lf, %lf\n",
		        densityR(param, h),
		        densityM(param, h),
		        densityO(param, h));
	}

	fclose(f);

	return 1;
}

static double atmo_deg2rad(double x)
{
	return x*M_PI/180.0;
}

static int
atmo_solver_plotAvgT(atmo_solverParam_t* param)
{
	ASSERT(param);

	FILE* f = fopen("plot_avgT.dat", "w");
	if(f == NULL)
	{
		LOGE("fopen failed");
		return 0;
	}

	// output a 3D plot average transmittance
	int    i;
	int    j;
	int    nh   = 40;
	int    nphi = 721;
	double Ha   = param->Ra - param->Rp;
	for(i = 0; i < nh; ++i)
	{
		double h;
		h = atmo_clampd(Ha*((double) i)/((double) (nh - 1)),
		                10.0, Ha - 10.0);

		for(j = 0; j < nphi; ++j)
		{
			double phi_deg;
			phi_deg = 180.0*((double) j)/
			                 ((double) (nphi - 1));

			double phi_rad;
			phi_rad = atmo_deg2rad(phi_deg);

			// canonical form of the scattering intensity
			// parameterization for P0 and V
			cc_vec3d_t P0 =
			{
				.z = h + param->Rp,
			};
			cc_vec3d_t V =
			{
				.x = sin(phi_rad),
				.z = cos(phi_rad),
			};

			// compute ray-sphere intersection
			// include a ray offset for the viewing vector
			// to ensure that the ray does not intersect at P0
			double Ro = 10.0;
			cc_vec3d_t Pa;
			cc_vec3d_t Pb;
			computePaPb(param, &P0, &V, Ro, &Pa, &Pb);

			cc_vec3d_t Vba;
			cc_vec3d_subv_copy(&Pb, &Pa, &Vba);
			double dist = cc_vec3d_mag(&Vba);

			cc_vec3d_t t;
			opticalDepth(param, &P0, &Pb, &t);

			cc_vec3d_t T;
			T.x = exp(-t.x);
			T.y = exp(-t.y);
			T.z = exp(-t.z);

			double avgT = (T.x + T.y + T.z)/3.0;

			LOGI("h=%lf, phi=%lf, avgT=%lf, dist=%lf, T=%lf,%lf,%lf",
			     h, phi_deg, avgT, dist,
			     T.x, T.y, T.z);

			if(j == 0)
			{
				fprintf(f, "%lf", avgT);
			}
			else
			{
				fprintf(f, " %lf", avgT);
			}
		}
		fprintf(f, "\n");
	}

	fclose(f);

	return 1;
}

static int atmo_solver_plotSpectralIrradiance(void)
{
	FILE* f = fopen("plot_spectralIrradiance.dat", "w");
	if(f == NULL)
	{
		LOGE("invalid");
		return 0;
	}

	int min = ATMO_SPECTRAL_IRRADIANCE_MIN;
	int max = ATMO_SPECTRAL_IRRADIANCE_MAX;
	int i;
	double lambda;
	for(i = min; i <= max; ++i)
	{
		lambda = (double) i;
		fprintf(f, "%lf %lf\n",
		        lambda, atmo_spectralIrradiance_get(lambda));
	}

	fclose(f);

	return 1;
}

static int atmo_solver_plotSpectralToRGB(void)
{
	FILE* f = fopen("plot_spectralToRGB.dat", "w");
	if(f == NULL)
	{
		LOGE("fopen failed");
		return 0;
	}

	int i;
	cc_vec3d_t rgb_none;
	cc_vec3d_t rgb_sum;
	cc_vec3d_t rgb_peak;
	for(i = ATMO_SPECTRAL_TO_RGB_MIN;
	    i <= ATMO_SPECTRAL_TO_RGB_MAX; ++i)
	{
		atmo_spectralToRGB_getRGB(ATMO_SPECTRAL_TO_RGB_NORMALIZE_NONE,
		                          i, &rgb_none);
		atmo_spectralToRGB_getRGB(ATMO_SPECTRAL_TO_RGB_NORMALIZE_SUM,
		                          i, &rgb_sum);
		atmo_spectralToRGB_getRGB(ATMO_SPECTRAL_TO_RGB_NORMALIZE_PEAK,
		                          i, &rgb_peak);
		fprintf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		        (double) i,
		        rgb_none.x, rgb_none.y, rgb_none.z,
		        rgb_sum.x,  rgb_sum.y,  rgb_sum.z,
		        rgb_peak.x, rgb_peak.y, rgb_peak.z);
	}

	fclose(f);

	return 1;
}

static int
atmo_solver_debugData(atmo_solver_t* self, cc_vec4d_t* data)
{
	ASSERT(self);
	ASSERT(data);

	atmo_solverParam_t* param = &self->param;

	// debug slices by height
	int texw = param->texture_width*param->texture_depth;
	int texh = 2*param->texture_height;

	texgz_tex_t* tex;
	tex = texgz_tex_new(texw, texh, texw, texh,
	                    TEXGZ_UNSIGNED_BYTE, TEXGZ_RGB,
	                    NULL);
	if(tex == NULL)
	{
		return 0;
	}

	// spectral intensity of incident light from the Sun
	double exposure = -2.5;
	cc_vec3d_t II =
	{
		.r = pow(2.0, exposure)*param->II_r,
		.g = pow(2.0, exposure)*param->II_g,
		.b = pow(2.0, exposure)*param->II_b,
	};

	LOGI("II: r=%lf, g=%lf, b=%lf", II.r, II.g, II.b);

	double   phi;
	double   delta;
	double   cos_phi;
	double   cos_delta;
	double   cos_theta;
	double   h;
	double   u;
	double   v;
	double   w;
	double   FR;
	double   FM;
	uint32_t k;
	uint32_t x;
	uint32_t y;
	uint32_t z;
	char fname[256];
	unsigned char pixel[4] = { 0, 0, 0, 255 };
	cc_vec4d_t fis;
	cc_vec3d_t is;
	cc_vec3d_t color;
	for(k = 1; k <= param->k; ++k)
	{
		for(x = 0; x < param->texture_width; ++x)
		{
			u = ((double) x)/((double) (param->texture_width - 1));
			h = getHeightU(param, u);

			for(z = 0; z < param->texture_depth; ++z)
			{
				w = ((double) z)/((double) (param->texture_depth - 1));
				cos_delta = getCosDeltaW(w);
				delta = acos(cos_delta);

				for(y = 0; y < param->texture_height; ++y)
				{
					v = ((double) y)/((double) (param->texture_height - 1));
					cos_phi = getCosPhiV(param, h, u, v);
					phi = acos(cos_phi);

					getData(param, k, x, y, z, data, &fis);

					cos_theta = cos(delta - phi);
					FR = atmo_phaseR(param, cos_theta);
					FM = atmo_phaseM(param, cos_theta);

					is.r = II.r*(FR*fis.r + FM*fis.a);
					is.g = II.g*(FR*fis.g + FM*fis.a);
					is.b = II.b*(FR*fis.b + FM*fis.a);

					// tone mapping and gamma correction
					cc_vec3d_copy(&is, &color);
					atmo_uncharted2Filmic(&color);
					atmo_gamma(&color);

					// clamp output
					color.r = atmo_clampd(color.r, 0.0, 1.0);
					color.g = atmo_clampd(color.g, 0.0, 1.0);
					color.b = atmo_clampd(color.b, 0.0, 1.0);

					// set output
					pixel[0] = (unsigned char) (255.0*color.r);
					pixel[1] = (unsigned char) (255.0*color.g);
					pixel[2] = (unsigned char) (255.0*color.b);
					texgz_tex_setPixel(tex, x*param->texture_depth + z,
					                   y, pixel);

					// set u,v,w
					pixel[0] = (unsigned char) (255.0*u);
					pixel[1] = (unsigned char) (255.0*v);
					pixel[2] = (unsigned char) (255.0*w);
					texgz_tex_setPixel(tex, x*param->texture_depth + z,
					                   param->texture_height + y, pixel);
				}
			}
		}
		snprintf(fname, 256, "atmo-slice-k%u.png", k);
		texgz_png_export(tex, fname);
	}
	texgz_tex_delete(&tex);

	atmo_solver_plotDensity(param);
	atmo_solver_plotAvgT(param);
	atmo_solver_plotSpectralIrradiance();
	atmo_solver_plotSpectralToRGB();

	return 1;
}
#else
static int
atmo_solver_debugData(atmo_solver_t* self, cc_vec4d_t* data)
{
	return 1;
}
#endif

static int
atmo_solver_paramValidate(atmo_solverParam_t* param)
{
	ASSERT(param);

	// perform minimalistic error checking

	if((param->Rp < 1.0) || (param->Ra < 1.0))
	{
		LOGE("invalid Rp=%lf, Ra=%lf",
		     param->Rp, param->Ra);
		return 0;
	}

	if((param->density_scale_height_rayleigh < 0.0) ||
	   (param->density_scale_height_mie      < 0.0))
	{
		LOGE("invalid Rp=%lf, Ra=%lf",
		     param->density_scale_height_rayleigh,
		     param->density_scale_height_mie);
		return 0;
	}

	if((param->phase_g_mie < -0.9999) ||
	   (param->phase_g_mie > 0.9999))
	{
		LOGE("invalid phase_g_mie=%lf",
		     param->phase_g_mie);
		return 0;
	}

	if((param->beta_s_r_rayleigh <= 0.0) ||
	   (param->beta_s_g_rayleigh <= 0.0) ||
	   (param->beta_s_b_rayleigh <= 0.0) ||
	   (param->beta_s_mie        <= 0.0) ||
	   (param->beta_a_mie        <= 0.0) ||
	   (param->beta_a_r_ozone    <= 0.0) ||
	   (param->beta_a_g_ozone    <= 0.0) ||
	   (param->beta_a_b_ozone    <= 0.0))
	{
		LOGE("invalid beta_s_r/g/b_rayleigh=%lf/%lf/%lf, beta_s_mie=%lf, beta_a_mie=%lf, beta_a_r/g/b_ozone=%lf/%lf/%lf",
		     param->beta_s_r_rayleigh,
		     param->beta_s_g_rayleigh,
		     param->beta_s_b_rayleigh,
		     param->beta_s_mie,
		     param->beta_a_mie,
		     param->beta_a_r_ozone,
		     param->beta_a_g_ozone,
		     param->beta_a_b_ozone);
		return 0;
	}

	if((param->II_r <= 0.0) ||
	   (param->II_g <= 0.0) ||
	   (param->II_b <= 0.0))
	{
		LOGE("invalid II_r/g/b=%lf/%lf/%lf",
		     param->II_r,
		     param->II_g,
		     param->II_b);
		return 0;
	}

	if((param->transmittance_steps < 1) ||
	   (param->transmittance_steps > 100))
	{
		LOGE("invalid transmittance_steps=%u",
		     param->transmittance_steps);
		return 0;
	}

	if((param->gather_m_steps < 1) ||
	   (param->gather_m_steps > 180))
	{
		LOGE("invalid gather_m_steps=%u",
		     param->gather_m_steps);
		return 0;
	}

	if((param->gather_n_steps < 1) ||
	   (param->gather_n_steps > 360))
	{
		LOGE("invalid gather_n_steps=%u",
		     param->gather_n_steps);
		return 0;
	}

	if((param->k < 1) || (param->k > 10))
	{
		LOGE("invalid k=%u", param->k);
		return 0;
	}

	if((param->texture_width  < 8) ||
	   (param->texture_height < 8) ||
	   (param->texture_depth  < 8) ||
	   (cc_find_pow2n(param->texture_width)  < 0) ||
	   (cc_find_pow2n(param->texture_height) < 0) ||
	   (cc_find_pow2n(param->texture_depth)  < 0))
	{
		LOGE("invalid texture_width=%u, texture_height=%u, texture_depth=%u",
		     param->texture_width,
		     param->texture_height,
		     param->texture_depth);
		return 0;
	}

	return 1;
}

static void
atmo_solver_step(atmo_solver_t* self, uint32_t k,
                 uint32_t x, uint32_t y, uint32_t z,
                 cc_vec4d_t* data)
{
	ASSERT(self);
	ASSERT(data);

	atmo_solverParam_t* param = &self->param;

	double width  = (double) param->texture_width;
	double height = (double) param->texture_height;
	double depth  = (double) param->texture_depth;

	double u = ((double) x)/(width - 1.0);
	double v = ((double) y)/(height - 1.0);
	double w = ((double) z)/(depth - 1.0);

	double h         = getHeightU(param, u);
	double cos_phi   = getCosPhiV(param, h, u, v);
	double cos_delta = getCosDeltaW(w);
	double phi       = acos(cos_phi);
	double delta     = acos(cos_delta);

	cc_vec4d_t fis;
	if(k == 1)
	{
		fIS1(param, h, phi, delta, &fis);
	}
	else
	{
		fISk(param, k, h, phi, delta, data, &fis);
	}

	setData(param, k, x, y, z, data, &fis);
}

static void
atmo_solver_finish(atmo_solver_t* self, cc_vec4d_t* data)
{
	ASSERT(self);
	ASSERT(data);

	atmo_solverParam_t* param = &self->param;

	uint32_t   x;
	uint32_t   y;
	uint32_t   z;
	double     mag0_r = 0.0;
	double     mag0_m = 0.0;
	cc_vec4d_t fis0;
	for(z = 0; z < param->texture_depth; ++z)
	{
		for(y = 0; y < param->texture_height; ++y)
		{
			for(x = 0; x < param->texture_width; ++x)
			{
				getData(param, 1, x, y, z, data, &fis0);
				mag0_r += cc_vec3d_mag((cc_vec3d_t*) &fis0);
				mag0_m += fis0.a;
			}
		}
	}
	LOGI("k=1, mag0_r=%lf, mag0_m=%lf", mag0_r, mag0_m);

	// compute total intensity for each k
	uint32_t   k;
	double     mag1_r = 0.0;
	double     mag1_m = 0.0;
	double     ATTEN_R = 0.97;
	double     ATTEN_M = 0.8;
	double     atten_r = ATTEN_R;
	double     atten_m = ATTEN_M;
	cc_vec4d_t fis1;
	for(k = 2; k <= param->k; ++k)
	{
		for(z = 0; z < param->texture_depth; ++z)
		{
			for(y = 0; y < param->texture_height; ++y)
			{
				for(x = 0; x < param->texture_width; ++x)
				{
					getData(param, k, x, y, z, data, &fis1);
					mag1_r += cc_vec3d_mag((cc_vec3d_t*) &fis1);
					mag1_m += fis1.a;
				}
			}
		}
		LOGI("k=%u, mag1_r=%lf, mag1_m=%lf, atten_r=%lf (%lf), atten_m=%lf (%lf)",
		     k, mag1_r, mag1_m, mag1_r/mag0_r, atten_r, mag1_m/mag0_m, atten_m);

		cc_vec4d_t fis;
		for(z = 0; z < param->texture_depth; ++z)
		{
			for(y = 0; y < param->texture_height; ++y)
			{
				for(x = 0; x < param->texture_width; ++x)
				{
					getData(param, k - 1, x, y, z, data, &fis0);
					getData(param, k,     x, y, z, data, &fis1);

					// optionally adjust the attenuation for each k
					#ifdef ATMO_ADJUST_ATTENUATION
					fis1.r *= atten_r*mag0_r/mag1_r;
					fis1.g *= atten_r*mag0_r/mag1_r;
					fis1.b *= atten_r*mag0_r/mag1_r;
					fis1.a *= atten_m*mag0_m/mag1_m;
					#endif

					cc_vec4d_addv_copy(&fis0, &fis1, &fis);
					setData(param, k, x, y, z, data, &fis);
				}
			}
		}

		mag1_r = 0.0;
		mag1_m = 0.0;
		atten_r *= ATTEN_R;
		if(k <= 3)
		{
			atten_m *= ATTEN_M;
		}
	}
}

static void atmo_solver_run(int tid, void* owner, void* task)
{
	ASSERT(owner);
	ASSERT(task);

	atmo_solver_t*      self  = (atmo_solver_t*) task;
	atmo_solverParam_t* param = &self->param;

	size_t count = param->k*
	               param->texture_width*
	               param->texture_height*
	               param->texture_depth;

	cc_vec4d_t* data;
	data = (cc_vec4d_t*)
	       CALLOC(count, sizeof(cc_vec4d_t));
	if(data == NULL)
	{
		LOGE("CALLOC failed");
		return;
	}

	uint32_t k; // k is base-1
	uint32_t x;
	uint32_t y;
	uint32_t z;
	uint32_t step  = 1;
	uint32_t steps = param->k*param->texture_depth;
	for(k = 1; k <= param->k; ++k)
	{
		for(z = 0; z < param->texture_depth; ++z)
		{
			for(y = 0; y < param->texture_height; ++y)
			{
				for(x = 0; x < param->texture_width; ++x)
				{
					atmo_solver_step(self, k, x, y, z, data);
				}
			}

			// update progress and check status
			pthread_mutex_lock(&self->mutex);
			self->progress = ((float) step)/((float) steps);
			LOGI("progress=%f", self->progress);
			if(self->status == ATMO_SOLVER_STATUS_STOPPING)
			{
				FREE(data);
				return;
			}
			pthread_mutex_unlock(&self->mutex);

			++step;
		}
	}

	atmo_solver_finish(self, data);
	atmo_solver_newImages(self, data);
	atmo_solver_exportData(self, data);
	atmo_solver_debugData(self, data);

	FREE(data);
}

/***********************************************************
* public                                                   *
***********************************************************/

atmo_solver_t* atmo_solver_new(vkk_engine_t* engine)
{
	ASSERT(engine);

	atmo_solver_t* self;
	self = (atmo_solver_t*)
	       CALLOC(1, sizeof(atmo_solver_t));
	if(self == NULL)
	{
		LOGE("CALLOC failed");
		return NULL;
	}

	self->engine = engine;

	atmo_solver_defaultParam(self, &self->param);

	// PTHREAD_MUTEX_DEFAULT is not re-entrant
	if(pthread_mutex_init(&self->mutex, NULL) != 0)
	{
		LOGE("pthread_mutex_init failed");
		goto fail_mutex;
	}

	self->jobq = cc_jobq_new(self, 1,
	                         CC_JOBQ_THREAD_PRIORITY_HIGH,
	                         atmo_solver_run);
	if(self->jobq == NULL)
	{
		goto fail_jobq;
	}

	atmo_solver_solve(self, &self->param);

	// success
	return self;

	// failure
	fail_jobq:
		pthread_mutex_destroy(&self->mutex);
	fail_mutex:
		FREE(self);
	return NULL;
}

void atmo_solver_delete(atmo_solver_t** _self)
{
	ASSERT(_self);

	atmo_solver_t* self = *_self;
	if(self)
	{
		atmo_solver_stop(self);
		cc_jobq_delete(&self->jobq);
		pthread_mutex_destroy(&self->mutex);
		atmo_solver_deleteImages(self);
		FREE(self);
		*_self = NULL;
	}
}

void atmo_solver_defaultParam(atmo_solver_t* self,
                              atmo_solverParam_t* param)
{
	ASSERT(self);
	ASSERT(param);

	atmo_solverParam_t default_param =
	{
		.Rp = ATMO_RP,
		.Ra = ATMO_RA,

		.density_scale_height_rayleigh = ATMO_DENSITY_SCALE_HEIGHT_RAYLEIGH,
		.density_scale_height_mie      = ATMO_DENSITY_SCALE_HEIGHT_MIE,

		.phase_g_mie = ATMO_PHASE_G_MIE,

		.beta_s_r_rayleigh = ATMO_BETA_S_R_RAYLEIGH,
		.beta_s_g_rayleigh = ATMO_BETA_S_G_RAYLEIGH,
		.beta_s_b_rayleigh = ATMO_BETA_S_B_RAYLEIGH,

		.beta_s_mie = ATMO_BETA_S_MIE,
		.beta_a_mie = ATMO_BETA_A_MIE,

		.beta_a_r_ozone = ATMO_BETA_A_R_OZONE,
		.beta_a_g_ozone = ATMO_BETA_A_G_OZONE,
		.beta_a_b_ozone = ATMO_BETA_A_B_OZONE,

		.transmittance_steps = ATMO_TRANSMITTANCE_STEPS,

		.gather_m_steps = ATMO_GATHER_M_STEPS,
		.gather_n_steps = ATMO_GATHER_N_STEPS,

		.k = ATMO_K,

		.texture_width  = ATMO_TEXTURE_WIDTH,
		.texture_height = ATMO_TEXTURE_HEIGHT,
		.texture_depth  = ATMO_TEXTURE_DEPTH,
	};

	atmo_solver_computeII(&default_param);

	memcpy(param, &default_param, sizeof(atmo_solverParam_t));
}

void atmo_solver_currentParam(atmo_solver_t* self,
                              atmo_solverParam_t* param)
{
	ASSERT(self);
	ASSERT(param);

	memcpy(param, &self->param, sizeof(atmo_solverParam_t));
}

vkk_image_t*
atmo_solver_image(atmo_solver_t* self, uint32_t k)
{
	ASSERT(self);

	float progress = 0.0f;
	atmo_solverStatus_e status;
	status = atmo_solver_status(self, &progress);

	// synchronization not required when solver is stopped
	if((status == ATMO_SOLVER_STATUS_STOPPED) &&
	   (self->image_array) &&
	   (k > 0) && (k <= self->param.k))
	{
		// k is base-1
		return self->image_array[k - 1];
	}
	return NULL;
}

atmo_solverStatus_e
atmo_solver_status(atmo_solver_t* self, float* _progress)
{
	ASSERT(self);
	ASSERT(_progress);

	atmo_solverStatus_e status;

	int pending = cc_jobq_pending(self->jobq);

	pthread_mutex_lock(&self->mutex);

	// check progress
	*_progress = self->progress;

	// update jobq status
	if(pending == 0)
	{
		self->status = ATMO_SOLVER_STATUS_STOPPED;
	}
	status = self->status;

	pthread_mutex_unlock(&self->mutex);

	return status;
}

int atmo_solver_solve(atmo_solver_t* self,
                      atmo_solverParam_t* param)
{
	ASSERT(self);
	ASSERT(param);

	if(atmo_solver_paramValidate(param) == 0)
	{
		return 0;
	}

	float progress = 0.0f;
	if(atmo_solver_status(self, &progress) !=
	   ATMO_SOLVER_STATUS_STOPPED)
	{
		return 0;
	}

	// synchronization not required when solver is stopped
	self->progress = 0.0f;
	atmo_solver_deleteImages(self);
	if(&self->param != param)
	{
		// optionally copy param
		memcpy(&self->param, param,
		       sizeof(atmo_solverParam_t));
	}

	self->status = ATMO_SOLVER_STATUS_RUNNING;
	if(cc_jobq_run(self->jobq, self))
	{
		return 1;
	}
	self->status = ATMO_SOLVER_STATUS_STOPPED;

	return 0;
}

void atmo_solver_stop(atmo_solver_t* self)
{
	ASSERT(self);

	pthread_mutex_lock(&self->mutex);
	if(self->status == ATMO_SOLVER_STATUS_RUNNING)
	{
		self->status = ATMO_SOLVER_STATUS_STOPPING;
	}
	pthread_mutex_unlock(&self->mutex);
}
