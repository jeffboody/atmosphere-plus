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
#include "libcc/math/cc_float.h"
#include "libcc/math/cc_pow2n.h"
#include "libcc/math/cc_ray3f.h"
#include "libcc/math/cc_sphere.h"
#include "libcc/math/cc_vec3d.h"
#include "libcc/math/cc_vec3f.h"
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
#define ATMO_DENSITY_SCALE_HEIGHT_RAYLEIGH 8000.0f
#define ATMO_DENSITY_SCALE_HEIGHT_MIE      1200.0f

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
#define ATMO_PHASE_G_MIE 0.8f

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
#define ATMO_BETA_S_R_RAYLEIGH 5.802e-6f
#define ATMO_BETA_S_G_RAYLEIGH 13.558e-6f
#define ATMO_BETA_S_B_RAYLEIGH 33.1e-6f
#define ATMO_BETA_S_MIE        3.996e-6f
#define ATMO_BETA_A_MIE        4.40e-6f
#define ATMO_BETA_A_R_OZONE    0.65e-6f
#define ATMO_BETA_A_G_OZONE    1.881e-6f
#define ATMO_BETA_A_B_OZONE    0.085e-6f

// transmittance numerical integration steps
#define ATMO_TRANSMITTANCE_STEPS 30

#define ATMO_STEP_THRESH 10.0f

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
#define ATMO_PARAM_PHI_WEIGHTED_POWER_PU  2.0f
#define ATMO_PARAM_PHI_WEIGHTED_POWER_PL  1.0f
#define ATMO_PARAM_PHI_WEIGHTED_POWER_PS  2.0f
#define ATMO_PARAM_PHI_WEIGHTED_POWER_WL1 (20.0f/32.0f)
#define ATMO_PARAM_PHI_WEIGHTED_POWER_WS0 (4.0f/32.0f)
#define ATMO_PARAM_PHI_WEIGHTED_POWER_WS1 (4.0f/32.0f)

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

	// integrate spectral irradiance and xyz
	int i;
	double si;
	cc_vec3d_t xyz  = { 0 };
	cc_vec3d_t xyz0;
	cc_vec3d_t xyz1;
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

	cc_vec3d_t II;
	cc_mat3d_mulv_copy(&Minv, &xyz, &II);

	param->II_r = (float) II.x;
	param->II_g = (float) II.y;
	param->II_b = (float) II.z;

	LOGI("II: r=%f, g=%f, b=%f",
	     param->II_r, param->II_g, param->II_b);
}

static uint32_t
atmo_clamp(uint32_t v, uint32_t min, uint32_t max)
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

static float getUHeight(atmo_solverParam_t* param, float h)
{
	ASSERT(param);

	float u;
	#if ATMO_PARAM_HEIGHT == ATMO_PARAM_HEIGHT_POWER
	u = powf(h/(param->Ra - param->Rp), 1.0f/2.0f);
	#else
	u = h/(param->Ra - param->Rp);
	#endif

	return cc_clamp(u, 0.0f, 1.0f);
}

static float getHeightU(atmo_solverParam_t* param, float u)
{
	ASSERT(param);

	float h;
	#if ATMO_PARAM_HEIGHT == ATMO_PARAM_HEIGHT_POWER
	h = powf(u, 2.0f)*(param->Ra - param->Rp);
	#else
	h = u*(param->Ra - param->Rp);
	#endif

	float Ha = param->Ra - param->Rp;

	return cc_clamp(h, 0.0f, Ha);
}

// compute height using double precision since the magnitude
// of points in the atmosphere produces very large numbers
static float
getHeightP(atmo_solverParam_t* param, cc_vec3f_t* P)
{
	ASSERT(P);

	cc_vec3d_t d =
	{
		.x = P->x,
		.y = P->y,
		.z = P->z,
	};

	float h  = (float) (cc_vec3d_mag(&d) - param->Rp);
	float Ha = param->Ra - param->Rp;

	return cc_clamp(h, 0.0f, Ha);
}

// compute Zenith using double precision since the magnitude
// of points in the atmosphere produces very large numbers
static void
getZenithP(const cc_vec3f_t* P, cc_vec3f_t* Zenith)
{
	ASSERT(P);

	cc_vec3d_t up =
	{
		.x = P->x,
		.y = P->y,
		.z = P->z,
	};
	cc_vec3d_normalize(&up);

	Zenith->x = (float) up.x;
	Zenith->y = (float) up.y;
	Zenith->z = (float) up.z;
}

static float
getCosPhiV(atmo_solverParam_t* param, float h, float u,
           float v)
{
	ASSERT(param);

	float Rp = param->Rp;

	float cos_phi;
	#if ATMO_PARAM_PHI == ATMO_PARAM_PHI_POWER
	float vv = 2.0f*v - 1.0f;
	cos_phi = cc_sign(vv)*powf(fabs(vv), 3.0f);
	#elif ATMO_PARAM_PHI == ATMO_PARAM_PHI_BODARE
	float ch = -sqrtf(h*(2.0f*Rp + h))/(Rp + h);
	if(v > 0.5f)
	{
		cos_phi = ch + powf(v - 0.5f, 5.0f)*(1.0f - ch);
	}
	else
	{
		cos_phi = ch - powf(v, 5.0f)*(1.0f + ch);
	}
	#elif ATMO_PARAM_PHI == ATMO_PARAM_PHI_WEIGHTED_POWER
	float hypH      = Rp + h;
	float oppH      = Rp;
	float adjH      = sqrt(hypH*hypH - oppH*oppH);
	float cos_phi_H = -adjH/hypH;

	float PU  = ATMO_PARAM_PHI_WEIGHTED_POWER_PU;
	float PL  = ATMO_PARAM_PHI_WEIGHTED_POWER_PL;
	float PS  = ATMO_PARAM_PHI_WEIGHTED_POWER_PS;
	float WL1 = ATMO_PARAM_PHI_WEIGHTED_POWER_WL1;
	float WS0 = ATMO_PARAM_PHI_WEIGHTED_POWER_WS0;
	float WS1 = ATMO_PARAM_PHI_WEIGHTED_POWER_WS1;
	float WS  = WS1*u + WS0;
	float WL  = WL1*u;
	float WU  = 1.0f - WL - WS;
	float epsilon = 0.00001f;
	if(v >= WL + WS)
	{
		cos_phi = powf((v - (WL + WS))/WU, PU);
	}
	else if(v >= WS)
	{
		cos_phi = -cos_phi_H*powf((v - WS)/
		                          cc_max(WL, epsilon), PL) +
		           cos_phi_H;
	}
	else
	{
		cos_phi = (-1.0f - cos_phi_H)*powf(1.0f - v/WS, PS) +
		          cos_phi_H;
	}
	#else
	cos_phi = 2.0f*v - 1.0f;
	#endif

	return cc_clamp(cos_phi, -1.0f, 1.0f);
}

static float
getVCosPhi(atmo_solverParam_t* param, float h,
           float cos_phi, float u)
{
	ASSERT(param);

	float Rp = param->Rp;

	float v;
	#if ATMO_PARAM_PHI == ATMO_PARAM_PHI_POWER
	v = 0.5f*(1.0f + cc_sign(cos_phi)*
	                 powf(fabsf(cos_phi), 1.0f/3.0f));
	#elif ATMO_PARAM_PHI == ATMO_PARAM_PHI_BODARE
	float ch = -sqrtf(h*(2.0f*Rp + h))/(Rp + h);
	if(cos_phi > ch)
	{
		v = 0.5f*powf((cos_phi - ch)/(1.0f - ch), 0.2f) + 0.5f;
	}
	else
	{
		v = 0.5f*powf((ch - cos_phi)/(1.0f + ch), 0.2f);
	}
	#elif ATMO_PARAM_PHI == ATMO_PARAM_PHI_WEIGHTED_POWER
	float hypH      = Rp + h;
	float oppH      = Rp;
	float adjH      = sqrt(hypH*hypH - oppH*oppH);
	float cos_phi_H = -adjH/hypH;

	float PU      = ATMO_PARAM_PHI_WEIGHTED_POWER_PU;
	float PL      = ATMO_PARAM_PHI_WEIGHTED_POWER_PL;
	float PS      = ATMO_PARAM_PHI_WEIGHTED_POWER_PS;
	float WL1     = ATMO_PARAM_PHI_WEIGHTED_POWER_WL1;
	float WS0     = ATMO_PARAM_PHI_WEIGHTED_POWER_WS0;
	float WS1     = ATMO_PARAM_PHI_WEIGHTED_POWER_WS1;
	float WS      = WS1*u + WS0;
	float WL      = WL1*u;
	float WU      = 1.0f - WL - WS;
	float epsilon = 0.00001f;
	if(cos_phi >= 0.0f)
	{
		v = WU*powf(cos_phi, 1.0f/PU) + (WL + WS);
	}
	else if(cos_phi >= cos_phi_H)
	{
		v = WL*powf((cos_phi - cos_phi_H)/
		            cc_max(-cos_phi_H, epsilon),
		            1.0f/PL) + WS;
	}
	else
	{
		v = WS*(1.0f - powf((cos_phi - cos_phi_H)/
		                    (-1.0f - cos_phi_H), 1.0f/PS));
	}
	#else
	v = (cos_phi + 1.0f)/2.0f;
	#endif

	return cc_clamp(v, 0.0f, 1.0f);
}

static float getCosDeltaW(float w)
{
	float cos_delta;
	#if ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_POWER
	float ww = 2.0f*w - 1.0f;
	cos_delta = cc_sign(ww)*powf(fabs(ww), 3.0f);
	#elif ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_BODARE
	cos_delta = tanf((2.0f*w - 1.0f + 0.26f)*0.75f)/
	            tanf(1.26f*0.75f);
	#else
	cos_delta = 2.0f*w - 1.0f;
	#endif

	return cc_clamp(cos_delta, -1.0f, 1.0f);
}

static float getWCosDelta(float cos_delta)
{
	float w;
	#if ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_POWER
	w = 0.5f*(1.0f + cc_sign(cos_delta)*
	                 powf(fabsf(cos_delta), 1.0f/3.0f));
	#elif ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_BODARE
	w =  0.5f*(atanf(cc_max(cos_delta, -0.1975f)*
	                 tanf(1.26f*1.1f))/1.1f +
	           (1.0f - 0.26f));
	#else
	w = (cos_delta + 1.0f)/2.0f;
	#endif

	return cc_clamp(w, 0.0f, 1.0f);
}

static cc_vec4f_t*
getDataK(atmo_solverParam_t* param,
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
        cc_vec4f_t* data, cc_vec4f_t* val)
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
        cc_vec4f_t* data, cc_vec4f_t* val)
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
static float densityR(atmo_solverParam_t* param, float h)
{
	ASSERT(param);

	return expf(-h/param->density_scale_height_rayleigh);
}

// Mie density function
static float densityM(atmo_solverParam_t* param, float h)
{
	ASSERT(param);

	return expf(-h/param->density_scale_height_mie);
}

// Ozone density function
static float densityO(atmo_solverParam_t* param, double h)
{
	ASSERT(param);

	float p = 1.0f - fabsf(h - 25000.0f)/15000.0f;
	if(p < 0.0)
	{
		return 0.0f;
	}

	return p;
}

// Rayleigh/Mie transmittance
static void
transmittance(atmo_solverParam_t* param, cc_vec3f_t* P1,
              cc_vec3f_t* P2, cc_vec3f_t* out)
{
	ASSERT(param);
	ASSERT(P1);
	ASSERT(P2);
	ASSERT(out);

	// initialize integration
	cc_vec3f_t P;
	cc_vec3f_t step;
	cc_vec3f_copy(P1, &P);
	cc_vec3f_subv_copy(P2, P1, &step);
	cc_vec3f_muls(&step, 1.0f/param->transmittance_steps);
	float ds  = cc_vec3f_mag(&step);
	float h   = getHeightP(param, &P);
	float pR0 = densityR(param, h);
	float pM0 = densityM(param, h);
	float pO0 = densityO(param, h);

	if(ds < ATMO_STEP_THRESH)
	{
		out->x = 0.0f;
		out->y = 0.0f;
		out->z = 0.0f;
		return;
	}

	// integrate transmittance
	int   i;
	float pR1;
	float pM1;
	float pO1;
	float tR = 0.0f;
	float tM = 0.0f;
	float tO = 0.0f;
	for(i = 0; i < param->transmittance_steps; ++i)
	{
		cc_vec3f_addv(&P, &step);

		h   = getHeightP(param, &P);
		pR1 = densityR(param, h);
		pM1 = densityM(param, h);
		pO1 = densityO(param, h);

		// apply trapesoidal rule
		tR += 0.5f*(pR0 + pR1)*ds;
		tM += 0.5f*(pM0 + pM1)*ds;
		tO += 0.5f*(pO0 + pO1)*ds;

		pR0 = pR1;
		pM0 = pM1;
		pO0 = pO1;
	}

	// extinction coefficients (scattering + absorption)
	cc_vec3f_t bR =
	{
		.x = param->beta_s_r_rayleigh,
		.y = param->beta_s_g_rayleigh,
		.z = param->beta_s_b_rayleigh,
	};

	float bM = param->beta_s_mie + param->beta_a_mie;

	cc_vec3f_t bO =
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
computePaPb(atmo_solverParam_t* param, cc_vec3f_t* P0,
            cc_vec3f_t* V, double Ro, cc_vec3f_t* Pa,
            cc_vec3f_t* Pb)
{
	ASSERT(param);
	ASSERT(P0);
	ASSERT(V);
	ASSERT(Pa);
	ASSERT(Pb);

	// initialize spheres
	// optionally include a radius offset to ensure that the
	// viewing vector does not intersect at P0
	cc_sphere_t sphereP;
	cc_sphere_t sphereA;
	cc_sphere_load(&sphereP, 0.0f, 0.0f, 0.0f, param->Rp - Ro);
	cc_sphere_load(&sphereA, 0.0f, 0.0f, 0.0f, param->Ra + Ro);

	// initialize ray
	cc_ray3f_t ray;
	cc_ray3f_load(&ray, P0->x, P0->y, P0->z, V->x, V->y, V->z);

	// intersect ray
	int   countA;
	int   countP;
	float nearA = 0.0f;
	float farA  = 0.0f;
	float nearP = 0.0f;
	float farP  = 0.0f;
	countA = cc_ray3f_intersect(&ray, &sphereA,
	                            &nearA, &farA);
	countP = cc_ray3f_intersect(&ray, &sphereP,
	                            &nearP, &farP);

	// check if ray intersects atmosphere
	if(countA)
	{
		cc_ray3f_getpoint(&ray, nearA, Pa);

		// check if ray intersects planet
		if(countP)
		{
			cc_ray3f_getpoint(&ray, nearP, Pb);
			return 0;
		}

		cc_ray3f_getpoint(&ray, farA, Pb);
		return 1;
	}

	cc_vec3f_copy(P0, Pa);
	cc_vec3f_copy(P0, Pb);
	return 0;
}

// compute the point Pc on the light vector
static int
computePc(atmo_solverParam_t* param, cc_vec3f_t* P,
          cc_vec3f_t* L, cc_vec3f_t* Pc)
{
	ASSERT(param);
	ASSERT(P);
	ASSERT(L);
	ASSERT(Pc);

	cc_vec3f_t Pa;
	cc_vec3f_t Sun;
	cc_vec3f_muls_copy(L, -1.0f, &Sun);
	double Ro = 10.0;
	return computePaPb(param, P, &Sun, Ro, &Pa, Pc);
}

// modified Rayleigh phase function
static float
atmo_phaseR(atmo_solverParam_t* param, float cos_theta)
{
	ASSERT(param);

	return 0.8f*(1.4f + 0.5f*cos_theta*cos_theta);
}

// Mie phase function
static float
atmo_phaseM(atmo_solverParam_t* param, float cos_theta)
{
	ASSERT(param);

	float g  = param->phase_g_mie;
	float g2 = g*g;
	float n1 = 3.0f*(1.0f - g2);
	float n2 = 1.0f + cos_theta*cos_theta;
	float d1 = 2.0f*(2.0f + g2);
	float d2 = powf(1.0f + g2 - 2.0f*g*cos_theta, 1.5f);
	return (n1/d1)*(n2/d2);
}

// factored single-scattered intensity
static void
fIS1(atmo_solverParam_t* param, float h, float phi,
     float delta, cc_vec4f_t* fis1)
{
	ASSERT(param);
	ASSERT(fis1);

	// initialize fis1
	fis1->r = 0.0f;
	fis1->g = 0.0f;
	fis1->b = 0.0f;
	fis1->a = 0.0f;

	// canonical form of the scattering intensity
	// parameterization for P0, V and L
	cc_vec3f_t P0 =
	{
		.z = h + param->Rp,
	};
	cc_vec3f_t V =
	{
		.x = sinf(phi),
		.z = cosf(phi),
	};
	cc_vec3f_t L =
	{
		.x = -sinf(delta),
		.z = -cosf(delta),
	};
	cc_vec3f_normalize(&V);
	cc_vec3f_normalize(&L);

	// compute ray-sphere intersection
	// include a ray offset for the viewing vector
	// to ensure that the ray does not intersect at P0
	double Ro = 10.0;
	cc_vec3f_t Pa;
	cc_vec3f_t Pb;
	cc_vec3f_t Pc;
	computePaPb(param, &P0, &V, Ro, &Pa, &Pb);
	if(computePc(param, &Pa, &L, &Pc) == 0)
	{
		// ray intersects planet
		return;
	}

	// initialize integration
	cc_vec3f_t P;
	cc_vec3f_t step;
	cc_vec4f_t fx0;
	cc_vec4f_t fx1;
	cc_vec3f_t tPPc;
	cc_vec3f_t tPaP;
	cc_vec3f_copy(&Pa, &P);
	cc_vec3f_subv_copy(&Pb, &Pa, &step);
	cc_vec3f_muls(&step, 1.0f/param->transmittance_steps);
	float ds = cc_vec3f_mag(&step);
	float pR = densityR(param, h);
	float pM = densityM(param, h);
	transmittance(param, &P, &Pc, &tPPc);
	transmittance(param, &Pa, &P, &tPaP);
	fx0.r = pR*expf(-tPPc.x -tPaP.x);
	fx0.g = pR*expf(-tPPc.y -tPaP.y);
	fx0.b = pR*expf(-tPPc.z -tPaP.z);
	fx0.a = pM*expf(-(tPPc.x + tPPc.y + tPPc.z)/3.0f
	                -(tPaP.y + tPaP.y + tPaP.z)/3.0f);

	if(ds < ATMO_STEP_THRESH)
	{
		return;
	}

	// integrate factored single-scattered intensity
	int i;
	for(i = 0; i < param->transmittance_steps; ++i)
	{
		cc_vec3f_addv(&P, &step);

		if(computePc(param, &P, &L, &Pc) == 0)
		{
			// PPc intersects planet
			fx0.r = 0.0f;
			fx0.g = 0.0f;
			fx0.b = 0.0f;
			fx0.a = 0.0f;
			continue;
		}

		h  = getHeightP(param, &P);
		pR = densityR(param, h);
		pM = densityM(param, h);
		transmittance(param, &P, &Pc, &tPPc);
		transmittance(param, &Pa, &P, &tPaP);
		fx1.r = pR*expf(-tPPc.x -tPaP.x);
		fx1.g = pR*expf(-tPPc.y -tPaP.y);
		fx1.b = pR*expf(-tPPc.z -tPaP.z);
		fx1.a = pM*expf(-(tPPc.x + tPPc.y + tPPc.z)/3.0f
		                -(tPaP.x + tPaP.y + tPaP.z)/3.0f);

		// apply trapesoidal rule
		fis1->r += 0.5f*(fx1.r + fx0.r)*ds;
		fis1->g += 0.5f*(fx1.g + fx0.g)*ds;
		fis1->b += 0.5f*(fx1.b + fx0.b)*ds;
		fis1->a += 0.5f*(fx1.a + fx0.a)*ds;

		cc_vec4f_copy(&fx1, &fx0);
	}

	// apply Rayleigh/Mie scattering coefficient
	fis1->r *= param->beta_s_r_rayleigh/(4.0*M_PI);
	fis1->g *= param->beta_s_g_rayleigh/(4.0*M_PI);
	fis1->b *= param->beta_s_b_rayleigh/(4.0*M_PI);
	fis1->a *= param->beta_s_mie/(4.0*M_PI);
}

static void
fISk_sample(atmo_solverParam_t* param, uint32_t k,
            cc_vec3f_t* P, cc_vec3f_t* V, cc_vec3f_t* L,
            cc_vec4f_t* data, cc_vec4f_t* fisk)
{
	ASSERT(param);
	ASSERT(P);
	ASSERT(V);
	ASSERT(L);
	ASSERT(data);
	ASSERT(fisk);

	// initialize fisk
	fisk->r = 0.0f;
	fisk->g = 0.0f;
	fisk->b = 0.0f;
	fisk->a = 0.0f;

	cc_vec3f_t Zenith;
	getZenithP(P, &Zenith);

	cc_vec3f_t Sun;
	cc_vec3f_muls_copy(L, -1.0f, &Sun);

	// compute u,v,w
	float h         = getHeightP(param, P);
	float cos_phi   = cc_vec3f_dot(&Zenith, V);
	float cos_delta = cc_vec3f_dot(&Zenith, &Sun);
	float u         = getUHeight(param, h);
	float v         = getVCosPhi(param, h, cos_phi, u);
	float w         = getWCosDelta(cos_delta);
	float width1    = (float) (param->texture_width  - 1);
	float height1   = (float) (param->texture_height - 1);
	float depth1    = (float) (param->texture_depth  - 1);

	// compute x0,y0,z0
	uint32_t x0;
	uint32_t y0;
	uint32_t z0;
	x0 = atmo_clamp((uint32_t) (u*width1),
	                0, param->texture_width  - 1);
	y0 = atmo_clamp((uint32_t) (v*height1),
	                0, param->texture_height - 1);
	z0 = atmo_clamp((uint32_t) (w*depth1),
	                0, param->texture_depth  - 1);

	#if ATMO_SAMPLE_MODE == ATMO_SAMPLE_MODE_LINEAR
	// compute x1,y1,z1
	uint32_t x1;
	uint32_t y1;
	uint32_t z1;
	x1 = atmo_clamp(x0 + 1, 0, param->texture_width  - 1);
	y1 = atmo_clamp(y0 + 1, 0, param->texture_height - 1);
	z1 = atmo_clamp(z0 + 1, 0, param->texture_depth  - 1);

	// sampling coordinates
	float uu = u*width1  - ((float) x0);
	float vv = v*height1 - ((float) y0);
	float ww = w*depth1  - ((float) z0);

	// sample corners
	cc_vec4f_t fisk000;
	cc_vec4f_t fisk001;
	cc_vec4f_t fisk010;
	cc_vec4f_t fisk011;
	cc_vec4f_t fisk100;
	cc_vec4f_t fisk101;
	cc_vec4f_t fisk110;
	cc_vec4f_t fisk111;
	getData(param, k, x0, y0, z0, data, &fisk000);
	getData(param, k, x0, y0, z1, data, &fisk001);
	getData(param, k, x0, y1, z0, data, &fisk010);
	getData(param, k, x0, y1, z1, data, &fisk011);
	getData(param, k, x1, y0, z0, data, &fisk100);
	getData(param, k, x1, y0, z1, data, &fisk101);
	getData(param, k, x1, y1, z0, data, &fisk110);
	getData(param, k, x1, y1, z1, data, &fisk111);

	// interpolate x
	cc_vec4f_t fiskx00;
	cc_vec4f_t fiskx01;
	cc_vec4f_t fiskx10;
	cc_vec4f_t fiskx11;
	cc_vec4f_lerp(&fisk000, &fisk100, uu, &fiskx00);
	cc_vec4f_lerp(&fisk001, &fisk101, uu, &fiskx01);
	cc_vec4f_lerp(&fisk010, &fisk110, uu, &fiskx10);
	cc_vec4f_lerp(&fisk011, &fisk111, uu, &fiskx11);

	// interpolate y
	cc_vec4f_t fiskxy0;
	cc_vec4f_t fiskxy1;
	cc_vec4f_lerp(&fiskx00, &fiskx10, vv, &fiskxy0);
	cc_vec4f_lerp(&fiskx01, &fiskx11, vv, &fiskxy1);

	// interpolate z
	cc_vec4f_lerp(&fiskxy0, &fiskxy1, ww, fisk);

	#else
	getData(param, k, x0, y0, z0, data, fisk);
	#endif
}

// factored multiple-scattered gathered intensity step
static void
fGk_step(atmo_solverParam_t* param, uint32_t k,
         cc_vec3f_t* P, cc_vec3f_t* V, cc_vec3f_t* L,
         cc_vec4f_t* data, float s, float xj, float yi,
         cc_vec4f_t* fgk)
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
	cc_vec3f_t omega =
	{
		.x = sinf(xj)*cosf(yi),
		.y = sinf(xj)*sinf(yi),
		.z = cosf(xj),
	};
	cc_vec3f_normalize(&omega);

	// compute phase
	float cos_theta;
	float FR;
	float FM;
	cc_vec3f_t minus_omega;
	cc_vec3f_t minus_V;
	cc_vec3f_muls_copy(&omega, -1.0f, &minus_omega);
	cc_vec3f_muls_copy(V,      -1.0f, &minus_V);
	cos_theta = cc_vec3f_dot(&minus_omega, &minus_V);
	FR = atmo_phaseR(param, cos_theta);
	FM = atmo_phaseM(param, cos_theta);

	// sample fisk
	cc_vec4f_t fisk;
	fISk_sample(param, k, P, &omega, L, data, &fisk);

	// compute sin_theta for domega
	// xj => omega spherical angle theta in (0, pi)
	float domega_sin_xj = sinf(xj);

	// add factored multiple-scattered gathered intensity
	fgk->r += FR*fisk.r*domega_sin_xj;
	fgk->g += FR*fisk.g*domega_sin_xj;
	fgk->b += FR*fisk.b*domega_sin_xj;
	fgk->a += FM*fisk.a*domega_sin_xj;
}

// factored multiple-scattered gathered intensity
static void
fGk(atmo_solverParam_t* param, uint32_t k, cc_vec3f_t* P,
    cc_vec3f_t* V, cc_vec3f_t* L, cc_vec4f_t* data,
    cc_vec4f_t* fgk)
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
	float x0 = 0.0f;
	float xn = (float) M_PI;
	float y0 = 0.0f;
	float ym = (float) (2.0*M_PI);

	// compute 2D trapezoidal rule step size
	int   m  = param->gather_m_steps;
	int   n  = param->gather_n_steps;
	float dx = (xn - x0)/((float) n);
	float dy = (ym - y0)/((float) m);

	// initialize fgk
	fgk->r = 0.0f;
	fgk->g = 0.0f;
	fgk->b = 0.0f;
	fgk->a = 0.0f;

	// apply 2D trapezoidal rule for corners
	fGk_step(param, k, P, V, L, data, 1.0f, x0, y0, fgk);
	fGk_step(param, k, P, V, L, data, 1.0f, xn, y0, fgk);
	fGk_step(param, k, P, V, L, data, 1.0f, x0, ym, fgk);
	fGk_step(param, k, P, V, L, data, 1.0f, xn, ym, fgk);

	// apply 2D trapezoidal rule for edges
	int   i;
	int   j;
	float yi;
	float xj;
	for(j = 1; j < n; ++j)
	{
		xj = ((float) j)*dx;
		fGk_step(param, k, P, V, L, data, 2.0f, xj, y0, fgk);
		fGk_step(param, k, P, V, L, data, 2.0f, xj, ym, fgk);
	}
	for(i = 1; i < m; ++i)
	{
		yi = ((float) i)*dy;
		fGk_step(param, k, P, V, L, data, 2.0f, x0, yi, fgk);
		fGk_step(param, k, P, V, L, data, 2.0f, xn, yi, fgk);
	}

	// apply 2D trapezoidal rule for center
	for(i = 1; i < m; ++i)
	{
		yi = ((float) i)*dy;
		for(j = 1; j < n; ++j)
		{
			xj = ((float) j)*dx;
			fGk_step(param, k, P, V, L, data, 4.0f, xj, yi, fgk);
		}
	}

	// apply 2D trapezoidal rule scale
	cc_vec4f_muls(fgk, 4.0f*dx*dy);
}

// factored multiple-scattered intensity
static void
fISk(atmo_solverParam_t* param, uint32_t k,
     float h, float phi, float delta, cc_vec4f_t* data,
     cc_vec4f_t* fisk)
{
	ASSERT(param);
	ASSERT(data);
	ASSERT(fisk);

	// initialize fisk
	fisk->r = 0.0f;
	fisk->g = 0.0f;
	fisk->b = 0.0f;
	fisk->a = 0.0f;

	// canonical form of the scattering intensity
	// parameterization for P0, V and L
	cc_vec3f_t P0 =
	{
		.z = h + param->Rp,
	};
	cc_vec3f_t V =
	{
		.x = sinf(phi),
		.z = cosf(phi),
	};
	cc_vec3f_t L =
	{
		.x = -sinf(delta),
		.z = -cosf(delta),
	};
	cc_vec3f_normalize(&V);
	cc_vec3f_normalize(&L);

	// compute ray-sphere intersection
	// include a ray offset for the viewing vector
	// to ensure that the ray does not intersect at P0
	double Ro = 10.0;
	cc_vec3f_t Pa;
	cc_vec3f_t Pb;
	computePaPb(param, &P0, &V, Ro, &Pa, &Pb);

	// initialize integration
	cc_vec3f_t P;
	cc_vec3f_t step;
	cc_vec4f_t fx0;
	cc_vec4f_t fx1;
	cc_vec3f_t tPaP;
	cc_vec4f_t fgk;
	cc_vec3f_copy(&Pa, &P);
	cc_vec3f_subv_copy(&Pb, &Pa, &step);
	cc_vec3f_muls(&step, 1.0f/param->transmittance_steps);
	float ds = cc_vec3f_mag(&step);
	float pR = densityR(param, h);
	float pM = densityM(param, h);
	fGk(param, k - 1, &P0, &V, &L, data, &fgk);
	transmittance(param, &Pa, &P, &tPaP);
	fx0.r = fgk.r*pR*expf(-tPaP.x);
	fx0.g = fgk.g*pR*expf(-tPaP.y);
	fx0.b = fgk.b*pR*expf(-tPaP.z);
	fx0.a = fgk.a*pM*expf(-(tPaP.x + tPaP.y + tPaP.z)/3.0f);

	if(ds < ATMO_STEP_THRESH)
	{
		return;
	}

	// integrate factored multiple-scattered intensity
	int i;
	for(i = 0; i < param->transmittance_steps; ++i)
	{
		cc_vec3f_addv(&P, &step);

		h  = getHeightP(param, &P);
		pR = densityR(param, h);
		pM = densityM(param, h);
		fGk(param, k - 1, &P, &V, &L, data, &fgk);
		transmittance(param, &Pa, &P, &tPaP);
		fx1.r = fgk.r*pR*expf(-tPaP.x);
		fx1.g = fgk.g*pR*expf(-tPaP.y);
		fx1.b = fgk.b*pR*expf(-tPaP.z);
		fx1.a = fgk.a*pM*expf(-(tPaP.x + tPaP.y + tPaP.z)/3.0f);

		// apply trapesoidal rule
		fisk->r += 0.5f*(fx1.r + fx0.r)*ds;
		fisk->g += 0.5f*(fx1.g + fx0.g)*ds;
		fisk->b += 0.5f*(fx1.b + fx0.b)*ds;
		fisk->a += 0.5f*(fx1.a + fx0.a)*ds;

		cc_vec4f_copy(&fx1, &fx0);
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
                      cc_vec4f_t* data)
{
	ASSERT(self);
	ASSERT(data);

	vkk_image_t* img;
	cc_vec4f_t*  datak;

	atmo_solverParam_t* param = &self->param;

	self->image_array = (vkk_image_t**)
	                    CALLOC(param->k,
	                           sizeof(vkk_image_t*));
	if(self->image_array == NULL)
	{
		LOGE("CALLOC failed");
		return 0;
	}

	int i;
	for(i = 0; i < param->k; ++i)
	{
		// k is base-1
		datak = getDataK(&self->param, i + 1, data);

		img = vkk_image_new(self->engine,
		                    param->texture_width,
		                    param->texture_height,
		                    param->texture_depth,
		                    VKK_IMAGE_FORMAT_RGBAF16,
		                    0, VKK_STAGE_FS,
		                    (const void*) datak);
		if(img == NULL)
		{
			goto failure;
		}

		self->image_array[i] = img;
	}

	// success
	return 1;

	// failure
	failure:
	{
		for(i = 0; i < param->k; ++i)
		{
			vkk_image_delete(&self->image_array[i]);
		}
		FREE(self->image_array);
		self->image_array = NULL;
	}
	return 0;
}

static int
atmo_solver_exportData(atmo_solver_t* self,
                       cc_vec4f_t* data)
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
	cc_jsmnStream_float(jsmn, param->Rp);
	cc_jsmnStream_key(jsmn, "%s", "Ra");
	cc_jsmnStream_float(jsmn, param->Ra);
	cc_jsmnStream_key(jsmn, "%s", "density_scale_height_rayleigh");
	cc_jsmnStream_float(jsmn, param->density_scale_height_rayleigh);
	cc_jsmnStream_key(jsmn, "%s", "density_scale_height_mie");
	cc_jsmnStream_float(jsmn, param->density_scale_height_mie);
	cc_jsmnStream_key(jsmn, "%s", "phase_g_mie");
	cc_jsmnStream_float(jsmn, param->phase_g_mie);
	cc_jsmnStream_key(jsmn, "%s", "beta_s_r_rayleigh");
	cc_jsmnStream_float(jsmn, param->beta_s_r_rayleigh);
	cc_jsmnStream_key(jsmn, "%s", "beta_s_g_rayleigh");
	cc_jsmnStream_float(jsmn, param->beta_s_g_rayleigh);
	cc_jsmnStream_key(jsmn, "%s", "beta_s_b_rayleigh");
	cc_jsmnStream_float(jsmn, param->beta_s_b_rayleigh);
	cc_jsmnStream_key(jsmn, "%s", "beta_s_mie");
	cc_jsmnStream_float(jsmn, param->beta_s_mie);
	cc_jsmnStream_key(jsmn, "%s", "beta_a_mie");
	cc_jsmnStream_float(jsmn, param->beta_a_mie);
	cc_jsmnStream_key(jsmn, "%s", "beta_a_r_ozone");
	cc_jsmnStream_float(jsmn, param->beta_a_r_ozone);
	cc_jsmnStream_key(jsmn, "%s", "beta_a_g_ozone");
	cc_jsmnStream_float(jsmn, param->beta_a_g_ozone);
	cc_jsmnStream_key(jsmn, "%s", "beta_a_b_ozone");
	cc_jsmnStream_float(jsmn, param->beta_a_b_ozone);
	cc_jsmnStream_key(jsmn, "%s", "II_r");
	cc_jsmnStream_float(jsmn, param->II_r);
	cc_jsmnStream_key(jsmn, "%s", "II_g");
	cc_jsmnStream_float(jsmn, param->II_g);
	cc_jsmnStream_key(jsmn, "%s", "II_b");
	cc_jsmnStream_float(jsmn, param->II_b);
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
	cc_vec4f_t* datak;
	char fname[256];
	FILE* f;
	size_t size = sizeof(cc_vec4f_t)*param->texture_width*
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
static float atmo_solver_luminance(cc_vec4f_t* v)
{
	cc_vec4f_t w =
	{
		.r = 0.2126f,
		.g = 0.7152f,
		.b = 0.0722f,
	};
	return cc_vec4f_dot(v, &w);
}

// https://64.github.io/tonemapping/
static float atmo_solver_reinhard(float x)
{
	return x / (1.0f + x);
}

// https://64.github.io/tonemapping/
static void
atmo_solver_reinhardLuminance(cc_vec4f_t* color)
{
	ASSERT(color);

	float l0 = atmo_solver_luminance(color);
	float l1 = atmo_solver_reinhard(l0);
	cc_vec4f_muls(color, l1/l0);
}

static void atmo_solver_gamma(cc_vec4f_t* color)
{
	ASSERT(color);

	color->r = powf(color->r, 1.0f/2.2f);
	color->g = powf(color->g, 1.0f/2.2f);
	color->b = powf(color->b, 1.0f/2.2f);
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
	float Ha = param->Ra - param->Rp;
	for(i = 0; i < nh; ++i)
	{
		float h;
		h = Ha*((float) i)/((float) (nh - 1));

		fprintf(f, "%f, %f, %f\n",
		        densityR(param, h),
		        densityM(param, h),
		        densityO(param, h));
	}

	fclose(f);

	return 1;
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
	int   i;
	int   j;
	int   nh   = 40;
	int   nphi = 721;
	float Ha   = param->Ra - param->Rp;
	for(i = 0; i < nh; ++i)
	{
		float h;
		h = cc_clamp(Ha*((float) i)/((float) (nh - 1)),
		             10.0f, Ha - 10.0f);

		for(j = 0; j < nphi; ++j)
		{
			float phi_deg;
			phi_deg = 180.0f*((float) j)/
			                 ((float) (nphi - 1));

			float phi_rad;
			phi_rad = cc_deg2rad(phi_deg);

			// canonical form of the scattering intensity
			// parameterization for P0 and V
			cc_vec3f_t P0 =
			{
				.z = h + param->Rp,
			};
			cc_vec3f_t V =
			{
				.x = sinf(phi_rad),
				.z = cosf(phi_rad),
			};
			cc_vec3f_normalize(&V);

			// compute ray-sphere intersection
			// include a ray offset for the viewing vector
			// to ensure that the ray does not intersect at P0
			double Ro = 10.0;
			cc_vec3f_t Pa;
			cc_vec3f_t Pb;
			computePaPb(param, &P0, &V, Ro, &Pa, &Pb);

			cc_vec3f_t Vba;
			cc_vec3f_subv_copy(&Pb, &Pa, &Vba);
			float dist = cc_vec3f_mag(&Vba);

			cc_vec3f_t t;
			transmittance(param, &P0, &Pb, &t);

			cc_vec3f_t T;
			T.x = exp(-t.x);
			T.y = exp(-t.y);
			T.z = exp(-t.z);

			float avgT = (T.x + T.y + T.z)/3.0f;

			LOGI("h=%f, phi=%f, avgT=%f, dist=%f, T=%f,%f,%f",
			     h, phi_deg, avgT, dist,
			     T.x, T.y, T.z);

			if(j == 0)
			{
				fprintf(f, "%f", avgT);
			}
			else
			{
				fprintf(f, " %f", avgT);
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
		fprintf(f, "%f %f %f %f %f %f %f %f %f %f\n",
		        (float) i,
		        rgb_none.x, rgb_none.y, rgb_none.z,
		        rgb_sum.x,  rgb_sum.y,  rgb_sum.z,
		        rgb_peak.x, rgb_peak.y, rgb_peak.z);
	}

	fclose(f);

	return 1;
}

static int
atmo_solver_debugData(atmo_solver_t* self, cc_vec4f_t* data)
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

	// spectral intensity of of incident light from the Sun
	float exposure = 0.0f;
	cc_vec4f_t II =
	{
		.r = powf(2.0f, exposure)*param->II_r,
		.g = powf(2.0f, exposure)*param->II_g,
		.b = powf(2.0f, exposure)*param->II_b,
	};

	LOGI("II: r=%f, g=%f, b=%f", II.r, II.g, II.b);

	float phi;
	float delta;
	float cos_phi;
	float cos_delta;
	float cos_theta;
	float h;
	float u;
	float v;
	float w;
	float FR;
	float FM;
	uint32_t k;
	uint32_t x;
	uint32_t y;
	uint32_t z;
	char fname[256];
	unsigned char pixel[4] = { 0, 0, 0, 255 };
	cc_vec4f_t fis;
	cc_vec4f_t is;
	cc_vec4f_t color;
	float      luminance;
	float      white_point;
	for(k = 1; k <= param->k; ++k)
	{
		white_point = 0.0f;
		for(x = 0; x < param->texture_width; ++x)
		{
			u = ((float) x)/((float) (param->texture_width - 1));
			h = getHeightU(param, u);

			for(z = 0; z < param->texture_depth; ++z)
			{
				w = ((float) z)/((float) (param->texture_depth - 1));
				cos_delta = getCosDeltaW(w);
				delta = acosf(cos_delta);

				for(y = 0; y < param->texture_height; ++y)
				{
					v = ((float) y)/((float) (param->texture_height - 1));
					cos_phi = getCosPhiV(param, h, u, v);
					phi = acosf(cos_phi);

					getData(param, k, x, y, z, data, &fis);

					cos_theta = cosf(delta - phi);
					FR = atmo_phaseR(param, cos_theta);
					FM = atmo_phaseM(param, cos_theta);

					is.r = II.r*(FR*fis.r + FM*fis.a);
					is.g = II.g*(FR*fis.g + FM*fis.a);
					is.b = II.b*(FR*fis.b + FM*fis.a);

					// find white point
					luminance = atmo_solver_luminance(&is);
					if(luminance > white_point)
					{
						white_point = luminance;
					}

					// tone mapping and gamma correction
					cc_vec4f_copy(&is, &color);
					atmo_solver_reinhardLuminance(&color);
					atmo_solver_gamma(&color);

					// clamp output
					color.r = cc_clamp(color.r, 0.0f, 1.0f);
					color.g = cc_clamp(color.g, 0.0f, 1.0f);
					color.b = cc_clamp(color.b, 0.0f, 1.0f);

					// set output
					pixel[0] = (unsigned char) (255.0f*color.r);
					pixel[1] = (unsigned char) (255.0f*color.g);
					pixel[2] = (unsigned char) (255.0f*color.b);
					texgz_tex_setPixel(tex, x*param->texture_depth + z,
					                   y, pixel);

					// set u,v,w
					pixel[0] = (unsigned char) (255.0f*u);
					pixel[1] = (unsigned char) (255.0f*v);
					pixel[2] = (unsigned char) (255.0f*w);
					texgz_tex_setPixel(tex, x*param->texture_depth + z,
					                   param->texture_height + y, pixel);
				}
			}
		}
		snprintf(fname, 256, "atmo-slice-k%u.png", k);
		texgz_png_export(tex, fname);

		LOGI("k=%u, white_point=%f", k, white_point);
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
atmo_solver_debugData(atmo_solver_t* self, cc_vec4f_t* data)
{
	return 1;
}
#endif

static int
atmo_solver_paramValidate(atmo_solverParam_t* param)
{
	ASSERT(param);

	// perform minimalistic error checking

	if((param->Rp < 1.0f) || (param->Ra < 1.0f))
	{
		LOGE("invalid Rp=%f, Ra=%f",
		     param->Rp, param->Ra);
		return 0;
	}

	if((param->density_scale_height_rayleigh < 0.0f) ||
	   (param->density_scale_height_mie      < 0.0f))
	{
		LOGE("invalid Rp=%f, Ra=%f",
		     param->density_scale_height_rayleigh,
		     param->density_scale_height_mie);
		return 0;
	}

	if((param->phase_g_mie < -0.9999f) ||
	   (param->phase_g_mie > 0.9999f))
	{
		LOGE("invalid phase_g_mie=%f",
		     param->phase_g_mie);
		return 0;
	}

	if((param->beta_s_r_rayleigh <= 0.0f) ||
	   (param->beta_s_g_rayleigh <= 0.0f) ||
	   (param->beta_s_b_rayleigh <= 0.0f) ||
	   (param->beta_s_mie        <= 0.0f) ||
	   (param->beta_a_mie        <= 0.0f) ||
	   (param->beta_a_r_ozone    <= 0.0f) ||
	   (param->beta_a_g_ozone    <= 0.0f) ||
	   (param->beta_a_b_ozone    <= 0.0f))
	{
		LOGE("invalid beta_s_r/g/b_rayleigh=%f/%f/%f, beta_s_mie=%f, beta_a_mie=%f, beta_a_r/g/b_ozone=%f/%f/%f",
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

	if((param->II_r <= 0.0f) ||
	   (param->II_g <= 0.0f) ||
	   (param->II_b <= 0.0f))
	{
		LOGE("invalid II_r/g/b=%f/%f/%f",
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
                 cc_vec4f_t* data)
{
	ASSERT(self);
	ASSERT(data);

	atmo_solverParam_t* param = &self->param;

	float width  = (float) param->texture_width;
	float height = (float) param->texture_height;
	float depth  = (float) param->texture_depth;

	float u = ((float) x)/(width - 1.0f);
	float v = ((float) y)/(height - 1.0f);
	float w = ((float) z)/(depth - 1.0f);

	float h         = getHeightU(param, u);
	float cos_phi   = getCosPhiV(param, h, u, v);
	float cos_delta = getCosDeltaW(w);
	float phi       = acos(cos_phi);
	float delta     = acos(cos_delta);

	cc_vec4f_t fis;
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
atmo_solver_finish(atmo_solver_t* self, cc_vec4f_t* data)
{
	ASSERT(self);
	ASSERT(data);

	atmo_solverParam_t* param = &self->param;

	// compute total intensity for each k
	uint32_t k;
	uint32_t x;
	uint32_t y;
	uint32_t z;
	cc_vec4f_t fis0;
	cc_vec4f_t fis1;
	cc_vec4f_t fis;
	for(k = 2; k <= param->k; ++k)
	{
		for(z = 0; z < param->texture_depth; ++z)
		{
			for(y = 0; y < param->texture_height; ++y)
			{
				for(x = 0; x < param->texture_width; ++x)
				{
					getData(param, k - 1, x, y, z, data, &fis0);
					getData(param, k,     x, y, z, data, &fis1);
					cc_vec4f_addv_copy(&fis0, &fis1, &fis);
					setData(param, k, x, y, z, data, &fis);
				}
			}
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

	cc_vec4f_t* data;
	data = (cc_vec4f_t*)
	       CALLOC(count, sizeof(cc_vec4f_t));
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
