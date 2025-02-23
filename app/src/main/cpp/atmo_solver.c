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

// radius of the planet and atmospheric boundary
// use doubles for radius for numerical stability
#define ATMO_RP 6371000.0
#define ATMO_RA 6471000.0

#define ATMO_DENSITY_SCALE_HEIGHT_RAYLEIGH 8000.0f
#define ATMO_DENSITY_SCALE_HEIGHT_MIE      1200.0f

#define ATMO_PHASE_G_MIE -0.85f

#define ATMO_BETA_R_RAYLEIGH 6.55e-6f
#define ATMO_BETA_G_RAYLEIGH 1.73e-5f
#define ATMO_BETA_B_RAYLEIGH 2.30e-5f
#define ATMO_BETA_MIE        2e-6f

#define ATMO_SPECTRAL_IRRADIANCE_R 0.1526f
#define ATMO_SPECTRAL_IRRADIANCE_G 0.191f
#define ATMO_SPECTRAL_IRRADIANCE_B 0.208f

#define ATMO_EXPOSURE 1.0f

#define ATMO_SPECTRAL_TO_RGB_R 133.3209f
#define ATMO_SPECTRAL_TO_RGB_G 88.51855f
#define ATMO_SPECTRAL_TO_RGB_B 112.7552f

#define ATMO_INTEGRATION_STEPS 30

#define ATMO_K 1

#define ATMO_TEXTURE_WIDTH  32
#define ATMO_TEXTURE_HEIGHT 256
#define ATMO_TEXTURE_DEPTH  32

// set defines to use non-linear parameterization
// requires corresponding change in sky_atmo.frag
#define ATMO_PARAM_NONLINEAR_H
// ATMO_PARAM_NONLINEAR_PHI seems buggy
//#define ATMO_PARAM_NONLINEAR_PHI
#define ATMO_PARAM_NONLINEAR_DELTA

/***********************************************************
* private - compute                                        *
***********************************************************/

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

// modified Rayleigh phase function
static float
phaseR(atmo_solverParam_t* param, float cos_theta)
{
	ASSERT(param);

	return 0.8f*(1.4f + 0.5f*cos_theta*cos_theta);
}

// Mie phase function
static float
phaseM(atmo_solverParam_t* param, float cos_theta)
{
	ASSERT(param);

	float g  = param->phase_g_mie;
	float g2 = g*g;
	float n1 = 3.0f*(1.0f - g2);
	float n2 = 1.0f + cos_theta*cos_theta;
	float d1 = 2.0f*(2.0f + g2);
	float d2 = powf(1.0f + g2 + 2.0f*g*cos_theta, 1.5f);
	return (n1/d1)*(n2/d2);
}

// compute height using double precision since the magnitude
// of points in the atmosphere produces very large numbers
static float
height(atmo_solverParam_t* param, cc_vec3f_t* P)
{
	ASSERT(P);

	cc_vec3d_t d =
	{
		.x = P->x,
		.y = P->y,
		.z = P->z,
	};

	return (float) (cc_vec3d_mag(&d) - param->Rp);
}

// Rayleigh/Mie transmittance
static void
transmittance(atmo_solverParam_t* param, cc_vec3f_t* P1,
              cc_vec3f_t* P2, cc_vec4f_t* out)
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
	cc_vec3f_muls(&step, 1.0f/param->integration_steps);
	float ds  = cc_vec3f_mag(&step);
	float h   = height(param, &P);
	float pR0 = densityR(param, h);
	float pM0 = densityM(param, h);

	// integrate transmittance
	int   i;
	float pR1;
	float pM1;
	float tR = 0.0f;
	float tM = 0.0f;
	for(i = 0; i < param->integration_steps; ++i)
	{
		cc_vec3f_addv(&P, &step);

		h   = height(param, &P);
		pR1 = densityR(param, h);
		pM1 = densityM(param, h);

		// apply trapesoidal rule
		tR += 0.5f*(pR0 + pR1)*ds;
		tM += 0.5f*(pM0 + pM1)*ds;

		pR0 = pR1;
		pM0 = pM1;
	}

	// apply Rayleigh/Mie scattering coefficient
	out->r = param->beta_r_rayleigh*tR;
	out->g = param->beta_g_rayleigh*tR;
	out->b = param->beta_b_rayleigh*tR;
	out->a = param->beta_mie*tM;
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

// Rayleigh/Mie factored single-scattered intensity
static void
fIS1(atmo_solverParam_t* param, float h, float phi,
     float delta, cc_vec4f_t* out)
{
	ASSERT(param);
	ASSERT(out);

	// initialize out
	out->r = 0.0f;
	out->g = 0.0f;
	out->b = 0.0f;
	out->a = 0.0f;

	// canonical form of the scattering intensity
	// parameterization for P0, V and L
	cc_vec3f_t P0 =
	{
		.z = h + param->Rp,
	};
	cc_vec3f_t V =
	{
		.x = sin(phi),
		.z = cos(phi),
	};
	cc_vec3f_t L =
	{
		.x = -sin(delta),
		.z = -cos(delta),
	};

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
	cc_vec4f_t tPPc;
	cc_vec4f_t tPaP;
	cc_vec3f_copy(&Pa, &P);
	cc_vec3f_subv_copy(&Pb, &Pa, &step);
	cc_vec3f_muls(&step, 1.0f/param->integration_steps);
	float ds = cc_vec3f_mag(&step);
	float pR = densityR(param, h);
	float pM = densityM(param, h);
	transmittance(param, &P, &Pc, &tPPc);
	transmittance(param, &Pa, &P, &tPaP);
	fx0.r = pR*exp(-tPPc.r -tPaP.r);
	fx0.g = pR*exp(-tPPc.g -tPaP.g);
	fx0.b = pR*exp(-tPPc.b -tPaP.b);
	fx0.a = pM*exp(-tPPc.a -tPaP.a);

	// integrate factored single-scattered intensity
	int i;
	for(i = 0; i < param->integration_steps; ++i)
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

		h  = height(param, &P);
		pR = densityR(param, h);
		pM = densityM(param, h);
		transmittance(param, &P, &Pc, &tPPc);
		transmittance(param, &Pa, &P, &tPaP);
		fx1.r = pR*exp(-tPPc.r -tPaP.r);
		fx1.g = pR*exp(-tPPc.g -tPaP.g);
		fx1.b = pR*exp(-tPPc.b -tPaP.b);
		fx1.a = pM*exp(-tPPc.a -tPaP.a);

		// apply trapesoidal rule
		out->r += 0.5f*(fx1.r + fx0.r)*ds;
		out->g += 0.5f*(fx1.g + fx0.g)*ds;
		out->b += 0.5f*(fx1.b + fx0.b)*ds;
		out->a += 0.5f*(fx1.a + fx0.a)*ds;

		cc_vec4f_copy(&fx1, &fx0);
	}

	// apply Rayleigh/Mie scattering coefficient
	out->r *= param->beta_r_rayleigh/(4.0*M_PI);
	out->g *= param->beta_g_rayleigh/(4.0*M_PI);
	out->b *= param->beta_b_rayleigh/(4.0*M_PI);
	out->a *= param->beta_mie/(4.0*M_PI);
}

/***********************************************************
* private - utility                                        *
***********************************************************/

static cc_vec4f_t*
atmo_getDataK(atmo_solverParam_t* param,
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
atmo_getData(atmo_solverParam_t* param,
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
atmo_setData(atmo_solverParam_t* param,
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
		datak = atmo_getDataK(&self->param, i + 1, data);

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
	cc_jsmnStream_key(jsmn, "%s", "beta_r_rayleigh");
	cc_jsmnStream_float(jsmn, param->beta_r_rayleigh);
	cc_jsmnStream_key(jsmn, "%s", "beta_g_rayleigh");
	cc_jsmnStream_float(jsmn, param->beta_g_rayleigh);
	cc_jsmnStream_key(jsmn, "%s", "beta_b_rayleigh");
	cc_jsmnStream_float(jsmn, param->beta_b_rayleigh);
	cc_jsmnStream_key(jsmn, "%s", "beta_mie");
	cc_jsmnStream_float(jsmn, param->beta_mie);
	cc_jsmnStream_key(jsmn, "%s", "spectral_irradiance_r");
	cc_jsmnStream_float(jsmn, param->spectral_irradiance_r);
	cc_jsmnStream_key(jsmn, "%s", "spectral_irradiance_g");
	cc_jsmnStream_float(jsmn, param->spectral_irradiance_g);
	cc_jsmnStream_key(jsmn, "%s", "spectral_irradiance_b");
	cc_jsmnStream_float(jsmn, param->spectral_irradiance_b);
	cc_jsmnStream_key(jsmn, "%s", "exposure");
	cc_jsmnStream_float(jsmn, param->exposure);
	cc_jsmnStream_key(jsmn, "%s", "spectral_to_rgb_r");
	cc_jsmnStream_float(jsmn, param->spectral_to_rgb_r);
	cc_jsmnStream_key(jsmn, "%s", "spectral_to_rgb_g");
	cc_jsmnStream_float(jsmn, param->spectral_to_rgb_g);
	cc_jsmnStream_key(jsmn, "%s", "spectral_to_rgb_b");
	cc_jsmnStream_float(jsmn, param->spectral_to_rgb_b);
	cc_jsmnStream_key(jsmn, "%s", "integration_steps");
	cc_jsmnStream_int(jsmn, param->integration_steps);
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
		datak = atmo_getDataK(param, k, data);
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

	if((param->beta_r_rayleigh <= 0.0f) ||
	   (param->beta_g_rayleigh <= 0.0f) ||
	   (param->beta_b_rayleigh <= 0.0f) ||
	   (param->beta_mie        <= 0.0f))
	{
		LOGE("invalid beta_r/g/b_rayleigh=%f/%f/%f, beta_mie=%f",
		     param->beta_r_rayleigh,
		     param->beta_g_rayleigh,
		     param->beta_b_rayleigh,
		     param->beta_mie);
		return 0;
	}

	if((param->spectral_irradiance_r <= 0.0f) ||
	   (param->spectral_irradiance_g <= 0.0f) ||
	   (param->spectral_irradiance_b <= 0.0f))
	{
		LOGE("invalid spectral_irradiance_r/g/b=%f/%f/%f",
		     param->spectral_irradiance_r,
		     param->spectral_irradiance_g,
		     param->spectral_irradiance_b);
		return 0;
	}

	if((param->exposure <= 0.0f) ||
	   (param->exposure > 1000.0f))
	{
		LOGE("invalid exposure=%f", param->exposure);
		return 0;
	}

	if((param->spectral_to_rgb_r <= 0.0f) ||
	   (param->spectral_to_rgb_g <= 0.0f) ||
	   (param->spectral_to_rgb_b <= 0.0f))
	{
		LOGE("invalid spectral_to_rgb_r/g/b=%f/%f/%f",
		     param->spectral_to_rgb_r,
		     param->spectral_to_rgb_g,
		     param->spectral_to_rgb_b);
		return 0;
	}

	if((param->integration_steps < 1) ||
	   (param->integration_steps > 100))
	{
		LOGE("invalid integration_steps=%u",
		     param->integration_steps);
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


	float Ra = param->Ra;
	float Rp = param->Rp;

	float h;
	#ifdef ATMO_PARAM_NONLINEAR_H
	h = u*u*(Ra - Rp);
	#else
	h = u*(Ra - Rp);
	#endif

	float phi;
	float cos_phi;
	#ifdef ATMO_PARAM_NONLINEAR_PHI
	float ch = -sqrtf(h*(2.0f*Rp + h))/(Rp + h);
	if(v > 0.5f)
	{
		cos_phi = ch + powf(v - 0.5f, 5.0f)*(1.0f - ch);
	}
	else
	{
		cos_phi = ch - powf(v, 5.0f)*(1.0f + ch);
	}
	#else
	cos_phi = 2.0f*v - 1.0f;
	#endif
	phi = acos(cos_phi);

	float delta;
	float cos_delta;
	#ifdef ATMO_PARAM_NONLINEAR_DELTA
	cos_delta = tan((2.0f*w - 1.0f + 0.26f)*0.75)/
	            tan(1.26f*0.75f);
	#else
	cos_delta = 2.0f*w - 1.0f;
	#endif
	delta = acos(cos_delta);

	cc_vec4f_t val;
	fIS1(param, h, phi, delta, &val);

	atmo_setData(param, k, x, y, z, data, &val);
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

	atmo_solver_newImages(self, data);
	atmo_solver_exportData(self, data);

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

		.beta_r_rayleigh = ATMO_BETA_R_RAYLEIGH,
		.beta_g_rayleigh = ATMO_BETA_G_RAYLEIGH,
		.beta_b_rayleigh = ATMO_BETA_B_RAYLEIGH,

		.beta_mie = ATMO_BETA_MIE,

		.spectral_irradiance_r = ATMO_SPECTRAL_IRRADIANCE_R,
		.spectral_irradiance_g = ATMO_SPECTRAL_IRRADIANCE_G,
		.spectral_irradiance_b = ATMO_SPECTRAL_IRRADIANCE_B,

		.exposure = ATMO_EXPOSURE,

		.spectral_to_rgb_r = ATMO_SPECTRAL_TO_RGB_R,
		.spectral_to_rgb_g = ATMO_SPECTRAL_TO_RGB_G,
		.spectral_to_rgb_b = ATMO_SPECTRAL_TO_RGB_B,

		.integration_steps = ATMO_INTEGRATION_STEPS,

		.k = ATMO_K,

		.texture_width  = ATMO_TEXTURE_WIDTH,
		.texture_height = ATMO_TEXTURE_HEIGHT,
		.texture_depth  = ATMO_TEXTURE_DEPTH,
	};

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
