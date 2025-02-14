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

#define LOG_TAG "mkatmo"
#include "libcc/math/cc_float.h"
#include "libcc/math/cc_sphere.h"
#include "libcc/math/cc_ray3f.h"
#include "libcc/math/cc_vec3d.h"
#include "libcc/math/cc_vec3f.h"
#include "libcc/math/cc_vec4f.h"
#include "libcc/cc_log.h"
#include "texgz/texgz_png.h"
#include "texgz/texgz_tex.h"

// radius of the planet and atmospheric boundary
// use doubles for radius for numerical stability
#define Rp 6371000.0
#define Ra 6471000.0

#define DENSITY_SCALE_HEIGHT_RAYLEIGH 8000.0f
#define DENSITY_SCALE_HEIGHT_MIE      1200.0f

#define PHASE_G_MIE -0.85f

#define BETA_R_RAYLEIGH 6.55e-6f
#define BETA_G_RAYLEIGH 1.73e-5f
#define BETA_B_RAYLEIGH 2.30e-5f
#define BETA_MIE        2e-6f

#define SPECTRAL_IRRADIANCE_R 0.1526f
#define SPECTRAL_IRRADIANCE_G 0.191f
#define SPECTRAL_IRRADIANCE_B 0.208f

#define EXPOSURE 1.0f

#define SPECTRAL_TO_RGB_R 133.3209f
#define SPECTRAL_TO_RGB_G 88.51855f
#define SPECTRAL_TO_RGB_B 112.7552f

#define INTEGRATION_STEPS 30

#define TEXTURE_SIZE 64

/***********************************************************
* private                                                  *
***********************************************************/

// Rayleigh density function
static float densityR(float h)
{
	return expf(-h/DENSITY_SCALE_HEIGHT_RAYLEIGH);
}

// Mie density function
static float densityM(float h)
{
	return expf(-h/DENSITY_SCALE_HEIGHT_MIE);
}

// modified Rayleigh phase function
static float phaseR(float cos_theta)
{
	return 0.8f*(1.4f + 0.5f*cos_theta*cos_theta);
}

// Mie phase function
static float phaseM(float cos_theta)
{
	float g  = PHASE_G_MIE;
	float g2 = g*g;
	float n1 = 3.0f*(1.0f - g2);
	float n2 = 1.0f + cos_theta*cos_theta;
	float d1 = 2.0f*(2.0f + g2);
	float d2 = powf(1.0f + g2 + 2.0f*g*cos_theta, 1.5f);
	return (n1/d1)*(n2/d2);
}

// compute height using double precision since the magnitude
// of points in the atmosphere produces very large numbers
static float height(cc_vec3f_t* P)
{
	ASSERT(P);

	cc_vec3d_t d =
	{
		.x = P->x,
		.y = P->y,
		.z = P->z,
	};

	return (float) (cc_vec3d_mag(&d) - Rp);
}

// Rayleigh/Mie transmittance
static void
transmittance(cc_vec3f_t* P1, cc_vec3f_t* P2, cc_vec4f_t* out)
{
	ASSERT(P1);
	ASSERT(P2);
	ASSERT(out);

	// initialize integration
	cc_vec3f_t P;
	cc_vec3f_t step;
	cc_vec3f_copy(P1, &P);
	cc_vec3f_subv_copy(P2, P1, &step);
	cc_vec3f_muls(&step, 1.0f/INTEGRATION_STEPS);
	float ds  = cc_vec3f_mag(&step);
	float h   = height(&P);
	float pR0 = densityR(h);
	float pM0 = densityM(h);

	// integrate transmittance
	int   i;
	float pR1;
	float pM1;
	float tR = 0.0f;
	float tM = 0.0f;
	for(i = 0; i < INTEGRATION_STEPS; ++i)
	{
		cc_vec3f_addv(&P, &step);

		h   = height(&P);
		pR1 = densityR(h);
		pM1 = densityM(h);

		// apply trapesoidal rule
		tR += 0.5f*(pR0 + pR1)*ds;
		tM += 0.5f*(pM0 + pM1)*ds;

		pR0 = pR1;
		pM0 = pM1;
	}

	// apply Rayleigh/Mie scattering coefficient
	out->r = BETA_R_RAYLEIGH*tR;
	out->g = BETA_G_RAYLEIGH*tR;
	out->b = BETA_B_RAYLEIGH*tR;
	out->a = BETA_MIE*tM;
}

// compute the points Pa and Pb on the viewing vector
static int
computePaPb(cc_vec3f_t* P0, cc_vec3f_t* V, double Ro,
            cc_vec3f_t* Pa, cc_vec3f_t* Pb)
{
	ASSERT(P0);
	ASSERT(V);
	ASSERT(Pa);
	ASSERT(Pb);

	// initialize spheres
	// optionally include a radius offset to ensure that the
	// viewing vector does not intersect at P0
	cc_sphere_t sphereP;
	cc_sphere_t sphereA;
	cc_sphere_load(&sphereP, 0.0f, 0.0f, 0.0f, Rp - Ro);
	cc_sphere_load(&sphereA, 0.0f, 0.0f, 0.0f, Ra + Ro);

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
computePc(cc_vec3f_t* P, cc_vec3f_t* L, cc_vec3f_t* Pc)
{
	ASSERT(P);
	ASSERT(L);
	ASSERT(Pc);

	cc_vec3f_t Pa;
	cc_vec3f_t Sun;
	cc_vec3f_muls_copy(L, -1.0f, &Sun);
	return computePaPb(P, &Sun, 0.0, &Pa, Pc);
}

// Rayleigh/Mie factored single-scattered intensity
static void
fIS1(float h, float phi, float delta, cc_vec4f_t* out)
{
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
		.z = h + Rp,
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
	computePaPb(&P0, &V, Ro, &Pa, &Pb);
	if(computePc(&Pa, &L, &Pc) == 0)
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
	cc_vec3f_muls(&step, 1.0f/INTEGRATION_STEPS);
	float ds = cc_vec3f_mag(&step);
	float pR = densityR(h);
	float pM = densityM(h);
	transmittance(&P, &Pc, &tPPc);
	transmittance(&Pa, &P, &tPaP);
	fx0.r = pR*exp(-tPPc.r -tPaP.r);
	fx0.g = pR*exp(-tPPc.g -tPaP.g);
	fx0.b = pR*exp(-tPPc.b -tPaP.b);
	fx0.a = pM*exp(-tPPc.a -tPaP.a);

	// integrate factored single-scattered intensity
	int i;
	for(i = 0; i < INTEGRATION_STEPS; ++i)
	{
		cc_vec3f_addv(&P, &step);

		if(computePc(&P, &L, &Pc) == 0)
		{
			// PPc intersects planet
			fx0.r = 0.0f;
			fx0.g = 0.0f;
			fx0.b = 0.0f;
			fx0.a = 0.0f;
			continue;
		}

		h  = height(&P);
		pR = densityR(h);
		pM = densityM(h);
		transmittance(&P, &Pc, &tPPc);
		transmittance(&Pa, &P, &tPaP);
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
	out->r *= BETA_R_RAYLEIGH/(4.0*M_PI);
	out->g *= BETA_G_RAYLEIGH/(4.0*M_PI);
	out->b *= BETA_B_RAYLEIGH/(4.0*M_PI);
	out->a *= BETA_MIE/(4.0*M_PI);
}

/***********************************************************
* public                                                   *
***********************************************************/

int main(int argc, const char** argv)
{
	texgz_tex_t* tex;
	tex = texgz_tex_new(TEXTURE_SIZE, TEXTURE_SIZE,
	                    TEXTURE_SIZE, TEXTURE_SIZE,
	                    TEXGZ_UNSIGNED_BYTE,
	                    TEXGZ_RGBA, NULL);
	if(tex == NULL)
	{
		return EXIT_FAILURE;
	}

	// spectral intensity of of incident light from the Sun
	cc_vec4f_t II =
	{
		.r = SPECTRAL_IRRADIANCE_R*EXPOSURE*SPECTRAL_TO_RGB_R,
		.g = SPECTRAL_IRRADIANCE_G*EXPOSURE*SPECTRAL_TO_RGB_G,
		.b = SPECTRAL_IRRADIANCE_B*EXPOSURE*SPECTRAL_TO_RGB_B,
	};

	int i;
	int j;
	for(i = 0; i < TEXTURE_SIZE; ++i)
	{
		for(j = 0; j < TEXTURE_SIZE; ++j)
		{
			float fi = (float) i;
			float fj = (float) j;
			float n1 = (float) (TEXTURE_SIZE - 1);
			float cos_phi   = cos((M_PI/2.0f)*fi/n1);
			float sin_phi   = sin((M_PI/2.0f)*fi/n1);
			float cos_delta = cos((M_PI/2.0f)*fj/n1);
			float sin_delta = sin((M_PI/2.0f)*fj/n1);
			cc_vec3f_t Sun    = { sin_delta, 0.0f, cos_delta };
			cc_vec3f_t V      = { sin_phi,   0.0f, cos_phi   };
			cc_vec3f_t Zenith = { 0.0f, 0.0f, 1.0f };
			cc_vec3f_t minusV;
			cc_vec3f_t L;
			cc_vec3f_normalize(&L);
			cc_vec3f_normalize(&V);
			cc_vec3f_normalize(&Zenith);
			cc_vec3f_muls_copy(&V, -1.0f, &minusV);
			cc_vec3f_muls_copy(&Sun, -1.0f, &L);

			float h         = 1.8f;
			float phi       = acos(cc_vec3f_dot(&V, &Zenith));
			float delta     = acos(cc_vec3f_dot(&Sun, &Zenith));
			float cos_theta = cc_vec3f_dot(&L, &minusV);

			cc_vec4f_t out;
			fIS1(h, phi, delta, &out);

			float FR = phaseR(cos_theta);
			float FM = phaseM(cos_theta);

			cc_vec4f_t IS1 =
			{
				.r = II.r*(FR*out.r + FM*out.a),
				.g = II.g*(FR*out.g + FM*out.a),
				.b = II.b*(FR*out.b + FM*out.a),
			};

			IS1.r = powf(1.0f - expf(-IS1.r), 1.0f/2.2f);
			IS1.g = powf(1.0f - expf(-IS1.g), 1.0f/2.2f);
			IS1.b = powf(1.0f - expf(-IS1.b), 1.0f/2.2f);
			IS1.r = cc_clamp(IS1.r, 0.0f, 1.0f);
			IS1.g = cc_clamp(IS1.g, 0.0f, 1.0f);
			IS1.b = cc_clamp(IS1.b, 0.0f, 1.0f);

			unsigned char pixel[4];
			pixel[0] = (unsigned char) (255.0*IS1.r);
			pixel[1] = (unsigned char) (255.0*IS1.g);
			pixel[2] = (unsigned char) (255.0*IS1.b);
			pixel[3] = 255;
			texgz_tex_setPixel(tex, j, i, pixel);

			LOGI("h=%f, phi=%f, delta=%f, cos_theta=%f, IS1: r=%f, g=%f, b=%f",
			     h, phi, delta, cos_theta,
			     IS1.r, IS1.g, IS1.b);
		}
	}

	texgz_png_export(tex, "atmo.png");
	texgz_tex_delete(&tex);

	return EXIT_SUCCESS;
}
