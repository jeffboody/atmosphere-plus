#version 450

// height parameterization
#define ATMO_PARAM_HEIGHT_LINEAR 0
#define ATMO_PARAM_HEIGHT_POWER  1

// view-zenith angle parameterization
#define ATMO_PARAM_PHI_LINEAR         0
#define ATMO_PARAM_PHI_POWER          1
#define ATMO_PARAM_PHI_WEIGHTED_POWER 3

// weighted power parameters
#define ATMO_PARAM_PHI_WEIGHTED_POWER_PU  2.0
#define ATMO_PARAM_PHI_WEIGHTED_POWER_PL  2.0
#define ATMO_PARAM_PHI_WEIGHTED_POWER_PS  2.0
#define ATMO_PARAM_PHI_WEIGHTED_POWER_WL0 (0.0/32.0)
#define ATMO_PARAM_PHI_WEIGHTED_POWER_WL1 (20.0/32.0)
#define ATMO_PARAM_PHI_WEIGHTED_POWER_WS0 (4.0/32.0)
#define ATMO_PARAM_PHI_WEIGHTED_POWER_WS1 (8.0/32.0)

// sun-zenith angle parameterization
#define ATMO_PARAM_DELTA_LINEAR 0
#define ATMO_PARAM_DELTA_POWER  1
#define ATMO_PARAM_DELTA_BODARE 2

// select parameterization
// requires corresponding change in atmo_solver.c
#define ATMO_PARAM_HEIGHT ATMO_PARAM_HEIGHT_POWER
#define ATMO_PARAM_PHI    ATMO_PARAM_PHI_WEIGHTED_POWER
#define ATMO_PARAM_DELTA  ATMO_PARAM_DELTA_POWER

// tone mapping
#define ATMO_TONE_MAPPING_REINHARD          0
#define ATMO_TONE_MAPPING_REINHARD_EXTENDED 1
#define ATMO_TONE_MAPPING_UNCHARTED2        2

// select tone mapping
#define ATMO_TONE_MAPPING ATMO_TONE_MAPPING_UNCHARTED2

layout(location=0) in vec3 varying_V;

layout(location=0) out float fragColor;

layout(std140, set=0, binding=1) uniform ub001_RaRp
{
	vec2 RaRp; // Ra, Rp
};

layout(std140, set=0, binding=2) uniform ub002_L4
{
	vec4 L4;
};

layout(std140, set=0, binding=3) uniform ub003_P0H
{
	vec4 P0H; // P0, h
};

layout(std140, set=1, binding=1) uniform ub101_Zenith4
{
	vec4 Zenith4;
};

layout(std140, set=1, binding=2) uniform ub102_IIE
{
	vec4 IIE; // II, Exposure
};

layout(set=1, binding=4) uniform sampler3D sampler104_fIS;
layout(set=1, binding=5) uniform sampler2D sampler105_T;

int intersect_planet(vec3 P0, vec3 V, out vec3 normal)
{
	// https://www.perplexity.ai
	// libcc/docs/ray_sphere_intersect.c
	// libcc/math/cc_ray3f.c

	// normalize for a planet radius of 1
	float Rp = RaRp[1];
	vec3  p0 = P0/Rp;
	float r  = 1.0;
	float a  = V.x*V.x + V.y*V.y + V.z*V.z;
	float b  = 2.0*(p0.x*V.x + p0.y*V.y + p0.z*V.z);
	float c  = p0.x*p0.x + p0.y*p0.y + p0.z*p0.z - r*r;

	float discriminant = b*b - 4*a*c;

	// check for ray-sphere intersection
	if (discriminant < 0)
	{
		return 0;
	}

	// compute the intersection distance
	float t1 = (-b - sqrt(discriminant))/(2.0*a);
	float t2 = (-b + sqrt(discriminant))/(2.0*a);

	// determine the nearest intersection
	if (t1 > t2)
	{
		float temp = t1;
		t1 = t2;
		t2 = temp;
	}

	// check the ray origin
	if (t1 < 0.0)
	{
		if (t2 < 0.0)
		{
			// ray origin is outside of the sphere but
			// ray does not have forward intersections
			return 0;
		}

		// ray origin is inside of the sphere and
		// ray has one forward intersection
		normal = normalize(p0);
		return 1;
	}

	// ray origin is outside of the sphere and
	// ray has two forward intersections
	normal = normalize(p0 + t1*V);
	return 2;
}

vec3 exposure(vec3 color)
{
	return pow(2.0, IIE.a)*color;
}

// https://64.github.io/tonemapping/
float luminance(vec3 v)
{
	return dot(v, vec3(0.2126, 0.7152, 0.0722));
}

// https://64.github.io/tonemapping/
vec3 change_luminance(vec3 c_in, float l_out)
{
	float l_in = luminance(c_in);
	return c_in * (l_out / l_in);
}

// https://64.github.io/tonemapping/
float reinhard(float x)
{
	return x / (1.0 + x);
}

// https://64.github.io/tonemapping/
vec3 reinhard_extended_luminance(vec3 v, float max_white_l)
{
	float l_old = luminance(v);
	float numerator = l_old * (1.0 + (l_old / (max_white_l * max_white_l)));
	float l_new = numerator / (1.0 + l_old);
	return change_luminance(v, l_new);
}

// https://64.github.io/tonemapping/
vec3 reinhard_luminance(vec3 v)
{
	float l0 = luminance(v);
	float l1 = reinhard(l0);
	return v*l1/l0;
}

// https://64.github.io/tonemapping/
vec3 uncharted2_tonemap_partial(vec3 x)
{
	float A = 0.15;
	float B = 0.50;
	float C = 0.10;
	float D = 0.20;
	float E = 0.02;
	float F = 0.30;
	return ((x*(A*x+C*B)+D*E)/(x*(A*x+B)+D*F))-E/F;
}

// https://64.github.io/tonemapping/
vec3 uncharted2_filmic(vec3 v)
{
    float exposure_bias = 2.0;
    vec3 curr = uncharted2_tonemap_partial(v * exposure_bias);

    vec3 W = vec3(11.2);
    vec3 white_scale = vec3(1.0) / uncharted2_tonemap_partial(W);
    return curr * white_scale;
}

vec3 gamma(vec3 v)
{
	return pow(v, vec3(1.0/2.2));
}

float getUHeight(float Ra, float Rp, float h)
{
	float u;
	#if ATMO_PARAM_HEIGHT == ATMO_PARAM_HEIGHT_POWER
	u = pow(h/(Ra - Rp), 1.0/2.0);
	#else
	u = h/(Ra - Rp);
	#endif
	return u;
}

float getVCosPhiPlanet(float Rp, float h, float cos_phi, float u)
{
	float v;
	#if ATMO_PARAM_PHI == ATMO_PARAM_PHI_POWER
	v = 0.5*(1.0 + sign(cos_phi)*pow(abs(cos_phi), 1.0/3.0));
	#elif ATMO_PARAM_PHI == ATMO_PARAM_PHI_WEIGHTED_POWER
	float hypH      = Rp + h;
	float oppH      = Rp;
	float adjH      = sqrt(hypH*hypH - oppH*oppH);
	float cos_phi_H = -adjH/hypH;

	float PU      = ATMO_PARAM_PHI_WEIGHTED_POWER_PU;
	float PL      = ATMO_PARAM_PHI_WEIGHTED_POWER_PL;
	float PS      = ATMO_PARAM_PHI_WEIGHTED_POWER_PS;
	float WL0     = ATMO_PARAM_PHI_WEIGHTED_POWER_WL0;
	float WL1     = ATMO_PARAM_PHI_WEIGHTED_POWER_WL1;
	float WS0     = ATMO_PARAM_PHI_WEIGHTED_POWER_WS0;
	float WS1     = ATMO_PARAM_PHI_WEIGHTED_POWER_WS1;
	float WS      = WS1*u + WS0*(1.0 - u);
	float WL      = WL1*u + WL0*(1.0 - u);
	float WU      = 1.0 - WL - WS;
	float epsilon = 0.00001;
	if(cos_phi >= 0.0)
	{
		v = WU*pow(cos_phi, 1.0/PU) + (WL + WS);
	}
	else if(cos_phi >= cos_phi_H)
	{
		v = WL*pow((cos_phi - cos_phi_H)/
		           max(-cos_phi_H, epsilon),
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
	return v;
}

float getVCosPhiSky(float Rp, float h, float cos_phi,
                    float u, vec3 V, out vec3 color)
{
	float v;
	#if ATMO_PARAM_PHI == ATMO_PARAM_PHI_POWER
	v = 0.5*(1.0 + sign(cos_phi)*pow(abs(cos_phi), 1.0/3.0));
	#elif ATMO_PARAM_PHI == ATMO_PARAM_PHI_WEIGHTED_POWER
	float hypH      = Rp + h;
	float oppH      = Rp;
	float adjH      = sqrt(hypH*hypH - oppH*oppH);
	float cos_phi_H = -adjH/hypH;

	float PU      = ATMO_PARAM_PHI_WEIGHTED_POWER_PU;
	float PL      = ATMO_PARAM_PHI_WEIGHTED_POWER_PL;
	float PS      = ATMO_PARAM_PHI_WEIGHTED_POWER_PS;
	float WL0     = ATMO_PARAM_PHI_WEIGHTED_POWER_WL0;
	float WL1     = ATMO_PARAM_PHI_WEIGHTED_POWER_WL1;
	float WS0     = ATMO_PARAM_PHI_WEIGHTED_POWER_WS0;
	float WS1     = ATMO_PARAM_PHI_WEIGHTED_POWER_WS1;
	float WS      = WS1*u + WS0*(1.0 - u);
	float WL      = WL1*u + WL0*(1.0 - u);
	float WU      = 1.0 - WL - WS;
	float epsilon = 0.00001;
	if(cos_phi >= 0.0)
	{
		v = WU*pow(cos_phi, 1.0/PU) + (WL + WS);

		#ifdef DEBUG_COS_PHI_SKY
		if(V.y < 0.0)
		{
			color = vec3(0.0, 0.0, 1.0);
			return -1.0;
		}
		#endif
	}
	else if(cos_phi >= cos_phi_H)
	{
		v = WL*pow((cos_phi - cos_phi_H)/
		           max(-cos_phi_H, epsilon),
		           1.0/PL) + WS;

		#ifdef DEBUG_COS_PHI_SKY
		if(V.y < 0.0)
		{
			color = vec3(1.0, 0.0, 0.0);
			return -1.0;
		}
		#endif
	}
	else
	{
		v = WS*(1.0 - pow((cos_phi - cos_phi_H)/
		                  (-1.0 - cos_phi_H), 1.0/PS));

		#ifdef DEBUG_COS_PHI_SKY
		if(V.y < 0.0)
		{
			color = vec3(0.0, 1.0, 0.0);
			return -1.0;
		}
		#endif
	}
	#else
	v = (cos_phi + 1.0)/2.0;
	#endif
	return v;
}

void main()
{
	vec3  V         = normalize(varying_V);
	vec3  L         = vec3(L4);
	vec3  Sun       = -vec3(L4);
	vec3  Zenith    = vec3(Zenith4);
	float cos_phi   = dot(V, Zenith);
	float cos_delta = dot(-L, Zenith);
	float Ra        = RaRp[0];
	float Rp        = RaRp[1];
	vec3  P0        = vec3(P0H);
	float h         = P0H[3];
	float PI        = 3.14159265358979323846;
	vec3  II        = vec3(IIE);

	// check if ray intersects planet
	vec3 color  = vec3(0.0, 0.0, 0.0);
	vec3 normal = vec3(0.0, 0.0, 1.0);
	vec3 debug  = vec3(0.0, 0.0, 0.0);
	if(intersect_planet(P0, V, normal) > 0)
	{
		float ndotsun = dot(normal, Sun);
		vec3  ambient = vec3(0.0, 0.0, 0.0);
		vec3  diffuse = vec3(0.0, 0.0, 0.0);
		float upc     = getUHeight(Ra, Rp, 0.0);
		float vpc     = getVCosPhiPlanet(Rp, 0.0, dot(normal, Sun), upc);
		float uap     = getUHeight(Ra, Rp, h);
		float vap     = getVCosPhiPlanet(Rp, h, cos_phi, uap);
		vec3  Tap     = vec3(texture(sampler105_T, vec2(uap, vap)));

		// apply diffuse lighting contribution
		if(ndotsun > 0.0)
		{
			vec3 Tpc    = vec3(texture(sampler105_T, vec2(upc, vpc)));
			vec3 albedo = vec3(0.3, 0.3, 0.3);
			diffuse = (albedo/PI)*clamp(0.0, 1.0, ndotsun)*Tpc;
		}

		// apply surface contribution
		color = (ambient + diffuse)*Tap;
	}

	// scattering texture coordinates
	float u = getUHeight(Ra, Rp, h);
	float v = getVCosPhiSky(Rp, h, cos_phi, u, V, debug);
	float w;
	#if ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_POWER
	w = 0.5*(1.0 + sign(cos_delta)*pow(abs(cos_delta), 1.0/3.0));
	#elif ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_BODARE
	w = 0.5*(atan(max(cos_delta, -0.1975)*tan(1.26*1.1))/1.1 +
	         (1.0 - 0.26));
	#else
	w = (cos_delta + 1.0)/2.0;
	#endif

	// apply sky and aerial perspective
	color += vec3(texture(sampler104_fIS, vec3(u, v, w)));

	// apply spectral irradiance
	color *= II;

	// convert to luminance
	fragColor = 0.2125*color.r + 0.7154*color.g + 0.0721*color.b;
}
