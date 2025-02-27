#version 450

// height parameterization
// LINEAR: u = h/Ha
// NONLINEAR: u = sqrt(h/Ha)
#define ATMO_PARAM_HEIGHT_LINEAR    0
#define ATMO_PARAM_HEIGHT_NONLINEAR 1

// view-azmuth angle parameterization
// LINEAR: v = (cos(phi) + 1)/2
// SHERVHEIM: v = 0.5*(1.0 + sign(cos_phi)*
//                     pow(abs(cos_phi), 1.0/3.0))
// BODARE: Efficient and Dynamic Atmospheric Scattering
//         WARNING: ATMO_PARAM_PHI_BODARE is buggy
#define ATMO_PARAM_PHI_LINEAR    0
#define ATMO_PARAM_PHI_SHERVHEIM 1
#define ATMO_PARAM_PHI_BODARE    2

// sun-azmuth angle parameterization
// linear:    w = (cos(delta) + 1)/2
// SHERVHEIM: w = 0.5*(1.0 + sign(cos_delta)*
//                     pow(abs(cos_delta), 1.0/3.0))
// BODARE: Efficient and Dynamic Atmospheric Scattering
#define ATMO_PARAM_DELTA_LINEAR    0
#define ATMO_PARAM_DELTA_SHERVHEIM 1
#define ATMO_PARAM_DELTA_BODARE    2

// select parameterization
// requires corresponding change in atmo_solver.c
#define ATMO_PARAM_HEIGHT ATMO_PARAM_HEIGHT_NONLINEAR
#define ATMO_PARAM_PHI    ATMO_PARAM_PHI_SHERVHEIM
#define ATMO_PARAM_DELTA  ATMO_PARAM_DELTA_SHERVHEIM

layout(location=0) in vec3 varying_V;

layout(location=0) out vec4 fragColor;

layout(std140, set=0, binding=1) uniform ub001_RaRp
{
	vec2 RaRp; // Ra, Rp
};

layout(std140, set=0, binding=2) uniform ub002_L4
{
	vec4 L4;
};

layout(std140, set=1, binding=0) uniform ub100_P0H
{
	vec4 P0H; // P0, h
};

layout(std140, set=1, binding=1) uniform ub101_Zenith4
{
	vec4 Zenith4;
};

layout(std140, set=1, binding=2) uniform ub102_II4
{
	vec4 II4;
};

layout(std140, set=1, binding=3) uniform ub103_phase_g_mie
{
	float phase_g_mie;
};

layout(set=1, binding=4) uniform sampler3D sampler104_fIS;

// modified Rayleigh phase function
float phaseR(float cos_theta)
{
	return 0.8*(1.4 + 0.5*cos_theta*cos_theta);
}

// Mie phase function
float phaseM(float cos_theta)
{
	float g  = phase_g_mie;
	float g2 = g*g;
	float n1 = 3.0*(1.0 - g2);
	float n2 = 1.0 + cos_theta*cos_theta;
	float d1 = 2.0*(2.0 + g2);
	float d2 = pow(1.0 + g2 + 2.0*g*cos_theta, 1.5);
	return (n1/d1)*(n2/d2);
}

void main()
{
	vec3  V         = normalize(varying_V);
	vec3  L         = vec3(L4);
	vec3  Zenith    = vec3(Zenith4);
	float cos_phi   = dot(V, Zenith);
	float cos_delta = dot(-L, Zenith);
	float cos_theta = dot(L, -V);
	float FR        = phaseR(cos_theta);
	float FM        = phaseM(cos_theta);
	float Ra        = RaRp[0];
	float Rp        = RaRp[1];
	float h         = P0H[3];

	float u;
	#if ATMO_PARAM_HEIGHT == ATMO_PARAM_HEIGHT_NONLINEAR
	u  = sqrt(h/(Ra - Rp));
	#else
	u  = h/(Ra - Rp);
	#endif

	float v;
	#if ATMO_PARAM_PHI == ATMO_PARAM_PHI_SHERVHEIM
	v = 0.5*(1.0 + sign(cos_phi)*pow(abs(cos_phi), 1.0/3.0));
	#elif ATMO_PARAM_PHI == ATMO_PARAM_PHI_BODARE
	float ch = -sqrt(h*(2.0*Rp + h))/(Rp + h);
	if(cos_phi > ch)
	{
		v = 0.5*pow((cos_phi - ch)/(1.0 - ch), 0.2) + 0.5;
	}
	else
	{
		v = 0.5*pow((ch - cos_phi)/(1.0 + ch), 0.2);
	}
	#else
	v = (cos_phi + 1.0)/2.0;
	#endif

	float w;
	#if ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_SHERVHEIM
	w = 0.5*(1.0 + sign(cos_delta)*pow(abs(cos_delta), 1.0/3.0));
	#elif ATMO_PARAM_DELTA == ATMO_PARAM_DELTA_BODARE
	w = 0.5*(atan(max(cos_delta, -0.1975)*tan(1.26*1.1))/1.1 +
	         (1.0 - 0.26));
	#else
	w = (cos_delta + 1.0)/2.0;
	#endif

	// sample fIS texture
	vec4 fIS = texture(sampler104_fIS, vec3(u, v, w));

	// apply constant phase function and
	// spectral intensity of incident light
	vec3 II = vec3(II4);
	vec3 IS = vec3(II.r*(FR*fIS.r + FM*fIS.a),
	               II.g*(FR*fIS.g + FM*fIS.a),
	               II.b*(FR*fIS.b + FM*fIS.a));

	// apply tone mapping and gamma correction
	vec3 color = pow(vec3(1.0) - exp(-IS), vec3(1.0/2.2));

	fragColor = vec4(color, 1.0);
}
