#version 450

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

layout(std140, set=0, binding=3) uniform ub003_P0H
{
	vec4 P0H; // P0, h
};

int intersect_planet(vec3 P0, vec3 V, out vec3 normal)
{
	// https://www.perplexity.ai
	// libcc/docs/ray_sphere_intersect.c
	// libcc/math/cc_ray3f.c

	// normalize for a planet radius of 1
	float Rp = RaRp[1];
	vec3  p0 = P0/Rp;
	float a  = V.x*V.x + V.y*V.y + V.z*V.z;
	float b  = 2.0*(p0.x*V.x + p0.y*V.y + p0.z*V.z);
	float c  = p0.x*p0.x + p0.y*p0.y + p0.z*p0.z - 1.0;

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

void main()
{
	vec3 V  = normalize(varying_V);
	vec3 P0 = vec3(P0H);

	// check if ray intersects planet
	vec3 normal = vec3(0.0, 0.0, 1.0);
	if(intersect_planet(P0, V, normal) > 0)
	{
		vec3 Sun    = -vec3(L4);
		vec3 color  = vec3(0.565, 0.502, 0.439)*dot(normal, Sun);
		fragColor   = vec4(color, 1.0);
		return;
	}

	fragColor = vec4(0.53, 0.81, 0.98, 1.0);
}
