#version 450

layout(location=0) in vec3 vertex;

layout(std140, set=0, binding=0) uniform ub000_mvp
{
	mat4 mvp;
};

layout(std140, set=0, binding=1) uniform ub001_RaRp
{
	vec2 RaRp; // Ra, Rp
};

void main()
{
	float Ra = RaRp[0];
	gl_Position = mvp*vec4(Ra*vertex, 1.0);
}
