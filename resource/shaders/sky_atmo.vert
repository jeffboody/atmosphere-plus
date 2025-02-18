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

layout(std140, set=1, binding=0) uniform ub100_P0H
{
	vec4 P0H; // P0, h
};

layout(location=0) out vec3 varying_V;

void main()
{
	float Ra    = RaRp[0];
	vec3  P0    = vec3(P0H);
	varying_V   = vertex - P0/Ra;
	gl_Position = mvp*vec4(Ra*vertex, 1.0);
}
