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

layout(location=0) out vec3  varying_vertex;

void main()
{
	varying_vertex = vertex;
	float Rp       = RaRp[1];
	gl_Position    = mvp*vec4(Rp*vertex, 1.0);
}
