#version 450

layout(location=0) in vec3 vertex;
layout(location=1) in vec3 normal;

layout(std140, set=0, binding=0) uniform ub000_mvp
{
	mat4 mvp;
};

layout(location=0) out vec3  varying_normal;

void main()
{
	varying_normal = normal;

	float Rp = 6371000.0;
	gl_Position = mvp*vec4(Rp*vertex, 1.0);
}
