#version 450

layout(location=0) in  vec3 vertex;
layout(location=1) in  vec3 V;
layout(location=0) out vec3 varying_V;

layout(std140, set=0, binding=0) uniform ub000_mvp
{
	mat4 mvp;
};

void main()
{
	varying_V   = V;
	gl_Position = mvp*vec4(vertex, 1.0);
}
