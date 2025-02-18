#version 450

layout(location=0) in vec3 varying_vertex;

layout(std140, set=0, binding=2) uniform ub002_L4
{
	vec4 L4;
};

layout(location=0) out vec4 fragColor;

void main()
{
	vec3 normal = normalize(varying_vertex);
	vec3 Sun    = -vec3(L4);
	vec3 color  = vec3(0.11, 0.57, 1.0)*dot(normal, Sun);
	fragColor   = vec4(color, 1.0);
}
