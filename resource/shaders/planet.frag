#version 450

layout(location=0) in vec3 varying_normal;

layout(std140, set=0, binding=1) uniform ub001_L
{
	vec4 L;
};

layout(location=0) out vec4 fragColor;

void main()
{
	vec3 normal = normalize(varying_normal);
	vec3 sun    = -vec3(L);
	vec3 color  = vec3(0.11, 0.57, 1.0)*dot(normal, sun);
	fragColor = vec4(color, 1.0);
}
