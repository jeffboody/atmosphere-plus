#version 450

layout(location=0) in vec3 varying_normal;

layout(location=0) out vec4 fragColor;

void main()
{
	vec3 normal = normalize(varying_normal);
	vec3 up     = vec3(0.0, 0.0, 1.0);
	vec3 color  = vec3(0.11f, 0.57f, 1.0f)*dot(normal, up);
	fragColor = vec4(color, 1.0f);
}
