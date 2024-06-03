#version 450

layout(location = 0) in vec3 position;

// ray positions and directions in world coordinates 
// to be passed to fragment shader
layout(location=0) out vec3 p;
layout(location=1) out vec3 d;

//push constants block
layout( push_constant ) uniform constants
{
	mat4 invView;
    	vec4 proj;
	float time;

} pushConstants;

void main() {
    
    float near = pushConstants.proj.x;
    float far = pushConstants.proj.y; 
    float aspect = pushConstants.proj.z;
    float fov = radians(0.5*pushConstants.proj.w);
    
    // position and direction of ray in world space for ray tracing
    // flip y for Vulkan
    vec3 ray = vec3(position.x*near*aspect*tan(fov), -position.y*near*tan(fov), -1.0*near);
    vec4 result = pushConstants.invView * vec4(ray, 1.0);
    p = result.xyz / result.w;

    result = pushConstants.invView * vec4(ray, 0.0);
    d = result.xyz;

    // trick the depth buffer
    gl_Position =  vec4( position, 1.0 );
}
