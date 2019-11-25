#version 430 core
#define IRRADIANCE_MAP_SAMPLES 1024
#define PI 3.14159265359
//output
layout (location = 0) out vec3 color;
layout (location = 1) out vec3 irradiance;

in vec3 worldPos;
uniform sampler2D u_enver;

float RadicalInverse_VdC(uint bits) 
{
    bits = (bits << 16u) | (bits >> 16u);
    bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
    bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
    bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
    bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
    return float(bits) * 2.3283064365386963e-10; // / 0x100000000
}

vec2 Hammersley(uint i, uint N)
{
    return vec2(float(i)/float(N), RadicalInverse_VdC(i));
}

const vec2 invAtan = vec2(0.1591, 0.3183);
vec2 sampleErMap(vec3 v)
{
    vec3 nv = normalize(v);
    vec2 uv = vec2(atan(nv.z, nv.x), asin(nv.y));
    uv *= invAtan;
    uv += 0.5;
    return uv;
}

vec3 fetchEnvColor(vec3 dir)
{   
    return texture(u_enver, sampleErMap(dir)).rgb;    
}

void main()
{		
    color = fetchEnvColor(worldPos);
    vec3 viewdir = normalize(worldPos);

    vec3 up = vec3(0.0, 1.0, 0.0);
    vec3 right = normalize(cross(up, viewdir));
    up = normalize(cross(viewdir, right));

    vec3 irrad = vec3(0.0, 0.0, 0.0);
    for(int i = 0; i < IRRADIANCE_MAP_SAMPLES; ++i)
    {
        vec2 usamp = Hammersley(i, IRRADIANCE_MAP_SAMPLES);
        // to spherical coordinates
        vec3 ssamp = vec3(
            cos(2 * PI * usamp.y) * 2.0 * sqrt(1.0 - usamp.x * usamp.x),
            sin(2 * PI * usamp.y) * 2.0 * sqrt(1.0 - usamp.x * usamp.x),
            usamp.x
        );

        vec3 wsamp = normalize(ssamp.x * right + ssamp.z * up + ssamp.y * viewdir); 

        // integrate
        irrad += texture(u_enver, sampleErMap(wsamp)).rgb * max(dot(wsamp, viewdir), 0.0);
    }

    irradiance = (irrad * PI) / (float(IRRADIANCE_MAP_SAMPLES));
}