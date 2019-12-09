#version 430 core

//output
layout (location = 0) out vec4 color;

in vec3 cmCoords;

uniform samplerCube u_skybox;

uniform float exposure;

vec3 TM_Exposure(vec3 l);
vec3 gammaCorrect(vec3 c);
float getLuminance(vec3 linearRGB);
vec3 inverseGammaCorrect(vec3 c);

void main()
{
    color = vec4(gammaCorrect(TM_Exposure(texture(u_skybox, cmCoords).rgb)), 1.0f); 
}

vec3 TM_Exposure(vec3 l)
{
    float il = getLuminance(l);
    float ol = 1.0 - exp(- (il * exposure));
    return l * (ol / il);
}
    
vec3 gammaCorrect(vec3 c)
{
    return pow(c, vec3(1.0f / 2.2f));
}

vec3 inverseGammaCorrect(vec3 c)
{
    return pow(c, vec3(2.2f));
}
    
float getLuminance(vec3 linearRGB)
{
    return dot(vec3(0.2126729,  0.7151522,  0.0721750), gammaCorrect(linearRGB));
}