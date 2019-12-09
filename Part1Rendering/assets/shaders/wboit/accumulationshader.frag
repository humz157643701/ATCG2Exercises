#version 440 core
     
//define some constants
#define MAX_SHININESS 4000
#define INSULATOR_BASE_REFLECTIVITY 0.04f
#define PI 3.14159265359f
#define MAX_DIR_LIGHTS 2
#define MAX_POINT_LIGHTS 8
#define MAX_AMBIENT_LIGHTS 2
    
    
//some structs for cleaner code ----------------------------------------------------------------------------------------
struct Camera
{
    mat4 viewMatrix;
    mat4 projectionMatrix;
    vec3 worldPos;
    float renderWidth;
    float renderHeight;
    float nearPlane;
    float farPlane;
};
    
struct DirectionalLight
{
    float lumint;
    vec3 color;
    vec3 direction;
};
    
struct PointLight
{
    float lumint;
    vec3 color;
    vec3 position;
};

struct AmbientLight
{
    float lumint;
    vec3 color;
};
    
struct Material
{
    bool constDiff;
    sampler2D diff;
    vec4 cdiff;
    bool constSpec;
    sampler2D spec;
    vec4 cspec;
    bool constShin;
    sampler2D shin;
    float cshin;
    bool normalMap;
    sampler2D norm;    
    bool nmisbm;
};

in struct VertexData
{
    vec3 viewspacePosition;
    vec2 uv;
    vec3 normal;
    mat3 tbn;
} vertexData;
    
//uniforms -------------------------------------------------------------------------------------------------------------
//camera
uniform Camera camera;
    
//lights
uniform DirectionalLight dirlights[MAX_DIR_LIGHTS];
uniform int dirlightcount;

uniform PointLight pointlights[MAX_POINT_LIGHTS];
uniform int pointlightcount;

uniform AmbientLight ambientlights[MAX_AMBIENT_LIGHTS];
uniform int ambientlightcount;
    
//material
uniform Material material;
uniform int normalsfromtex;

uniform int usedepthfunc;

//particles special case
uniform bool isParticle = false;
uniform float particleExponent = 1.0;
in float particleLife;
    
//fragment shader output -----------------------------------------------------------------------------------------------
layout(location = 0) out vec4 weightedcolor;
layout(location = 1) out float alpha;
    
//forward declaration of some helper functions -------------------------------------------------------------------------

vec3 readNormal(vec2 tc, vec3 vnorm);

//retrieve light values
vec3 calculateDirectionalLightIntensity(int i);
vec3 calculatePointLightIntensity(vec3 P, int i);
vec3 calcAmbientLightIntensity(int i);
    
//fresnel approximation
vec3 fresnelSchlick(vec3 L, vec3 H, vec3 F0);
    
//functions applying the shading model
vec3 shade(vec3 N, vec3 L, vec3 V, vec3 diffc, float shn, vec3 F0, float metallic);
vec3 ambientShade(vec3 diffc, vec3 F0, float metallic);
float specToShininess(float spec);
    
//tone mapping and gamma correction
vec3 inverseGammaCorrect(vec3 c);

//depth weighting function
float depthweight(float z, float a, float s, float e);
    
//main function --------------------------------------------------------------------------------------------------------
void main()
{    
    //prepare vectors
    vec3 P = vertexData.viewspacePosition;    
    vec3 V = normalize(-P);
    vec2 tc = vertexData.uv;

    vec3 outcol = vec3(0.0f);  
   
    //read material properties
    vec4 matdiffdata = material.constDiff ? material.cdiff : texture(material.diff, tc);
    vec3 mat_diffuse = inverseGammaCorrect(matdiffdata.rgb);
    if(isParticle)
    {
        mat_diffuse *= material.cdiff.rgb;
    }
    float mat_alpha = matdiffdata.a;
    if(isParticle)
    {
        mat_alpha = pow(particleLife, particleExponent) * mat_alpha;
    }
    float mat_shininess = material.constShin ? specToShininess(material.cshin) : specToShininess(texture(material.shin, tc).r);
    vec3 N = material.normalMap ? readNormal(tc, normalize(vertexData.normal)) : normalize(vertexData.normal); 
    //color = vec4(N, 1.0f); return;
    vec3 F0 = material.constSpec ? material.cspec.rgb : texture(material.spec, tc).rgb;
    if(isParticle)
        F0 = vec3(0.0f, 0.0f, 0.0f);
    float mat_metallic = 1.0f - (material.constSpec ? material.cspec.a : texture(material.spec, tc).a);
    if(isParticle)
        mat_metallic = 0.0f;

    //directional lights
    for(int i = 0; i < dirlightcount; ++i)
    {
        vec3 L = normalize(-dirlights[i].direction);
        if(isParticle)
            N = L;
        vec3 intensity = calculateDirectionalLightIntensity(i);
        outcol += shade(N, L, V, mat_diffuse, mat_shininess, F0, mat_metallic) * intensity;
    }

    //point lights
    for(int i = 0; i < pointlightcount; ++i)
    {
        vec3 L = normalize(pointlights[i].position - P);
        if(isParticle)
            N = L;
        vec3 intensity = calculatePointLightIntensity(P, i);
        outcol += shade(N, L, V, mat_diffuse, mat_shininess, F0, mat_metallic) * intensity;
    }

    //ambient lights
    for(int i = 0; i < ambientlightcount; ++i)
    {
        outcol += ambientShade(mat_diffuse, F0, mat_metallic) * calcAmbientLightIntensity(i);
    }
    float weight = mat_alpha * depthweight(P.z, mat_alpha, camera.nearPlane, camera.farPlane);
    weightedcolor = vec4(outcol * weight, weight);
    alpha = mat_alpha;
}
    
//retrieve light values
vec3 calculateDirectionalLightIntensity(int i)
{
    return dirlights[i].lumint * dirlights[i].color;
}
    
vec3 calculatePointLightIntensity(vec3 fpos, int i)
{
    vec3 d = fpos - pointlights[i].position;
    return (pointlights[i].lumint * pointlights[i].color) * (1.0f / dot(d, d));
}
    
vec3 calcAmbientLightIntensity(int i)
{
    return ambientlights[i].lumint * ambientlights[i].color;
}

vec3 fresnelSchlick(vec3 L, vec3 H, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - max(dot(L, H), 0.0f), 5.0);
}
 
vec3 shade(vec3 N, vec3 L, vec3 V, vec3 diffc, float shn, vec3 F0, float metallic)
{
    vec3 H = normalize(V + L);
    float NdotL = max(dot(N, L), 0.0f);
    float NdotH = max(dot(N, H), 0.0f);
    vec3 F;
    if(!isParticle)
        F = fresnelSchlick(L, H, vec3(F0));
    else
        F = vec3(0.0f, 0.0f, 0.0f);
    return ((vec3(1.0f) - F) * (diffc / PI) * (1.0f - metallic) + F * pow(NdotH, shn) * ((shn + 8) / (8 * PI))) * NdotL;
}
 
vec3 ambientShade(vec3 diffc, vec3 F0, float metallic)
{
    vec3 ambientColor = diffc * vec3(1.0f - F0) * (1.0f - metallic);
    return ambientColor;
}
    
float specToShininess(float spec)
{
    const float f = sqrt(MAX_SHININESS - 1.0f);
    return  1.0f + (f * spec) * (f * spec);
}
    
vec3 inverseGammaCorrect(vec3 c)
{
    return pow(c, vec3(2.2f));
}

vec3 readNormal(vec2 tc, vec3 vnorm)
{
    if(material.nmisbm)
    {
        float height = texture(material.norm, tc).r;
        vec2 hgrad = vec2(dFdx(height), dFdy(height));
        vec3 normal = vec3(-hgrad, 1.0f);
        return normalize(vertexData.tbn * normal);
    }
    else
    {
        return normalize(vertexData.tbn * (texture(material.norm, tc).xyz * 2.0f - 1.0f));
    }
}

float depthweight(float z, float a, float s, float e)
{
    //rescale w.r.t s and e
    float absz = (abs(z) - abs(s)) * (500.0 / (abs(e) - abs(s)));
    return usedepthfunc == 1 ? a * max(0.01, min(3000.0, 10.0 / (0.00001 + (absz / 10.0) * (absz / 10.0) * (absz / 10.0) + (absz / 200.0) * (absz / 200.0) * (absz / 200.0) * (absz / 200.0) * (absz / 200.0) * (absz / 200.0)))) : 1.0f;
    //return usedepthfunc == 1 ? a * max(0.01, min(3000.0, 3000.0 - (3000.0 / 500.0) * absz)) : 1.0f;
}